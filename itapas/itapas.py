"""iTAPAS (incremental Traffic Assignment by Paired Alternative Segments)

References:
    Bar-Gera, H. (2010). Traffic assignment by paired alternative segments.
        Transportation Research Part B, 44(8-9), 1022-1046.
    Xie, J. & Xie, C. (2016). New insights and improvements of using paired
        alternative segments for traffic assignment.
        Transportation Research Part B, 93, 406-424.
"""

import time
import random
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra as sp_dijkstra


class ITAPAS:
    """iTAPASによる静的利用者均衡配分ソルバー"""

    def __init__(self, network, cost_func, cost_deriv, cost_integral):
        self.net = network
        self.cost_func = cost_func
        self.cost_deriv = cost_deriv
        self.cost_integral = cost_integral

        E = network.num_links
        N = network.num_nodes

        # リンクインデックス検索キャッシュ
        self._link_index_cache = {}
        for e in range(E):
            self._link_index_cache[(network.link_from[e], network.link_to[e])] = e

        # CSR構造の事前構築 (data だけ差し替えて再利用)
        dummy = np.ones(E)
        self._graph_template = csr_matrix(
            (dummy, (network.link_from, network.link_to)), shape=(N, N)
        )

        # 起点別需要の事前グループ化
        self._demand_by_origin = {}
        for r in network.origins:
            dests, vols = [], []
            for (o, d), vol in network.demand.items():
                if o == r:
                    dests.append(d)
                    vols.append(vol)
            self._demand_by_origin[r] = (np.array(dests, dtype=int),
                                          np.array(vols, dtype=float))

        # 入リンクを numpy 配列化 (MFS高速化)
        # in_link_ptr[node] ~ in_link_ptr[node+1] が node への入リンク範囲
        in_counts = np.zeros(N, dtype=int)
        for e in range(E):
            in_counts[network.link_to[e]] += 1
        self._in_ptr = np.zeros(N + 1, dtype=int)
        np.cumsum(in_counts, out=self._in_ptr[1:])
        self._in_edges = np.empty(E, dtype=int)
        pos = self._in_ptr[:-1].copy()
        for e in range(E):
            to = network.link_to[e]
            self._in_edges[pos[to]] = e
            pos[to] += 1

        # PAS セグメントを numpy 配列に変換するキャッシュ
        self._seg_arrays = {}

    def _get_seg_array(self, seg):
        """タプルを numpy 配列に変換 (キャッシュ付き)"""
        arr = self._seg_arrays.get(seg)
        if arr is None:
            arr = np.array(seg, dtype=int)
            self._seg_arrays[seg] = arr
        return arr

    def _build_graph(self, link_cost):
        """コスト付き疎行列を構築 (構造キャッシュ使用)"""
        g = self._graph_template.copy()
        g.data[:] = np.maximum(link_cost[self._graph_template.data != 0] if False else link_cost, 1e-15)
        # 直接 data を上書き (構造が同一なので安全)
        # ただし csr_matrix のコンストラクタ経由で data 順序が変わりうるので再構築
        net = self.net
        data = np.maximum(link_cost, 1e-15)
        g = csr_matrix((data, (net.link_from, net.link_to)),
                        shape=(net.num_nodes, net.num_nodes))
        return g

    def _all_shortest_paths(self, link_cost):
        """全起点からの最短経路を一括計算"""
        g = self._build_graph(link_cost)
        origins = np.array(self.net.origins, dtype=int)
        dist, pred = sp_dijkstra(g, directed=True, indices=origins,
                                  return_predecessors=True)
        # dist[i] = origin_list[i] からの距離, pred[i] = predecessor
        return {r: (dist[i], pred[i]) for i, r in enumerate(origins)}

    def solve(self, max_iter=50, gap_threshold=1e-6, time_limit=1200,
              epsilon=1e-14, theta=1e-14, mu=1e-3,
              inner_iterations=20, verbose=True):
        """iTAPASアルゴリズムによる均衡配分の実行"""
        net = self.net
        E = net.num_links
        start_time = time.time()

        # === Step 0: 初期化 (All-or-Nothing) ===
        link_flow = np.zeros(E)
        link_flow_by_origin = {}

        link_cost = self._compute_cost(link_flow)
        sp_data = self._all_shortest_paths(link_cost)

        for r in net.origins:
            dist, pred = sp_data[r]
            local_flow = np.zeros(E)
            dests, vols = self._demand_by_origin[r]
            for idx in range(len(dests)):
                d, vol = int(dests[idx]), vols[idx]
                node = d
                while node != r:
                    p = pred[node]
                    if p < 0:
                        break
                    e = self._link_index_cache.get((p, node))
                    if e is not None:
                        local_flow[e] += vol
                    node = p
            link_flow_by_origin[r] = local_flow
            link_flow += local_flow

        link_cost = self._compute_cost(link_flow)
        link_cost_de = self._compute_cost_deriv(link_flow)

        total_cost = float(np.dot(link_flow, link_cost))
        sp_data = self._all_shortest_paths(link_cost)
        gap = self._compute_duality_gap_from_sp(link_flow, link_cost, sp_data)
        rel_gap = gap / total_cost if total_cost > 0 else float('inf')

        log = [{'iter': 0, 'gap': gap, 'rel_gap': rel_gap,
                'time': time.time() - start_time, 'num_pas': 0}]
        if verbose:
            print(f"Iter  0: rel_gap={rel_gap:.6e}, "
                  f"gap={gap:.2f}, time={log[-1]['time']:.2f}s")

        if rel_gap < gap_threshold:
            return link_flow, log

        # PASセット: {(seg1, seg2): origin}
        pas_set = {}

        # ポテンシャルリンク検出用ベクトル
        link_from = net.link_from
        link_to = net.link_to

        # === メインループ ===
        for iteration in range(1, max_iter + 1):
            if time.time() - start_time > time_limit:
                if verbose:
                    print("Time limit reached.")
                break

            # --- Step 1: PAS識別 + 即座のフローシフト ---
            cs = max(rel_gap, 1e-10)
            new_pas_count = 0
            eps_cs = epsilon * cs
            theta_cs = theta * cs

            sp_data = self._all_shortest_paths(link_cost)

            for r in net.origins:
                dist, pred = sp_data[r]
                local_flow = link_flow_by_origin[r]

                # ベクトル化: ポテンシャルリンク検出
                rc = dist[link_from] - dist[link_to] + link_cost
                mask = (local_flow > eps_cs) & (rc > theta_cs)
                pot_indices = np.where(mask)[0]

                if len(pot_indices) == 0:
                    continue

                # 縮約コスト降順ソート
                pot_rc = rc[pot_indices]
                order = np.argsort(-pot_rc)
                pot_indices = pot_indices[order]

                for e0 in pot_indices:
                    if local_flow[e0] <= eps_cs:
                        continue

                    pas = self._mfs(r, int(e0), local_flow, dist, pred,
                                    epsilon, link_flow)
                    if pas is None:
                        continue

                    seg1, seg2 = pas
                    key = (seg1, seg2)
                    pas_set[key] = r
                    new_pas_count += 1

                    # 即座にフローシフト
                    self._shift(key, link_flow, local_flow,
                                link_cost, link_cost_de, epsilon, mu,
                                init_flag=True)

            # ランダムサブセットに追加シフト
            pas_keys = list(pas_set.keys())
            if len(pas_keys) > 100:
                pas_keys = random.sample(pas_keys, 100)
            for key in pas_keys:
                r = pas_set[key]
                self._shift(key, link_flow, link_flow_by_origin[r],
                            link_cost, link_cost_de, epsilon, mu,
                            init_flag=True)

            # --- Step 2: フローシフトフェーズ (単一起点, 原著準拠) ---
            for _ in range(inner_iterations):
                to_remove = []
                for key in list(pas_set.keys()):
                    r = pas_set[key]
                    result = self._shift(key, link_flow, link_flow_by_origin[r],
                                         link_cost, link_cost_de, epsilon, mu,
                                         init_flag=False)
                    if result is None:
                        to_remove.append(key)
                for key in to_remove:
                    pas_set.pop(key, None)

            # --- Step 3: 収束判定 ---
            link_cost = self._compute_cost(link_flow)
            link_cost_de = self._compute_cost_deriv(link_flow)
            total_cost = float(np.dot(link_flow, link_cost))
            sp_data = self._all_shortest_paths(link_cost)
            gap = self._compute_duality_gap_from_sp(link_flow, link_cost, sp_data)
            rel_gap = gap / total_cost if total_cost > 0 else float('inf')

            elapsed = time.time() - start_time
            log.append({'iter': iteration, 'gap': gap, 'rel_gap': rel_gap,
                        'time': elapsed, 'num_pas': len(pas_set)})

            if verbose:
                print(f"Iter {iteration:3d}: rel_gap={rel_gap:.6e}, "
                      f"gap={gap:.2f}, PAS={len(pas_set)}, "
                      f"new={new_pas_count}, time={elapsed:.2f}s")

            if rel_gap < gap_threshold:
                if verbose:
                    print("Converged.")
                break

        return link_flow, log

    def _compute_cost(self, flow):
        net = self.net
        return self.cost_func(flow, net.link_fftt, net.link_capacity,
                              net.link_alpha, net.link_beta)

    def _compute_cost_deriv(self, flow):
        net = self.net
        return self.cost_deriv(flow, net.link_fftt, net.link_capacity,
                               net.link_alpha, net.link_beta)

    def _shortest_path(self, link_cost, origin):
        g = self._build_graph(link_cost)
        dist, pred = sp_dijkstra(g, directed=True, indices=origin,
                                  return_predecessors=True)
        return dist, pred

    def _get_spt_path(self, pred, origin, target):
        """SPT上のorigin→targetの経路をリンクインデックスのタプルで返す"""
        path = []
        node = target
        while node != origin and node >= 0:
            p = pred[node]
            if p < 0:
                return None
            e = self._link_index_cache.get((p, node))
            if e is None:
                return None
            path.append(e)
            node = p
        path.reverse()
        return tuple(path)

    def _mfs(self, origin, e0, local_flow, dist, pred, epsilon, link_flow):
        """Maximum Flow Search によるPAS識別"""
        net = self.net
        i = net.link_from[e0]
        j = net.link_to[e0]

        # SPT上のorigin→jパスのノード集合
        spt_nodes = set()
        node = j
        while node != origin and node >= 0:
            spt_nodes.add(node)
            node = pred[node]
        if node == origin:
            spt_nodes.add(origin)

        # 逆方向トレース (配列ベースの入リンク参照)
        visited = {i: 0}
        backward = [e0]
        current = i
        in_ptr = self._in_ptr
        in_edges = self._in_edges

        for _ in range(net.num_nodes * 2):
            # current への最大フロー入リンク (配列参照)
            start = in_ptr[current]
            end = in_ptr[current + 1]
            if start == end:
                return None

            best_e = -1
            best_f = -1.0
            for idx in range(start, end):
                e = in_edges[idx]
                f = local_flow[e]
                if f > best_f:
                    best_f = f
                    best_e = e

            if best_e < 0 or best_f < epsilon:
                return None

            prev = net.link_from[best_e]

            if prev in spt_nodes:
                backward.append(best_e)
                seg2 = tuple(reversed(backward))
                seg1 = self._get_spt_path(pred, prev, j)
                if seg1 is None or len(seg1) == 0:
                    return None
                return (seg1, seg2)

            if prev in visited:
                cpos = visited[prev]
                cycle = list(backward[cpos + 1:]) + [best_e]
                if cycle:
                    cf = min(local_flow[e] for e in cycle)
                    if cf > epsilon:
                        for e in cycle:
                            local_flow[e] -= cf
                            link_flow[e] -= cf
                if local_flow[e0] < epsilon:
                    return None
                visited = {i: 0}
                backward = [e0]
                current = i
                continue

            visited[prev] = len(backward)
            backward.append(best_e)
            current = prev

        return None

    def _shift(self, pas_key, link_flow, local_flow, link_cost, link_cost_de,
               epsilon, mu, init_flag):
        """PAS上のNewtonステップによるフローシフト (numpy 配列インデックス版)"""
        seg1, seg2 = pas_key
        s1 = self._get_seg_array(seg1)
        s2 = self._get_seg_array(seg2)

        # ボトルネックフロー
        f1 = float(local_flow[s1].min()) if len(s1) > 0 else 0.0
        f2 = float(local_flow[s2].min()) if len(s2) > 0 else 0.0

        if f1 < epsilon and f2 < epsilon:
            return None

        # セグメントコスト
        t1 = float(link_cost[s1].sum())
        t2 = float(link_cost[s2].sum())

        if not init_flag:
            if abs(t2 - t1) < mu * (t2 + t1):
                return None

        # Newton step
        dt1 = float(link_cost_de[s1].sum())
        dt2 = float(link_cost_de[s2].sum())
        denom = dt1 + dt2
        if denom < 1e-20:
            return (0.0,)

        delta = (t2 - t1) / denom

        if delta > 0:
            delta = min(delta, f2)
        else:
            delta = max(delta, -f1)

        if abs(delta) < epsilon:
            return (0.0,)

        # フローシフト (numpy 一括操作)
        link_flow[s1] += delta
        local_flow[s1] += delta
        link_flow[s2] -= delta
        local_flow[s2] -= delta

        # 影響リンクのコスト更新
        affected = np.union1d(s1, s2)
        net = self.net
        fa = np.maximum(link_flow[affected], 0.0)
        link_cost[affected] = self.cost_func(
            fa, net.link_fftt[affected], net.link_capacity[affected],
            net.link_alpha[affected], net.link_beta[affected])
        link_cost_de[affected] = self.cost_deriv(
            fa, net.link_fftt[affected], net.link_capacity[affected],
            net.link_alpha[affected], net.link_beta[affected])

        return (delta,)

    def _compute_duality_gap_from_sp(self, link_flow, link_cost, sp_data):
        """事前計算済み最短経路データから双対ギャップを計算"""
        total_tt = float(np.dot(link_flow, link_cost))
        sp_tt = 0.0
        for r in self.net.origins:
            dist, _ = sp_data[r]
            dests, vols = self._demand_by_origin[r]
            sp_tt += float(np.dot(vols, dist[dests]))
        return total_tt - sp_tt

    def _compute_duality_gap(self, link_flow, link_cost):
        """双対ギャップ = Σ(flow*cost) - Σ(demand*shortest_path_cost)"""
        sp_data = self._all_shortest_paths(link_cost)
        return self._compute_duality_gap_from_sp(link_flow, link_cost, sp_data)
