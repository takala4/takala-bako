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

        # CSR 構造キャッシュ: data のマッピングを事前計算
        # coo → csr 変換時のインデックス並び替えを記録
        from scipy.sparse import coo_matrix
        coo = coo_matrix((np.arange(E, dtype=float),
                          (network.link_from, network.link_to)),
                         shape=(N, N))
        csr = coo.tocsr()
        # csr.data には coo の data が並び替えられて入っている
        # csr.data[k] = link_index → self._csr_data_map[k] = link_index
        self._csr_data_map = csr.data.astype(int)
        self._csr_indices = csr.indices.copy()
        self._csr_indptr = csr.indptr.copy()
        self._csr_shape = (N, N)

    def _build_graph(self, link_cost):
        """コスト付き疎行列を構築 (構造キャッシュ使用, data のみ差し替え)"""
        data = np.maximum(link_cost[self._csr_data_map], 1e-15)
        return csr_matrix((data, self._csr_indices, self._csr_indptr),
                          shape=self._csr_shape, copy=False)

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
        """PAS上のNewtonステップによるフローシフト

        小配列 (2-10要素) に対して numpy 呼出オーバーヘッドを回避し、
        Python ループで直接計算する。
        """
        seg1, seg2 = pas_key

        # ボトルネックフロー・セグメントコスト・導関数を一括計算 (Python loop)
        f1 = float('inf')
        f2 = float('inf')
        t1 = 0.0
        t2 = 0.0
        dt1 = 0.0
        dt2 = 0.0
        for e in seg1:
            lf = local_flow[e]
            if lf < f1:
                f1 = lf
            t1 += link_cost[e]
            dt1 += link_cost_de[e]
        for e in seg2:
            lf = local_flow[e]
            if lf < f2:
                f2 = lf
            t2 += link_cost[e]
            dt2 += link_cost_de[e]

        if f1 < epsilon and f2 < epsilon:
            return None

        if not init_flag:
            if abs(t2 - t1) < mu * (t2 + t1):
                return None

        denom = dt1 + dt2
        if denom < 1e-20:
            return 0.0  # PAS を維持

        delta = (t2 - t1) / denom

        if delta > 0.0:
            if delta > f2:
                delta = f2
        else:
            if delta < -f1:
                delta = -f1

        if abs(delta) < epsilon:
            return 0.0  # PAS を維持

        # フローシフト + コスト更新 (Python loop, 関数呼出を最小化)
        net = self.net
        fftt = net.link_fftt
        cap = net.link_capacity
        alpha = net.link_alpha
        beta = net.link_beta
        for e in seg1:
            link_flow[e] += delta
            local_flow[e] += delta
            f = link_flow[e] if link_flow[e] > 0.0 else 0.0
            ratio = f / cap[e]
            ratio_b = ratio ** beta[e]
            link_cost[e] = fftt[e] * (1.0 + alpha[e] * ratio_b)
            link_cost_de[e] = fftt[e] * alpha[e] * beta[e] * ratio_b / f if f > 1e-20 else 0.0
        for e in seg2:
            link_flow[e] -= delta
            local_flow[e] -= delta
            f = link_flow[e] if link_flow[e] > 0.0 else 0.0
            ratio = f / cap[e]
            ratio_b = ratio ** beta[e]
            link_cost[e] = fftt[e] * (1.0 + alpha[e] * ratio_b)
            link_cost_de[e] = fftt[e] * alpha[e] * beta[e] * ratio_b / f if f > 1e-20 else 0.0

        return delta

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
