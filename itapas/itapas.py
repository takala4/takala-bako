"""iTAPAS (incremental Traffic Assignment by Paired Alternative Segments)

hanqiu92/itapas の参考実装に忠実な再実装 (graph-tool 非依存)。

References:
    Bar-Gera, H. (2010). Traffic assignment by paired alternative segments.
        Transportation Research Part B, 44(8-9), 1022-1046.
    Xie, J. & Xie, C. (2016). New insights and improvements of using paired
        alternative segments for traffic assignment.
        Transportation Research Part B, 93, 406-424.
"""

import time
import numpy as np
from scipy.sparse import csr_matrix, coo_matrix
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

        self._link_index_cache = {}
        for e in range(E):
            self._link_index_cache[(network.link_from[e], network.link_to[e])] = e

        # CSR 構造キャッシュ
        coo = coo_matrix((np.arange(E, dtype=float),
                          (network.link_from, network.link_to)), shape=(N, N))
        csr = coo.tocsr()
        self._csr_data_map = csr.data.astype(int)
        self._csr_indices = csr.indices.copy()
        self._csr_indptr = csr.indptr.copy()
        self._csr_shape = (N, N)

        # 起点別需要
        self._demand_by_origin = {}
        for r in network.origins:
            dests, vols = [], []
            for (o, d), vol in network.demand.items():
                if o == r:
                    dests.append(d)
                    vols.append(vol)
            self._demand_by_origin[r] = (np.array(dests, dtype=int),
                                          np.array(vols, dtype=float))

        # 入リンク配列 (MFS用)
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

    def _build_graph(self, link_cost):
        data = np.maximum(link_cost[self._csr_data_map], 1e-15)
        return csr_matrix((data, self._csr_indices, self._csr_indptr),
                          shape=self._csr_shape, copy=False)

    def _all_shortest_paths(self, link_cost):
        g = self._build_graph(link_cost)
        origins = np.array(self.net.origins, dtype=int)
        dist, pred = sp_dijkstra(g, directed=True, indices=origins,
                                  return_predecessors=True)
        return {r: (dist[i], pred[i]) for i, r in enumerate(origins)}

    def _shortest_path(self, link_cost, origin):
        g = self._build_graph(link_cost)
        dist, pred = sp_dijkstra(g, directed=True, indices=origin,
                                  return_predecessors=True)
        return dist, pred

    def solve(self, max_iter=50, gap_threshold=1e-6, time_limit=1200,
              epsilon=1e-12, theta=1e-12, mu=1e-3,
              inner_iterations=20, verbose=True):
        """iTAPASアルゴリズムによる均衡配分の実行 (参考実装準拠)"""
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
        gap = self._compute_duality_gap(link_flow, link_cost)
        rel_gap = gap / total_cost if total_cost > 0 else float('inf')

        log = [{'iter': 0, 'gap': gap, 'rel_gap': rel_gap,
                'time': time.time() - start_time, 'num_pas': 0}]
        if verbose:
            print(f"Iter  0: rel_gap={rel_gap:.6e}, "
                  f"gap={gap:.2f}, time={log[-1]['time']:.2f}s")

        if rel_gap < gap_threshold:
            return link_flow, log

        # PASセット: {(seg1, seg2): origin} — 参考実装と同じ単一起点
        pas_set = {}

        link_from = net.link_from
        link_to = net.link_to

        # === メインループ ===
        for iteration in range(1, max_iter + 1):
            if time.time() - start_time > time_limit:
                if verbose:
                    print("Time limit reached.")
                break

            new_pas_count = 0

            # --- Step 1: 各起点についてPAS識別 ---
            for r in net.origins:
                dist, pred = self._shortest_path(link_cost, r)
                local_flow = link_flow_by_origin[r]

                # 起点別補完スラック (参考実装準拠)
                rc = dist[link_from] - dist[link_to] + link_cost
                comp_slack = float(np.dot(local_flow, rc))

                if comp_slack > 1e-8:
                    mask = (local_flow > epsilon * comp_slack) & \
                           (rc > theta * comp_slack)
                    pot_indices = np.where(mask)[0]

                    # ポテンシャルリンクを処理 (即座シフトなし、PAS登録のみ)
                    for e0 in pot_indices:
                        if local_flow[e0] <= epsilon:
                            continue
                        e0_tuple = (int(link_from[e0]), int(link_to[e0]))

                        pas = self._mfs(link_flow, local_flow, link_cost,
                                        link_cost_de, pred, r, e0_tuple,
                                        epsilon)
                        if pas is not None:
                            if pas not in pas_set:
                                pas_set[pas] = r
                                new_pas_count += 1

                # ランダムシフト (起点ループ内、参考実装準拠)
                if pas_set:
                    pas_keys = list(pas_set.keys())
                    if len(pas_keys) > 100:
                        import random
                        pas_keys = [pas_keys[i] for i in
                                    np.random.choice(len(pas_keys), 100, replace=False)]
                    for pas in pas_keys:
                        self._shift(link_flow, link_flow_by_origin[pas_set[pas]],
                                    link_cost, link_cost_de, pas, epsilon, mu, True)

            # --- Step 2: フローシフト (内部イテレーション) ---
            for _ in range(inner_iterations):
                to_remove = []
                for pas in list(pas_set.keys()):
                    result = self._shift(link_flow, link_flow_by_origin[pas_set[pas]],
                                         link_cost, link_cost_de, pas,
                                         epsilon, mu, False)
                    if not result:
                        to_remove.append(pas)
                for pas in to_remove:
                    pas_set.pop(pas, None)

            # --- Step 3: 収束判定 ---
            link_cost = self._compute_cost(link_flow)
            link_cost_de = self._compute_cost_deriv(link_flow)
            total_cost = float(np.dot(link_flow, link_cost))
            gap = self._compute_duality_gap(link_flow, link_cost)
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

    def _mfs(self, flow_map, local_flow_map, cost_map, cost_de_map,
             prev_map, r, e0, epsilon):
        """Maximum Flow Search (参考実装準拠)

        Parameters
        ----------
        flow_map : ndarray — 集計リンクフロー
        local_flow_map : ndarray — 起点別リンクフロー
        cost_map, cost_de_map : ndarray — リンクコスト・導関数
        prev_map : ndarray — SPT predecessor
        r : int — 起点ノード
        e0 : tuple (i, j) — ポテンシャルリンクのノードペア
        epsilon : float — 最小フロー閾値
        """
        net = self.net
        i, j = e0
        e0_idx = self._link_index_cache.get(e0)

        # 特殊ケース: e0の始点が起点の場合
        if i == r:
            seg1 = self._get_spt_path(prev_map, r, j)
            if seg1 is None:
                return None
            return (seg1, (e0_idx,))

        # SPTパス上のノード集合と predecessor
        spt_nodes = set()
        node = j
        while node != r and node >= 0:
            spt_nodes.add(node)
            node = prev_map[node]
        if node == r:
            spt_nodes.add(r)

        in_ptr = self._in_ptr
        in_edges = self._in_edges

        while True:
            # prev_new: 逆方向トレースの predecessor マップ
            prev_new = {}
            status = {}  # -1: SPTパス上, 0: 未訪問, 1: トレース済み

            # SPTパス上のノードに status=-1 を設定
            for v in spt_nodes:
                status[v] = -1

            # 初期設定
            prev_new[j] = i
            status[j] = 1
            status[i] = 1

            current = i
            found = False

            while True:
                # current への最大フロー入リンク
                start = in_ptr[current]
                end = in_ptr[current + 1]
                best_e = -1
                best_f = -1.0
                best_m = -1
                for idx in range(start, end):
                    e = in_edges[idx]
                    f = local_flow_map[e]
                    if f > best_f:
                        best_f = f
                        best_e = e
                        best_m = net.link_from[e]

                if best_e < 0:
                    return None

                m = best_m
                prev_new[current] = m

                m_status = status.get(m, 0)

                if m_status == -1:
                    # SPTパスに到達 → PAS発見
                    # cheap segment: SPT上で m → j
                    seg1 = self._get_spt_path(prev_map, m, j)
                    if seg1 is None:
                        return None

                    # expensive segment: prev_new を辿って m → j
                    seg2_links = []
                    node = j
                    while node != m:
                        p = prev_new.get(node)
                        if p is None:
                            return None
                        e = self._link_index_cache.get((p, node))
                        if e is None:
                            return None
                        seg2_links.append(e)
                        node = p
                    seg2_links.reverse()
                    seg2 = tuple(seg2_links)

                    return (seg1, seg2)

                elif m_status == 1:
                    # サイクル検出 → サイクルフロー除去
                    cycle_links = []
                    curr_m = prev_new[m]
                    cycle_links.append(self._link_index_cache.get((curr_m, m)))
                    prev_m = prev_new.get(curr_m)
                    while curr_m != m:
                        cycle_links.append(self._link_index_cache.get((prev_m, curr_m)))
                        curr_m = prev_m
                        prev_m = prev_new.get(curr_m)

                    cycle_links = [e for e in cycle_links if e is not None]
                    if cycle_links:
                        delta = min(local_flow_map[e] for e in cycle_links)
                        for e in cycle_links:
                            local_flow_map[e] -= delta
                            flow_map[e] -= delta
                        # コスト更新
                        n = self.net
                        for e in cycle_links:
                            f = max(flow_map[e], 0.0)
                            ratio = f / n.link_capacity[e]
                            rb = ratio ** n.link_beta[e]
                            cost_map[e] = n.link_fftt[e] * (1.0 + n.link_alpha[e] * rb)
                            cost_de_map[e] = n.link_fftt[e] * n.link_alpha[e] * n.link_beta[e] * rb / f if f > 1e-20 else 0.0

                    if local_flow_map[e0_idx] < epsilon:
                        return None

                    # 外側ループに戻ってリトライ
                    break
                else:
                    current = m
                    status[m] = 1

    def _shift(self, link_flow, local_flow, link_cost, link_cost_de,
               pas, epsilon, mu, init_flag):
        """PAS上のNewtonステップによるフローシフト (参考実装準拠)

        Returns
        -------
        tuple or None
            tuple: PAS維持 (フローシフト実行 or 両セグメント空)
            None: PAS削除 (均衡化済み or 片側フロー枯渇)
        """
        seg1, seg2 = pas

        # ボトルネックフロー・コスト (Python loop, 小配列高速)
        f1 = float('inf')
        f2 = float('inf')
        t1 = 0.0
        t2 = 0.0
        for e in seg1:
            lf = local_flow[e]
            if lf < f1:
                f1 = lf
            t1 += link_cost[e]
        for e in seg2:
            lf = local_flow[e]
            if lf < f2:
                f2 = lf
            t2 += link_cost[e]

        if f1 < epsilon and f2 < epsilon:
            return (0.0, 0.0)  # 両セグメント空 → PAS維持

        # 参考実装準拠: f1<ε OR f2<ε OR 均衡 → PAS削除 (init時は除く)
        if not init_flag:
            if f1 < epsilon or f2 < epsilon or abs(t2 - t1) < mu * (t2 + t1):
                return None

        # Newton step
        dt1 = 0.0
        dt2 = 0.0
        for e in seg1:
            dt1 += link_cost_de[e]
        for e in seg2:
            dt2 += link_cost_de[e]

        denom = dt1 + dt2
        if denom < 1e-20:
            return (f1, f2)

        delta = (t2 - t1) / denom

        if delta > 0.0:
            if delta > f2:
                delta = f2
        else:
            if delta < -f1:
                delta = -f1

        # フローシフト + コスト更新
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
            rb = ratio ** beta[e]
            link_cost[e] = fftt[e] * (1.0 + alpha[e] * rb)
            link_cost_de[e] = fftt[e] * alpha[e] * beta[e] * rb / f if f > 1e-20 else 0.0
        for e in seg2:
            link_flow[e] -= delta
            local_flow[e] -= delta
            f = link_flow[e] if link_flow[e] > 0.0 else 0.0
            ratio = f / cap[e]
            rb = ratio ** beta[e]
            link_cost[e] = fftt[e] * (1.0 + alpha[e] * rb)
            link_cost_de[e] = fftt[e] * alpha[e] * beta[e] * rb / f if f > 1e-20 else 0.0

        return (f1 + delta, f2 - delta)

    def _compute_duality_gap(self, link_flow, link_cost):
        """双対ギャップ = Σ(flow*cost) - Σ(demand*shortest_path_cost)"""
        net = self.net
        total_tt = float(np.dot(link_flow, link_cost))
        sp_tt = 0.0
        sp_data = self._all_shortest_paths(link_cost)
        for r in net.origins:
            dist, _ = sp_data[r]
            dests, vols = self._demand_by_origin[r]
            sp_tt += float(np.dot(vols, dist[dests]))
        return total_tt - sp_tt
