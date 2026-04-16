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
        self._link_index_cache = {}
        for e in range(E):
            self._link_index_cache[(network.link_from[e], network.link_to[e])] = e

        # 疎行列の構造をキャッシュ (data だけ差し替えて再利用)
        N = network.num_nodes
        self._graph_row = network.link_from
        self._graph_col = network.link_to
        self._graph_shape = (N, N)

        # 起点別需要の事前グループ化
        self._demand_by_origin = {}
        for r in network.origins:
            self._demand_by_origin[r] = [
                (d, vol) for (o, d), vol in network.demand.items() if o == r
            ]

    def _build_graph(self, link_cost):
        """コスト付き疎行列を構築 (構造キャッシュ使用)"""
        data = np.maximum(link_cost, 1e-15)
        return csr_matrix((data, (self._graph_row, self._graph_col)),
                          shape=self._graph_shape)

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

        for r in net.origins:
            dist, pred = self._shortest_path(link_cost, r)
            local_flow = np.zeros(E)
            for d, vol in self._demand_by_origin[r]:
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

        total_cost = np.sum(link_flow * link_cost)
        gap = self._compute_duality_gap(link_flow, link_cost)
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

        # === メインループ ===
        for iteration in range(1, max_iter + 1):
            if time.time() - start_time > time_limit:
                if verbose:
                    print("Time limit reached.")
                break

            # --- Step 1: PAS識別 + 即座のフローシフト ---
            cs = max(rel_gap, 1e-10)
            new_pas_count = 0

            for r in net.origins:
                dist, pred = self._shortest_path(link_cost, r)
                local_flow = link_flow_by_origin[r]

                # ポテンシャルリンク (縮約コスト降順)
                potential = []
                eps_cs = epsilon * cs
                theta_cs = theta * cs
                for e in range(E):
                    if local_flow[e] <= eps_cs:
                        continue
                    rc = dist[net.link_from[e]] - dist[net.link_to[e]] + link_cost[e]
                    if rc > theta_cs:
                        potential.append((rc, e))

                potential.sort(reverse=True)

                for _, e0 in potential:
                    if local_flow[e0] <= eps_cs:
                        continue

                    pas = self._mfs(r, e0, local_flow, dist, pred,
                                    epsilon, link_flow)
                    if pas is None:
                        continue

                    seg1, seg2 = pas
                    key = (seg1, seg2)

                    # PAS登録 (既存なら起点を更新)
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
            total_cost = np.sum(link_flow * link_cost)
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

    def _shortest_path(self, link_cost, origin):
        graph = self._build_graph(link_cost)
        dist, pred = sp_dijkstra(graph, directed=True, indices=origin,
                                  return_predecessors=True)
        return dist, pred

    def _find_link(self, from_node, to_node):
        return self._link_index_cache.get((from_node, to_node))

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

        # 逆方向トレース
        visited = {i: 0}
        backward = [e0]
        current = i

        for _ in range(net.num_nodes * 2):
            # current への最大フロー入リンク
            best_e = -1
            best_f = -1.0
            for e in net.in_links.get(current, []):
                if local_flow[e] > best_f:
                    best_f = local_flow[e]
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
                # サイクル除去
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
        """PAS上のNewtonステップによるフローシフト"""
        seg1, seg2 = pas_key

        f1 = min((local_flow[e] for e in seg1), default=0.0)
        f2 = min((local_flow[e] for e in seg2), default=0.0)

        if f1 < epsilon and f2 < epsilon:
            return None

        t1 = sum(link_cost[e] for e in seg1)
        t2 = sum(link_cost[e] for e in seg2)

        if not init_flag:
            if abs(t2 - t1) < mu * (t2 + t1):
                return None

        dt1 = sum(link_cost_de[e] for e in seg1)
        dt2 = sum(link_cost_de[e] for e in seg2)
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

        for e in seg1:
            link_flow[e] += delta
            local_flow[e] += delta
        for e in seg2:
            link_flow[e] -= delta
            local_flow[e] -= delta

        net = self.net
        for e in set(seg1) | set(seg2):
            f = max(link_flow[e], 0.0)
            link_cost[e] = self.cost_func(
                f, net.link_fftt[e], net.link_capacity[e],
                net.link_alpha[e], net.link_beta[e])
            link_cost_de[e] = self.cost_deriv(
                f, net.link_fftt[e], net.link_capacity[e],
                net.link_alpha[e], net.link_beta[e])

        return (delta,)

    def _compute_duality_gap(self, link_flow, link_cost):
        """双対ギャップ = Σ(flow*cost) - Σ(demand*shortest_path_cost)"""
        net = self.net
        total_tt = np.sum(link_flow * link_cost)
        sp_tt = 0.0
        for r in net.origins:
            dist, _ = self._shortest_path(link_cost, r)
            for d, vol in self._demand_by_origin[r]:
                sp_tt += vol * dist[d]
        return total_tt - sp_tt
