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

        self._link_index_cache = {}
        for e in range(network.num_links):
            self._link_index_cache[(network.link_from[e], network.link_to[e])] = e

    def solve(self, max_iter=50, gap_threshold=1e-6, time_limit=1200,
              epsilon=1e-14, theta=1e-14, mu=1e-3,
              inner_iterations=20, verbose=True):
        """iTAPASアルゴリズムによる均衡配分の実行

        Returns
        -------
        link_flow : ndarray
            均衡リンクフロー
        log : list of dict
            各イテレーションの収束情報
        """
        net = self.net
        E = net.num_links
        N = net.num_nodes
        start_time = time.time()

        # === Step 0: 初期化 (All-or-Nothing) ===
        link_flow = np.zeros(E)
        link_flow_by_origin = {}

        link_cost = self._compute_cost(link_flow)

        for r in net.origins:
            dist, predecessors = self._shortest_path(link_cost, r)
            local_flow = np.zeros(E)
            for (o, d), vol in net.demand.items():
                if o != r:
                    continue
                node = d
                while node != r:
                    pred = predecessors[node]
                    if pred < 0:
                        break
                    e = self._find_link(pred, node)
                    if e is not None:
                        local_flow[e] += vol
                    node = pred
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
                  f"gap={gap:.2f}, time={log[-1]['time']:.1f}s")

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
                dist, predecessors = self._shortest_path(link_cost, r)
                local_flow = link_flow_by_origin[r]

                # ポテンシャルリンク特定 (縮約コスト基準)
                potential_links = []
                for e in range(E):
                    if local_flow[e] <= epsilon * cs:
                        continue
                    i_node = net.link_from[e]
                    j_node = net.link_to[e]
                    rc = dist[i_node] - dist[j_node] + link_cost[e]
                    if rc > theta * cs:
                        potential_links.append((rc, e))

                # 縮約コスト降順でソート (最も非均衡なリンクから処理)
                potential_links.sort(reverse=True)

                for _, e0 in potential_links:
                    if local_flow[e0] <= epsilon * cs:
                        continue

                    pas = self._mfs(r, e0, local_flow, dist, predecessors,
                                    epsilon, link_flow)
                    if pas is None:
                        continue

                    seg1, seg2 = pas
                    key = (seg1, seg2)

                    if key not in pas_set:
                        pas_set[key] = r
                        new_pas_count += 1

                    # 即座にフローシフト (参考実装に準拠)
                    self._shift(key, link_flow, link_flow_by_origin[pas_set[key]],
                                link_cost, link_cost_de, epsilon, mu, init_flag=True)

            # ランダムサブセットにも追加シフト
            if len(pas_set) > 100:
                for key in random.sample(list(pas_set.keys()), 100):
                    r = pas_set[key]
                    self._shift(key, link_flow, link_flow_by_origin[r],
                                link_cost, link_cost_de, epsilon, mu, init_flag=True)

            # --- Step 2: フローシフトフェーズ (全起点でシフト) ---
            for _ in range(inner_iterations):
                keys_to_remove = []
                for key in list(pas_set.keys()):
                    any_shifted = False
                    for r in net.origins:
                        lf = link_flow_by_origin[r]
                        # この起点がPASセグメント上にフローを持つか確認
                        seg1, seg2 = key
                        has_flow = any(lf[e] > epsilon for e in seg2) or \
                                   any(lf[e] > epsilon for e in seg1)
                        if not has_flow:
                            continue
                        result = self._shift(key, link_flow, lf,
                                             link_cost, link_cost_de,
                                             epsilon, mu, init_flag=False)
                        if result is not None:
                            any_shifted = True
                    if not any_shifted:
                        keys_to_remove.append(key)
                for key in keys_to_remove:
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
                      f"new={new_pas_count}, time={elapsed:.1f}s")

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
        net = self.net
        N = net.num_nodes
        data = np.maximum(link_cost, 1e-15)
        graph = csr_matrix((data, (net.link_from, net.link_to)), shape=(N, N))
        dist, predecessors = sp_dijkstra(graph, directed=True, indices=origin,
                                         return_predecessors=True)
        return dist, predecessors

    def _find_link(self, from_node, to_node):
        return self._link_index_cache.get((from_node, to_node))

    def _get_spt_path(self, predecessors, origin, target):
        """SPT上のorigin→targetの経路をリンクインデックスのタプルで返す"""
        path_links = []
        node = target
        while node != origin and node >= 0:
            pred = predecessors[node]
            if pred < 0:
                return None
            e = self._find_link(pred, node)
            if e is None:
                return None
            path_links.append(e)
            node = pred
        path_links.reverse()
        return tuple(path_links)

    def _mfs(self, origin, e0, local_flow, dist, predecessors,
             epsilon, link_flow):
        """Maximum Flow Search によるPAS識別

        ポテンシャルリンク e0 の tail から逆方向に最大フロー入リンクを辿り、
        SPTパス上のノードに到達したらPASとして返す。
        """
        net = self.net
        i = net.link_from[e0]  # tail
        j = net.link_to[e0]    # head

        # SPT上のorigin→jパスのノード集合
        spt_path_nodes = set()
        node = j
        while node != origin and node >= 0:
            spt_path_nodes.add(node)
            node = predecessors[node]
        if node == origin:
            spt_path_nodes.add(origin)

        # 逆方向トレース
        max_steps = net.num_nodes * 2
        visited = {i: 0}
        # backward_edges: 参考実装と同じ構成
        # backward_edges[0] = e0 (i→j)
        # backward_edges[k] = link(prev_k → prev_{k-1}) (逆方向トレースで発見)
        backward_edges = [e0]
        current = i

        for _ in range(max_steps):
            # current への最大フロー入リンク
            best_link = -1
            best_flow = -1.0
            for e in net.in_links.get(current, []):
                if local_flow[e] > best_flow:
                    best_flow = local_flow[e]
                    best_link = e

            if best_link < 0 or best_flow < epsilon:
                return None

            prev_node = net.link_from[best_link]

            if prev_node in spt_path_nodes:
                # SPTパスに到達 → PAS発見
                backward_edges.append(best_link)
                # expensive segment: backward_edges を逆順にする
                seg2 = tuple(reversed(backward_edges))

                # cheap segment: SPT上でprev_node→j
                seg1 = self._get_spt_path(predecessors, prev_node, j)
                if seg1 is None or len(seg1) == 0:
                    return None

                return (seg1, seg2)

            if prev_node in visited:
                # サイクル検出 → サイクルフロー除去
                cycle_pos = visited[prev_node]
                cycle_links = list(backward_edges[cycle_pos + 1:]) + [best_link]

                if cycle_links:
                    cycle_flow = min(local_flow[e] for e in cycle_links)
                    if cycle_flow > epsilon:
                        for e in cycle_links:
                            local_flow[e] -= cycle_flow
                            link_flow[e] -= cycle_flow

                if local_flow[e0] < epsilon:
                    return None

                # リトライ
                visited = {i: 0}
                backward_edges = [e0]
                current = i
                continue

            visited[prev_node] = len(backward_edges)
            backward_edges.append(best_link)
            current = prev_node

        return None

    def _shift(self, pas_key, link_flow, local_flow, link_cost, link_cost_de,
               epsilon, mu, init_flag):
        """PAS上のNewtonステップによるフローシフト"""
        seg1, seg2 = pas_key

        # セグメント上のボトルネックフロー (local = origin-specific)
        f1 = min((local_flow[e] for e in seg1), default=0.0)
        f2 = min((local_flow[e] for e in seg2), default=0.0)

        if f1 < epsilon and f2 < epsilon:
            return None

        # セグメントコスト
        t1 = sum(link_cost[e] for e in seg1)
        t2 = sum(link_cost[e] for e in seg2)

        # 均衡判定 (init_flagがFalseの場合のみ)
        if not init_flag:
            if abs(t2 - t1) < mu * (t2 + t1):
                return None

        # Newton step
        dt1 = sum(link_cost_de[e] for e in seg1)
        dt2 = sum(link_cost_de[e] for e in seg2)
        denom = dt1 + dt2
        if denom < 1e-20:
            return (0.0,)

        delta = (t2 - t1) / denom

        # フロー非負制約
        if delta > 0:
            delta = min(delta, f2)
        else:
            delta = max(delta, -f1)

        if abs(delta) < epsilon:
            return (0.0,)

        # フローシフト実行
        for e in seg1:
            link_flow[e] += delta
            local_flow[e] += delta
        for e in seg2:
            link_flow[e] -= delta
            local_flow[e] -= delta

        # 影響リンクのコスト更新
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
            for (o, d), vol in net.demand.items():
                if o == r:
                    sp_tt += vol * dist[d]
        return total_tt - sp_tt
