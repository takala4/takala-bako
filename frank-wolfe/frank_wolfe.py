"""Frank-Wolfe 法による静的利用者均衡配分

Wardrop の第一原則 (利用者均衡) を Beckmann の等価最小化問題として定式化し、
Frank-Wolfe (Conditional Gradient) 法で解く。

References:
    LeBlanc, L. J., Morlok, E. K., & Piercey, W. P. (1975).
        An efficient approach to solving the road network equilibrium
        traffic assignment problem.
        Transportation Research, 9(5), 309-318.
    Sheffi, Y. (1985). Urban Transportation Networks.
        Prentice-Hall, Englewood Cliffs, NJ.
"""

import time
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra as sp_dijkstra


class FrankWolfe:
    """Frank-Wolfe 法による利用者均衡配分ソルバー

    Parameters
    ----------
    network : Network
        ネットワークデータ
    cost_func : callable
        リンクコスト関数 cost(flow, fftt, capacity, alpha, beta)
    cost_deriv : callable
        リンクコスト導関数
    cost_integral : callable
        リンクコスト積分 (Beckmann 目的関数用)
    """

    def __init__(self, network, cost_func, cost_deriv, cost_integral):
        self.net = network
        self.cost_func = cost_func
        self.cost_deriv = cost_deriv
        self.cost_integral = cost_integral

        self._link_index_cache = {}
        for e in range(network.num_links):
            self._link_index_cache[(network.link_from[e], network.link_to[e])] = e

    def solve(self, max_iter=100, gap_threshold=1e-4, verbose=True):
        """Frank-Wolfe 法による均衡配分の実行

        Parameters
        ----------
        max_iter : int
            最大イテレーション数
        gap_threshold : float
            相対双対ギャップの収束閾値
        verbose : bool
            進捗表示

        Returns
        -------
        link_flow : ndarray
            均衡リンクフロー
        log : list of dict
            各イテレーションの収束情報
        """
        net = self.net
        E = net.num_links
        start_time = time.time()

        # === Step 0: 初期化 (All-or-Nothing) ===
        link_cost = self._compute_cost(np.zeros(E))
        link_flow = self._aon_assign(link_cost)

        link_cost = self._compute_cost(link_flow)
        total_cost = np.sum(link_flow * link_cost)
        gap = self._compute_duality_gap(link_flow, link_cost)
        rel_gap = gap / total_cost if total_cost > 0 else float('inf')

        log = [{'iter': 0, 'gap': gap, 'rel_gap': rel_gap,
                'time': time.time() - start_time}]
        if verbose:
            print(f"Iter  0: rel_gap={rel_gap:.6e}, gap={gap:.2f}, "
                  f"time={log[-1]['time']:.2f}s")

        # === メインループ ===
        for iteration in range(1, max_iter + 1):
            # Step 1: 現在コストでの AON 配分 → 補助フロー
            link_cost = self._compute_cost(link_flow)
            aux_flow = self._aon_assign(link_cost)

            # Step 2: 探索方向
            direction = aux_flow - link_flow

            # Step 3: 最適ステップサイズの探索 (二分法)
            step_size = self._line_search(link_flow, direction)

            # Step 4: フロー更新
            link_flow = link_flow + step_size * direction

            # Step 5: 収束判定
            link_cost = self._compute_cost(link_flow)
            total_cost = np.sum(link_flow * link_cost)
            gap = self._compute_duality_gap(link_flow, link_cost)
            rel_gap = gap / total_cost if total_cost > 0 else float('inf')

            elapsed = time.time() - start_time
            log.append({'iter': iteration, 'gap': gap, 'rel_gap': rel_gap,
                        'time': elapsed, 'step_size': step_size})

            if verbose:
                print(f"Iter {iteration:3d}: rel_gap={rel_gap:.6e}, "
                      f"gap={gap:.2f}, step={step_size:.6f}, "
                      f"time={elapsed:.2f}s")

            if rel_gap < gap_threshold:
                if verbose:
                    print("Converged.")
                break

        return link_flow, log

    def _compute_cost(self, flow):
        net = self.net
        return self.cost_func(flow, net.link_fftt, net.link_capacity,
                              net.link_alpha, net.link_beta)

    def _compute_cost_integral(self, flow):
        net = self.net
        return self.cost_integral(flow, net.link_fftt, net.link_capacity,
                                  net.link_alpha, net.link_beta)

    def _beckmann_objective(self, flow):
        """Beckmann 目的関数 Z(f) = Σ ∫_0^{f_a} t_a(x) dx"""
        return np.sum(self._compute_cost_integral(flow))

    def _shortest_path(self, link_cost, origin):
        net = self.net
        N = net.num_nodes
        data = np.maximum(link_cost, 1e-15)
        graph = csr_matrix((data, (net.link_from, net.link_to)), shape=(N, N))
        dist, predecessors = sp_dijkstra(graph, directed=True, indices=origin,
                                         return_predecessors=True)
        return dist, predecessors

    def _aon_assign(self, link_cost):
        """All-or-Nothing 配分: 全 OD 需要を最短経路に割り当て"""
        net = self.net
        E = net.num_links
        aon_flow = np.zeros(E)

        for r in net.origins:
            dist, predecessors = self._shortest_path(link_cost, r)
            for (o, d), vol in net.demand.items():
                if o != r:
                    continue
                node = d
                while node != r:
                    pred = predecessors[node]
                    if pred < 0:
                        break
                    e = self._link_index_cache.get((pred, node))
                    if e is not None:
                        aon_flow[e] += vol
                    node = pred

        return aon_flow

    def _line_search(self, flow, direction, tol=1e-8, max_bisect=50):
        """二分法による最適ステップサイズの探索

        min_{λ∈[0,1]} Z(flow + λ * direction)
        一階条件: Σ t_a(f + λd) * d_a = 0
        """
        lo, hi = 0.0, 1.0

        for _ in range(max_bisect):
            if hi - lo < tol:
                break
            mid = (lo + hi) / 2.0
            trial_flow = flow + mid * direction
            trial_cost = self._compute_cost(trial_flow)
            # 勾配: Σ t_a(f + λd) * d_a
            grad = np.sum(trial_cost * direction)
            if grad < 0:
                lo = mid
            else:
                hi = mid

        return (lo + hi) / 2.0

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
