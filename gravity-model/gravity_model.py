"""重力モデルによるOD需要推計とパラメータ推定

交通工学における古典的な重力モデルのキャリブレーション。
観測トリップ長分布 (TLFD) または平均トリップコストに合致するよう
抵抗関数のパラメータを推定する。

重力モデル (二重制約):
    T_ij = a_i * b_j * O_i * D_j * f(c_ij)
    a_i, b_j は行和・列和制約を満たすバランシング係数 (Furness法)

推定手法:
    Hyman法: 平均トリップコストの一致による逐次更新
    TLFD法:  トリップ長頻度分布の適合 (最小二乗 or カイ二乗)

References:
    Hyman, G. M. (1969). The calibration of trip distribution models.
        Environment and Planning, 1, 105-112.
    Ortúzar, J. de D. & Willumsen, L. G. (2011).
        Modelling Transport (4th ed.). Wiley.
    Wilson, A. G. (1970). Entropy in Urban and Regional Modelling. Pion.
"""

import numpy as np
from scipy.optimize import minimize_scalar


def impedance_exp(cost, beta):
    """指数型抵抗関数: f(c) = exp(-β * c)"""
    return np.exp(-beta * cost)


def impedance_power(cost, gamma):
    """べき乗型抵抗関数: f(c) = c^(-γ)"""
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(cost > 0, cost ** (-gamma), 0.0)


def impedance_combined(cost, beta, gamma):
    """複合型抵抗関数: f(c) = c^(-γ) * exp(-β * c)"""
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(cost > 0, cost ** (-gamma) * np.exp(-beta * cost), 0.0)


class GravityModel:
    """交通工学における重力モデル

    Parameters
    ----------
    O : array-like (n_zones,)
        生成交通量
    D : array-like (n_zones,)
        集中交通量
    cost_matrix : array-like (n_zones, n_zones)
        ゾーン間コスト行列 (時間, 距離等)
    impedance : str or callable
        抵抗関数 ('exp', 'power', 'combined', or callable)
    """

    def __init__(self, O, D, cost_matrix, impedance='exp'):
        self.O = np.asarray(O, dtype=float)
        self.D = np.asarray(D, dtype=float)
        self.cost = np.asarray(cost_matrix, dtype=float)
        self.n_zones = len(self.O)

        if impedance == 'exp':
            self.impedance_func = impedance_exp
        elif impedance == 'power':
            self.impedance_func = impedance_power
        elif impedance == 'combined':
            self.impedance_func = impedance_combined
        elif callable(impedance):
            self.impedance_func = impedance
        else:
            raise ValueError(f"Unknown impedance type: {impedance}")
        self.impedance_type = impedance

    def furness(self, F, max_iter=100, tol=1e-6):
        """Furness法 (二重制約バランシング)

        Parameters
        ----------
        F : ndarray (n_zones, n_zones)
            抵抗関数値の行列 (対角要素は0)

        Returns
        -------
        T : ndarray (n_zones, n_zones)
            バランシング済みOD行列
        """
        T = np.outer(self.O, self.D) * F
        for _ in range(max_iter):
            # 行和制約
            row_sum = T.sum(axis=1)
            row_factor = np.where(row_sum > 0, self.O / row_sum, 0.0)
            T = T * row_factor[:, None]

            # 列和制約
            col_sum = T.sum(axis=0)
            col_factor = np.where(col_sum > 0, self.D / col_sum, 0.0)
            T = T * col_factor[None, :]

            # 収束判定
            row_err = np.max(np.abs(T.sum(axis=1) - self.O))
            col_err = np.max(np.abs(T.sum(axis=0) - self.D))
            if row_err < tol and col_err < tol:
                break

        return T

    def predict(self, param):
        """パラメータを指定してOD行列を推計 (二重制約)

        Parameters
        ----------
        param : float or tuple
            抵抗関数のパラメータ

        Returns
        -------
        T : ndarray (n_zones, n_zones)
        """
        if isinstance(param, (list, tuple, np.ndarray)):
            F = self.impedance_func(self.cost, *param)
        else:
            F = self.impedance_func(self.cost, param)
        np.fill_diagonal(F, 0.0)
        return self.furness(F)

    def mean_trip_cost(self, T):
        """OD行列の平均トリップコスト"""
        total = T.sum()
        if total > 0:
            return np.sum(T * self.cost) / total
        return 0.0

    def trip_length_distribution(self, T, bins):
        """OD行列のトリップ長頻度分布 (TLFD)

        Parameters
        ----------
        T : ndarray (n_zones, n_zones)
            OD行列
        bins : array-like
            コストのビン境界

        Returns
        -------
        freq : ndarray (len(bins)-1,)
            各ビンのトリップ数
        """
        mask = np.ones_like(T, dtype=bool)
        np.fill_diagonal(mask, False)
        costs_flat = self.cost[mask]
        trips_flat = T[mask]
        freq, _ = np.histogram(costs_flat, bins=bins, weights=trips_flat)
        return freq

    def calibrate_hyman(self, observed_mean_cost, param0=None,
                        max_iter=50, tol=1e-4, verbose=False):
        """Hyman法によるパラメータ推定

        観測平均トリップコストに一致するようパラメータを逐次更新する。

        β_{k+1} = β_k * c̄_obs / c̄_model(β_k)

        Parameters
        ----------
        observed_mean_cost : float
            観測平均トリップコスト
        param0 : float, optional
            初期パラメータ値
        max_iter : int
            最大イテレーション数
        tol : float
            収束判定閾値 (相対誤差)
        verbose : bool
            進捗表示

        Returns
        -------
        result : dict
        """
        c_obs = observed_mean_cost

        if param0 is None:
            param0 = 1.0 / c_obs

        param = param0
        history = []

        for i in range(max_iter):
            T = self.predict(param)
            c_model = self.mean_trip_cost(T)

            rel_err = abs(c_model - c_obs) / c_obs if c_obs > 0 else float('inf')
            history.append({'iter': i, 'param': param,
                            'mean_cost_model': c_model, 'rel_error': rel_err})

            if verbose:
                print(f"Iter {i:3d}: param={param:.6f}, "
                      f"c_model={c_model:.4f}, c_obs={c_obs:.4f}, "
                      f"rel_err={rel_err:.6f}")

            if rel_err < tol:
                break

            # Hyman 更新: c_model > c_obs → β が小さすぎ → 増加
            if c_model > 0:
                param = param * c_model / c_obs

        T_final = self.predict(param)
        return {
            'param': param,
            'predicted': T_final,
            'mean_cost_model': self.mean_trip_cost(T_final),
            'mean_cost_observed': c_obs,
            'rel_error': rel_err,
            'history': history,
        }

    def calibrate_tlfd(self, observed_tlfd, bins, param_range=(0.001, 5.0),
                       method='chi2'):
        """トリップ長頻度分布 (TLFD) への適合によるパラメータ推定

        Parameters
        ----------
        observed_tlfd : array-like (len(bins)-1,)
            観測トリップ長頻度分布
        bins : array-like
            コストのビン境界
        param_range : tuple
            パラメータ探索範囲
        method : str
            適合度指標 ('chi2', 'rmse', 'rmsn')

        Returns
        -------
        result : dict
        """
        obs = np.asarray(observed_tlfd, dtype=float)
        obs_total = obs.sum()

        if method == 'chi2':
            # カイ二乗統計量
            def objective(param):
                T = self.predict(param)
                pred = self.trip_length_distribution(T, bins)
                # 正規化して比較
                pred_norm = pred / pred.sum() * obs_total if pred.sum() > 0 else pred
                safe = np.maximum(pred_norm, 1e-10)
                chi2 = np.sum((obs - pred_norm) ** 2 / safe)
                return chi2

        elif method == 'rmse':
            # RMSE
            def objective(param):
                T = self.predict(param)
                pred = self.trip_length_distribution(T, bins)
                pred_norm = pred / pred.sum() * obs_total if pred.sum() > 0 else pred
                return np.sqrt(np.mean((obs - pred_norm) ** 2))

        elif method == 'rmsn':
            # 正規化RMSE (比率で比較)
            def objective(param):
                T = self.predict(param)
                pred = self.trip_length_distribution(T, bins)
                obs_prop = obs / obs_total if obs_total > 0 else obs
                pred_prop = pred / pred.sum() if pred.sum() > 0 else pred
                return np.sqrt(np.mean((obs_prop - pred_prop) ** 2))

        else:
            raise ValueError(f"Unknown method: {method}")

        res = minimize_scalar(objective, bounds=param_range, method='bounded')
        best_param = res.x

        T_final = self.predict(best_param)
        pred_tlfd = self.trip_length_distribution(T_final, bins)

        return {
            'param': best_param,
            'predicted': T_final,
            'predicted_tlfd': pred_tlfd,
            'observed_tlfd': obs,
            'objective': res.fun,
            'bins': np.asarray(bins),
            'mean_cost_model': self.mean_trip_cost(T_final),
        }

    @staticmethod
    def mean_cost_from_od(T_obs, cost_matrix):
        """観測OD行列とコスト行列から平均トリップコストを計算"""
        total = T_obs.sum()
        if total > 0:
            return np.sum(T_obs * cost_matrix) / total
        return 0.0

    @staticmethod
    def tlfd_from_od(T_obs, cost_matrix, bins):
        """観測OD行列からトリップ長頻度分布を計算"""
        mask = np.ones_like(T_obs, dtype=bool)
        np.fill_diagonal(mask, False)
        freq, _ = np.histogram(cost_matrix[mask], bins=bins,
                                weights=T_obs[mask])
        return freq
