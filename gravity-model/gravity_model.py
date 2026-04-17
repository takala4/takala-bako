"""重力モデルによるOD需要推計とパラメータ推定

重力モデル:
    T_ij = O_i * D_j * f(c_ij) / Σ_j D_j * f(c_ij)   (生成量制約)
    T_ij = a_i * b_j * O_i * D_j * f(c_ij)             (二重制約)

抵抗関数 f(c):
    指数型:  f(c) = exp(-β * c)
    べき乗型: f(c) = c^(-γ)

パラメータ推定:
    最尤推定 (Poisson回帰) または最小二乗法

References:
    Ortúzar, J. de D. & Willumsen, L. G. (2011).
        Modelling Transport (4th ed.). Wiley.
    Wilson, A. G. (1970). Entropy in Urban and Regional Modelling. Pion.
"""

import numpy as np
from scipy.optimize import minimize_scalar, minimize


def impedance_exp(cost, beta):
    """指数型抵抗関数: f(c) = exp(-β * c)"""
    return np.exp(-beta * cost)


def impedance_power(cost, gamma):
    """べき乗型抵抗関数: f(c) = c^(-γ)"""
    return np.where(cost > 0, cost ** (-gamma), 0.0)


class GravityModel:
    """重力モデル

    Parameters
    ----------
    origins : array-like (n_zones,)
        生成交通量 O_i
    destinations : array-like (n_zones,)
        集中交通量 D_j
    cost_matrix : array-like (n_zones, n_zones)
        ゾーン間のコスト行列 c_ij
    impedance : str
        抵抗関数の種類 ('exp' or 'power')
    """

    def __init__(self, origins, destinations, cost_matrix, impedance='exp'):
        self.O = np.asarray(origins, dtype=float)
        self.D = np.asarray(destinations, dtype=float)
        self.cost = np.asarray(cost_matrix, dtype=float)
        self.n_zones = len(self.O)

        if impedance == 'exp':
            self.impedance_func = impedance_exp
        elif impedance == 'power':
            self.impedance_func = impedance_power
        else:
            raise ValueError(f"Unknown impedance type: {impedance}")
        self.impedance_type = impedance

    def predict_singly_constrained(self, param):
        """生成量制約モデルによるOD行列の推計

        Parameters
        ----------
        param : float
            抵抗関数のパラメータ (β or γ)

        Returns
        -------
        T : ndarray (n_zones, n_zones)
            推計OD行列
        """
        F = self.impedance_func(self.cost, param)
        # 対角要素 (ゾーン内) を除外
        np.fill_diagonal(F, 0.0)

        denom = F @ self.D  # Σ_j D_j * f(c_ij)
        denom = np.where(denom > 0, denom, 1e-20)

        T = np.outer(self.O, self.D) * F / denom[:, None]
        return T

    def predict_doubly_constrained(self, param, max_iter=100, tol=1e-6):
        """二重制約モデル (Furness法) によるOD行列の推計

        Parameters
        ----------
        param : float
            抵抗関数のパラメータ
        max_iter : int
            バランシングの最大イテレーション数
        tol : float
            収束判定閾値

        Returns
        -------
        T : ndarray (n_zones, n_zones)
            推計OD行列
        """
        F = self.impedance_func(self.cost, param)
        np.fill_diagonal(F, 0.0)

        a = np.ones(self.n_zones)
        b = np.ones(self.n_zones)

        for _ in range(max_iter):
            # b_j = D_j / Σ_i a_i * O_i * f(c_ij)
            col_sum = (a * self.O) @ F
            b = np.where(col_sum > 0, self.D / col_sum, 0.0)

            # a_i = O_i / Σ_j b_j * D_j * f(c_ij)
            # 実際には a_i = 1 / Σ_j b_j * D_j * f(c_ij) (O_i は外で掛ける)
            row_sum = F @ (b * self.D)
            a_new = np.where(row_sum > 0, 1.0 / row_sum, 0.0)

            if np.max(np.abs(a_new - a)) < tol:
                a = a_new
                break
            a = a_new

        T = np.outer(a * self.O, b * self.D) * F
        return T

    def estimate(self, observed, constraint='doubly', method='mle',
                 param_range=(0.01, 5.0)):
        """観測OD行列からパラメータを推定

        Parameters
        ----------
        observed : array-like (n_zones, n_zones)
            観測OD行列
        constraint : str
            制約の種類 ('singly' or 'doubly')
        method : str
            推定手法 ('mle' or 'lsq')
        param_range : tuple
            パラメータの探索範囲

        Returns
        -------
        result : dict
            推定結果 {'param': float, 'objective': float, 'predicted': ndarray}
        """
        T_obs = np.asarray(observed, dtype=float)
        mask = T_obs > 0  # 正の需要がある OD ペアのみ

        if constraint == 'singly':
            predict_func = self.predict_singly_constrained
        else:
            predict_func = self.predict_doubly_constrained

        if method == 'mle':
            # Poisson 最尤推定: max Σ T_obs * log(T_pred) - T_pred
            def neg_loglik(param):
                T_pred = predict_func(param)
                T_pred = np.maximum(T_pred, 1e-20)
                ll = np.sum(T_obs[mask] * np.log(T_pred[mask]) - T_pred[mask])
                return -ll

            res = minimize_scalar(neg_loglik, bounds=param_range, method='bounded')
            best_param = res.x
            best_obj = -res.fun

        elif method == 'lsq':
            # 最小二乗法: min Σ (T_obs - T_pred)^2
            def sse(param):
                T_pred = predict_func(param)
                return np.sum((T_obs - T_pred) ** 2)

            res = minimize_scalar(sse, bounds=param_range, method='bounded')
            best_param = res.x
            best_obj = res.fun

        else:
            raise ValueError(f"Unknown method: {method}")

        T_pred = predict_func(best_param)

        # 適合度指標
        rmse = np.sqrt(np.mean((T_obs - T_pred) ** 2))
        r_squared = 1.0 - np.sum((T_obs - T_pred) ** 2) / np.sum((T_obs - T_obs.mean()) ** 2)
        total_obs = T_obs.sum()
        total_pred = T_pred.sum()

        return {
            'param': best_param,
            'objective': best_obj,
            'predicted': T_pred,
            'rmse': rmse,
            'r_squared': r_squared,
            'total_observed': total_obs,
            'total_predicted': total_pred,
        }

    def estimate_multivariate(self, observed, impedance_func_custom,
                               param0, constraint='doubly', method='mle'):
        """複数パラメータの抵抗関数に対する推定

        Parameters
        ----------
        observed : array-like
            観測OD行列
        impedance_func_custom : callable
            f(cost, *params) → ndarray のカスタム抵抗関数
        param0 : array-like
            パラメータ初期値
        constraint : str
        method : str

        Returns
        -------
        result : dict
        """
        T_obs = np.asarray(observed, dtype=float)
        mask = T_obs > 0
        original_func = self.impedance_func

        def predict(params):
            self.impedance_func = lambda c, p=None: impedance_func_custom(c, *params)
            if constraint == 'singly':
                T = self.predict_singly_constrained(0)  # param は impedance_func 内で使用
            else:
                T = self.predict_doubly_constrained(0)
            return T

        if method == 'mle':
            def neg_loglik(params):
                T_pred = predict(params)
                T_pred = np.maximum(T_pred, 1e-20)
                return -np.sum(T_obs[mask] * np.log(T_pred[mask]) - T_pred[mask])
            res = minimize(neg_loglik, param0, method='Nelder-Mead')
        else:
            def sse(params):
                T_pred = predict(params)
                return np.sum((T_obs - T_pred) ** 2)
            res = minimize(sse, param0, method='Nelder-Mead')

        best_params = res.x
        T_pred = predict(best_params)
        self.impedance_func = original_func

        rmse = np.sqrt(np.mean((T_obs - T_pred) ** 2))
        r_squared = 1.0 - np.sum((T_obs - T_pred) ** 2) / np.sum((T_obs - T_obs.mean()) ** 2)

        return {
            'params': best_params,
            'objective': res.fun,
            'predicted': T_pred,
            'rmse': rmse,
            'r_squared': r_squared,
        }
