# 重力モデル パラメータ推計

交通需要の重力モデルによるOD行列推計とパラメータ推定。

## モデル

$$T_{ij} = a_i \cdot b_j \cdot O_i \cdot D_j \cdot f(c_{ij})$$

- $O_i$: ゾーン $i$ の生成交通量
- $D_j$: ゾーン $j$ の集中交通量
- $c_{ij}$: ゾーン間コスト
- $f(c)$: 抵抗関数

## 抵抗関数

| 型 | 数式 | パラメータ |
|---|---|---|
| 指数型 | $f(c) = \exp(-\beta c)$ | $\beta$ |
| べき乗型 | $f(c) = c^{-\gamma}$ | $\gamma$ |

## 推定手法

- **Poisson 最尤推定 (MLE)**: $\max \sum T_{ij}^{obs} \log T_{ij}^{pred} - T_{ij}^{pred}$
- **最小二乗法 (LSQ)**: $\min \sum (T_{ij}^{obs} - T_{ij}^{pred})^2$

## 制約条件

- **生成量制約**: 行和が $O_i$ に一致
- **二重制約 (Furness法)**: 行和・列和がそれぞれ $O_i$, $D_j$ に一致

## 参考文献

- Ortúzar, J. de D. & Willumsen, L. G. (2011). *Modelling Transport* (4th ed.). Wiley.
- Wilson, A. G. (1970). *Entropy in Urban and Regional Modelling*. Pion.

## ファイル構成

| ファイル | 内容 |
|---|---|
| `gravity_model.py` | 重力モデルクラス (推計・推定) |
| `demo.ipynb` | 合成データでのパラメータ復元実験 |

## 依存ライブラリ

- NumPy
- SciPy

## 使用例

```python
from gravity_model import GravityModel

gm = GravityModel(O, D, cost_matrix, impedance='exp')
result = gm.estimate(T_observed, constraint='doubly', method='mle')
print(f"beta = {result['param']:.4f}, R2 = {result['r_squared']:.4f}")
```
