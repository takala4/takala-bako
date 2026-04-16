# Frank-Wolfe 法による利用者均衡配分

Wardrop の第一原則 (利用者均衡) を Beckmann の等価最小化問題として定式化し、
Frank-Wolfe (Conditional Gradient) 法で解く Python 実装。

## アルゴリズム概要

1. **初期化**: 自由走行時間での All-or-Nothing 配分
2. **補助フロー計算**: 現在コストでの All-or-Nothing 配分
3. **探索方向**: 補助フロー - 現在フロー
4. **ステップサイズ決定**: 二分法による Beckmann 目的関数の最小化
5. **フロー更新**: f = f + λ(y - f)
6. **収束判定**: 相対双対ギャップが閾値を下回れば終了

## 参考文献

- LeBlanc, L. J., Morlok, E. K., & Piercey, W. P. (1975). An efficient approach to solving the road network equilibrium traffic assignment problem. *Transportation Research*, 9(5), 309-318.
- Sheffi, Y. (1985). *Urban Transportation Networks*. Prentice-Hall.

## ファイル構成

| ファイル | 内容 |
|---|---|
| `frank_wolfe.py` | Frank-Wolfe ソルバー |
| `network.py` | ネットワークデータ構造 (GMNS Plus 形式読み込み) |
| `cost_function.py` | BPR リンクコスト関数 |
| `demo.ipynb` | Sioux Falls ネットワークでの使用例 |

## 入力データフォーマット

[GMNS Plus Dataset](https://github.com/HanZhengIntelliTransport/GMNS_Plus_Dataset) 形式に準拠 (`itapas/data/` のデータを共用)。

## 依存ライブラリ

- NumPy
- SciPy
- Pandas
- matplotlib (ノートブック可視化用)

## 使用例

```python
from network import Network
from cost_function import bpr_cost, bpr_cost_derivative, bpr_cost_integral
from frank_wolfe import FrankWolfe

net = Network("node.csv", "link.csv", "demand.csv")
solver = FrankWolfe(net, bpr_cost, bpr_cost_derivative, bpr_cost_integral)
link_flow, log = solver.solve(max_iter=200, gap_threshold=1e-3)
```
