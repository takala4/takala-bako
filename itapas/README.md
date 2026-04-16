# iTAPAS (incremental Traffic Assignment by Paired Alternative Segments)

静的利用者均衡 (User Equilibrium) 交通量配分の Python 実装。

## アルゴリズム概要

iTAPAS は、Paired Alternative Segments (PAS) を用いた均衡配分アルゴリズムである。
全経路を列挙せず、同一始点・終点を持つ経路セグメントのペア上でフローを均衡化する。

主要ステップ:
1. **初期化**: All-or-Nothing 配分
2. **PAS 識別**: Maximum Flow Search (MFS) によりポテンシャルリンクから PAS を特定
3. **フローシフト**: Newton ステップにより PAS 上のフローを均衡化
4. **収束判定**: 相対双対ギャップが閾値を下回れば終了

## 参考文献

- Bar-Gera, H. (2010). Traffic assignment by paired alternative segments. *Transportation Research Part B*, 44(8-9), 1022-1046.
- Xie, J. & Xie, C. (2016). New insights and improvements of using paired alternative segments for traffic assignment. *Transportation Research Part B*, 93, 406-424.

## ファイル構成

| ファイル | 内容 |
|---|---|
| `itapas.py` | iTAPAS ソルバー本体 |
| `network.py` | ネットワークデータ構造 (GMNS Plus 形式読み込み) |
| `cost_function.py` | BPR リンクコスト関数 |
| `demo.ipynb` | Sioux Falls ネットワークでの使用例 |
| `data/SiouxFalls/` | サンプルデータ (GMNS Plus 形式) |

## 入力データフォーマット

[GMNS Plus Dataset](https://github.com/HanZhengIntelliTransport/GMNS_Plus_Dataset) 形式に準拠。

必要ファイル:
- `node.csv`: ノード定義 (node_id, zone_id, x_coord, y_coord)
- `link.csv`: リンク定義 (link_id, from_node_id, to_node_id, capacity, lanes, vdf_fftt, vdf_alpha, vdf_beta)
- `demand.csv`: OD 需要 (o_zone_id, d_zone_id, volume)

## 依存ライブラリ

- NumPy
- SciPy
- Pandas
- matplotlib (ノートブック可視化用)

## 使用例

```python
from network import Network
from cost_function import bpr_cost, bpr_cost_derivative, bpr_cost_integral
from itapas import ITAPAS

net = Network("data/SiouxFalls/node.csv",
              "data/SiouxFalls/link.csv",
              "data/SiouxFalls/demand.csv")

solver = ITAPAS(net, bpr_cost, bpr_cost_derivative, bpr_cost_integral)
link_flow, log = solver.solve(max_iter=50, gap_threshold=1e-3)
```
