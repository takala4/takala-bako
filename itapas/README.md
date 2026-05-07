# iTAPAS (incremental Traffic Assignment by Paired Alternative Segments)

静的利用者均衡 (User Equilibrium, UE) 交通量配分の Python 実装。
オプションとして **リンクトール (与件)** と **システム最適 (System Optimal, SO)** にも対応。

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
| `tutorial_toll_so.ipynb` | トール (与件) / SO / Pigovianトール のチュートリアル |
| `data/SiouxFalls/` | サンプルデータ (GMNS Plus 形式) |

## 入力データフォーマット

[GMNS Plus Dataset](https://github.com/HanZhengIntelliTransport/GMNS_Plus_Dataset) 形式に準拠。

必要ファイル:
- `node.csv`: ノード定義 (node_id, zone_id, x_coord, y_coord)
- `link.csv`: リンク定義 (link_id, from_node_id, to_node_id, capacity, lanes, vdf_fftt, vdf_alpha, vdf_beta, [vdf_toll])
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

# 標準UE配分
link_flow, log = solver.solve(max_iter=50, gap_threshold=1e-3)

# リンクトールを与件としたUE配分 (一般化費用 = t(f) + toll)
# tolls省略時は link.csv の vdf_toll 列が使われる (無ければゼロ)
import numpy as np
tolls = np.zeros(net.num_links)
tolls[10:20] = 5.0
link_flow_ue, _ = solver.solve(max_iter=50, gap_threshold=1e-3, tolls=tolls)

# システム最適 (SO) 配分
# マージナルコスト mc(f) = t(f) + f·t'(f) を用いた均衡として求解
# SOはコスト勾配が急なため、mu (PAS削除閾値) はUEより小さめを推奨
link_flow_so, _ = solver.solve(max_iter=50, gap_threshold=1e-3,
                                inner_iterations=20, mu=1e-7, mode='SO')

# Pigovian (一次最適) トール τ* = f · t'(f) は solve() 後に自動計算
# これを与件トールとしてUE配分するとSOフローを再現する
optimal_toll = solver.optimal_toll
link_flow_check, _ = solver.solve(mode='UE', tolls=optimal_toll)
# → link_flow_check ≈ link_flow_so

# SO + トール (一般化費用ベースのSO均衡) も同様に組合せ可能
link_flow_so_t, _ = solver.solve(mode='SO', tolls=tolls,
                                  inner_iterations=20, mu=1e-7)
```

### `solve()` の主要パラメータ

| 引数 | 既定値 | 説明 |
|---|---|---|
| `mode` | `'UE'` | `'UE'` (Wardrop第一原理) または `'SO'` (Wardrop第二原理) |
| `tolls` | `None` | リンクトール配列 (長さ`num_links`)。`None`なら`network.link_toll`を使用 |
| `gap_threshold` | `1e-6` | 相対双対ギャップの収束閾値 |
| `inner_iterations` | `20` | 各外側反復でのフローシフト反復回数 |
| `mu` | `1e-3` | PAS均衡判定閾値 (SOでは`1e-7`程度を推奨) |
