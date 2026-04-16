# takala-bako

個人的によく使うコードをストックするリポジトリです。

## プログラム一覧

### 数理最適化

| ディレクトリ | 内容 | 言語 |
|---|---|---|
| [LCP](LCP/) | 線形相補性問題 (Linear Complementarity Problem) ソルバー群 | Python |
| [hitchcock](hitchcock/) | ヒッチコック輸送問題 | Python |
| [projection](projection/) | Capped Simplex への射影 | Python |

### ネットワーク・交通工学

| ディレクトリ | 内容 | 言語 |
|---|---|---|
| [dijkstra](dijkstra/) | ダイクストラ法 (優先度付きキュー使用) | C |
| [itapas](itapas/) | iTAPAS 利用者均衡配分 (Xie & Xie, 2016) | Python |
| [UE](UE/) | 利用者均衡配分 (User Equilibrium) | Python |
| [input_network_data](input_network_data/) | 道路ネットワークデータ読み込みモジュール | C |
| [make_network_data](make_network_data/) | 格子状ネットワーク作成ツール | Python |

### 数値計算

| ディレクトリ | 内容 | 言語 |
|---|---|---|
| [gauss-seidel](gauss-seidel/) | ガウス＝ザイデル法 (連立一次方程式の反復解法) | C |
| [SIR](SIR/) | SIR感染症モデルソルバー | Python |
| [tropical](tropical/) | トロピカル線形代数 (max-plus代数) による最短経路問題 | Python |

## ドキュメント

* [コーディング規約](document/coding_conventions.md)
* [データフォーマット](document/file_format.md)
