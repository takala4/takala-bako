# LCP (Linear Complementarity Problem) ソルバー

線形相補性問題に対する各種最適化アルゴリズムの実装。

## ソルバー一覧

| ファイル | アルゴリズム |
|---|---|
| `LCP_AVE.py` | Averaging method |
| `LCP_DCA.py` | Difference of Convex functions Algorithm |
| `LCP_DCAE.py` | DCA with Envelope method |
| `LCP_DCAT.py` | DCA with Trust region |
| `LCP_FISTA.py` | Fast Iterative Shrinkage-Thresholding Algorithm |
| `LCP_FS.py` | Facchinei-Solodov algorithm (1995, 1997) |
| `LCP_FW.py` | Frank-Wolfe algorithm |
| `LCP_NAG.py` | Nesterov Accelerated Gradient descent |
| `LCP_Luca.py` | Fischer-Burmeister function based method |

## 入力データ

`Lemke100/` ディレクトリにサンプルデータ (M.csv, b.csv) があります。

## 依存ライブラリ

* NumPy
* Pandas
* Gurobi (gurobipy)

## Notebooks

* `LCP_Test.ipynb` : ソルバーのテスト実行
* `Merit.ipynb` : メリット関数の可視化
* `NAG.ipynb` : NAGアルゴリズムの実験
* `AdaHyGrad_test.ipynb` : AdaHyGrad テスト
* `Optimal_Stop.ipynb` : 最適停止条件の検証
