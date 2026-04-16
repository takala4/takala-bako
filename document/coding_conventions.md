# コーディング規約

本リポジトリにおけるコーディング規約を定める。


## 全般

### ディレクトリ構成

* 各アルゴリズム・ツールごとに独立したディレクトリを作成する
* 各ディレクトリには `README.md` を配置し、内容・使い方を記載する
* データファイルは各ディレクトリ内の `data/` サブディレクトリに配置する

### ファイル

* 文字コード: UTF-8
* 改行コード: CRLF
* Git で追跡しないファイルは `.gitignore` に記載する
  * `__pycache__/`, `*.pyc`, `.ipynb_checkpoints/`, `.DS_Store` 等

### コメント

* コードの意図が自明でない箇所にコメントを付与する
* デバッグ用の出力 (`printf("dbg...")`, `print("debug")` 等) はコミットに含めない


## Python

### バージョン

* Python 3 を対象とする

### 命名規則

| 対象 | 規則 | 例 |
|---|---|---|
| クラス名 | CapWords (PascalCase) | `Model`, `SIRSolver` |
| 関数・メソッド名 | snake_case | `calc_step_size`, `input_data` |
| 変数名 | snake_case | `num_nodes`, `step_size` |
| 定数 | UPPER_SNAKE_CASE | `MAX_ITER`, `TOLERANCE` |
| ファイル名 | 内容を表す名前 (CapWords または snake_case) | `LCP_DCA.py`, `SIR_Solver.py` |

### インポート

* 標準ライブラリ、サードパーティ、ローカルの順にグループ化する
* 使用するモジュールは全てファイル先頭で `import` する

### データ入出力

* CSV/DAT ファイルの読み書きは `file_format.md` の仕様に従う


## C

### 命名規則

| 対象 | 規則 | 例 |
|---|---|---|
| 構造体名 | CapWords + `_t` サフィックス | `Link_t`, `Node_t` |
| 関数名 | snake_case | `input_parameter_data`, `add_heap` |
| 変数名 | snake_case | `num_nodes`, `link_id` |
| マクロ・定数 | UPPER_SNAKE_CASE | `MAX_SIZE` |

### メモリ管理

* `malloc` で確保したメモリは使用後に必ず `free` で解放する
* ファイルポインタは使用後に必ず `fclose` で閉じる

### エラー処理

* ファイルオープン失敗時は `stderr` にメッセージを出力し、`exit(1)` で終了する


## Jupyter Notebook

* ノートブックのファイル名は内容を表す簡潔な名前とする
* 実行結果 (output) をコミットに含めてよいが、大きな画像出力は避ける
