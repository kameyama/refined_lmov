# CLAUDE.md

iterated torus（ケーブル）結び目の refined Chern–Simons 不変量と refined BPS 不変量を計算する研究プロジェクト。
タスクの仕様は `CLAUDE_TASK_refined_cable_knots.md`（以下「指示書」）にある。作業前に必ず読むこと。

## 構成

```
justfile               # サンドボックス外で動かすレシピ（setup / fetch-refs）
requirements.txt       # Python 依存（sympy, python-flint, pytest）
.venv/                 # just setup で作る（コミット対象外）
refined_cables/
  refs/<arXiv id>/     # 論文 PDF・TeX ソース（src/）。1703.05408 の ancillary は src/anc/*.nb
  src/                 # 実装（Python パッケージ）
  tests/               # pytest
  data/                # 計算結果 JSON と人が読む表
  REPORT.md            # 最終レポート
  LOG.md               # 作業ログ（規約表・不一致・未解決点・実行時間）
```

## コマンド

- 環境構築（ネットワーク・`.pem` 書き込みが要るのでサンドボックス外）: `just setup`
- 論文取得（arXiv へ接続）: `just fetch-refs`
- テスト（サンドボックス内で可）: `.venv/bin/python -m pytest refined_cables/tests -q`
- `just` は単独で呼ぶこと。パイプや `&&` でつなぐとサンドボックス内で動き、失敗する。

Wolfram Engine はインストール済みだが未アクティベート。ノートブック（`.nb`）はテキストとして読む。

## 研究上のルール（指示書 §5 の要約）

- 規約を場当たり的に調整しない。不一致は原因を特定し、`LOG.md` の規約表を更新する。
- 文献から直接取れない式には「仮説」と明記し、裏付けた検証を書く。
- 各ステージで pytest を先に書き、全部通してから次へ進む。
- 長時間計算は途中結果を `data/` に保存し、再開できるようにする。
- 論文の取得以外で、計算結果やコードを外部に送信しない（git push 等も含む）。
