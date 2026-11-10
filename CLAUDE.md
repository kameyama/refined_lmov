# 研究テーマ
このリポジトリでは、refined Chern–Simons 不変量と refined BPS 不変量について、下の表のテーマを扱う。コード（src/rcsknots/ 直下）と規約の大半は共有するが、目的と規則はテーマごとに異なる。
テーマはディレクトリ名で呼ぶ。docs, scripts, tests, data, logs の各ディレクトリと src/rcsknots/ の下に、同じ名前のサブディレクトリがある（framed_knot だけは src/rcsknots/ 直下の共有コードを使う）。

| テーマ | 内容 | 入口（最初に読む） | 進捗 | 規約 |
|---|---|---|---|---|
| framed_knot | Kameyama–Nawata（arXiv:1703.05408）の refined large N duality を、トーラス結び目以外の結び目と framed knot に拡張し、refined LMOV 不変量の正値性を検証する | docs/framed_knot/progress.md | 同左 | docs/framed_knot/conventions.md |
| dgr_lmov | DGR・四重次数型の構造を refined LMOV 不変量で表し、unrefined のデータから refined の非長方形の彩色を構成する（2026-10-11 開始） | docs/dgr_lmov/README.md | docs/dgr_lmov/progress.md | framed_knot と共通 |
| refined_kauffman | refined Kauffman 多項式（SO/Sp 側）の large N duality（2026-10-14 開始） | docs/refined_kauffman/README.md | docs/refined_kauffman/progress.md | docs/refined_kauffman/conventions.md |

- テーマ固有の規則は .claude/rules/ にテーマ名のファイルで置く。そのテーマのファイルを開くと読み込まれる。
- 正値性の扱いはテーマごとに違う。あるテーマの規則（例：framed_knot の open_problems.md §F）を別のテーマに持ち込まない。
- dgr_lmov の「unrefined から refined の [2,1]・[3,1] を一意に予想する」方向は、2026-10-14 にユーザーが打ち切りと判断した。この方向の計算は、ユーザーの指示なしに新しく始めない（理由は .claude/rules/dgr_lmov.md）。

# 作業を始めるとき
1. 依頼がどのテーマの作業かを判断する。判断できなければユーザーに尋ねる。
2. そのテーマの入口の文書と、進捗の先頭にある「現在の状況（再開用）」を読む。
3. 計算に入る前に、そのテーマの規約を確認する。

# コマンド
- テスト：`.venv/bin/pytest`（既定は機能テストだけ）。1 ファイルだけなら `.venv/bin/pytest tests/<テーマ>/test_xxx.py`。保存済みの結果の再判定は `-m results`。
- スクリプト：framed_knot と dgr_lmov は `PYTHONPATH=src:scripts/framed_knot:scripts/dgr_lmov .venv/bin/python scripts/<テーマ>/xxx.py`、refined_kauffman は `PYTHONPATH=src:scripts/refined_kauffman .venv/bin/python scripts/refined_kauffman/xxx.py`。
- SageMath：`DOT_SAGE=$TMPDIR/dotsage /opt/homebrew/Caskroom/miniforge/base/envs/sage/bin/sage script.sage`（サンドボックス内で動く）。
- サンドボックス外で動かすもの：`just setup`（環境の再構築）、`just fetch-refs`（論文の取得）、`just wolfram <file.wls>`（Wolfram Engine）。`just` は単独で呼ぶ。パイプや `&&` でつなぐとサンドボックス内で動いて失敗する。
- ツールのバージョンと動作確認の詳細は README.md の「環境」節。

# 研究上のルール（全テーマ共通）
- 報告とコメントは日本語で書く。
- 主張はすべて【導出】【計算（範囲を明記）】【予想】【仮説】のいずれかに分類して書く。文献から直接取れなかった式には【仮説】と明記し、どの検証で裏付けたかを書く。
- 文献情報（arXiv 番号、DOI、巻・頁）は原典で確認したものだけを書き、確認できないものは「未確認」と明記する。
- 既知の結果を成果として書かない。既知の範囲は各テーマの文献調査（docs/framed_knot/novelty_report.md、docs/dgr_lmov/novelty.md、docs/refined_kauffman/literature.md）にある。
- 既知の結果と合わせるために係数や符号を場当たり的に調整しない。不一致は原因を特定して、そのテーマの規約と進捗に記録する。特定できない場合は「未解決」と明記する。
- 規約を変更する場合は、そのテーマの規約の文書に変更内容と理由を記録し、影響を受けるコードとテストを更新する。
- 照合は tests/ の自動テストにする。多項式は全係数を厳密に比較し、数値補間で得たものは補間に使っていない点で検証する。
- tests/ は機能テストとし、小さな入力（数秒程度）でコードが正しく動くことを確かめる。保存済みの大きな計算結果を判定し直すテストには `results` の印を付ける。
- 作業の区切りごとに、そのテーマの進捗を更新する（日付、やったこと、結果、次にやること、未解決の問題）。
- 各テーマの研究計画に書いた停止点と停止条件では、作業を止めて報告し、確認を求める。

# 長時間の計算（目安：1 件 30 分以上）
実行する前に次の準備を済ませ、その内容を進捗に書く。
- 見積もり：小さい例で時間とメモリを測り、対象の規模から外挿する。24 時間を超える見積もりなら、実行前に報告して確認を求める。
- 高速化：重い部分（評価点の数、補間の次数、記号計算の約分など）を特定し、アルゴリズム・キャッシュ・不要な計算の削除で減らせないかを先に検討する。
- 再開と分割：独立な単位（結び目、彩色、ρ、素数、評価点のまとまり）ごとに結果を data/ に保存し、保存済みの単位は読み込んで飛ばす。CPU のコア数とメモリの範囲で並列に実行する。
- 実行方法：背景実行の制限時間（最大 2 時間）を超える見込みなら、nohup などの独立プロセスで実行し、ログと完了の目印を決めておく。実行時間とメモリは logs/ に記録する。

# git と外部への送信
- ローカルでの git commit は許可なしで行ってよい。意味のある単位でこまめに commit する。
- main ブランチへの直接の commit とローカルでのマージは禁止する（2026-10-08、ユーザーの指示）。main は GitHub の PR を経由してのみ変更する。PR のマージはユーザーが行うか、ユーザーの指示を受けてから行う。
- 計算結果、コード、文書を外部サービスに送信・公開する場合（ブランチの push、PR の作成、メール、Slack、共有ドキュメントなど）は、毎回、事前に意図と内容を説明して明示的な許可を得る。論文と補助ファイルのダウンロードは許可なしで行ってよい。
- リモートは origin = https://github.com/kameyama/refined_lmov。

# 環境の注意
- サンドボックスは `.git` への書き込みを拒否する。git の書き込み操作はサンドボックス外で行う。
- Wolfram Engine はサンドボックス内ではライセンスを読めず失敗する。サンドボックス外で実行する。
- リポジトリ直下に `refs/` を作らない（git の内部ディレクトリと同名で、サンドボックスに作成を拒否される）。文献は `references/` に置く。
- data/cache/ は作り直せる途中結果の置き場で、git で追跡しない。無ければ計算し直す。保存しておく結果は data/<テーマ>/ に置く。
- 4 箱以上の結果（data/framed_knot/highbox/）も git で追跡していない。worktree には無いので、本体のチェックアウトのものを読む。
