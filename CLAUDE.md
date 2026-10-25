# プロジェクト概要
Kameyama–Nawata "Refined large N duality for knots"（arXiv:1703.05408）は、主にトーラス結び目・トーラス絡み目について refined Chern–Simons（CS）不変量を resolved conifold 上の refined BPS 不変量 N̂_{R,g,β,J_r} に書き換え、正値性予想を検証した。本研究では、これをトーラス結び目以外の結び目クラスへ拡張する。
現在の対象：iterated torus 結び目（トーラス結び目のケーブル）。最優先は [(2,3),(2,13)]（代数的）と [(2,3),(2,1)]（非代数的）。
必要な彩色：|λ| ≤ 3 の全ヤング図形 [1], [2], [1,1], [3], [2,1], [1,1,1]。
主眼：(1) 1703.05408 の結果の再現による規約の固定、(2) refined ケーブル公式（Rosso–Jones 型の Schur→Macdonald 持ち上げ）と DAHA（Cherednik–Danilenko）の照合、(3) 新しい結び目の refined BPS 不変量の正値性・整数性の判定（代数的 vs 非代数的、長方形 vs [2,1]）。

# リポジトリの構成（2026-10-14 に整理）
docs, src/rcsknots, scripts, tests, data, logs の各ディレクトリの中で、本研究（framed knot への拡張と正値性の検証）のファイルは `framed_knot/`、別の研究（DGR・四重次数型の構造と refined LMOV）のファイルは `dgr_lmov/` に置く（一覧は README.md の「構成」）。
src/rcsknots/ 直下は両方の研究で共有する計算コード、src/rcsknots/dgr_lmov/ は別の研究だけが使うモジュール。data/cache/ は共有のキャッシュ（追跡しない）。
スクリプトは `PYTHONPATH=src:scripts/framed_knot:scripts/dgr_lmov` で実行する（dgr_lmov のスクリプトは framed_knot のスクリプトも import する）。

# セッション開始時に必ず行うこと
1. docs/framed_knot/progress.md を読み、現在のステージと未完了の作業を把握する。
2. 作業内容に応じて docs/framed_knot/research_plan.md（ステージ 0〜7 の計算計画）と docs/framed_knot/novelty_report.md（既知の結果と新規性の範囲）を参照する。
3. 計算に入る前に docs/framed_knot/conventions.md の規約を確認する。
4. 別の研究「DGR・四重次数型の構造と refined LMOV」（2026-10-11 開始）の作業では、docs/dgr_lmov/README.md を読み、そこに書いた規則に従う。正値性の扱いが本研究と異なり、この研究では open_problems.md §F は適用しない。その研究の進捗は docs/dgr_lmov/progress.md に書く。
5. 第 3 の研究「refined Kauffman 多項式の large N duality」（2026-10-14 開始）の作業では、docs/refined_kauffman/README.md を読み、そこに書いた規則に従う。規約は docs/refined_kauffman/conventions.md、進捗は docs/refined_kauffman/progress.md に書く。過去の試行錯誤（Dropbox の Mathematica ノートブック）は読むだけで書き換えない。 この研究のファイルは docs/refined_kauffman/、scripts/refined_kauffman/、tests/refined_kauffman/、data/refined_kauffman/、logs/refined_kauffman/、src/rcsknots/refined_kauffman/（この研究だけが使うモジュール：Koornwinder・D_N/C_N の DAHA・flint の有理関数・Kauffman の各計算）に置く。スクリプトは `PYTHONPATH=src:scripts/refined_kauffman` で実行する。

# 規約
- 変数 (a, q, t) の対応（AS・DGR・1703.05408）、Macdonald 多項式の正規化、S・T 行列と γ 因子、framing（ケーブルの q は companion の Seifert framing で測る）、reduced/unreduced と U(N)/SU(N)、プレシズムの規約、refined BPS 抽出式は docs/framed_knot/conventions.md に従う。
- 規約を変更する場合は、conventions.md に変更内容と理由を記録し、影響を受けるコードとテストを更新する。

# 研究上のルール
- 報告とコメントは日本語で書く。
- 主張はすべて「解析的に導出」「計算機で確認（範囲を明記）」「予想」のいずれかに分類して書く。文献から直接取れなかった式（特に refined ケーブル公式）には「仮説」と明記し、どの検証で裏付けたかを書く。
- 文献情報（arXiv 番号、DOI、巻・頁）は原典で確認したものだけを書き、確認できないものは「未確認」と明記する。
- novelty_report.md で既知とされた結果（トーラス結び目の refined 不変量、4_1 の [1]/[2]/[1,1] の BPS 不変量、DAHA による iterated torus の定義など）を、本研究の成果として書かない。
- 既知の結果と合わせるために係数や符号を場当たり的に調整しない。不一致は原因を特定して conventions.md と progress.md に記録し、特定できない場合は「未解決」と明記する。
- 照合は tests/ の自動テストにする。対象は、1703.05408 の値の再現、t=q での colored HOMFLY との一致、ブレイドから計算した HOMFLY との一致、DAHA の明示例との一致、R↔Rᵀ の対称性。多項式は全係数を厳密に比較し、数値補間で得たものは補間に使っていない点で検証する。
- tests/ は機能テストとする。コードが正しく動くことを、小さな入力（数秒程度）で確かめる。主要な結果は、コードとテストが信頼できることによって正しいと判断する。
  - 保存済みの大きな計算結果を判定し直すテストには `results` の印を付け、既定の実行から外す（pyproject.toml）。内容や理論に問題が出たときだけ `-m results` で実行する。
  - 主要な結果は、計算時の判定とログを一次の記録として docs/framed_knot/RESULTS.md に整理する。
- 正値性の判定で重視するのは、\hat f_ρ から得る refined LMOV 不変量の、G 展開の前（太字の係数 \hat𝐍_{ρ;i,j,k}）と後（\hat N_{ρ,g,β,J_r}）の両方での正整数係数性である。f^q, f^{\bar t} の係数は refined LMOV 不変量ではなく、補足にとどめる。
- ヤング図 [1] は特別に扱う。
  - [1] では \hat f_{[1]} は原理的に rCS_{[1]} と同一になる（\hat f_{[1]} = (a^{1/2} − a^{−1/2}) rCS_{[1]}）。
  - したがって [1] の場合に限り、homological な不変量（ホモロジーの Poincaré 多項式として解釈される量）はすべて正整数係数でなければならない。
  - KN のアプローチに肯定的な立場をとるなら、G 展開ができるかどうか（展開が閉じ、正整数係数になるか）も [1] で確認できなければならない。
  - 新しい結び目では、まず [1] で G 展開の前と後の両方の正整数係数性を確認し、その結果を彩色付きの判定より先に報告する。
- 長時間の計算は途中結果を data/ に保存して再開できるようにし、実行時間とメモリを logs/ に記録する。
- 長時間（目安：1 件 30 分以上）が予想される計算は、実行する前に次の準備を済ませる。準備の内容は実行前に progress.md に書く。
  - 見積もり：小さい例で時間とメモリを測り、対象の規模から外挿する。見積もりが 24 時間を超えるなら、実行前に報告して確認を求める。
  - 高速化：重い部分（評価点の数、補間の次数、記号計算の約分など）を特定し、アルゴリズム・キャッシュ・不要な計算の削除で減らせないかを先に検討する。
  - 再開：停止しても途中から再開できるようにする。彩色ごと・ρ ごと・素数ごと・評価点のまとまりごとなど細かい単位で結果を保存し、保存済みの単位は読み込んで飛ばす。
  - 分割と並列：独立な単位（結び目、彩色、ρ、素数、評価点のまとまり）に分け、CPU のコア数とメモリの範囲で並列に実行する。
  - 実行方法：ツールの背景実行の制限時間（最大 2 時間）を超える見込みなら、nohup などの独立プロセスで実行し、ログと完了の目印（最後の行など）を決めておく。
- 非代数的なケーブルについて、DAHA の値を入力とする正値性の計算（3 箱への拡張や破れ方の分類を含む）は計画に入れない。理由は docs/framed_knot/open_problems.md D5。
- 作業の区切りごとに docs/framed_knot/progress.md を更新する（日付、やったこと、結果、次にやること、未解決の問題）。
- 停止点では、必ず作業を止めて報告し、確認を求める。停止点は次のとおり。
  - ステージ 0 の後（規約表の確定）
  - ステージ 2 の後（トーラス結び目での 1703.05408 の再現）
  - ステージ 5 の後（refined ケーブル公式と DAHA の照合。本計算に入る前）
  - research_plan.md の停止条件に当たったとき

# 外部への送信
- 計算結果、コード、文書を外部サービス（GitHub などへの push、メール、Slack、共有ドキュメントなど）に送信・公開する場合は、事前に意図と内容を説明して明示的な許可を得る。論文と補助ファイルのダウンロードは許可なしで行ってよい。
- ローカルでの git commit は許可なしで行ってよい。意味のある単位でこまめに commit する。
- **main ブランチへの直接の commit とローカルでのマージは禁止する。main は GitHub の PR を経由してのみ変更する**（2026-10-08、ユーザーの指示）。
  作業は main 以外のブランチで行い、ブランチを push して PR を作る。PR のマージはユーザーが行うか、ユーザーの指示を受けてから行う。ブランチの push と PR の作成も外部への送信なので、事前に許可を得る。

# 環境
- Python：.venv（Python 3.14.8、rcsknots は editable インストール）。テストは `.venv/bin/pytest` で実行する。
- 環境の再構築は `just setup`、論文の取得は `just fetch-refs`、Wolfram スクリプトの実行は `just wolfram <file.wls>`（どちらもネットワークが要るのでサンドボックス外で動く just レシピ）。`just` は単独で呼ぶこと。パイプや `&&` でつなぐとサンドボックス内で動いて失敗する。
- SageMath 10.9：conda 環境 `sage`。`/opt/homebrew/Caskroom/miniforge/base/envs/sage/bin/sage script.sage` で実行する（`DOT_SAGE=$TMPDIR/dotsage` を付ける）。サンドボックス内で動く。
- Wolfram Engine 15.0.0：`wolframscript` で使える。ただしサンドボックス内ではライセンスを読めず失敗するので、サンドボックス外で実行する。1703.05408 のノートブックは `NotebookImport[file, "Input" -> "HeldExpression"]` で入力セルを取り出して実行できる。読むだけなら `scripts/framed_knot/nb2txt.py` でテキスト化してもよい。
- サンドボックスは `.git` への書き込みを拒否するので、git の書き込み操作はサンドボックス外で行う（`git -C <repo> ...`）。リモートは GitHub の https://github.com/kameyama/refined_lmov（origin、2026-10-07 にユーザーの指示で設定）。push は外部への送信なので、毎回事前に許可を得る。
- リポジトリ直下に `refs/`（git の内部ディレクトリと同名）を作ると、サンドボックスの保護機能に作成を拒否され、後で削除される。文献は `references/` に置く。
- ツールのバージョンと確認内容の詳細は README.md の「環境」節を参照。
