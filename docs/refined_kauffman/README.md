# 第 3 の研究：refined Kauffman 多項式の large N duality（2026-10-14 開始）

このディレクトリは、本研究（KN 1703.05408 の検証、docs/framed_knot/）とも、第 2 の研究（DGR・四重次数型の構造、docs/dgr_lmov/）とも**別の研究**の文書である。
ユーザーの指示（2026-10-14）：「3 つめの研究として Refined Kauffman の large N duality を研究せよ。ディレクトリは新しく分けよ。
Dropbox の『Mathematica/Refined Kauffman』に過去の試行錯誤の結果がある。これはうまくいかなかったが同じ方向性で研究せよ。
まずは unrefined の Kauffman の large N duality を調べ、そののち過去の試行錯誤にミスがないか調査せよ。」

## 目的
- SO(N)（または Sp(N)）Chern–Simons の彩色 Kauffman 不変量の large N duality（Mariño 0904.1088：SO(N) の自由エネルギーから U(N) の composite 表現の寄与の 1/2 を引いた残りが、crosscap c = 1, 2 の整数 BPS 不変量で書ける）を、refined の設定（Macdonald 多項式、Aganagic–Schaeffer 1202.4456 の orientifold の refined topological string と q-brane / t-brane）に持ち上げ、refined Kauffman 不変量から refined BPS 不変量を抽出する式を作り、整数性と正値性を調べる。
- 過去の試行錯誤（2018 年 7〜10 月、Mathematica ノートブック 42 本。場所：`/Users/masaya/Library/CloudStorage/Dropbox/Mathematica/Refined Kauffman/`、git の外）と同じ方向で進める。まず unrefined の構造を固め、次に過去の計算のミスの有無を調べる。

## 規則
共通の規則（CLAUDE.md の研究上のルール）はこの研究にも適用する：主張の分類（【導出】【計算（範囲）】【予想】【仮説】）、文献情報の原典確認、既知の結果を成果として書かないこと、場当たり的な調整の禁止、自動テスト、外部への送信の事前許可、停止点での報告。

この研究に固有の点：
- 対象は Kauffman（SO/Sp）側なので、本研究の docs/framed_knot/open_problems.md §F（rCS の正値性の扱い）や docs/framed_knot/conventions.md の U(N) 固有の規約はそのままは適用しない。規約はこのディレクトリの conventions.md に別に定める（unrefined は Mariño の規約を基準にする）。
- unrefined の Kauffman の BPS 不変量（Mariño の \hat N_{R,g,c,β}）は整数だが正とは限らない。refined 版で何に正値性を課すかは、unrefined の構造を固めてから決める（研究計画 K2 以降）。
- refined では正値性を期待する（2026-10-19、ユーザーの指示）。unrefined のクロスキャップの BPS 数の符号は β ごとに変わる（ε_β）ので、refined の正値性は「ブロックごとに ε_β と同じ符号にそろう」こととして判定する（report.md §8 の予想）。unrefined で符号が混在することを、refined で正値性を期待しない理由にしない。
- 過去のノートブックは読むだけで書き換えない。再現はこのリポジトリの Python（src/rcsknots/）か Sage で行い、tests/ で照合する。

## ファイル
| ファイル | 内容 |
|---|---|
| research_plan.md | 計画（段階 K0〜、停止点） |
| conventions.md | Kauffman 側の規約（変数、unknot、composite 表現、符号、framing） |
| literature.md | 文献調査（unrefined の Kauffman large N duality、refined の候補、既知の結果の範囲） |
| past_attempt.md | 2018 年の試行錯誤の整理（ノートブックごとの内容、仮説の式、検査の結果、ミスの候補と検証結果） |
| report.md | 本文（K0〜K2' の結果と導出、現時点の判断） |
| literature_deep_research_2026-10-14.md | 外部の deep research の報告（未検証の要約） |
| deep_research_prompt_AB.md | 正値性の障害（A：凍結した状態の多重項、B：E と C の分離）についての deep research への指示（閉じた形） |
| literature_deep_research_2026-10-21_AB.md | その報告（外部の deep research の出力、未検証の要約。要点は literature.md §5） |
| progress.md | 進捗（先頭に再開用の「現在の状況」、その下に日付ごとの記録。記録は古いものが上で、新しいものは末尾に足す） |
| theory_orientifold_index.md | orientifold の開弦の refined 指数の理論（全 d の構造、E と C の分離、閉弦との照合、捻ったトレース、§1〜§12。§12.23〜§12.31 は n = 1 の E を状態の空間の同定から作る原理と、n = 2 の背景・結び目の部分） |
| tex/report.tex, report.pdf | 計算結果の報告書（TeX、生成した式は tex/gen/） |
| paper/ | 論文の草稿。ja_v2/（第 2 版、2026-11-08 から本体。main.tex、sections/、refs.bib、main.pdf）、ja/（第 1 版、記録として残す）、claims.md（主張の格付け）、outline.md（第 1 版の構成案） |
| daha_fast.md | DAHA-Jones の高速版の設計と確認 |

未解決の問題は progress.md の「現在の状況」§5 と、論文の第 2 版の §10.3 にまとめてある（open_problems.md は作っていない）。

## 作業を再開するとき
1. progress.md の先頭「現在の状況（再開用）」を読む。§0 に別の PC・別のセッションで再開する手順（clone するブランチ、Python と TeX の環境、データ、文献、ユーザーとの取り決め）、§1 に論文の段階、§2 に研究の到達点と次の候補、§4 に試して失敗した方法がある。
2. 論文の作業なら paper/ja_v2/（本体）と paper/claims.md（主張と根拠）を読む。書き方の規則は .claude/rules/refined_kauffman_paper.md。
3. research_plan.md の状態の表と停止点は 2026-10-18 時点のもので、それ以降の進み方は progress.md に記録してある。
4. 規約は conventions.md（この研究固有）。本研究の docs/framed_knot/conventions.md は U(N) 側の参照用。
