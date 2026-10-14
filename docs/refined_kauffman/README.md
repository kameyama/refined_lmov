# 第 3 の研究：refined Kauffman 多項式の large N duality（2026-10-14 開始）

このディレクトリは、本研究（KN 1703.05408 の検証、docs/ 直下）とも、第 2 の研究（DGR・四重次数型の構造、docs/dgr_lmov/）とも**別の研究**の文書である。
ユーザーの指示（2026-10-14）：「3 つめの研究として Refined Kauffman の large N duality を研究せよ。ディレクトリは新しく分けよ。
Dropbox の『Mathematica/Refined Kauffman』に過去の試行錯誤の結果がある。これはうまくいかなかったが同じ方向性で研究せよ。
まずは unrefined の Kauffman の large N duality を調べ、そののち過去の試行錯誤にミスがないか調査せよ。」

## 目的
- SO(N)（または Sp(N)）Chern–Simons の彩色 Kauffman 不変量の large N duality（Mariño 0904.1088：SO(N) の自由エネルギーから U(N) の composite 表現の寄与の 1/2 を引いた残りが、crosscap c = 1, 2 の整数 BPS 不変量で書ける）を、refined の設定（Macdonald 多項式、Aganagic–Schaeffer 1202.4456 の orientifold の refined topological string と q-brane / t-brane）に持ち上げ、refined Kauffman 不変量から refined BPS 不変量を抽出する式を作り、整数性と正値性を調べる。
- 過去の試行錯誤（2018 年 7〜10 月、Mathematica ノートブック 42 本。場所：`/Users/masaya/Library/CloudStorage/Dropbox/Mathematica/Refined Kauffman/`、git の外）と同じ方向で進める。まず unrefined の構造を固め、次に過去の計算のミスの有無を調べる。

## 規則
共通の規則（CLAUDE.md の研究上のルール）はこの研究にも適用する：主張の分類（【導出】【計算（範囲）】【予想】【仮説】）、文献情報の原典確認、既知の結果を成果として書かないこと、場当たり的な調整の禁止、自動テスト、外部への送信の事前許可、停止点での報告。

この研究に固有の点：
- 対象は Kauffman（SO/Sp）側なので、本研究の docs/open_problems.md §F（rCS の正値性の扱い）や docs/conventions.md の U(N) 固有の規約はそのままは適用しない。規約はこのディレクトリの conventions.md に別に定める（unrefined は Mariño の規約を基準にする）。
- unrefined の Kauffman の BPS 不変量（Mariño の \hat N_{R,g,c,β}）は整数だが正とは限らない。refined 版で何に正値性を課すかは、unrefined の構造を固めてから決める（研究計画 K2 以降）。
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
| progress.md | 進捗（日付ごと、新しいものが上） |
| open_problems.md | 未解決 |

## 作業を再開するとき
1. progress.md の先頭「現在の状況（再開用）」を読む。
2. research_plan.md の状態の表と停止点を確認する。
3. 規約は conventions.md（この研究固有）。本研究の docs/conventions.md は U(N) 側の参照用。
