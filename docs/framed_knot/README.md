# framed_knot の文書

Kameyama–Nawata（arXiv:1703.05408）の refined large N duality を、トーラス結び目以外の結び目と framed knot に拡張するテーマの文書。
作業を再開するときは progress.md の先頭「現在の状況（再開用）」から読む。規則は .claude/rules/framed_knot.md。

## 直下（入口になる文書）

| ファイル | 中身 |
|---|---|
| progress.md | 進捗の記録。先頭に「現在の状況（再開用）」 |
| conventions.md | 規約の正本（変数、正規化、framing、refined BPS の抽出式） |
| open_problems.md | 未解決の問題 |
| novelty_report.md | 既知の結果と新規性の範囲（文献調査） |

## ディレクトリ

| ディレクトリ | 中身 |
|---|---|
| reports/ | 結果の報告。RESULTS.md（主要な結果の一覧）、REPORT_overview.md（全体の整理）、REPORT.md（第 1 期：ケーブル）、REPORT_knot_table.md（第 2 期：結び目表）、REPORT_phase3.md〜REPORT_phase5.md（第 3〜5 期） |
| plans/ | 研究計画と停止条件。research_plan.md（第 1 期）、research_plan_phase2.md（不採用のケーブル案）、research_plan_knot_table.md（第 2 期）、research_plan_phase3.md、research_plan_phase4.md、research_plan_21.md（別研究 [2,1]） |
| categorification/ | 圏化の文献地図（categorification_literature_map.md）と、非代数的なケーブルの [1] の圏化の状況（categorification_status.md） |
| literature/ | 文献ごとのノートと文献調査（KN、AS・GN、DAHA、GGS、[2,1] の調査 w0_21_*、framing の考察 framing_large_n.md） |
| paper/ | 論文。en/（英語版）、ja/（日本語の準備稿）、phase1_draft/（第 1 期だけの旧草稿）、generated/（スクリプトが生成する表）、references.bib、準備の文書（outline.md、claims.md、notation.md など）。詳細は paper/README.md |

文書の中では、このディレクトリからの相対パス（例：reports/RESULTS.md）か、リポジトリ直下からのパス（docs/framed_knot/reports/RESULTS.md）で他の文書を指す。
