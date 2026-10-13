# 論文の準備（docs/framed_knot/paper/）

論文（framed knot への KN の拡張と、KN のデータを超えた G 展開の前後の正整数性の検証）を書く前の準備の文書。
準備の段階は日本語、最終版は英語。投稿先は未定（arXiv には上げる）。2026-10-09 時点。

| ファイル | 中身 |
|---|---|
| outline.md | 章立てと、執筆前にすることの一覧（§3） |
| claims.md | 主張の格付け表（定理・命題・計算結果・予想・仮説・既知）と、書く前に解決した点（§7） |
| notation.md | 記号表と表記の決めごと |
| references.bib | 文献データベース（arXiv API で書誌を確認、69 件） |
| data_availability.md | データとコードの公開の準備（公開は未実施、許可が要る） |
| generated/summary_table.tex | 主結果の表の原型（scripts/framed_knot/paper_summary.py が data/framed_knot/bps/ から生成） |
| generated/*.tex（その他） | 第 1 期の草稿の表（scripts/framed_knot/paper_tables.py） |
| ja/main.tex, ja/sections/*.tex | **本文の日本語の準備稿（2026-10-09 に執筆開始）**。`cd docs/framed_knot/paper/ja && latexmk -lualatex main.tex`（サンドボックスの中ではフォントのキャッシュを書けないので外で実行）。各主張の直前に claims.md の ID を `% [ID]` で書く。`% TODO` は執筆中の未了の箇所 |
| main.tex | 第 1 期だけを扱った旧草稿。準備・計算法・BPS の抽出・検証の節は再利用する。軸が変わったので構成は outline.md に従って書き直す |

関連：docs/framed_knot/literature/framing_large_n.md（framed の refined large N の物理的な考察）、docs/framed_knot/novelty_report.md §10（新規性）、docs/framed_knot/conventions.md（規約の正本）。
