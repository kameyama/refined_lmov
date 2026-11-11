# 論文の準備（docs/framed_knot/paper/）

論文（framed knot への KN の拡張と、KN のデータを超えた G 展開の前後の正整数性の検証）を書く前の準備の文書。
準備の段階は日本語、最終版は英語。投稿先は未定（arXiv には上げる）。2026-11-11 更新。

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
| en/main.tex, en/sections/*.tex, en/tables/*.tex | **英語版の本文（2026-11-11 に新しく書き起こした初稿、29 ページ）**。`cd docs/framed_knot/paper/en && latexmk -pdf main.tex`（pdflatex、サンドボックスの中で動く）。内容の根拠は ja/ と claims.md。構成は ja/ から変えてある（命題を使う節へ移動、framing の走査を §7.2 に集約、予想 A・B を §7.3 にまとめた）。en/tables/ は generated/ の表の見出しを英語にした写しで、数値は同じ |
| phase1_draft/main.tex | 第 1 期だけを扱った旧草稿（日本語）。参考として残してある |

関連：docs/framed_knot/literature/framing_large_n.md（framed の refined large N の物理的な考察）、docs/framed_knot/novelty_report.md §10（新規性）、docs/framed_knot/conventions.md（規約の正本）。
