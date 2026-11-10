---
paths:
  - "docs/refined_kauffman/**"
  - "scripts/refined_kauffman/**"
  - "tests/refined_kauffman/**"
  - "data/refined_kauffman/**"
  - "logs/refined_kauffman/**"
  - "src/rcsknots/refined_kauffman/**"
---
# refined_kauffman の規則

目的・文書の一覧・再開の手順は docs/refined_kauffman/README.md。論文の書き方は .claude/rules/refined_kauffman_paper.md。

- 対象は Kauffman（SO/Sp）側なので、framed_knot の docs/framed_knot/open_problems.md §F（rCS の正値性の扱い）や docs/framed_knot/conventions.md の U(N) 固有の規約はそのままは適用しない。規約は docs/refined_kauffman/conventions.md に別に定める（unrefined は Mariño の規約を基準にする）。
- unrefined の Kauffman の BPS 不変量（Mariño の \hat N_{R,g,c,β}）は整数だが正とは限らない。refined 版で何に正値性を課すかは、unrefined の構造を固めてから決める（研究計画 K2 以降）。
- refined では正値性を期待する（2026-10-19、ユーザーの指示）。unrefined のクロスキャップの BPS 数の符号は β ごとに変わる（ε_β）ので、refined の正値性は「ブロックごとに ε_β と同じ符号にそろう」こととして判定する（docs/refined_kauffman/report.md §8 の予想）。unrefined で符号が混在することを、refined で正値性を期待しない理由にしない。
- 過去の試行錯誤（Dropbox の Mathematica ノートブック。場所は docs/refined_kauffman/README.md）は読むだけで書き換えない。再現はこのリポジトリの Python（src/rcsknots/）か Sage で行い、tests/ で照合する。
