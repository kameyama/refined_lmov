---
paths:
  - "docs/refined_kauffman/paper/**"
  - "docs/refined_kauffman/tex/**"
---
# refined_kauffman の論文の書き方

## 版
- paper/ja_v2 が日本語の本体（第 2 版、2026-11-08 作成。肯定的な結果を本体にし、多箱の否定的な結果は 1 章、brane 2・image 2 の探索は付録）。日本語の推敲は ja_v2 で行う。paper/ja は第 1 版で、記録として残す。
- paper/en は英語版（2026-11-10 作成）。ja_v2 の翻訳ではなく、hep-th の論文の形に構成を変えたもの。refs.bib は ja_v2 と共有する。査読・校正の指摘と判定は paper/en/review_20261110.md と review_20261110_round2.md。
- 生成ファイルは tex/gen（ja と ja_v2 で共有）と tex/gen_en（英語版）。tex/gen_en は `PYTHONPATH=src .venv/bin/python scripts/refined_kauffman/rk_report_tex_en.py` が tex/gen の見出しを英語に置き換えて作る（数式は同じ）。
- 主張と根拠は paper/claims.md。

## 引用と記法
- refs.bib のエントリは skill `inspire-bibtex` で作る（INSPIRE の BibTeX をそのまま使い、キーと shorthand だけ替える。INSPIRE に無い情報を補ったエントリには出典をコメントで書く）。
- 引用は \cite{キー} だけで書く。参考文献は refs.bib と biblatex（style=alphabetic、backend=bibtex）で組み、各エントリに `shorthand = {キー}` を付けてラベルを bib のキーと同じにする（[KN17]、[AgSc12]、[Marino09] など）。
- 本文に著者の略号や名前を書かない。「KN [KN17]」「AgSc の refined unknot」は誤りで、「[KN17]」「[AgSc12] の refined unknot」と書く。見出しで引用するときは \texorpdfstring{\cite{KN17}}{[KN17]} を使う。
- 表現（ヤング図）を文中に書くときは「表現 [1]」「R = [1]」のように、何の [1] かが分かる語を添える。裸の「[1]」は引用ラベルと区別できない。
- 数式の中の表現はヤング図で書く（main.tex の \yone, \ytwo, \yoneone, \ythree, \ytwoone, \yoneoneone。caption や見出しでも使える）。例：$g^{\mathrm{true}}_{\yone} = (t/q)^{\delta/2} g_{\yone}$。

## 用語（日本語版）
- 日本語の専門文献で定着している訳語（規格化、被覆、商、鏡像、多重項、二重項、母関数、核、整数計画、Möbius の帯、法束など）だけ日本語にし、それ以外は原語かカタカナで書く（タワー、ハンドル、ツイスト、ストランド、reduced、composite、brane、image、framing、hook など）。
- 「塔」「取っ手」「撚り」「簡約」「走査」「準位」「枠」のような訳語や、「移り目」のような造語は使わない。論文だけの用語は初出で定義する。

## コンパイルと確認
- `TEXMFVAR=$TMPDIR/texmfvar lualatex main.tex` → `bibtex main` → `lualatex main.tex` を 2 回（ja_v2 と en で同じ）。
- コンパイル後は、ログに `!` のエラーと「undefined」がないことを確認し、変更した箇所のページを画像にして見る（`gs -sDEVICE=png16m`）。
