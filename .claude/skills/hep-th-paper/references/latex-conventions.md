# LaTeX の慣習

## クラスと体裁

- 既存の原稿は、今のプリアンブルと文献スタイルのまま直す（テーマごとの .claude/rules/ の規則が優先する）。
- 新しい原稿は投稿先を先に決める。
  - JHEP 系（hep-th の標準）：`\documentclass[a4paper,11pt]{article}` と `\usepackage{jheppub}`。`\title`、`\author`、`\affiliation`、`\emailAdd`、`\abstract{...}`、`\keywords{...}`、`\preprint{...}` を書くと、`\maketitle` がタイトルページと目次を作る。templates/main.tex がこの骨格である。jheppub.sty は TeX Live に入っていないので、JHEP のサイト（https://jhep.sissa.it/jhep/help/JHEP_TeXclass.jsp）から取って原稿のディレクトリに置く（手元では references/1703.05408/src/jheppub.sty にもある）。
  - math-ph・math.QA 寄り（CMP、LMP など）：`amsart`。
  - PRD・PRL：REVTeX 4-2。
- 文献スタイル：JHEP 体裁なら `JHEP.bst`（または `utphys.bst`）。.bib は INSPIRE の BibTeX をそのまま使う（skill `inspire-bibtex`）。

## 数式

- 本文から参照する式にだけ番号を付ける。参照は `\eqref{}`。式番号は節ごと（`\numberwithin{equation}{section}`。jheppub は既定でそうなる）。
- `eqnarray` を使わず、`align`、`equation` と `split`、`multline` を使う。
- ディスプレイ式は文の一部である。式の末尾にカンマかピリオドを付け、式のあとの "where" は字下げしない（式の直後に空行を入れない）。
- 記法は初出で 1 回だけ定義する。同じ量に 2 つの記号を使わない。
- 長い多項式の表は本文に並べず、付録か ancillary file に置く。
- 節・式・表・図の参照の書き方（"section 2" か "Section 2" か、"eq. (2.3)" か "(2.3)" か）は原稿の中で統一する。JHEP の慣習は小文字の "section 2"、"table 1"、"figure 1"、"appendix A"。

## 定理環境

- hep-th の論文では定理環境は控えめに使う。物理的な主張は散文の "We conjecture that ..." で書くことが多い。
- 数学的に厳密で証明のある主張だけを Proposition / Theorem にする。証明がなければ Conjecture。計算で確かめただけの主張を Proposition にしない。
- Remark を多用しない。本文に入れられる注意は散文か脚注にする。5 個を超えたら見直す。
- このリポジトリの原稿には `hypothesis`、`result`（Computational Result）の環境がある。使うかどうかは claims.md の格付けとユーザーの判断に従い、勝手に Proposition に変えない。

## 引用

- `\cite{a,b,c}` でまとめて引用する。数字のスタイルなら [1–5] のように範囲にまとめる（`cite` パッケージか、natbib・biblatex の sort&compress）。
- 文の主語に引用を使うときは著者名を添える（"Aganagic and Shakirov [3] formulated ..."）。引用ラベルを文の主語にしない。abstract には `\cite` を書かず、著者名だけを書く。
- "see e.g. [12–15] and references therein"、"as shown in [7]" は普通の言い回しである。
- 引用する層：原典、レビュー、最近の関連研究。研究論文の引用は数十本が普通で、20 本に満たないときは、先行研究の扱いが薄くないかを確かめる。水増しはしない。
- 予想・公式・データを借りた箇所には、その場で出典を書く。

## 強調と脚注

- 太字で強調しない。定義する語を初出で `\emph{}` にする程度にとどめる。
- 補足的な技術上の注意は脚注に書いてよい。
- 表は booktabs（`\toprule`、`\midrule`、`\bottomrule`）で組み、caption に範囲と記号の意味を書く。caption は完結した文にする。

## arXiv に出すとき

- 主カテゴリ hep-th、クロスリストは math-ph、math.QA、math.GT など内容に応じて。
- Comments にページ数と図・表の数、ancillary file の有無を書く。
- 投稿はユーザーが行う。
