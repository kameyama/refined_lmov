# 構成の規則

hep-th の慣習は、JHEP の公式マニュアルと分野の実務を合わせたものである。確立した規範ではないので、迷ったら exemplars.md の模範論文と投稿先の author guidelines を見る。

## 全体の並び

title・author・abstract → 目次 → Introduction → 準備（review）→ 主要結果 → 例と検証 → Discussion → Acknowledgments（番号なし）→ Appendix → 参考文献。

節の数は本文で 5〜8 程度を目安にする。2〜3 段落ごとに `\subsection` を切らない。1 ページに満たない節は、隣の節に統合できないかを先に考える。`\paragraph{...}` を小見出しの代わりに乱発しない。

## Abstract

- 1 段落、5〜10 文程度。引用、ディスプレイ式、脚注を入れない。
- "We formulate ...", "We conjecture that ...", "We show that, for ..., ..." のように、定式化したこと・予想したこと・確かめた範囲を書き分ける。確かめた範囲（結び目、彩色、箱の数など）を具体的に書く。
- 背景の説明は 1 文まで。"In recent years, ... has attracted much attention" のような前置きで始めない。
- "novel"、"comprehensive framework"、"deep connection" のような自己評価を書かない。

## Introduction

1. 分野の文脈：何が知られているか。原典を引用しながら、この論文に必要な流れだけを書く。
2. 未解決点：既存の結果がどこで止まっているか。
3. この論文の寄与：何をしたかを散文で書く。主結果の式を 1〜2 本、ここで示してよい。ML 論文のような "Our contributions are:" の箇条書きにしない。
4. 既知の結果との境界：どこまでが既知で、どこからが新しいか。既知の結果を成果のように書かない（各テーマの文献調査を見る）。
5. 最後の段落："The paper is organized as follows. In section 2 we review ..."。

長さは 2〜4 ページが普通で、背景の教科書的な復習で膨らませない。

## 準備の節

- 標準的な対象（Macdonald 多項式、HOMFLY-PT 多項式、LMOV の再定式化など）の定義は、引用と、必要な式だけにする。教科書的な再説明を長く書かない。長い復習が要るなら付録に置く。
- この論文に固有の記法と規約（変数の対応、正規化、framing の取り方など）は、省かずにここで 1 回だけ定義する。読者が既存の文献と比べられるように、文献の規約との対応を書く。

## 主要結果と検証

- 主張の格付けに合わせて書き分ける。
  - 証明がある：Proposition / Theorem と証明。
  - 証明がない一般的な主張：Conjecture。散文で "We conjecture that ..." と書くのが hep-th では普通で、定理環境にするのは主要な予想だけでよい。出典がある予想は出典を示す。
  - 計算で確かめた事実："We have checked (3.5) for all knots in table 2 with up to three boxes." のように、範囲を明記した 1〜2 文と、付録の表への参照。
- 検証の記録（計算ログ、PASS/FAIL の表、スクリプト名、実行時間）は本文に書かない。結果の表は付録に置く。計算方法の説明が要るなら付録に "Details of the computation" として散文で書く。コードを公開する場合は、付録か脚注で 1 回だけ触れる。
- 例は、一般的な主張のあとに具体的な式で示す。式を並べるだけにせず、何が読み取れるかを文で書く。

## Discussion

- 散文で書く。まとめ直しは 1 段落までにし、残りは未解決の問題と今後の方向に使う。
- 未解決問題を箇条書きにするなら 3〜5 項目までで、各項目は 1〜3 文の完結した文にする。

## Acknowledgments と AI の開示

- 番号なしの節。議論の相手、グラント、計算資源を書く。
- AI の使用の開示もここに書く。arXiv は生成 AI の significant use の報告を求め、AI を著者にすることを認めず、内容の責任はすべて著者にあるとしている。開示文は著者が確定する。案："The author used Claude (Anthropic) for computations, literature search and drafting. All results, references and text were checked by the author, who takes full responsibility for the content."

## Appendix

長い計算、規約の一覧、対称関数の公式、数値データの表を置く。本文から必ず参照する。参照されない付録は削る。
