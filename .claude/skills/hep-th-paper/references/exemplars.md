# 模範論文の構造

構成と語り口をまねるための要約である。本文を写さない（arXiv の論文は再配布できるライセンスではないので、このスキルにも本文を入れていない）。数は、references/<arXiv 番号>/src/ の TeX ソースを 2026-11-11 に数えたもの。原文を見たいときは、そのソースか PDF を読む（worktree に無ければ本体のチェックアウトの references/ にある）。

## Kameyama–Nawata, "Refined large N duality for knots"（arXiv:1703.05408）

ユーザー自身の論文で、今の原稿の直接の先行研究でもある。共著者が見慣れている文体はこれなので、基準はまずここに置く。

- 体裁：article 11pt a4paper と jheppub。`\preprint`、`\affiliation`、`\emailAdd`、`\abstract{}` を使う。文献スタイルは halpha（アルファベットのラベル）。
- Abstract：1 段落、6 文。引用も数式の表示もない。流れは "We formulate ..." → "By studying ..., we provide ..." → "This form enables us to relate ..." → "Assuming that ..., the duality predicts ..." → "Thus, this formulation can be interpreted as a positivity conjecture ..." → "We also discuss ..."。仮定（"Assuming that ..."）と、予想であること（"positivity conjecture"）を abstract の中で明示している。
- 本文の節は 7 つ：Introduction / Deformed conifold side / Resolved conifold side（小節 2 つ）/ Large N duality for torus knots / Large N duality for torus links / Extension to non-torus knots / Discussion。小節は論文全体で 2 つだけ。Acknowledgements は番号なし。
- 付録は 3 つ：Symmetric functions（番号なしの小見出し 3 つ）/ Explicit formulas of refined reformulated invariants / Tables of BPS degeneracies。表は 33 個で、ほぼすべて付録 C にある。本文は式と議論、データは付録、という分け方である。
- Introduction：約 80 行のソース。't Hooft → Gopakumar–Vafa → Ooguri–Vafa・LMOV → 結び目ホモロジーと Gukov–Schwarz–Vafa → Aganagic–Shakirov、と原典を引きながら文脈を 2 段落で述べ、3 段落目を "The goal of this paper is to formulate ..." で始める。続いて M 理論の brane の配置を式で示し、最後の段落は "The organization of this paper is as follows. In §2, we review ..."。寄与の箇条書きはない。
- 引用：57 本。`\cite{a,b,c}` でまとめる。文の主語にするときは著者名を書く（"Gukov, Schwarz and Vafa \cite{...} proposed that ..."）。
- 定理環境：使っていない。主要な予想は "Now let us formulate the conjecture of refined large N duality for torus knots, which is the main claim of this paper." と散文で導入し、その 1 か所だけ枠（tcolorbox）で囲んでいる。links の節では `\emph{}` で予想の文を示す。つまり、枠は論文全体で主要な予想 1 つにだけ使い、まとめや要点の囲みには使っていない。
- 箇条書き：3 か所。どれも直前の式の性質や解釈を、式や完結した文で 2〜3 項目並べるもので、寄与や要点の一覧ではない。
- 太字のラベル、`\paragraph`、Remark：ない。脚注は 4 個。
- 語り口："we shall ..."、"let us ..."、"To this end, ..."、"It is known that ..."、"we conjecture that ..."。確かめた範囲は "we present some examples of ... for the Hopf link and the T_{2,4} link" のように対象を具体的に書く。
- Discussion：3 段落ほど。1 段落目で何をしたかと仮定を 1 回だけまとめ、残りは未解決の問題（数学的な定義、他のゲージ群への拡張）を散文で書く。箇条書きはない。

## Aganagic–Shakirov, "Knot Homology from Refined Chern-Simons Theory"（arXiv:1105.5117）

分野の標準的な書き方の参照。組版は harvmac（plain TeX）なので、体裁ではなく構成と語り口を見る。

- Abstract：1 段落、9 文。"We formulate ..."、"We give an explicit solution ..."、"We conjecture that ..."、"For knots in S^3 colored by fundamental representation, the theory ends up computing ..."、"As a byproduct, we show that ..."、"We also provide ..."。定式化したこと、予想、示したこと、当てはまる範囲を文ごとに書き分けている。
- 本文の節は 8 つ（Introduction、review の節が 2 つ、M 理論からの導出、理論の定義、large N duality、knot homology と例、Some Future Directions）と Acknowledgments。付録は 3 つ（Macdonald Polynomials / The refined Chern-Simons matrix model / Verlinde Formula）。標準的な対称関数の復習は付録に置いている。
- Introduction の最後は "The organization of the paper is as follows. In section 2 we review ..."。その中に "In all the cases that we have checked the knot polynomial we obtain agrees with the superpolynomial of [...]" とあり、確認が有限個の例であることをそのまま書いている。
- 例の節：Trefoil knot、the (2,5) torus knot、the (2,2m−1) torus knot、the (3,4) torus knot、"Prediction: colored superpolynomial of the trefoil" と、例ごとに小節を立て、既知の結果との一致と、新しい予言を区別している。
- 引用：68 本。脚注 26 個。"conjecture" という語は約 27 回現れ、定理環境は使わない。番号つきの列挙は 1 か所（superpolynomial を得る手順の 4 項目）。

同じ著者の短い解説版 arXiv:1202.2489 は、同じ内容を短くまとめるときの比較に使える（ここでは要約していない）。

## 今の原稿に当てはめるときの目安

- 節の数は 7〜8、小節は必要なところだけ。
- 主要な予想は 1 つか 2 つに絞り、散文で導入する。枠で囲むなら主要な予想 1 つまで。
- データの表は付録、本文は式と議論。
- 引用は 50 本前後が、この話題の論文の普通の規模である。
- 定理環境を多く使う今の原稿（Proposition、Hypothesis、Computational Result）は、数学寄りの読者には合うが、上の 2 本とは違う形である。どちらに寄せるかはユーザーが決めることなので、勝手に変えずに尋ねる。
