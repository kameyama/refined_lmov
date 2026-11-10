# 英語版（paper/en、続き 147 の後）の第 3 回査読・校正・レビュー（2026-11-10）

対象は paper/en の main.tex・sections/*.tex・tex/gen_en（49 ページ）。今回はサブエージェントに任せず、全文を通読して式を手で検算し、疑問点を src/・scripts/・tests/ と照合した。論文のファイルは変更していない。
第 1 回・第 2 回（review_20261110.md、review_20261110_round2.md）で判定済みの指摘は繰り返さない。判定と対応はまだ行っていない。

確認の方法の区別：
- 【検算】本文の式を手で展開して確かめた。
- 【コード照合】リポジトリのコード・テスト・生成ファイルと突き合わせた。
- 【未確認】原典に当たれていない。

---

## A. 誤り（修正が必要）

1. **(2.3) eq:Zexpand の符号 $(-1)^{|\mu|+Nm}$ は $(-1)^{|\mu|}$ が正しい。**【検算・コード照合】
   (2.1) で $x$ 側の分割を $\mu_x = m^N - \nu$ と置き直すと $(-1)^{|\nu|} = (-1)^{Nm-|\mu_x|}$ で、右辺の書き直しの $(-1)^{Nm}$ と打ち消し合う。$N = m = 1$ では $(1-yx)(1-y/x) = -y[K_{(1)}(x) - K_{(1)}(y)]$ で、$W_{(1)}$ の係数は $-y$、$W_\emptyset$ の係数は $+yK_{(1)}(y)$。本文の式は両方の符号が逆になる。コード（kauffman_hopf.py の `(-1)**k * Wk`）は $(-1)^{|\mu|}$ で、計算結果には影響しない。この $(-1)^{Nm}$ は第 1 回の A3 の対応で入れたもの。直後の「The sign $(-1)^{|\mu|}$ in (2.3) is the hallmark …」とも合わなくなっている。

2. **Proposition 2.6（unknot）の証明の核と積公式が (2.2) と合わない。**【検算】
   証明は $\prod_k (1 - t^{N-k}y_j)(1 - t^{N-k}y_j^{-1})$、命題は $\prod_j (y_j^{\pm};t)_\infty/(Uy_j^{\pm};t)_\infty$ と書くが、(2.2) の核は $(1-y x)(1-y x^{-1})$ なので $\prod_k (1 - t^{N-k}y_j)(1 - t^{k-N}y_j)$、積公式は $\prod_j \frac{(y_j;t)_\infty\,(ty_j/U;t)_\infty}{(Uy_j;t)_\infty\,(ty_j;t)_\infty}$。$f_{10} = (1-t^N)(t^N+t)/(t^N(t-1))$ は後者から出る（$-\sum_k (t^{N-k}+t^{k-N})$ を計算して確認）。$y^{-1}$ の形では $\log$ が $Y$ の対称関数にならない。

3. **§5 の $g^{\mathrm{true}}_{\yone}$ の符号が §2 の $g_{\yone}$ と逆。**【検算・コード照合】
   - Hypothesis 5.1 は $f_{\yone} = -W_{\SO} = \tfrac12 h_{\yone} + g_{\yone}$、$\tfrac12 h_{\yone} = +\mathcal P\,(\lambda-\lambda^{-1})/(\bq-\bq^{-1})$ と書く。unknot（$\mathcal P = 1$）では §4.2 の $\tfrac12 f^{\cov}_{10} = -(\lambda-\lambda^{-1})/z_t$ と符号が逆。
   - §5.2 は $g^{\mathrm{true}}_{\yone} = (t/q)^{\delta/2}g_{\yone}$ と書くが、(5.2) と (N)「$\SO(2)$ の点で $1$」は $W_{\SO} - \tfrac12 h$ の符号（コードの `g_family = Z_V·mono·unknot_SO − h_half`）で、§2 の $g_{\yone} = f_{\yone} - \tfrac12 f^{\cov}_{\yone}$（unknot で $-1$、Proposition 2.6）とは逆。
   - (8.2) の $\Delta = \mathrm{even}_{\ba}(g^{\SO}_{\yone}) + 1$ は、右辺を展開すると $1 - \mathrm{even}_{\ba}$[(5.2)] になる。つまり (8.2) の $g^{\SO}_{\yone}$ は §2 の符号で、§8 冒頭の「DAHA 側の規約 $g^{\mathrm{true}} = W_{\SO} - W^{\cov}$」とも (5.2) とも逆。
   - 第 2 回の A19 は §8 冒頭だけを直した。§5 の定義（Hypothesis 5.1 と $g^{\mathrm{true}}$ の定義式）に符号を入れ、(8.2) がどちらの規約かを明記する必要がある。Result 6.3(2) の $g_{\yone} = -\mirror(g^{\mathrm{Mar}}_{\yone})$ とは、$g^{\mathrm{true}} = -(t/q)^{\delta/2}g_{\yone}$ とすれば整合する。

4. **reduced の式で $JD$（規格化済み）をさらに $K_R(t^\rho)$ で割っている。**【コード照合】
   (6.1) は $JD = \{\hat\gamma(P_b)\}/P_b(q^{-\rho_k})$（reduced）と定義するのに、§6.1 と付録 A.2 は $\Rred_R = T_R^{-F}\mirror(JD)/K_R(t^\rho)$、§6.1 と Result A.1(2) は「$JD = T_R^{-F}\mirror(W_R)$（$W_R$ は unreduced）」と書く。コード（`daha_jones` は非 reduced の $Z = \{\hat\gamma(K_\lambda)\}$ を返す）と付録 D.2（「$Z$ is the unnormalized DAHA-Jones value」）が正しい。本文の 3 か所の $JD$ を $Z$（非規格化）に替えるか、$K_R(t^\rho)$ を消す。(6.2) の $W_b = \mirror(JD_b/M_b)$ も同じ。

5. **§10.3「each handle contributes $-(\ql+\ql^{-1})+(\qr+\qr^{-1}) = G$」は $-G$。**【検算】
   §2.1 の定義は $G = \ql+\ql^{-1}-\qr-\qr^{-1}$ で、同じ表に「[KN17] では $-G$」とある（conventions.md の「KN の G とは $(-1)^g$ だけ違う」とも一致）。

6. **Proposition 3.5 の証明の「a pair $\{\omega\},\{\omega t\}$ exchanged by $\tilde I$」は「電荷 $\omega$ の対と電荷 $\omega t$ の対の 2 組」。**【検算】
   $f^{\cov} = 2(\omega+\omega t)/(1-t)$、$f = 2\omega/(1-t)$ になるのは 2 組の対のときだけ。1 組の対は両方の状態が同じ電荷をもつ（Hypothesis 3.1）。

## B. 査読（内容・論理）

7. **Conjecture 6.6 の $|R|$ 偶数の場合が一意に読めない。**「the same holds for $2(1+t)\hat g_R$」は、$X_1[2(1+t)\hat g]$（$t = q$ の部分は $2(1+q)\hat g|_{t=q}$）なのか $2(1+t)X_1[\hat g]$ なのかが決まらない。両者の差は $2(t-q)\hat g|_{t=q}(\tlam)$ で、$\hat g|_{t=q}$ が $\lambda^2 = q$ で消えるときだけ同値。コードは 2 箱で $X_1[\hat g]$ の分子、4 箱で $X_1[(1+t)\hat g]$ を見ている【コード照合】。式で書き、同値になる理由（か、どちらを主張するか）を添える。括弧内の「quotient has the denominator $2(1+t)$」も、$\hat g|_{t=q}$ 由来の $(1+q)$ の分母に触れていない。

8. **(E2) の「$+1$」の位置付け。** $E|_{t=q} = -(1+q)\psi^2(g_{\yone}+1)$ の $+1$ は Hypothesis 3.9 と (3.4) からは出ない（そのまま当てはめると $-(1+q)\psi^2 g_{\yone}$）。unknot の $g^{\mathrm{odd}}_{\ytwo} = -\tfrac12$ を別扱いにする規約（Result 6.3(2) の「subtracting the unknot values $\mp\tfrac12$」）だと思うが、§3.4 では「removes the crosscap constant」としか書かれていない。仮説の一部なのか導出なのかを明記する。

9. **Proposition 9.3「$\Nfix$ is a Laurent polynomial if and only if …」。** (9.6) の第 1 項は $y$ の冪級数なので、主張は「各ブロック $n_k$ が $q, t$ の Laurent 多項式」のはず。$k$ について有限（$K = rs-3$）で切れることは別の主張で、命題からは出ない。

10. **Result 6.1(3) の $f_{(3,1)}, f_{(4)}, \dots$ は $g$ の冪和成分。** $f$ は §2 で $\SO$ 側の自由エネルギーの記号（$1/(1-t)$ の分母をもつ）。$g_\rho := \sum_R\chi_R(\rho)g_R$ と書く。

11. **付録 A の「Every Result in the text is checked by an automated test」と双対 Cauchy 恒等式の検証範囲。** 本文は $(n,m) = (1,1),(2,1),(1,2),(2,2),(2,3),(3,2)$ と書くが、tests/refined_kauffman/test_rk_koornwinder.py にあるのは $(2,1)$ と $(1,2)$ だけ【コード照合】。残りがスクリプトでの確認なら、テストに加えるか書き方を合わせる。

12. **表 3（補間）の区切りと本文が合わない。** $T(2,7)$・$T(2,9)$ の 2 箱（冪 15 個・19 個に対し節点 6〜7 個）と cover の adjoint の $T(2,7)$・$T(2,9)$（節点 6 個）は Newton 補間では決まらないので restricted-range のはずだが、罫線の上にあり、付録 A.2 の列挙（「three boxes of $T(2,5)$ and beyond, four boxes, composites with $m \ge 2$, two boxes of the cables」）にも入っていない。

13. **要旨の「for torus knots with up to four boxes」。** 4 箱は trefoil だけ、$T(2,5)$–$T(2,9)$ は 3 箱まで、$T(3,4)$ は基本表現だけ。第 2 回 A4(iii) と同じ型の強すぎる表現が残っている。

14. **Result 7.2(4) と表 4。** 本文は「two-box computation was done for $\Cab{13}$ only」だが、表 4 に $\Cab{15}$ の $(2)$ と adjoint の行がある。$(1,1)$ だけが無いことを Result に書く。

## C. 校正（記法・用語・英語）

15. **$p_2$, $q_2$ の衝突。** §7 のケーブルのパラメータ $(p_2, q_2)$、$F = p_2q_2$ が、冪和 $p_2$（「the $p_2$ direction」）と変数 $q$ に重なる。第 2 回 B36 で $(p,q)$ から替えた結果、今度は冪和と衝突した。$(k_2, l_2)$ などに。
16. **付録 B の $Q$ の二重定義。** (B.1) では $Q = t^{2N-1}$、その 5 行後から「$(Q, T)$ は [AM22] の変数 $= (q^{1/2}, t^{1/2})$」。さらに「$q/t - t/q$ in the variables of [AM22]」は直前の $(Q,T)$ の宣言と合わない。
17. **付録 C.1 Result C.3「write $T(r,s)$ with $s = p$」。** 本文は $T(3,4)$ を three-strand と呼び、$T(r,s)$ の $r$ が strand 数に読める。$\lfloor r/2\rfloor$、$r'$ の $r$ がどちらかを 1 か所で固定する。
18. **作業ノート調の残り。** 付録 C.3 Result C.8「The extra $(q-t)qt^2$ of an earlier interpretation … was the $E$ of the background」「outside the earlier family with $|e| \le |R|$」は、読者が知らない過去の試行を指す。削るか、何と比べているかを書く。
19. **§11 の「tripling of the closed-string Kähler variable」。** 本文で「tripling」が出るのは付録 B の [AM22] の説明だけ。Remark 6.7 への参照を付ける。
20. **Macdonald 次元の記号。** $\dim_R$（§2.2）、$P_\mu(t^\rho)$・$P_R(t^\rho)$（§2.2・(2.7)）、$K_R(t^\rho)$（§6.1・付録）、$P_{\yone}$（Prop. 9.3 の証明）が混在し、$D_N$ と $\U(n)$ のどちらかが式から分からない。(2.7) が $\SO$ 側と cover 側の両方に使う式なら、そう書く。
21. **§2.1 の表の $\delta$ の行。** 「for thin knots it equals the $\delta$-degree …」は本文で使われない。削ってよい。
22. **main.tex の日付** `November 10, 2026 (draft)` と付録 D 冒頭「at submission they may be moved to ancillary files」は投稿前に消す。
23. 軽微：Result 4.1(3)「up to the sign and the unknot factor」→ 符号を明記（$f_{11} = -f_{10}z_q(\dots)$ から決まる）。§9.1「happens to have the three terms of the reduced superpolynomial of the trefoil in suitable variables」は、直後の「$S = V$ for the trefoil」と重複。

## D. 参考文献

24. arXiv のページで著者・題名・関連 DOI を確認した【確認済み】：AgSh11、KN17（10.1142/S0218216520410011）、KNTZ19（10.1007/s11005-020-01318-5）、Cherednik11、LM01、CE15、Witten11、Young16、AM22（10.1007/JHEP11(2022)097）、GW05、SV00、AgSc12（10.1007/JHEP09(2012)084）。arXiv のページには Journal reference の欄が無く、巻・頁はそこからは確認できなかった。
25. KN17 と KNTZ19 は journal と `year = 2017/2019` だけで巻・頁が無く、「J. Knot Theory Ramifications (2017)」と出力される。出版年・巻・頁を入れる（数値は【未確認】。出版社のページで確認する）。
26. AgSh11、Cherednik11、LM01、CE15、Witten11 は @misc（arXiv のみ）だが、いずれも雑誌に掲載されているはず【未確認】。Young16、GW05 も掲載の有無を確認する。
27. 引用キーの過不足は無い（本文の 41 キーすべてが refs.bib にあり、未使用のエントリも無い）【コード照合】。

## E. コンパイル

28. lualatex → bibtex → lualatex ×2 で 49 ページ。`!` のエラー、未定義の参照、Overfull \hbox は無い（作業用ディレクトリで確認）。

## F. 検算して正しかった式（記録）

- §2.1 の変数の関係（$G = z_qz_t$、$\sigma$、$\tlam = -\ba\bt$、$\ba_{\mathrm{KN}} = \ba(t/q)^{1/4}$）、(2.5) の $M(q^{-1})$ と逆行列、(2.4) の Möbius 反転、(2.8) の $f_{10}$。
- Prop. 3.2 の証明、(3.4) の全 $d$ と奇数 $d$ の関係、$E|_{t=q} = -(1+q)\psi^2g^{\mathrm{odd}}$。
- §4 の $g_{11}$・$g_{1b}$ と $f_{1b}$ の整合、Result 4.3(2) と [Marino09] の対応、$(\tlam^2-q)(q-t)$ の太字変数での形。
- (5.2) の 11 項：$c = 2$ 部分の因数分解、$c = 1$ の 3 ブロックの $G$ 展開、(N)、図 2、Result 5.3 の単項式の指数 $k_1, k_2, k_3$、項数 $3m^2+7m+1$。
- Result 6.1(3) の Schur の組合せ 3 つが $p_4$・$p_2^2$ 方向を打ち消すこと（指標表）。
- §7 の genus（8, 9）、パターンと companion の語の第 1 列、6/9・8/11 の内訳。
- Conjecture 9.1：$F_{j'}$ の状態数 $\binom{2j'+3}{3}$、$X_{j'}$ の $(j'+1)^2$、$F = \psi^2X + 2P$ の恒等式。Prop. 9.3 の証明の恒等式と (9.6)、Result 9.4 の $n_{2\delta} = \Sym^{2\delta}V - (\rho_0-\rho_1)^2$。
- §10.4 の remark の mod 2 の同値、付録 B の Prop. B.2 の証明、$G = T_+T_-$。

## 総合評価

結果の中身（整数性、普遍因子、符号の性質の範囲）に関わる誤りは見つからなかった。見つかった誤りは A1〜A6 の式の書き方で、計算結果は変わらない。ただし A3（$g^{\mathrm{true}}$ の符号）は §2・§5・§8 をまたぐ規約の不整合で、読者が (5.2) を §2 の定義から再現しようとすると符号が合わない。A1・A2・A4 と合わせて、定義を追って再計算する読者が最初につまずく場所なので、投稿前に直す必要がある。

---

## G. 判定と対応（2026-11-10、続き 148）

ユーザーの指示で A1〜A6 を本文に反映した。B・C・D は未対応。

| 番号 | 対応 |
|---|---|
| A1 | setup.tex：(2.3) の符号を $(-1)^{|\mu|}$ にし、和の添字の置き直しで $(-1)^{Nm}$ が打ち消すことを 1 文で補った。 |
| A2 | setup.tex：Proposition 2.7（unknot）の積公式を $\prod_j (y_j;t)_\infty(ty_j/U;t)_\infty/((Uy_j;t)_\infty(ty_j;t)_\infty)$ に、証明の核を $(1-t^{N-k}y_j)(1-t^{k-N}y_j)$ にした。 |
| A3 | torus.tex：Hypothesis 5.1 を $f_{\yone} = -W_{\SO} = -\tfrac12 h_{\yone} + g_{\yone}$（$\tfrac12 f^{\cov}_{\yone} = -\tfrac12 h_{\yone}$）に、$g^{\mathrm{true}}_{\yone} = -(t/q)^{\delta/2}g_{\yone}$ にして、符号の取り方（$W_{\SO} - \tfrac12 h$、unknot で $+1$）を書いた。(5.2) と生成ファイルは変えていない。sp.tex：(8.2) を $\Delta = 1 - \mathrm{even}_{\ba}(g^{\SO}_{\yone})$ にし、$g^{\SO}_{\yone}$ が (5.2) の $g^{\mathrm{true}}$ であることと §2 の符号での形を注記した。 |
| A4 | generalN.tex：非規格化の値 $Z := \{\hat\gamma(P_b)\}$ を (6.1) の後で定義し、§6.1 の 3 か所と (6.2) の $JD$ を $Z$ にした。app_methods.tex：Result A.1(2) と A.2 冒頭も $Z$ に。 |
| A5 | multibox.tex：ハンドルの寄与は「our convention では $-G$」と書いた。 |
| A6 | index.tex：「two pairs exchanged by $\tilde I$, one of charge $\omega$ and one of charge $\omega t$」にした。 |

コンパイル：52 ページ、`!` のエラー・未定義の参照・Overfull なし。変更箇所（7, 9, 15, 17, 22 ページのうち 9 と 22）を画像で確認した。
日本語版（ja_v2）に同じ箇所があるかは調べていない。

---

## H. B・C・D の判定と対応（2026-11-10、続き 149）

| 番号 | 判定 | 対応 |
|---|---|---|
| B7 | 妥当 | trefoil の 2 箱で確かめた【計算】：$X_1[(1+t)\hat g_R]$ は普遍因子で割り切れ、商の分母は $2\times$ 単項式（$R = (2), (1,1)$）。$\hat g_R|_{t=q}$ の分母は $2(1+q)\times$ 単項式で、$\lambda^2 = q$ で $\hat g_R|_{t=q} = 0$。Conjecture 6.5 の偶数の場合を式 (6.5) で書き、Result 6.4 の形と同値になる条件を添えた。$T(2,5)$–$T(2,9)$・ケーブルの 2 箱で $\hat g_R|_{t=q}$ が $\lambda^2 = q$ で消えるかは調べていない。 |
| B8 | 妥当 | Hypothesis 3.8 に「unknot の定数 $-1$ は Adams 補正に入れず、$g^{\mathrm{odd}}$ の unknot の値（$p_2$ 方向で $-\tfrac12$）を $C$ から引く」ことを書き、(E2) と本文の「$+1$」の説明を合わせた。 |
| B9 | 妥当 | 「every block $n_k$ is a Laurent polynomial in $q, t$」に。 |
| B10 | 妥当 | $g_\rho := \sum_R\chi_R(\rho)g_R$ を定義して $f_\rho$ を $g_\rho$ に。 |
| B11 | 妥当 | $(1,1), (2,2), (2,3), (3,2)$ を実行して恒等式を確かめた【計算】（$(2,3)$・$(3,2)$ は各約 40 秒）。テストに $(1,1)$・$(2,2)$ を加え、本文は「4 組は自動テスト、$(2,3)$・$(3,2)$ は同じ検査の直接実行」と書いた。 |
| B12 | 妥当 | progress.md の続き 88 以降の記録どおり、$T(2,7)$・$T(2,9)$ の 2 箱と adjoint は restricted-range。表 4 の 4 行を罫線の下に移し、付録 A.2 の列挙に加えた。 |
| B13 | 妥当 | 要旨を「two-strand torus knots with up to three boxes (four for the trefoil), $T(3,4)$ in the fundamental representation」に。 |
| B14 | 妥当 | Result 7.2(4) に、$\Cab{15}$ の $(2)$ と adjoint は補間したが使っていないことを書いた。 |
| C15 | 妥当 | ケーブルのパラメータを $(\mathsf p, \mathsf q)$ に（cables.tex と表 2）。 |
| C16 | 妥当 | (B.1) の Kähler 変数を $\mathsf Q$ にし、[AM22] の $(Q, T)$ との区別を書いた。「$q/t - t/q$ in the variables of [AM22]」は削った（原典の記号は【未確認】）。 |
| C17 | 妥当 | Result C.3 に「$s = p < r$（本文の $T(3,4)$ は $r = 4$, $s = 3$）」と書いた。 |
| C18 | 妥当 | 2 文を、過去の試行を指さない形に書き直した。 |
| C19 | 妥当 | Remark 6.6 と付録 B への参照を付けた。 |
| C20 | 妥当 | §2.2 の unknot を $K_\mu(t^\rho) = \dim_\mu$ に、§2.4 に $P_R(t^\rho)$（$\U(n)$）と $K_R(t^\rho)$（$D_N$）の区別と $\Rred_R = W_R/K_R(t^\rho)$ を書いた。 |
| C21 | 見送り | 害が無く、第 1 回の対応（$\delta$ = 種数）で入れた説明なので残す。 |
| C22 | 見送り | 草稿の間は残す。投稿時に消す。 |
| C23 | 見送り | 軽微。 |
| D25 | 対応 | Crossref で確認【確認済み】：KNTZ19 は Lett. Math. Phys. 110 (2020) 2573–2583。KN17 は article no. 2041001、オンライン 2020-07-08。**KN17 の巻は Crossref では 32（issue 06、印刷 2023-05）と出るが、論文番号 2041001 の形と合わないので refs.bib には巻を入れていない。著者（ユーザー）の確認が要る。** |
| D26 | 対応 | Crossref で確認【確認済み】：Cherednik11 は IMRN 2013(23) 5366–5425、LM01 は J. Knot Theory Ramifications 11 (2002) 173–197、CE15 は Ann. Fac. Sci. Toulouse Math. 25 (2016) 433–471、Young16 は Commun. Math. Phys. 380 (2020) 273–322、Witten11 は Quantum Topol. 3, 1–137。AgSh11 は Commun. Math. Phys. 333 (2015) 187–228 だが、掲載版の題名は「Knot homology and refined Chern–Simons index」で arXiv の題名と違う（著者の一致からの同定で、arXiv 番号との対応は Crossref では確認できない）。Witten11 の年は Crossref では 2011-11-06（オンライン）で、refs.bib の 2012 は巻 3 の年としての記入【未確認】。GW05 は Crossref に掲載版が見つからず @misc のまま。 |

refs.bib は ja_v2 と共有なので、日本語版の文献表も同じに変わる。
コンパイル：53 ページ、`!` のエラー・未定義の参照・Overfull なし。21・33・52 ページを画像で確認した。tests/refined_kauffman/test_rk_koornwinder.py は 2 passed。

---

## I. 文献表を INSPIRE の BibTeX に置き換え（2026-11-10、続き 150、ユーザーの指示）

ユーザーの指示「文献は INSPIRE で調べてそこ記載の bibtex を使え」「調べ方と作り方を skill にしてそれを使え」により、.claude/skills/inspire-bibtex を作り、その手順で ja_v2/refs.bib（英語版と共有）の 37 件を作り直した。§H の D25・D26 の手入力はこれで置き換わった。

- INSPIRE の BibTeX をそのまま使ったもの（23 件）：Marino09、BFM04、BFM05、NRZ13、AgSh11、KN17、AgSc12、AM22、PU14、CD14、GV98、OV99、LM00、LMV00、PBR10、Stevan10、GS11、MV10、GN13、Young16、KNTZ19、SV00、GW05（最後の 2 件は arXiv のみで、型を @misc に替えて month を除いた）。
- INSPIRE の記録に掲載誌が無く、Crossref で補ったもの（6 件）：DGR05、Cherednik11、LM01、CE15、Witten11、ORS12。
- INSPIRE に記録が無く、手入力のエントリを Crossref（DOI）で確認して残したもの（6 件）：Mimachi01、Rains05、RW21、vD96、Sahi99、Koornwinder92。
- 書籍で未確認のまま残したもの（2 件）：Macdonald95、Cvitanovic08。

INSPIRE の記載で以前の refs.bib と違った点：
- KN17 は J. Knot Theor. Ramifications 32 (2023) 06, 2041001（Crossref と同じ値）。§H で保留にした巻は INSPIRE の記載に従った。
- SV00 の著者は「Sinha, Shishir」（以前は Sanjib と誤っていた）。
- AgSh11 の題名は掲載版の「Knot Homology and Refined Chern-Simons Index」。
- 著者名は INSPIRE の表記で、アクセント記号が無いものがある（Marino、Stosic など）。CD14 は INSPIRE の記録に arXiv 番号が無い。

コンパイル：英語版 53 ページ、日本語版（ja_v2）50 ページ。どちらも `!` のエラーと未定義の参照なし。
