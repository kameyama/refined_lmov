# 英語版（paper/en、続き 144 の修正後）の第 2 回査読・校正・レビュー（2026-11-10）

文脈を持たないサブエージェント 3 体による報告（前回の review_20261110.md は読ませていない）。対象は paper/en の main.tex・sections/*.tex・main.pdf（49 ページ）だけ。ファイルは変更していない。
判定と対応はまだ行っていない（§D に書く）。

---

## A. 査読（hep-th の査読者役）

### 【重大】

1. **再定式化行列 $M$ の定義が曖昧で、§8.1・付録 D の使い方と整合しない。** $M_{RS}=\sum_\mu\chi_R(\mu)\chi_S(\mu)z_\mu^{-1}\prod_i(t_M^{\mu_i}-t_M^{-\mu_i})/(t_M-t_M^{-1})$ は括弧の掛かり方が二通りに読める（各因子ごとに割ると $\mu=(2)$ で $z_q^+$、全体を一つの $(t_M-t_M^{-1})$ で割る形と §8.1 の「$p_2$ 成分を $q^{1/2}+q^{-1/2}$ で割る」は後者）。[Marino09] の $M$ には分母がなく、「$M(q^{-1})$」の引数の意味（[KN17] の $M(x)$ の $x$）も書かれていない。→ 一意に読める式、$|R|=2$ の $2\times2$ 行列の明示、引数の定義、[Marino09] との関係。
2. **$f_R$ の抽出手続き（Schur 展開と $d$ の Möbius 反転）が一般 $R$ について定義されていない。** $s_R(Y^d)$ の plethysm を冪和基底に移して反転する手順が無い。冪和成分 $f_\rho$ の定義、$f_R=\sum_\rho\chi_R(\rho)f_\rho/z_\rho$、(3.7) $g=g^{\mathrm{odd}}-\tfrac12 g^{\mathrm{odd}}(x^2)$ の導出（$\sum_d\frac1d g(x^d)=\sum_{d\,\mathrm{odd}}\frac1d g^{\mathrm{odd}}(x^d)$ の再編成）と「向き付け部分は全 $d$」の前提を §2.3 に置く。
3. **被覆の正規化の自由度が極の消去で固定され、同じ極の消去が被覆仮説の証拠として引用されている（循環の疑い）。** 枠付けの単項式、$(t/q)^{\delta|R|/2}$、$GL_N$ の $c$、ケーブルの $F=p_2q_2$ と genus 公式がデータ合わせ。→ 規格化の自由度と固定した条件の表を一か所に置き、「導出／$t=q$ で一致／極の消去で fit」を明記し、fit した自由度の数と独立に検証した条件の数を対比。ケーブルでは被覆の規格化が仮定であることを Result 6.8 と Table 1 に注記。
4. **Abstract・§1 の記述が本文より強い。** (i)「$G$-expansion with aligned signs」：$\lambda^{-4m},\dots$ のブロックは $(\bt^2-1)\times$ 正でそれ自体は aligned でない。(ii)「the obstruction is localized …」は本文では conjecture。(iii)「$T(2,5)$, $T(2,7)$, $T(2,9)$ up to three boxes」：$T(2,7),T(2,9)$ は 3 箱のみ。(iv)「positivity holds on the cover side」：brane–image に switching block の mixed がある。(v)「we derive」：Hypothesis 3.1 の下。→ 各文を本文の範囲に。
5. **Gukov–Walcher・NRZ13 の Kauffman superpolynomial との比較がない。** 技術的に容易なはずで「not examined [open]」で済ませるべきでない。§9.1 の「$\Rred_{\yone}$ の 9 項中 4 項が負」は [AgSc12] の refined $\SO(2N)$ 不変量が Kauffman ホモロジーの Poincaré 多項式でないことを示唆し、明示すべき。→ trefoil と $T(2,5)$ で $W^{\SO}_{\yone}$ と [GW05]・[NRZ13] を比べ、変数の対応がつくか Result として書く。
6. **[CE15] との重複と、Cherednik11 の手法に関する記述の確認。** §6.1 の composite の DAHA 計算は [CE15] と重なる。何を取り何が新しいか（一般 $N$ の補間、一様正規化、$|\alpha|=3$ の Pieri）を書く。「the shift-operator method of [Cherednik11] is limited to $t=1,q$」は同論文が一般 $t$ の hyperpolynomial を多数与えていることと整合しない。原典で確認し、誤りなら削除。

### 【中】

7. §8 冒頭「$f^{\cov}$ … does not enter the crosscap part $g$」は $g=f-\tfrac12 f^{\cov}$ と矛盾。意図は「brane/image の分解は $g$ の計算に不要」。
8. Prop 3.5：「boundaries」と「boxes」の混同（$\tilde I$ 不変の条件は巻き数の一致）、符号の不一致（$(1-q^{-1})q^{1/2}=z_q$ なので「$-z_q$」は誤記）、最後の文が fixed sector の定義とずれる。
9. Hopf link の fixed sector と $n$ の定義が不明（セクター $(2,2),(2,4)$ で何が brane/image の箱数で、なぜ $n=2$ か）。
10. 「criterion」の呼称と Conjecture 5.5 の条件の不備：2 ストランドに限った予想は「判定基準」でない。$(\bt^2-1)\times$ のブロックは $t=q$ で消えるので $\varepsilon_\beta$ が定義できない。$\lambda^{-4m},\dots,\lambda^{-(2m+4)}$ は $m=1$ で空。→ 性質 P(K) を定義し、Conjecture（全 $m$）と Result（$T(3,4)$・ケーブルで破れる）を分ける。
11. $c=1$／$c=2$ の分割の定義が複数（$\ba$ の偶奇、$\ql$ の整数冪、$z$ の偶奇、hatted では $\beta\equiv|R|-1$）。一つの定義と同値性を一か所に。
12. Conjecture 6.6 の分母：偶数セクターでは $\hat g|_{t=q}$ が $\tfrac12\psi^2 g^{\mathrm{odd}}$ を含み半整数。「in every color」は 4 箱までの検証。
13. 検証の弱い箇所の表示（$\Cab{13}$ の $(1,1)$ は $N=2,3$ だけ）が Abstract と Table 1 に無い。Table 2 に $T(2,7),T(2,9)$ の adjoint の行が無い（3 箱の $f^{\cov}$ に必要なはず。adjoint があれば 2 箱も計算できるはず）。
14. 補間法の説明が追えない（「$2Nj+$const」の論理、制限付き範囲法の仮定）。アルゴリズムとして付録 A に。
15. $V=1+t+qt$ を「trefoil の reduced superpolynomial」と呼ぶのは誤解を招く：lower series は $\kappa$ から来る結び目に依らない項。DGR 変数との対応「$a^2\to-t/q$」も曖昧。→ $V$ を $\kappa$ の Fock 空間（3 セル）として定義。
16. Result に検証事実と解釈が混在（Result 8.3 の「sign rule」「lower/upper series」が定義前、Result 8.4 の「$\tilde I$ acts as diag(1,1,−1)」は解釈、Result 9.5 の「attempts were made」）。
17. Result 8.10 と付録 F は圧縮されすぎ（$C_j$, $r'$, $D_{\delta-i}$, $\tilde D_i$, $a(m,i)$, $\Omega'$, fac, cells, background）。
18. Prop 8.11 の論理：(E2) が pair を見ないことから「(E1) だけでは固定されない」は従わない。「(E2) は制約しない［導出］、(E1) だけでは不足［計算］」に分ける。
19. §7「$g^{\mathrm{true}}_{\yone}=W_{\SO}-W^{\cov}_{\yone}$」は §5 の $f_{\yone}=-W_{\SO}=\tfrac12h+g$ と符号・$\tfrac12$ が合わない（$W^{\cov}_{\yone}$ 未定義）。Sp の被覆を SO と同じに取るのは仮説。
20. 整数計画法の定式化が一か所にない（変数、制約、目的関数、一意性の判定法）。付録に。
21. Hypothesis 3.1 は [AgSc12] (3.13) の閉弦の構造を開弦に外挿したもので結論をほぼ含む。出自を明記。Remark 3.3 は長い。
22. Remark 3.12 の「difference × sum は単一ハンドルの因子ではない」は根拠が無い。
23. Result 2.1（$D_N$ Macdonald ＝ Koornwinder の特殊化、既知）と Mimachi の恒等式の数値確認は付録 A へ。
24. 被覆の群の呼称：§1 の「$\U(N)$ invariant in composite representations」と §2 の $\U(2N-1)$。一文で説明。
25. 普遍因子 $(\ba^2-\bq^2)$：$\lambda^2=t$ は $N=1$（$\SO(2)$）で、アーベル群の refined CS は $t$ に依らないはず。解析的に示せる可能性。少なくとも unknot・Hopf・基本表現で確かめ「unexplained」を弱める。

### 【軽微】

26. 記号の多重使用（$c$, $R$, $N$, $S$, $P$, $\kappa$, $\Gamma$）。式 (6.1) の $c$ を別記号に。
27. §2.6 の証明の prefactor $t^{-(N-k)}$ は (2.3) の「no prefactor」と不整合。
28. §4.1 の $c_{2j}(n)$ の閉形式は Conjecture とするか付録へ。
29. Table 1 の「Signs」列が blocks と components で不統一、「0–2/18–22」が読みにくい。Table 2 の「e.g.」は全列挙に。
30. Result 8.9 の $T(r,s)$, $s=2,3,4$ の $s$ が strand 数。統一。
31. §1.3「Prop 8.8 holds for fourteen torus knots」：検証は Result 8.9。
32. Discussion (5) の Hilbert scheme の結果は本文に対応が無く唐突。
33. 文献の出版情報（AgSh11 CMP 333 (2015)、Cherednik11 IMRN 2013、Witten11 Quantum Topol. 3 (2012) など。要原典確認）。Paul–Borhade–Ramadevi（2010）、Stevan（2010）の追加を検討（要確認）。
34. 英語：「hallmark of a fermionic kernel」「fitted on」「Of the $n=1$ principle, the following hold also」「turns into the knot part」。
35. trefoil の $g^{\mathrm{true}}_{\yone}$（11 項）を §5 に明示。
36. 構成の重複と前方参照（Result 3.6 が §4, §8, §9 を前方参照；§9.7・§10.1・§1.3 の重複；Result 9.3 は 9.2 とほぼ重複；§4.1・§6.1 の計算技術は付録 A へ）。
37. コードとデータ：「available from the author」では再現性が担保されない。公開リポジトリと版を示す。

### 総合評価（査読者役）
定義は自然で、極の消去・整数性・普遍因子・$E/C$ 構造・$\Sp$ との関係・Prop 8.8 は実質的な貢献。Prop 3.2, 3.7, 8.8 は確認した限り正しい（Prop 8.8 の恒等式は項ごとに検算して一致）。一方で (i) $M$ と $f_R$ の抽出が再現できる形で定義されていない、(ii) 被覆の規格化の自由度と証拠が循環に見える、(iii) Abstract が本文より強い、(iv) GW05/NRZ13・CE15 との比較が欠ける、(v) §8.2 後半と付録 F が長い。major revision。1〜6 と 10〜13 が解消されれば十分な内容。

---

## B. 校正（copy editor 役）

### 【英語】
1. index.tex「A $\tilde I$-invariant curve」→「An $\tilde I$-invariant curve」。
2. hopf.tex L36 懸垂分詞「Splitting the probes …, the sector-wise cover gives …」→「When the probes are split …, the sector-wise cover gives the sum $P$ over the states in the fixed sector.」
3. hopf.tex L40「Besides it we impose two conditions, used throughout the paper:」→「In addition, we impose two conditions, which are used throughout the paper:」
4. generalN.tex L87 懸垂修飾「At fixed $N = 2, 3$, obtained through …」→「At fixed $N = 2, 3$, where the $\SO$ data are obtained from A-type data through …」
5. generalN.tex L80「(…; the seven shapes not available there by the $GL_5$ DAHA)」→「… were computed with the $GL_5$ DAHA」。
6. generalN.tex L126「\cite{CD14} define, …」→「the iterated DAHA-Jones polynomial of \cite{CD14} is defined … by」。
7. sp.tex L18「Relation \eqref{eq:sp1} and \eqref{eq:sp2}」→「Relations」。
8. sp.tex L27「satisfies $\Sp(\ba) = +\,\SO(-\ba)$」→「$X^{\Sp}(\ba) = +X^{\SO}(-\ba)$, whereas for the fundamental representation the sign is $-$」。
9. multibox.tex L79 Result 9.4 項目 3 は主動詞の無い断片 → 書き直し案あり。
10. cover.tex L93 Prop 8.12 の「Hence」が論理的に繋がらない →「Hence (E2) imposes no constraint on these signs; since (E1) alone does not fix them (computation), …」。
11. discussion.tex L6「has no $N$-dependent poles」→「is free of $N$-dependent poles」。
12. 「needs $N \le 13$ for the interpolation」→「needs the values for all $N \le 13$」。
13. intro.tex L24 →「determined uniquely under Hypotheses \ref{hyp:QP} (numbers of states) and \ref{hyp:E2} (condition (E2))」。
14. 数式記号で始まる文 4 か所（「$\dim_R$ is …」「$f_R$ closes …」「$(q-1)(t-1)$ is …」「$\tilde I = …$, and」）。
15. 「i.e.\ 」/「i.e.,」、「e.g.\ 」/「e.g.,」の混在 →「i.e.,」「e.g.,」。
16. 「iff」（2 か所）→「if and only if」。
17. 「double twist knots」→「double-twist knots」。「closed string free energy」等 →「closed-string」「open-string」。
18. Oxford comma の不統一。
19. torus.tex L41 Result 5.3 (c=1) の「$\lambda^{-4m}, \dots, \lambda^{-(2m+4)}$」は $m=1$ で空 →「(absent for $m = 1$)」。
20. intro.tex L75 \open は討論以外でも使う →「Throughout, \open\ marks an open point; …」。

### 【用語・記法の統一】
21. ヤング図の [2] と (2) が混在（Table 2、Result 6.2、Result 6.1）。[·] は引用ラベルと紛れるので (·) に統一し、§2.1 の宣言も直す。
22. 「$\ba$-blocks」/「$\ba$ blocks」/「$\lambda^\beta$ block」→「$\ba$ block」に統一し、$\ba=\lambda$ なので同じものと §2.1 で断る。
23. 「general $N$」と「large $N$」の使い分け：補間で厳密に得た量は「general $N$」、形容詞は「large-$N$」、題名の「large $N$ duality」だけハイフン無し。
24. 「switching block」/「switching term」/「the $n = 1$ switching」→「switching block」。
25. $S\otimes S$ / $S^2$ / $S_m^2$ → 統一。
26. 「(E1)(E2)」「(E1)+(E2)」「(E1) and (E2)」「(a)(b)」→「(E1) and (E2)」。「mod 2」と「\pmod 2」。
27. $M^{-1}(q^{-1})$（hopf.tex）と $M(q^{-1})^{-1}$ → 統一。
28. mirror の演算子表記 → \DeclareMathOperator{\mirror}{mirror}。
29. setup.tex「$([1],[1])$」→「$(1;\bar1)$」。
30. sp.tex「$C_n$」→「$C_N$」。$GL_N$ をローマン $\mathrm{GL}_N$ に。

### 【記法：多重定義・未定義】
31. 多重使用の改名：(a) 定数 $c$ → $\mathsf c_K$；(b) $z$（Mariño の $z$、Hopf の第 2 成分の probe 変数、母関数の形式変数）；(c) $S$（superpolynomial、Result 9.2 の $S$、$S$ 行列、$S_R$、$S_d(x)$）；(d) $U$（$t^N$、$U_n := t^n$、upper series $U_j$、多項式 $U_n$）；(e) $Q$；(f) $B$（$\Bim$、$B_{12}$・$E_B$、$B_N$）；(g) $\epsilon_R$ と $\varepsilon$ → 言葉で；(h) $w$；(i) $\rho$；(j) $g$ と $G^g$ の $g$；(k) $X$；(l) $T$。§2.1 に「節ごとに再定義する記号」の注記。
32. 未定義：$F^\varepsilon$、$I_X$（定義は $I$）、$\mathrm{fac}$、$c_1$（App. F）、Prop 8.8 証明の $q_n$。
33. hopf.tex「the Hopf link of $\U(N)$ refined Chern–Simons theory」：cover は $\U(N_A)$。
34. hopf.tex「virtual normalization $\tilde Z = y^N z^N Z$」が (eq:Z) の $\tilde Z$ と同じ記号で別物 → $Z^{\mathrm{virt}}$。
35. setup.tex L62「$\prod_i(1 - y_\alpha x_i)\cdots$」添字 $\alpha$ が自由 → $\prod_{\alpha,i}$。
36. generalN.tex L142「$(p,q)$-cable」の $q$ が変数と衝突 →「$(p_2,q_2)$-cable」。
37. cover.tex Result 8.7「the term $j = 2\delta$」→ ブロック添字は $k = 2\delta$。Result 8.9・8.10 の $T(r,s)$ の $s$ が strand 数で本文の「three-strand $T(3,4)$」と逆 → 明記か入れ替え。
38. Result 8.2 は $0\le j\le 2m-2$ で $R_j$ を置くが Conjecture 8.3・8.5 は $j = m + j' \ge m$ だけ。$j < m$ で $R_j = 0$ ならそう書く。
39. index.tex L43「normalize by $q^{1/2}$ to $-z_q$」：$q^{1/2}(1-q^{-1}) = +z_q$。符号を確認。
40. sp.tex L10「$g^{\mathrm{true}}_{\yone} = W_{\SO} - W^{\cov}_{\yone}$」と Hypothesis 5.1 の符号・規格化の対応を一文で補う。
41. Result 4.5「the ten sectors with $|R| \le 4$, $|S| \le 2$」：組は 33 あり、どの 10 か再構成できない。条件を明記。「4 of the 25 pairs」→「$|R| = |S| = 4$ (25 pairs)」。

### 【LaTeX・参照】
42. generalN.tex §6.4 見出し「at $t = q$」が \texorpdfstring 外（hyperref 警告 2 件の原因）。
43. 「$E_{\ytwo} = -E_{\yoneone}$ of Proposition~\ref{prop:alld}」（Result 6.3 項目 1、Result 6.11 項目 2）は Remark 3.3（rem:Esector、未参照）に。
44. intro.tex L29「(Proposition \ref{prop:rows}) holds for fourteen torus knots」→「(Proposition \ref{prop:rows}, Result \ref{cr:rows})」。
45. index.tex L73「(\S\ref{sec:positivity})」→ \S\ref{sec:EC}。hopf.tex L40・L59 も。
46. 未参照のラベル：app:largeN（本文から参照する）、hyp:cablecover（討論の仮説の列挙に加える）。
47. Table 1 の p 列が両端揃えで Underfull 多数 → `>{\raggedright\arraybackslash}p{…}`。
48. 生成ファイルの先頭が導入文と同じ段落に続く（付録 B・C・D・E・G）→ \input の前に空行。
49. Table 5 が $T(2,3)$ のブロック一覧の途中に割り込む → [t] か \FloatBarrier。
50. Result 4.2 が enumerate で始まる → 導入文。
51. 生成付録 C：行ラベルが未定義の $[A^{-5}]$, $[A^{-9/2}]: 0$ → $[\lambda^{-5}]$ 等、0 の半整数行は削るか注記。
52. 生成付録 E：$\Nfix$ を $N$、$\Bim$ を $x$ と書いている → 生成側で \Nfix・\Bim を出す。付録 D の $X$ は $X_1$。
53. Table 1 最終行「$W_R$: $(-1)^{|R|}$」→ 式で。Table 2 の cover $m = 2,3,4$ の行は $m$ ごとに分ける。
54. 表示式末尾の句読点（式 (11)、生成付録 D の末尾）。
55. hopf.tex L4 の段落が Underfull。

### 【数値・範囲の整合】
56. abstract「for $T(2,5)$, $T(2,7)$, $T(2,9)$ up to three boxes」→ $T(2,7)$・$T(2,9)$ は 3 箱のみ。
57. Result 6.10 項目 4「with the data for $N \le 13$ three degrees remained undetermined」と Table 3「nodes $N = 4$–14」：$N = 14$ で解決したと明記。
58. Result 6.11 項目 3 の内訳から 6/9・8/11 が出ることを一文で書く。他の数値は一致を確認済み。
59. §10.1 項目 3「[computation for $T(2,3)$, $T(2,5)$; …]」に $T(2,7)$ を加える。
60. intro「fitted on $T(2,3)$–$T(2,9)$ … predicted for $T(2,11)$–$T(2,15)$」と Result 8.2「for $T(2,7)$–$T(2,15)$ … as a prediction」を揃える。

### 【参考文献】
61. [AM22]「M. Y. Avetisyan and R. L. Mkrtchyan」と [NRZ13]「P. Ramadevi」だけイニシャル。Mkrtchyan は [MV10] でフルネーム。
62. DOI はあるが誌名・巻・頁が無い：[AgSc12]（JHEP 09 (2012) 084）、[AM22]（JHEP 11 (2022) 097）、[KN17]、[KNTZ19]。
63. arXiv のみだが掲載誌があると思われるもの（要原典確認）：[AgSh11]、[Cherednik11]、[CE15]、[Witten11]、[Young16]。
64. JHEP の巻表記が混在。
65. [Rains05] 題名「BCn-symmetric」の下付き → `{$BC_n$}`。
66. doi=true なのに DOI の無い掲載論文（[Marino09] ほか）。[Koornwinder92] は編者・出版社が無い。「Experimental Math.」→「Experiment. Math.」。
67. shorthand の年と掲載年のずれ（方針として可）。

### 総評（校正役）
文法上の重大な誤りは少ない。主な直し所は懸垂分詞 2 件、断片文、「A $\tilde I$」、引用ラベルを主語にした文、Prop 8.12 の「Hence」。記号の多重使用（$c$, $z$, $S$, $U$, $Q$, $X$, $T$, $\rho$, $\epsilon/\varepsilon$）が最も目立つ。相互参照は解決しているが、Remark 3.3 に帰すべき 2 か所、Result 8.9 に帰すべき序論、§9.5 を指すべき箇所、未参照の付録 D・Hypothesis 6.9 を直す。数値は整合（唯一 abstract の「up to three boxes」）。組版は Table 1 の両端揃え、生成付録の段落結合、Table 5 の割り込み、見出しの hyperref 警告、ヤング図の混在。参考文献は著者名形式と誌名の補完。

---

## C. レビュー（編集者役：構成・序論・環境・記法・慣習）

### 【構成】
1. 本文に「どう計算したか」が多く混ざり、構造的な主張が埋もれている（§4.1、§6.1 の規約、§6.2 と Table 2・3、Result 6.1、Result 7.1 項目 1）。Appendix A に移し、§6 は 5 頁→約 3 頁に。
2. (E1)(E2) と Hypothesis 4.3・4.4 が Hopf の節に埋まっている → §3.4 の直後に §3.5「Conditions fixing $E$」。Result 3.5（boundary factor の表）は §4 の末尾に。
3. §3.5 Closed strings（Result 3.8、Prop 3.9、Remark 3.10）は開弦の論旨から外れる → 一段落残して付録へ。Remark 3.3（trace form）も長い。
4. $T(3,4)$ の結果が三か所に散っている（§5.3、Result 6.7(4)、Table 1）→ §5.4「$T(3,4)$: failure of ($c{=}1$)」を新設。
5. §6.6 ケーブルは独立の節に昇格。
6. §7 Sp は小節が無く流れを切る → §6 の後に置いて本文を絞るか、§9 の後へ。p.23 の長い段落は Remark に。
7. §8 は最も密度が高い → 本筋（brane–brane = [KN17]、Result 8.2、Prop 8.9、Result 8.10）だけ残し、Result 8.7・8.8・8.10 後半・8.11・Prop 8.12 は付録「Further structure of the fixed sector」へ。§8 冒頭の導入は小節 8.1「Setup」に。
8. Result 8.2 と Conjecture 8.6 がほぼ同じ式を二度書いている → Conjecture を正式な主張とし、Result は検証だけに。
9. §9.5 の否定的な試行のリストは一段落と framed unknot の式と Remark だけ残し、詳細は付録へ。
10. 要約が四重（§1.3、§9.7、§10.1、§10.4）→ §9.7 削除、§10.1 は表か削除、§10.4 は §10.3 の前書きに吸収。
11. 付録 B〜E の自動生成の係数リスト 12 頁は ancillary files に。残すなら例だけ、二段組か表に。
12. 肯定・否定のバランスは妥当。要旨の否定的結論は一文に。
13. §2.1 が見出し直後にいきなり表。§9.1 の一般論が長い。

### 【序論・要旨】
14. 要旨が未定義の記号・用語で書かれている（$z_q^{n-1}(z_q^+)^n$、$\sigma$-symmetric、$\lambda^\beta$ block、$G$-expansion、frozen states、brane–image、$2(1+t)$）→ 150〜200 語に。
15. 節への道案内が二重（§1.2 末尾と §1.4 第 2 段落）。
16. 格付けの説明（§1.4）が Table 1 と §1.3 で Result 等を参照した後に来る → §1.3 の前へ。
17. Table 1 は独立に読めない：未定義語、「0–2/18–22」、Signs 列が Table 4 と重複、p.6 に単独。→ Object / Integer / Universal factor / Criterion の 4 列で記号化し、符号統計は Table 4 に一本化。
18. §1.3 の各項目が長い。Definition の項の「naive refinement fails」は動機。Hopf の項は主結果として弱い。
19. §1.1 に crosscap invariant とは何か、なぜ refinement で正値性が期待されるかの説明が無い。$\Cab{13}$ の記号は §1.3 で定義。
20. 図が一枚も無い → orientifold の模式図、$\lambda^\beta \times G^g$ の符号の格子、lower/switching/upper の模式図。

### 【環境と格付け】
21. Result の使い方が三通り（肯定的主張、実装のチェック、否定的な試行の列挙）→ 「検証範囲を明記した肯定的主張」に限定。
22. Remark 2.2 の「None of five choices …」は計算機による検証 → Result（negative）か [computation] の印。本文中の「(a conjecture)」「(an interpretation)」「(computation)」の括弧書きを統一。
23. Conjecture 5.6 の中の解説文を外へ。(c=2) の一般形と $T(2,2m+1)$ の場合を分ける。
24. (E1)（導出）と (E2)（仮説）に格付けを付記。
25. 長すぎる列挙（Result 6.3、6.7、6.8、8.2、F.1・F.3）→「主張 1〜2 文 ＋ 検証範囲 1 文」。
26. Prop 8.12 は Remark に。§2.4 の「被覆は A 型」を番号付き Hypothesis に（Hypothesis 5.1 はその基本表現版）。
27. §2.1 の直後に「Terminology」の小節（ブロック、$G$ 展開、aligned、refined part）。

### 【記法・用語】
28. 文字の重複（$c$, $P$, $Q$, $X$, $R$, $S$, $T$, $\beta$, $n$）→ 少なくとも $c$, $P$, $Q$, $X$ を改名、付録に記号表。
29. §8 の $r = q^{1/2}$, $s = t^{1/2}$ が $T(r,s)$ と同居 → $q^{1/2}, t^{1/2}$ で書く。
30. 変数系の多さ（$(q,t,\lambda)$、bold、$(\ql,\qr)$、$\tlam$、$a$、$\ba_{\rm KN}$、$U=T$、[AM22] の $(Q,T,X)$、§8 の $(r,s,\xi,y)$）→ 各節の冒頭で宣言。Result 5.3 の $\lambda^\beta$ と $\ba$ の混在。
31. jargon の初出定義（switching block、rows、background・cells、Lefschetz weight、hatted、uniform normalization）。
32. 表記の揺れ（$([1],[1])$、「Appendix A」と「§3」の混在、[2]/(2)/\ytwo）。
33. 具体例の不在：§5 に trefoil の $g^{\rm true}_{\yone}$ の分解と符号を例示。

### 【hep-th の慣習】
34. 式番号：参照される表示式の多くが番号なし → 全表示式に番号。
35. Acknowledgements、Data/code availability、著者の所属・メール、pdftitle/pdfauthor、「(draft)」。
36. 引用の様式：見出しに引用ラベルを置かない（「Dictionary with [Marino09]」「Relation to [KN17]」）。「the $\lambda$ of [NRZ13]」「(6.15) of [AgSc12]」「Theorem 2.4 of [CD14]」は人名でないと読めない。方針自体も一貫していない（Ooguri–Vafa、Koornwinder、Askey–Wilson、Cherednik–Danilenko は名前で書いている）。→ 固有名詞化した構成は名前＋\cite で書いてよいと規則を緩め、見出しからラベルを外す。
37. Table 5 の「sign per $\lambda^\beta$」列は読めない → $\beta$ ごとに列を立て記号で。Table 1 も記号化。Table 2 の「nodes」列は表現ごとに行を分ける。
38. 太字「[open]」は見慣れない →「remains open」か脚注。[derivation]/[computation] は一覧表に。
39. 節の粒度の不均衡（§7 小節なし、§9 は 7 小節）。
40. 文献の雑誌情報（校正と同件）。

### 提案する目次（編集者役）
1 Introduction（1.1 Background and questions、1.2 Grading、1.3 Main results）／2 Refined crosscap invariants（2.1 Variables and terminology、2.2 SO(2N) invariants and the dual Cauchy kernel、2.3 Definition and scheme（被覆は A 型の Hypothesis）、2.4 Uniform normalization、2.5 Unknot）／3 The orientifold index（3.1 Setup（Figure 1）、3.2 All-d structure and trace form、3.3 Boundary zero modes、3.4 E–C ambiguity and (E1)(E2)、3.5 Closed strings（一段落））／4 Unknot and Hopf link／5 Torus knots, fundamental: the criterion（5.4 T(3,4)）／6 Higher representations at general N／7 Algebraic cables of the trefoil／8 Refined Sp(2N) and the SO(−N) correspondence／9 Structure of the cover／10 Multi-box crosscap invariants: obstructions（4 小節）／11 Open problems／Acknowledgements, Data availability／A Computational methods／B Closed strings／C Further structure of the fixed sector／D Explicit invariants: trefoil／E Unrefined check

### 総合所見（編集者役）
筋は通っており格付けの姿勢は好ましい。最大の問題は密度と重複（計算法が本文に混在、要約が四重、生成データ 12 頁）で、再配置で本文 31→24〜25 頁、全体 35 頁前後にできる。要旨と Table 1 は非専門家向けに書き直す。図が無いこと、文字の重複、§8 の $r,s$ の衝突は費用が小さい。見出しから引用ラベルを外す、全表示式に番号、Acknowledgements とデータ公開の節、著者の所属は投稿前に必ず。

---

## D. 第 2 回の指摘の妥当性の判定と対応（2026-11-10、続き 145）

判定は「妥当」「部分的に妥当」「不当」「不採用（方針）」。対応は paper/en の全面的な再構成（下の目次）に含めた。

### 査読（A）

| # | 判定 | 対応 |
|---|---|---|
| A1 | 妥当 | $M(x)$ を分母を一つにした式 (2.7) で書き直し、冪和の基底で対角（$z_q$, $z_q^+$）であること、Schur 基底の 2×2 行列 (2.8)、[Marino09] の行列を 1 箱の因子で割ったものであることを明記。コードの実装（$M(q^{-1})^{-1} = [[r, 1/r],[1/r, r]]/(r^2 - r^{-2})$）と一致することを確認。ja_v2 の定義も同じ形に直した。 |
| A2 | 妥当 | §2.2 に冪和成分 $f_\rho$、$L_\rho$ からの Möbius 反転 (2.6)、$f_R = \sum_\rho \chi_R(\rho) f_\rho / z_\rho$ を書いた。(3.7)（odd-$d$ との関係）は §3.4 で「向き付け部分は全 $d$、クロスキャップは odd $d$ という [Marino09] の組合せを全 $d$ の展開と比べると偶数 $d$ の項が Adams 作用素で表せる」と説明。 |
| A3 | 妥当 | §2.4 に Table 2「規格化の自由度・固定した条件・その後の独立な検証」を置き、fit した量が少数で、整数性は固定に使っていないことを明記。ケーブルは Hypothesis 7.1 と経験的な定数 $\mathsf c_K$ に依ることを §7 末尾に明記。 |
| A4 | 妥当 | 要旨・§1.3 を本文の格付けに合わせた：「aligned up to $(\bt^2-1)$」、「we conjecture that the obstruction …」、「$T(2,7), T(2,9)$ with two and three boxes」（2 箱も計算して追加）、「on the cover side」の表現を brane–brane と brane–image に分けて記述、「derive … under a classification of one-particle states modelled on [AgSc12]」。 |
| A5 | 妥当だが新しい計算 | Gukov–Walcher の超多項式との比較は行っていない（別作業）。§5.3 と討論に未検討と明記し、$\Rred_{\yone}$ が bold 変数で正でないことから「この変数では Kauffman ホモロジーの Poincaré 多項式でない」と明示した。 |
| A6 | 妥当 | 「shift 作用素の方法は $t = 1, q$ に限る」の文は原典で確認できないので削除。§6.1 に [CE15] との関係（composite の DAHA 計算を一般 $N$ の補間・一様な規格化・$|\alpha| = 3$ に広げた）を書いた。 |
| A7 | 妥当 | §9.1 冒頭を「brane/image の分解は $g$ の計算に不要」に書き直した。 |
| A8 | 妥当 | Prop 3.5 を箱数（巻き数）の一致で述べ、証明を補い、$q^{1/2}(1-q^{-1}) = z_q$ に直した。 |
| A9 | 妥当 | §4.3 に Hopf の fixed sector（各成分に brane 1 箱・image 1 箱）と $n = 2$ の意味を書いた。 |
| A10 | 妥当 | 「criterion」を「sign property」（Definition 5.7）とし、Conjecture 5.8 は「全ての $T(2,2m+1)$ が sign property をもつ」、$T(3,4)$・ケーブルでの破れは Result（5.9、7.3）に。$(\bt^2-1)\times$ のブロックは $t = q$ で消えるので $\varepsilon_\beta$ を課さないと明記。$m = 1$ で空の範囲に「absent for $m = 1$」。 |
| A11 | 妥当 | Definition 5.4 で [1] の $c = 1/2$ を $\ba$ の冪の偶奇で定義し、(P) により $\ql$ の整数冪と同値、hatted の多箱では Result 10.4 の規則、と一か所にまとめた。§2.1 の用語の一覧にも入れた。 |
| A12 | 妥当 | Conjecture 6.6 を奇数・偶数のセクターに分け、偶数では $2(1+t)\hat g_R$ について述べ、$\hat g|_{t=q}$ の分母 2 を明記。検証範囲を添えた。 |
| A13 | 妥当 | Table 1 と §7 に $(1,1)$ の照合が $N = 2, 3$ だけであることを注記。$T(2,7)$・$T(2,9)$ の随伴と 2 箱は計算済みだったので表と本文に加えた（2 箱の和の符号 10/26、14/34 は新しい計算、logs/refined_kauffman/b_2box_structure_m{3,4}_flint.log）。 |
| A14 | 妥当 | 付録 A.2 に補間の手続き（傾きから冪の範囲、Newton 補間、制限付き範囲法の帯の仮定と一意性の判定）を書いた。 |
| A15 | 妥当 | $V$ を $\kappa$ から来る 3 セルの Fock 空間として定義し、trefoil の超多項式との一致は「適当な変数で同じ 3 項」にとどめ、「$a^2 \to -t/q$」の記述を消した。 |
| A16 | 妥当 | Result 8.3（sign rule）と 8.4 を、Conjecture 9.1（閉形式と $\tilde I$ の作用）＋ Result 9.2（検証範囲）に組み替え、「attempts」は付録 C.2 の negative Result に。 |
| A17 | 妥当 | Result 8.10 と付録 F を付録 C に移し、$C_j$, $r'$, $[n]_x$, $\mathcal U_j$, $\Omega^{(0)}$, cells, background を定義。 |
| A18 | 妥当 | Prop 8.11 を Remark C.5 にし、「(E2) は制約しない［導出］、(E1) だけでは不足［計算］」に分けた。 |
| A19 | 妥当 | §8 冒頭に DAHA 側の規約での $g^{\mathrm{true}} = W_{\SO} - W^{\cov}$ と $W^{\cov}_{\yone}$（$\tfrac12 h$ の鏡像）の関係、Sp の被覆を Hypothesis 8.1 に。 |
| A20 | 妥当 | 付録 A.4 に整数計画の定式化（未知数、制約 (i)〜(iv)、一意性の判定、対の最小数）を書いた。 |
| A21 | 妥当 | Hypothesis 3.1 に「[AgSc12] (3.13) の閉弦の分類を開弦に外挿」と明記。Remark 3.4 を短縮。 |
| A22 | 妥当 | Remark B.4 を「In our reading, …」と弱めた。 |
| A23 | 妥当 | $D_N$ Macdonald ＝ Koornwinder の特殊化を引用で述べ、数値確認は付録 A の実装の検証に移した。 |
| A24 | 妥当 | §2.3 に「[Marino09] の $\U(N)$ の composite 不変量を $\lambda$ をずらして評価したものが $\U(N_A)$ の被覆（Result 2.5 の行列式のずらしの不変性）」と一文。 |
| A25 | 妥当だが新しい導出 | 解析的な導出は行っていない。Result 6.4 の後と討論に「$N = 1$ への形式的な外挿として示せるかは未検討」と明記。 |
| A26 | 妥当 | 定数 $c$ → $\mathsf c_K$。他の重複は §2.1 の注記と改名（$P \to \Phi$、$Q \to \Phi_-$、$B \to \mathcal B$、$X$ の衝突の解消、$U_n$・$z$ の衝突の解消）。 |
| A27 | 妥当 | §2.5 の証明から前因子 $t^{-(N-k)}$ を外した。 |
| A28 | 妥当 | $c_{2j}(n)$ の閉形式を Conjecture A.3 として付録 A.3 に。 |
| A29 | 妥当 | Table 1 を記号化（Integer / Universal factor / Sign property）し、符号の数は Table 3 に一本化。「e.g.」を消した。 |
| A30 | 妥当 | $T(r,s)$ の $s$ をストランド数と明記（付録 C.1）、Result 9.4 は $p$ をストランド数として書いた。 |
| A31 | 妥当 | 序論の参照を Result 9.4 に。 |
| A32 | 妥当 | Hilbert scheme の結果を討論の 2 文に圧縮。 |
| A33 | 部分的に妥当 | arXiv で確認できたものだけ補った（CD14、AgSc12・AM22 は DOI から誌名・巻・番号、KN17・KNTZ19 は誌名、PBR10・Stevan10 を追加して序論で引用）。AgSh11・Cherednik11・Witten11 は arXiv に掲載情報が無いので未確認のまま。 |
| A34 | 妥当 | 語を直した（「hallmark」は残した）。 |
| A35 | 妥当 | §5.2 に trefoil の 11 項の式 (5.2) と $c = 1/2$ への分解、Figure 2 の符号の格子を置いた。 |
| A36 | 妥当 | 構成を変更（下の目次）。§9.7 の要約を削除、§10.1 を表に、§10.4 を討論の末尾に吸収、計算技術を付録 A に。 |
| A37 | 部分的に妥当 | 公開先は未定なので「Data and code availability」の節に「著者から入手可能、最終版で公開リポジトリを添える」と書いた。 |

### 校正（B）

- B1〜B20（英語）：すべて妥当。再構成の際に反映（B10 は Remark C.5 で論理を分けた、B19 は「absent for $m = 1$」、B20 は \open の説明を全体に）。
- B21〜B30（統一）：妥当。ヤング図は本文で $(2), (1,1)$ に統一、「$\ba$ block」に統一（$\lambda^\beta$ block と同じと §2.1 で断る）、general $N$／large-$N$ の使い分け、switching block、$S\otimes S$、(E1) and (E2)、$M(q^{-1})^{-1}$、\mirror 演算子、$(1;\bar1)$、$C_N$、$\GL_N$。
- B31〜B41（多重定義・未定義）：妥当。B31 の改名は上記。B32 の $F^\varepsilon$・$I_X$・fac・$c_1$・$q_n$ を言葉で置き換えた。B33 は $\U(N_A)$ に。B34 は $Z^{\mathrm{virt}}$。B35 は $\prod_{\alpha,i}$。B36 は $(p_2, q_2)$。B37 は $k = 2\delta$、ストランド数を $p$ に。B38 は「$j < m$ の $R_j$ は (9.3) で与えない」と明記。B39 は $+z_q$。B40 は §8 冒頭の文。B41 は「refined part が 0 でない 10 成分」（14 成分中 (2),(1,1) と (1,1),(2) など 4 つが 0）と「$|R| = |S| = 4$ の 25 組」。
- B42〜B55（LaTeX・参照）：妥当。見出しの \texorpdfstring、Remark 3.3 への参照、Result 9.4 への参照、§10.4 への参照、付録 D.3 と Hypothesis 7.1 の参照、Table 1 の \raggedright、生成ファイルの前の \medskip\noindent、Table 5 を [t] に、Result 4.2 の導入文、生成付録 C の $[\lambda^{-5}]$、生成付録 E の $\mathcal N$・$\mathcal B$・$X_1$（生成スクリプトの規則で対応）、Table 1 最終行を式で、cover の $m = 2,3,4$ の行の「e.g.」を消去、式 (8.3) の句読点。B55 は段落を書き直した。
- B56〜B60（数値）：妥当。要旨の「up to three boxes」、$N = 14$ で決まったことの明記、6/9・8/11 の内訳、§10.1 の $T(2,7)$、「fitted on $T(2,3)$–$T(2,9)$」の範囲を Result 9.2 で説明。
- B61〜B67（参考文献）：B61 は arXiv にフルネームが無いので未対応（AM22・NRZ13）。B62 は DOI から AgSc12・AM22 の誌名・巻・番号、KN17・KNTZ19 の誌名を補った（巻・頁は未確認）。B63 は arXiv に無く未確認のまま。B64 は JHEP の巻表記を統一。B65 は既に `{$BC_n$}`。B66 は Koornwinder92 に出版社、Experiment. Math. に。B67 は方針として可。

### レビュー（C）

- C1〜C13（構成）：妥当。目次を組み直した：1 Introduction（背景と問い・格付け・主結果・構成）／2 Refined crosscap invariants（変数と用語、SO(2N) と双対 Cauchy 核、定義と計算の手順、規格化（Table 2）、unknot）／3 The orientifold index（設定（Figure 1）、全 $d$、境界ゼロモード、$E$–$C$ の不定性と (E1)(E2)、閉弦（一段落））／4 Unknot and Hopf link（境界因子の表を含む）／5 Torus knots, fundamental: the sign property（5.4 $T(3,4)$）／6 Higher representations at general $N$／7 Algebraic cables of the trefoil／8 Refined Sp(2N) and the SO(−N) correspondence／9 Structure of the cover（Setup、分割、$n = 1$ の状態空間：Conjecture → Result、Prop rows、Figure 3）／10 Multi-box crosscap invariants: obstructions（4 小節＋composite）／11 Discussion and open problems／Data and code availability／A Computational methods／B Closed strings／C Further structure of the fixed sector／D Explicit invariants（生成データ、E の unrefined check を D.5 に）。C11 の ancillary file 化は投稿時の判断として付録 D 冒頭に注記（部分的に妥当）。
- C14〜C20（序論・要旨）：妥当。要旨を 200 語程度に、道案内を一本化、格付けを主結果の前に、Table 1 を記号化、§1.3 を (Q1)〜(Q5) の答えとして短縮、§1.1 に crosscap number と「refinement が $(-1)^F$ の相殺を解く」の説明、$\mathcal C_{13}$ の定義、図を 3 枚（orientifold の模式図、符号の格子、lower/switching/upper の模式図）。
- C21〜C27（環境と格付け）：妥当。Result は肯定的主張に限り、否定的な結果は「negative」と明記した Result に、実装の確認は付録 A へ。Remark 2.2 → Result 2.1（negative）。Conjecture 5.8 の解説文を外へ、(E1)(E2) に格付けを付記、長い列挙を分割、Prop 8.11 → Remark、§2.3 に一般形の Hypothesis 2.3（被覆は A 型）、§2.1 に Terminology。
- C28〜C33（記法・用語）：妥当。改名は上記。§9 の $r, s$ を $q^{1/2}, t^{1/2}$ に。各節の冒頭で使う変数系を宣言。jargon の初出定義。trefoil の例を §5.2 に。
- C34〜C40（慣習）：C34 妥当（全表示式に番号、\numberwithin）。C35 部分（Data and code availability の節、pdftitle/pdfauthor を追加。Acknowledgements は書く内容を著者が決めるので未追加、所属も未記入、「(draft)」は準備稿なので残す）。C36 部分：見出しから引用ラベルを外した（「Dictionary with the unrefined theory at $t = q$」「Relation to [KN17]」は段落題に）。本文で著者名を書かない方針は著者の規則（CLAUDE.md）なので維持し、構成の固有名詞（Ooguri–Vafa、Koornwinder、Askey–Wilson、DAHA-Jones）だけ名前で書く現状を続ける（不採用：方針）。C37 部分：Table 1 は記号化、Table 2（補間）は \resizebox で収めた。Table 5（生成）の $\beta$ ごとの列化は行の長さが結び目ごとに違うので見送り（不採用）。C38 妥当（\open を斜体の (open) に）。C39 妥当（再構成で解消）。C40 は B62 と同件。
