<!-- 出典：ユーザーが deep research（文献調査ツール）で得た報告（2026-10-14、~/Downloads/refined_kauffman.md）。本文は未検証の外部の要約で、書誌と主張は原典で確認してから使う。確認済みの点は literature.md と progress.md に記す。 -->

# refined Kauffman 多項式の large N duality：§1 の「極の不一致」を既存文献で解けるか（2026年10月時点の調査）

**結論：§1 の障害（D_N の unknot の値がもつ N 依存の分母と、composite 側の分母の位置が合わないこと）を直接扱った文献は見つからず、既存文献だけでは解決しない。** ただし、障害を回避または再定式化するのに使える部品が三つの系統にある。(i) Aganagic–Schaeffer の第二量子化（積表示）による orientifold の整数性、(ii) Kameyama–Nawata と Cheng–Sułkowski の refined 開弦展開（分母を Macdonald 基底・計量・refined holonomy に吸収する）、(iii) Rains の lifted/virtual Koornwinder 多項式と、その lifted Cauchy 恒等式。T（= t^n）に関する極が明示されている。

## TL;DR
- **Q1・Q6（refined かつ開弦の orientifold 整数性）：該当文献は無い。** Aganagic–Schaeffer（arXiv:1202.4456）は、refined の場合について「第一量子化の定義はまだ無い」と明言している。彼らが与えた refined の整数構造は閉弦の積表示 (3.13) だけである。2019〜2026\[1\] 年の SO/Sp 関連の refined 研究は、閉弦（Avetisyan–Mkrtchyan arXiv:2205.12832）と 5d ゲージ理論の頂点（ON 面・O 頂点）に限られる。refined Kauffman の LMOV 型整数性を定式化した文献も、その不可能性を論じた文献も見つからなかった。
- **Q2（分母を許す refined BPS 展開）：U(N) では確立している。** Kameyama–Nawata（arXiv:1703.05408）は、unknot の分母 Π(1 − q^{a}t^{l+1}) を Macdonald 基底と計量 g_λ に吸収した形で refined LMOV を定式化している。分母は彩色ごとに有限で、N（a）に依らない。Cheng–Sułkowski（arXiv:2104.00713）は同じ分母を\[2\] refined holonomy の置換 (q,t)_n/(t,t)_n に吸収している。本件の困難は、D_N\[3\] の分母が N（Kähler 変数）に依存する点にある。この型を扱った文献は無い。
- **最も有望な次の一手（推測）：** Mariño 型の母関数を Rains の lifted Cauchy 恒等式 Σ(−1)^{|λ|} \hat K_λ(x;t,q,T;t₀…t₃) \tilde K_{λ'}(y;q,t,T;t₀…t₃) = Π(1 − x_i y_j) の上で組み直す。そのうえで、T に関する極（Rains–Warnaar (2.6.6) に明示）を「許される分母」として扱う。これが Kauffman 側の自然な refined Ooguri–Vafa 演算子の候補である。補助線として Avetisyan–Mkrtchyan の「クロスキャップの数に依存する Kähler 変数のずれ」も使える。

## 変数の規約（論文ごとの対応）
| 文献 | q, t の定義 | Kähler 変数 | 備考 |
|---|---|---|---|
| Aganagic–Schaeffer 1202.4456 | q = exp(2πi/(k+βy))、t = exp(2πiβ/(k+βy))（y は双対コクセター数、SO(2N) では 2N−2）。t = q が unrefined | 閉弦では t^{(2N−1)d}。結び目変数は **q** = √t、**t** = −√(q/t)、**a** = t^{(2N−1)/2}（(6.15)） | Ω背景は z₁→qz₁、z₂→t^{−1}z₂\[1\] |
| Kameyama–Nawata 1703.05408 | q = e^{iε₁}、t = e^{iε₂}、Ω背景は (z₁,z₂)→(qz₁, t^{−1}z₂) | a := e^{−τ}√(q/t)、安定極限は a = t^N | 超多項式への変換は a = −**a**²**t**、q^{1/2} = −**qt**、t^{1/2} = **q** |
| Cheng–Sułkowski 2104.00713 | q = e^{−ε₂}、t = e^{ε₁}、(z₁,z₂)→(qz₁, t^{−1}z₂) | Q（閉）、z（開） | ブレーンの種類（q/t/q̄/t̄）で式が変わる |
| Avetisyan–Mkrtchyan 2205.12832 | **q = e^{g_s/2}、t = e^{g_s y/2}**（AS の q, t の平方根にあたる） | τ = a g_s（D 型では a = N−1） | 本文の (3.12)〜(3.16) はこの規約\[4\] |
| Rains–Warnaar 1506.02755、LRW 2007.03174 | Macdonald の標準的な (q,t) | T（T = t^n で n 変数の Koornwinder に戻る） | Koornwinder の境界パラメータは t₀…t₃\[5\]\[6\] |

## Key Findings（問いごと）

### Q1. orientifold 上の refined 開弦 BPS 整数性
**判定：直接の文献は無い。§1 の障害は解かない。**

- **Aganagic–Schaeffer, "Orientifolds and the Refined Topological String", arXiv:1202.4456, JHEP 09 (2012) 084。**
  - refined orientifold の定義は M 理論の二次量子化トレース (3.9) Z = Tr(−1)^F q^{S₁−S_R} t^{S_R−S₂} e^{−βH} である。開弦は「M5 を加えて同じトレースを計算する」とだけ述べている。\[1\]
  - 閉弦の積表示 (3.13)：Z_ref(X/I;q,t) = Π_{l₁,l₂≥1}(1 − q^{m₁+l₁} t^{m₂+2l₂} Q^β)^{M_β^{(m₁,m₂)}} · Π_{l₁≥1}(1 − q^{m₁+l₁} t^{m₂} Q)^{\tilde M_β^{(m₁,m₂)}}。\[1\]
    - 向き付けられた M2 は z₂ の偶または奇モードだけが残るので t^{2l₂} になる。\[1\]
    - 固定面に凍結された M2 は量子二重対数（t の昇冪なし）になる。\[1\]
  - §4 の強い整数性（mod 2 の条件 N^g_β ≡ N^{g/2,c=1}_{β/2} など）は unrefined である。開弦版 (4.17) N̂_{(□,□),g,β} ≡ N̂^{c=1}_{□,g/2,β/2}（g 偶）も unrefined で、unknot・trefoil・T(2,5) で検証されている。\[1\]
  - unrefined の開弦の第一量子化形 (4.5) は F = ½Σ_{R₁,R₂}Σ_d (1/d) f^{cov}_{R₁,R₂}(q^d,Q^d) Tr_{R₁⊗R₂}V^d − Σ_R Σ_{d 奇/偶}(1/d) f^{unor}_R(q^d,Q^d) Tr_R V^d である。\[1\]
  - refined については原文で "we do not yet have a first-quantized definition, but we still demand integrality from the second quantized M-theory index" と述べている。\[1\]
  - §6 では、SO(2N) refined CS のトーラス結び目不変量は一般に Kauffman 超多項式と変数変換で結べない（Hopf link だけは結べる）と報告している。\[1\]\[7\]
  - 閉弦 (7.8) の第三項は −Σ_d (1/2d) t^{(2N−1)d} t^{d/2} q^{−d/2}/((q^{d/2}−q^{−d/2})(t^{d/2}+t^{−d/2})) で、著者自身が「分母の正符号を第一量子化の側から理解するのは興味深い」として未解決のまま残している。\[1\]
- **Walcher の real GV、Krefl–Walcher、M 理論による解釈（arXiv:1405.6019、"M-theory interpretation of the real topological string"）。** いずれも unrefined である。refined\[8\] かつ開弦の版は見つからなかった。
- **refined 頂点と orientifold（2012 年以降）。**
  - "Refined topological vertex with ON-planes", arXiv:2201.12264, JHEP 08 (2022) 006 は、反射像の頂点 C^R を導入し、「refined では C と C^R の両方が必要」としている。対象は\[9\]\[10\] 5-brane web の閉弦（Nekrasov 分配関数）である。
  - Sp(N) と反対称物質（arXiv:1809.00629）、O 頂点と O7⁺（arXiv:2412.19655、unrefined）、対称物質の頂点（arXiv:2510.15662）も同じく閉弦・ゲージ理論に限られる。\[11\]\[12\]\[13\]
  - arXiv:2510.15662 の付録 B は、自分たちの "topological vertex with orientifold" が real topological vertex とは別物だと明記している。\[12\]
  - いずれも probe brane（結び目）の refined 整数性は扱っていない。

### Q2. 分母 (1 − q^k t^l a^m) を許す refined BPS 展開
**判定：N に依存しない分母については先例がある。N に依存する分母（本件）については先例が無い。条件付きで使える。**

- **Kameyama–Nawata, "Refined large N duality for knots", arXiv:1703.05408（v1 の題名は "...for torus knots"）, J. Knot Theory Ramif.（DOI 10.1142/S0218216520410011）。**\[14\]\[15\]\[16\]
  - unreduced の unknot は rCS‾_λ(○) = Π_{x∈λ}(t^{l'(x)} − a q^{a'(x)})/(1 − q^{a(x)} t^{l(x)+1}) である。分母は hook 型で有限、a に依らない。\[2\]
  - q-brane の双対性 (4.1)：Σ_λ rCS‾_λ(K;a,q,t) g_λ(q,t) P_λ(x;q,t) = exp(Σ_d Σ_μ (1/d) f^q_μ(a^d,q^d,t^d)/(q^{d/2}−q^{−d/2}) s_μ(x^d))。ここで (q/t)^{|λ|}/g_λ = ⟨P_λ,P_λ⟩_{q,t} である。\[2\]
  - f^q_μ = Σ_ρ M_{μρ}(t) \hat f_ρ と書く。\hat f_ρ は整数 N̂_{ρ,g,β,J_r,S_R} の母関数で、(q^{1/2}−q^{−1/2})^g (t^{−1/2}−t^{1/2})^g (q/t)^{J_r−S_R−β/2} a^β で展開される。\[2\]
  - 分母は「単一の (q^{d/2}−q^{−d/2})」と Macdonald の計量・基底に押し込まれ、整数性は N̂ の側に置かれる。
  - SO/Sp について §6 で「[Mn10]（Mariño）の予想を refine したときに正値性が見えるかは自然な問い」と述べ、未着手のままにしている。\[2\]
- **Shi Cheng, Piotr Sułkowski, "Refined open topological strings revisited", arXiv:2104.00713（Phys. Rev. D 104, 106012 (2021)）。**\[3\]\[17\]
  - t-brane の refined OV (2.40)：Z = ΠΠΠ_{n≥0}(1 − q^{−s+n+1/2} t^{r+1/2} z^d Q^β)^{(−1)^{2s}N^{(s,r)}_{(d,β)}} = exp(Σ(−1)^{2s}N q^{−ns}t^{n(r+1/2)}/(n(q^{n/2}−q^{−n/2})) z^{nd}Q^{nβ})。\[3\]
  - 分母は q の一変数だけである。代わりに holonomy を refine する。q-brane では Tr_ν V → (−√q z)^n (q,t)_n/(t,t)_n とする。\[3\]
  - strip 幾何では Z が箙の motivic 母関数に一致する。\[3\]\[18\]
  - 示唆：**分母を BPS 展開から外して「refined holonomy」に吸収する**方式は文献上確立している。
- **Cherednik の hyperpolynomial（Q4 も参照）**には分母の構造の議論が無い。arXiv:1509.08351 は「J 多項式を使うのは DAHA 超多項式に分母を出さないため」としている。\[19\]
- **Lee–Rains–Warnaar（arXiv:2007.03174）の一般的な観察**：無重複な分岐則の q,t 版の係数は、ほぼ常に 1 − q^k t^ℓ 型の二項の積の比に分解する。これは「分母が有限個の\[5\] hook 型因子である」ことの数学側の傍証だが、BPS の整数性とは結びつけられていない。

### Q3. Koornwinder／lifted Koornwinder の評価式の分母、T の極
**判定：極の位置そのものは明示されている。「安定化しない」ことを論じた文献は無い。条件付きで使える。**

- **Rains–Warnaar, "Bounded Littlewood identities", arXiv:1506.02755, Mem. Amer. Math. Soc. 270 (2021), no. 1317。**\[20\]\[21\]
  - lifted Koornwinder \tilde K_λ(x;q,t,T;t₀,…,t₃) の t_r に関する極 (2.6.6)：t₀t₁t₂t₃ = q^{1−â(s)} t^{\hat l(s)+1} T^{−2} = q^{2−λ_i−j} t^{i+λ'_j} T^{−2}（s = (i,j) ∈ λ）。â(s) = λ_i + j − 1、\hat l(s) = λ'_j + i − 1 である。\[6\]
  - 例として \tilde K_1 の係数の分母は (1−t)(1 − t₀t₁t₂t₃T²/t²) である。\[6\]
  - T = t^n では (2.6.7) により Koornwinder K_λ に戻る（l(λ) > n なら 0）。\[5\]\[6\]
  - virtual Koornwinder \hat K_λ は Q = q^m で有界な K_{m^n−λ} に特殊化する (2.6.8)。\[6\]
  - Q = 0 では Macdonald の P_λ と Littlewood 型の無限積の積になる（[108, Cor. 7.21] の引用）。\[6\]
  - lifted Cauchy (2.6.10)：Σ_λ(−1)^{|λ|} \hat K_λ(x;t,q,T;t₀…t₃) \tilde K_{λ'}(y;q,t,T;t₀…t₃) = Π_{i,j}(1 − x_i y_j)。\[6\]
  - Theorem 1.1（有界 Littlewood）は Σ_{λ₁≤m} P_λ Π_{l'(s) 偶}(1−q^{m−a'}t^{l'})/(1−q^{m−a'−1}t^{l'+1}) Π_{l(s) 偶}(1−q^{a}t^{l+1})/(1−q^{a+1}t^{l}) = (x₁⋯x_n)^{m/2} P^{(B_n,B_n)}_{(m/2)^n}(x;q,t,t) である。\[6\]
- **Lee–Rains–Warnaar, "An Elliptic Hypergeometric Function Approach to Branching Rules", arXiv:2007.03174, SIGMA 16 (2020) 142。**\[22\]\[23\]
  - (1.8b)：P_{(m^r)}(q,t) = Σ_{λ⊂(m^r), λ 偶} d_λ(q^{−m}, t^r/T; q,t) \tilde K_{(m^r)−λ}(q,t,T; 1,−1,t^{1/2},−t^{1/2})。\[5\]
  - (a) 境界パラメータの違いについて：
    - Lemma 4.1 は t = q で o_λ = \tilde K_λ(q,q,T;1,−1,q^{1/2},−q^{1/2}) を示す。t\[5\] = q では二つの型が一致する。
    - (1,−1,t^{1/2},−t^{1/2}) については「C⁺_λ(t^{2n−2}t₀t₁t₂t₃/q) = 0 となるので、t を特殊化する前に T = t^n とすることが重要」という注意がある。**この型では極と零点が衝突する**ことが文献上明記されている。\[5\]
    - (1.5a) より P^{(C_n,B_n)}(q,t,s) = K(q,t; s^{1/2}, −s^{1/2}, q^{1/2}, −q^{1/2}) なので、AS\[5\] の D_N 型 (±1, ±q^{1/2}) は s = 1 の場合にあたる。
    - Conjecture 1.4 は (±q^{1/2}, ±t^{1/2})（= (C_n,B_n) で s = t）について、2-core が空の分割の和で書かれる分岐則を予想している。\[5\]
    - (±1,±q^{1/2}) と (±1,±t^{1/2}) の q ≠ t での違いを正面から比較した文献は見つからなかった。
  - (c) Hoshino–Shiraishi との関係：postscript（pp. 48–50）は、Hoshino–Shiraishi の A_{2n−1} 型・C_n 型の漸近自由固有関数 \tildeφ^{(A_{2n−1})}、φ^{(C_n)}（x^λ\tildeφ^{(A_{2n−1})}(t^n q^{λ₁},…) = P_λ(x^±;q,t) など）の分岐則の予想（HS の Conjectures 9.12, 9.13）と (1.12a) の関係を論じている。結論部分は今回確認していない。\[5\]
- **Hoshino–Shiraishi：**
  - "Macdonald Polynomials of Type C_n with One-Column Diagrams and Deformed Catalan Numbers", arXiv:1801.09939, SIGMA（2018）。§9 に C_n 型の予想がある。\[24\]
  - "Branching Rules for Koornwinder Polynomials with One Column Diagrams and Matrix Inversions", arXiv:2002.02148, SIGMA 16 (2020) 084。BC_n→BC_{n−1} の一列分岐則と、B_n q-Toda の分岐の予想を扱う。\[25\]\[26\]
  - いずれも N（T）依存の分母の「安定化」は論じていない。
- **(b) T の極の意味：** 「Kähler 変数の極」としての物理的な解釈を与えた文献は無い。
- **推測（筆者による代入、要検証）：** AS 型の (1,−1,q^{1/2},−q^{1/2}) では t₀t₁t₂t₃ = q となり、(2.6.6) は T² に関する条件になる。T = t^N として箱ごとに代入すると次のようになる。
  - [1,1] の箱 (1,1) から q t^{2N−3} = 1 が出て、§1 の因子 (1 − q t^{2N−3}) と一致する。
  - [2] の箱 (1,1) からは q²t^{2N−2} = 1 が出る。これは (1 − q t^{N−1}) を因子に含む。
  - [3] の (1 − q²t^{2N−3}) は再現されない。
  - したがって「D_N の分母 ＝ lifted Koornwinder の T 極」は部分的にしか成り立たない。評価式（主特殊化）固有の分母が別にあると考えられる。

### Q4. 非 ADE 型の refined CS、DAHA の B/C 型 hyperpolynomial、composite DAHA
**判定：Koornwinder（C^∨C_n）型の refined CS は無い。composite DAHA の安定化は示されているが、large N duality には使われていない。**

- **Cherednik, "Jones polynomials of torus knots via DAHA", arXiv:1111.6195, IMRN 2013(23) 5366–5425。** B/C 型の hyper-polynomial を「Kauffman 不変量を一般化し、超多項式より一つ多いパラメータをもつ」ものとして定義した。後半は C^∨C₁ の DAHA-Jones を扱い、C^∨C_n への拡張の準備になっている（arXiv:1408.4348 の記述による）。\[27\]\[28\]\[29\]
- **arXiv:1509.08351（Cherednik–Danilenko, "DAHA approach to iterated torus links"）。** 「B, C, D 型の DAHA-hyperpolynomial の存在と対称性は未解決」と明記している。\[19\]
- **arXiv:2304.02200（Cherednik, "Zeta-polynomials, superpolynomials, DAHA and plane curve singularities"）。** B, C, D でも a-安定化が成り立つことを予想として掲げ、整数係数の H^C_{(r,s)}(λ;q,t,u,a) の存在を予想している。どちらも証明は無い。\[30\]
- **Cherednik–Elliot, "Refined composite invariants of torus knots via DAHA", arXiv:1503.01441, Ann. Fac. Sci. Toulouse Math. (6) 25 (2016) no. 2–3, 433–471。**
  - 頁は numdam の記録で 433–471（DOI 10.5802/afst.1501）である。
  - Theorem 2.3 で composite DAHA 超多項式の安定化と q = 1 での評価を示し、Theorem 2.4 で composite super-duality を示す。\[31\]
  - 本文に「composite の場合の安定化はより微妙」「composite DAHA 超多項式は正値性を欠く」とある。annulus の HOMFLY-PT ホモロジーとの関係が主要な未解決問題とされている。\[31\]\[32\]
  - unreduced 化で unknot の何を掛けるかは今回確認できず、未確認。large N duality や Mariño 予想に使った例は見つからなかった。
- **Elliot–Gukov, "Exceptional knot homology", arXiv:1505.01635, J. Knot Theory Ramif. (2016)。** DAHA を BPS 状態の物理に位置づけ、A 型・D 型の超多項式を中間段階で使う。orientifold\[33\] の整数性は扱っていない。
- **Elliot の学位論文：未確認**（検索で見つからなかった）。
- **AS の脚注 8：** 非単純レース型は (2,0) 理論の外部自己同型によるひねりで得られるはずだとし、(R,R^∨) 型の対称 Macdonald 多項式が必要だと述べる。Koornwinder\[1\] 多項式による refined CS の構成は文献に無い。

### Q5. (q,t)-Littlewood 恒等式と orientifold
**判定：直接の文献は無い。**

- Kawanaka 型の恒等式を orientifold の refined 分配関数や refined Ooguri–Vafa 演算子として解釈した文献は見つからなかった。
- 数学側で最も近いもの：
  - Rains–Warnaar の Q = 0 の virtual Koornwinder の積公式（Π_i Π_r(t_r x_i;q)_∞/(x_i²;q)_∞ · Π_{i<j}(t x_i x_j;q)_∞/(x_i x_j;q)_∞ が Macdonald の P_λ に掛かる形）。\[6\]
  - LRW の予想 Littlewood 恒等式 Σ_λ κ^{(1)}_λ(q,t)P_λ(x;q,t) = Π_i(tx_i²;q²)_∞/(x_i²;q²)_∞ Π_{i<j}(tx_ix_j;q)_∞/(x_ix_j;q)_∞（2-core が空の λ の和）。Schur の場合は Seamus P. Albion, "Proof of some Littlewood identities conjectured by Lee, Rains and Warnaar", Proc. Amer. Math. Soc. Ser. B 11 (2024) 133–146（arXiv:2108.12766）に証明がある（LRW Conjecture 9.5 の Schur 版）。
  - これらの因子 Π_{i<j}(t x_i x_j;q)_∞/(x_i x_j;q)_∞ は、AS (3.13) の「向き付けられた M2 の t^{2l₂} 型の積」と形は似ている。対応を示した文献は無い。
- AS 閉弦の (1+q^d)/(1+t^d) 構造との関係：
  - AS 自身が (7.8) の t^{d/2}+t^{−d/2} 分母を未解明としている。
  - **Avetisyan–Mkrtchyan, "On refined Chern-Simons / topological string duality for classical gauge groups", arXiv:2205.12832, JHEP 11 (2022) 097** は、AS の D 型の結果と自分たちの GV 部分が一致することを確認したうえで、次のように分解している（規約は q = e^{g_s/2}、t = e^{g_s y/2}）。\[4\]\[34\]
    - クロスキャップ一つ (3.14)：Σ_{n 奇}(1/n) e^{τ₁n/2}/(q^n − q^{−n})、τ₁ = g_s(1+ay−y)。\[4\]
    - クロスキャップ二つ（Klein 瓶）(3.16)：Σ (1/2n)(q^n/t^n − t^n/q^n)/((q^{2n}−q^{−2n})(t^n+t^{−n})) e^{τ₂n}、τ₂ = (1+2ay−2y)g_s/2。\[4\]
    - 「Kähler パラメータの三重化」（τ₀, τ₁, τ₂ がクロスキャップの数ごとに異なる量子シフトをもつ）を提案している。\[4\]\[34\]
    - Klein 瓶項は unrefined では消え、refined でのみ非零になる。\[4\]\[34\]
    - B 型のクロスキャップ一つ (3.22) の分母は (q^{n/2}−q^{−n/2})(t^{n/2}+t^{−n/2}) である。\[4\]
  - これは (1+t^d) 型の分母に「クロスキャップ二つ」という幾何的な意味を与える唯一の文献である。ただし閉弦のみで、著者自身も仮説と明記している。\[4\]

### Q6. 2019〜2026 年の新しい結果
**判定：refined Kauffman の LMOV 型整数性は、肯定（定式化）も否定（不可能性）も文献に無い。**

- refined、閉弦：Avetisyan–Mkrtchyan 2205.12832（上記）、Avetisyan–Mkrtchyan "On partition functions of refined Chern-Simons theories on S³", arXiv:2107.08679, JHEP 10 (2021) 033、refined\[4\] Vogel 普遍性 arXiv:2504.13831。\[34\]\[35\]\[36\]
- refined 頂点と O 面：2201.12264、2412.19655、2510.15662（いずれも 5d の閉弦）。\[9\]\[12\]\[13\]
- unrefined の Mariño 予想の検証：
  - "BPS Spectra of Complex Knots", arXiv:2410.10468。Singh–Chbili による。Manolescu の準交代結び目族で |R| ≤ 2 まで整数性を確認し、極値 refined BPS 整数についての予想も述べている。
  - "Checks of integrality properties in topological strings", arXiv:1702.06316, JHEP 08 (2017) 139。\[37\]\[38\]
  - L. Chen–Q. Chen, "Orthogonal Quantum Group Invariants of Links", arXiv:1007.1656。\[39\]\[40\]
  - Bouchard–Florea–Mariño, hep-th/0411227。\[41\]
- Kauffman ホモロジー：
  - Nawata–Ramadevi–Zodinmawia, "Colored Kauffman Homology and Super-A-polynomials", arXiv:1310.2240, JHEP 01 (2014) 126。「SO/Sp の refined CS 不変量は Kauffman ホモロジーの Poincaré 多項式と異なる」と明記している。[r]\[42\] 彩色の Kauffman ホモロジーが [r] 彩色の HOMFLY ホモロジーを含み、[r,r] 彩色の HOMFLY ホモロジーと構造が類似することも示した。付録 B で 8 の字結び目について Mariño 予想を確認している。\[42\]\[43\]\[44\]
  - "On knots, complements, and 6j-symbols", arXiv:2012.12008。薄い結び目について、[r²] 彩色の HOMFLY-PT ホモロジーから [r] 彩色の Kauffman ホモロジーへの次数の付け替え規則を与えている。\[45\]
  - Gukov–Walcher, "Matrix factorizations and Kauffman homology", hep-th/0512298（基本表現での refined BPS 縮退度）。\[7\]\[46\]
- composite／annular の圏化：
  - sutured annular KR（Trans. AMS 370 (2018)）。\[47\]
  - HOMFLY-PT の solid torus スケイン加群の圏化（Gorsky–Hogancamp–Wedrich, "Derived traces of Soergel categories", arXiv:2002.06110, IMRN。Soergel 両加群の導来水平トレース）。
  - "Evaluations of annular Khovanov–Rozansky homology"（Gorsky–Wedrich, arXiv:1904.04481, Math. Z. 303 (2023) 25）。
  - composite 表現そのものの refined／圏化不変量を large N に使った文献は無い。
- 関連（2026 年）：arXiv:2603.03628（指標による HOMFLY–Kauffman 関係、unrefined）。\[48\]

## 本件への含意と推奨（推測は推測と明記）

1. **整数性の基準を第一量子化から第二量子化に移す（文献に基づく）。** AS は refined の第一量子化形（½ と Σ_{d 奇}）を定義していない。§1\[1\] の「½F_cov ＋ unoriented」への分解は、refined では文献の裏づけが無い。AS (3.13) の積表示（向き付けられた部分は t^{2l₂}、固定面の部分は t の昇冪なし）を開弦に拡げ、「Kauffman 側の母関数が整数冪の積に書けるか」を判定基準にするのが、文献と整合する唯一の路線である。
2. **分母を holonomy の側に吸収する（KN/CS 方式の転用、推測）。** U(N) では Macdonald 基底・g_λ・refined holonomy (q,t)_n/(t,t)_n が分母を吸収している。D_N\[2\]\[3\] の分母は N に依存するので、吸収先は N（Kähler 変数）に依存する holonomy、すなわち lifted Koornwinder 基底になるはずである。具体的には B_λ(v) に A 型 Macdonald を使うのをやめ、Rains の (2.6.10) で x を probe brane の holonomy、T を Kähler 変数とする。こうすれば、§1 で見た「composite ([1],[1]) が p₁² 方向だけにある」という制約そのものが基底の取り替えで変わる可能性がある。
3. **境界パラメータの選択を検証する。** LRW は (1,−1,t^{1/2},−t^{1/2}) で極と零点が衝突すること（T = t^n を先に特殊化すべきこと）を明記している。§1\[5\] で「[1,1] に極が無い」と見えた現象は、この特殊化の順序に起因する可能性がある（推測）。AS 型 (±1,±q^{1/2}) との比較は文献に無いので、自前で行う必要がある。Conjecture 1.4 の (±q^{1/2},±t^{1/2}) も候補になる。
4. **Kähler 変数のずれをクロスキャップごとに変える（推測）。** Avetisyan–Mkrtchyan の三重化（τ₀ = yτ + (1−y)g_s/2、τ₁ = yτ + (1−y)g_s、τ₂ = yτ + (1−2y)g_s/2）を開弦に持ち込む案がある。具体的には、Kauffman\[4\] 側の [2] と [1,1] で Kähler 変数の t/q のずれが異なると仮定し、§1 の極の位置のずれ（t^{N−1} 対 t^{2N−3}）を吸収できるかを試す。
5. **最小の検証計算。** 上の Q3 の代入結果（[1,1] と [2] は lifted Koornwinder の T 極で説明でき、[3] は説明できない）を 3〜4 箱まで系統的に確かめる。評価式の分母のうち「T 極由来」と「主特殊化由来」の部分を切り分けることが、障害の正体を確定する最短の手順である。

## Caveats
- 検索回数に上限があり、INSPIRE を直接は照会していない。2019〜2026 年の論文のうち、検索に索引されていないもの（とくに SO(2N) refined CS の large N の BPS 解析）が存在する可能性は否定できない。
- Kameyama–Nawata の式番号は v1 による（v2 では一部が一つずれる）。予想 (4.5)〜(4.7) と §6 の SO/Sp への言及は v1 での確認である。
- Cheng–Sułkowski の掲載誌（Phys. Rev. D 104, 106012）は検索結果の DOI から判断した。Cherednik–Elliot の頁番号は numdam の記録（433–471）に従った。
- 「Q3 の推測」の代入は筆者が行ったもので、どの文献にも無い。Rains の (2.6.6) は lifted 多項式そのものの極であり、主特殊化の分母とは別の量である点に注意が必要である。
- Hoshino–Shiraishi の予想と LRW (1.12a) の関係について、postscript の結論は確認できていない。

## 出典

1. [Orientifolds and the Refined Topological String](https://arxiv.org/pdf/1202.4456)
2. <https://arxiv.org/pdf/1703.05408v2>
3. [Refined open topological strings revisited](https://arxiv.org/pdf/2104.00713)
4. <https://www.arxiv.org/pdf/2205.12832>
5. [An Elliptic Hypergeometric Function Approach to Branching Rules](https://arxiv.org/pdf/2007.03174)
6. <https://arxiv.org/pdf/1506.02755>
7. [Matrix Factorizations and Kauffman Homology](https://www.researchgate.net/publication/2065294_Matrix_Factorizations_and_Kauffman_Homology)
8. [M-theory interpretation of the real topological string](https://arxiv.org/pdf/1405.6019)
9. [Refined topological vertex with ON-planes](https://arxiv.org/pdf/2201.12264)
10. [JHEP08(2022)006 Published for SISSA by Springer Received: May 19, 2022](https://d-nb.info/1270133861/34)
11. [Refined topological vertex for a 5D 𝑆⁢𝑝⁢(𝑁) gauge theory with antisymmetric matter](https://arxiv.org/html/1809.00629)
12. [Topological Vertex for Symmetric matter](https://arxiv.org/html/2510.15662)
13. [O-vertex, O7+-plane, and Topological Vertex](https://arxiv.org/pdf/2412.19655)
14. [Refined large N duality for knots](https://arxiv.org/abs/1703.05408)
15. <https://arxiv.org/pdf/1703.05408v1>
16. [Refined large N duality for knots](https://www.worldscientific.com/doi/10.1142/S0218216520410011)
17. [Refined open topological strings revisited](https://doi.org/10.1103/physrevd.104.106012)
18. [\[2104.00713\] Refined open topological strings revisited](https://arxiv.org/abs/2104.00713)
19. [DAHA approach to iterated torus links](https://arxiv.org/html/1509.08351)
20. [\[1506.02755\] Bounded Littlewood identities](https://arxiv.org/abs/1506.02755)
21. [1 Introduction](https://arxiv.org/html/2507.00433)
22. [\[2007.03174\] An Elliptic Hypergeometric Function Approach to Branching Rules](https://arxiv.org/abs/2007.03174)
23. [An Elliptic Hypergeometric Function Approach to Branching Rules](https://sigma-journal.com/2020/142/)
24. [Macdonald Polynomials of Type 𝑪\_𝒏 with One-ColumnDiagrams and Deformed Catalan Numbers](https://arxiv.org/html/1801.09939)
25. [\[2002.02148\] Branching Rules for Koornwinder Polynomials with One Column Diagrams and Matrix Inversions](https://arxiv.org/abs/2002.02148)
26. [Branching Formula for 𝑞-Toda Functions of Type B](https://arxiv.org/html/2103.04378)
27. [Jones Polynomials of Torus Knots via DAHA](https://doi.org/10.1093/imrn/rns202)
28. [\[1111.6195\] Jones polynomials of torus knots via DAHA](https://arxiv.org/abs/1111.6195)
29. [DAHA and iterated torus knots](https://arxiv.org/pdf/1408.4348)
30. [Zeta-polynomials, superpolynomials, DAHA and plane curve singularities](https://arxiv.org/html/2304.02200)
31. [Refined composite invariants of torus knots via DAHA](https://arxiv.org/html/1503.01441)
32. [arxiv.org](https://arxiv.org/pdf/1503.01441)
33. [Exceptional Knot Homology](https://ar5iv.labs.arxiv.org/html/1505.01635)
34. [On refined Chern-Simons / topological string duality for classical gauge groups](https://arxiv.org/abs/2205.12832)
35. [link.springer.com](<https://link.springer.com/article/10.1007/JHEP10(2021)033>)
36. [On Refined Vogel’s universality](https://arxiv.org/html/2504.13831)
37. [Checks of integrality properties in topological strings](https://arxiv.org/pdf/1702.06316)
38. [Checks of integrality properties in topological strings](<https://link.springer.com/article/10.1007/JHEP08(2017)139>)
39. [Orthogonal Quantum Group Invariants of Links](https://arxiv.org/html/1007.1656)
40. [Orthogonal Quantum Group Invariants of Links](https://arxiv.org/pdf/1007.1656)
41. [Topological Open String Amplitudes On Orientifolds](https://arxiv.org/pdf/hep-th/0411227)
42. [Colored Kauffman Homology and Super-A-polynomials](https://arxiv.org/pdf/1310.2240)
43. [Colored Kauffman Homology and Super-A-polynomials](https://arxiv.org/pdf/1310.2240v3)
44. [11institutetext: NIKHEF theory group,](https://ar5iv.arxiv.org/html/1310.2240)
45. [On knots, complements, and 𝟔⁢𝒋-symbols](https://arxiv.org/html/2012.12008v2)
46. [Matrix Factorizations and Kauffman Homology](https://arxiv.org/pdf/hep-th/0512298)
47. [Sutured annular Khovanov-Rozansky homology](https://www.ams.org/journals/tran/2018-370-02/S0002-9947-2017-07117-8/S0002-9947-2017-07117-8.pdf)
48. [Keywords.](https://arxiv.org/html/2603.03628)
