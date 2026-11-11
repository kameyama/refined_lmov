# W0（C）：ホモロジーを作らない構成の調査（双曲結び目の [2,1] を目標に）（2026-10-08）

plans/research_plan_21.md の W0「ホモロジーを作らない構成」の項の成果物。計算はしていない。
本研究の主結果（ホモロジーではない DAHA/EHA の [2,1] で、rCS 自体は負の係数をもつが、BPS 不変量 \hat N の段で正になる。代数的な iterated torus 結び目）と同じ型の構成を、双曲結び目について探した。

- 方法：arXiv の TeX ソースを `references/<id>/src/` に取得し、本文を grep と通読で確認した（「本文」）。ソースを読んでいないものは arXiv の概要ページだけで確認した（「概要」）。
  行番号は `references/<id>/src/<ファイル>` の行番号。
- 主張の分類：【明記】＝原典にそう書いてある、【導出】＝原典の記述から本調査で導いた、【推測】＝本調査の見立て、【未確認】。
- 用語：conventions.md の約束に従い、superpolynomial は Poincaré 多項式の意味に限る。Cherednik などが「superpolynomial」と呼ぶ DAHA・インスタントン側の量は「（Cherednik の）𝔥^daha」「ℋ^inst」と書く。

## 結論（要約）

1. **双曲結び目で [2,1]（= sl_N の ω_1+ω_{N−1}、随伴型）を、a 依存性と t の精密化つきで定義し計算した、ホモロジーを作らない構成は見つからなかった。**
   最も近いのは次の二つで、どちらも [2,1] には届いていない。
   - Cherednik 2607.25666：双曲結び目（K12n242, K12n725 など、Coxeter 結び目）に DAHA 型の 𝔥^daha_λ と、F_q 点の数え上げによる ℋ^inst を与えた。a と t の両方をもつ。**ただし彩色は minuscule の ω_m（列 [1^m]）だけ**【明記】。非 minuscule への拡張は「motivic 側が知られていない」と前作 2410.14703 に明記。
   - Gruen–San Martín Suárez 2508.15171：sl_3 の F_K を**全ての既約表現**について、正のブレイド結び目（10_139, 10_152 などの双曲結び目を含む）で構成した【明記：概要】。ただし sl_3 の固定階数で a が無く、t も無い。計算例はトーラス結び目だけ【本文 grep】。
2. **3d–3d 対応・F_K 系**：双曲結び目で計算されているのは sl_2（4_1：GM、正のブレイド結び目・double twist：Park）と、対称表現の a 変形・(a,t) 変形（4_1, 5_2：EGGKPS、2110.13768）だけ。t 変形は superpolynomial（ホモロジーの Poincaré 多項式）を入力にして作るので、本研究の rCS（index）とは段が違う。
   KN が引く GPV16・GPPV17 の正値性は、U(1)_β の flavor 対称性がある **Seifert 多様体**（L(p,1) とその Hopf ファイバーの結び目）でのもの【明記】で、双曲結び目には直接は及ばない。
3. **BPS 不変量の段の正値性**と整合する主張は、文献には次の形でしかない。
   - knots–quivers（KRSS 1707.04017、Kucharski 2005.13394）：対称表現の LMOV 不変量を対称 quiver の motivic DT 不変量で書き、DT 不変量は**正整数**（Efimov）【明記】。双曲結び目 4_1, 5_2 などで quiver が見つかっている。ただし入力は対称彩色の HOMFLY-PT／superpolynomial で、正なのは quiver 変数での DT 不変量。KN の \hat N との対応（変数変換の符号）は本調査では確かめていない。
   - GPV16・GPPV17：refined の homological block（S 変換の後の量）が q と −t について正【明記】。Seifert 多様体に限る。
   - Cherednik 2607.25666：入力の段（ℋ^inst）での非負性の予想と、「𝔥^daha が正 ⇔ ℋ^inst と一致」という予想【明記】。BPS 段の話は無い。
   **「入力（rCS／DAHA）は負でも BPS の段で正」という本研究の型の主張は、双曲結び目について文献に見つからなかった。**
4. 本研究にとって新しい手がかり【導出・推測】：Cherednik の λ = (2,2) の Coxeter 結び目は双曲で、無彩色の 𝔥^daha に**負の係数**がある（ℋ^inst は正だが superdual でない）【明記：本文 3055–3080】。
   本文の規則（2430–2445 行）を当てはめるとブレイドは (T₂T₁²T₂)²T₁T₂（3 本・10 交点の正のブレイド）になる【導出】。どの結び目かは本文に無く【未確認】。
   本研究の CLAUDE.md の規則では [1] の homological な量は正でなければならないので、この 𝔥^daha_{(2,2)} は Poincaré 多項式ではない。BPS 段の正値性を [1] で試す材料にはなるが、[2,1] には届かない。

## 1. 構成ごとの表

列：(i) 双曲結び目での定義と計算例、(ii) 彩色と a 依存性（t の有無も）、(iii) 本研究の rCS（index）の規約との対応、(iv) 正値性の主張（どの段か）。

| 構成・文献 | 確認方法 | (i) 双曲結び目 | (ii) 彩色・a・t | (iii) rCS（index）との対応 | (iv) 正値性・整数性 |
|---|---|---|---|---|---|
| **GPV16** Gukov–Putrov–Vafa, *Fivebranes and 3-manifold homology*, 1602.05302 | 本文 `3manifolds-7.tex`：2560–2620（L(p,1) の refined CS）、2785–2830（L(3,1) の refined block と正値性） | なし（L(p,1)、Seifert 多様体、球面 Seifert 多様体） | 閉多様体（結び目なし）。SU(2)／U(N)。t は U(1)_β の fugacity y（t = −y²） | y ≠ 1 が AS の refined CS（2594–2600 行）。ただし Ẑ は CS 分配関数を S 変換した「homological block」で、rCS そのものではない | refined block の係数が全て正整数で Poincaré 多項式の下界と解釈【明記 2808】。Ẑ_a = Ẑ_a\|_{t=−1} + (1+t)P_{a,+}、P_{a,+} ∈ ℤ₊[[q,t]]【明記 2810–2814】。**入力（分配関数）の段を S 変換した量**での正値性 |
| **GPPV17** Gukov–Pei–Putrov–Vafa, *BPS spectra and 3-manifold invariants*, 1701.06567 | 本文 `3ManBPS-10.tex`：797–809（U(1)_β）、1823–1906（正値性）、1909–（refined CS との比較）、2818–2835（L(p,1) の結び目）、2876–（同） | なし。結び目は L(p,1) の Hopf ファイバーの unknot だけ | SU(2) の彩色 λ（Hopf ファイバーの unknot）。t は U(1)_β（Seifert 多様体にしかない【明記 799】） | refined CS の S, T 行列（t = q^β）と比較し、t = q^β への簡約を論じる（1909–、2876–） | ±Ẑ_a(q,−t) ∈ q^Δ ℤ₊[t][[q]]【明記 1825–1830】。「homological-flavor locking」で説明【明記 1898】。結び目を入れても正【明記 2827–2833】。KN はこれを「同じ現象」として引く（KN `refined-integrality.tex:872`） |
| **GM** Gukov–Manolescu, *A two-variable series for knot complements*, 1904.06057 | 本文 `surgeries.tex`：111（概要）、272–303、2921–（4_1）、1223 | 4_1：数値的に係数を求めた（再帰）。4_1 の −1/2 手術で双曲閉多様体の Ẑ を初めて計算と明記（298–303） | SU(2)。彩色は x = q^n の解析接続（対称表現 [n] の族）。a・t なし | colored Jones の解析接続（resurgence／量子 A 多項式）。index ではない | 整数係数を予想。正値性の主張なし。4_1 の手術の値は符号が混在（300 行の式） |
| **Park 2019** *Higher rank Ẑ and F_K*, 1909.13002 | 本文 `sigma20-044.tex`：44（概要）、352、356–377 | なし（F_K^G は**トーラス結び目だけ**）。352 行で、それまでの双曲の計算は 4_1（GM）と Park 2020 だけと記す | 任意のルート系 G、x は重みの多変数（全ての彩色を含む族）。a なし（対称表現への特殊化で a = q^N の規則性を観察） | index ではない | 整数係数（予想） |
| **Park 2020** *Large color R-matrix for knot complements and strange identities*, 2004.02087 | 本文 `FKexamples.tex`：58（概要）、104、119、147、535–543 | 正のブレイド結び目（10_139, 10_152 を含む）で定義と証明、double twist 結び目・ファイバー付き強擬正結び目・Whitehead 絡み目で実験的に計算 | sl_2 のみ。a・t なし | なし | F_K ∈ ℤ[q][[x^{-1}]]（正のブレイド）【明記 462】。正値性の主張なし |
| **Gruen–San Martín Suárez** *A large color R-matrix for sl_3*, 2508.15171 | 概要＋本文 grep（`main.tex`：98, 550–） | 任意の正のブレイド結び目で F_K^{sl_3} を構成（双曲を含む）。計算例はトーラス結び目だけ | sl_3 の**全ての既約表現**（[2,1] = 随伴を含む族）。a・t なし | なし | F_K^{sl_3} ∈ ℤ[q^{±1}][[x,y]]【明記：概要】。正値性の主張なし |
| **EGGKPS** Ekholm–Gruen–Gukov–Kucharski–Park–Sułkowski, *Ẑ at large N*, 2005.13349 | 本文 `HOMFLYPTFK.tex`：834–（a 変形）、1073–1164（4_1）、1165–1251（5_2）、1252–1346（t 変形）、1498–1520（t の意味）、1660 | 4_1：量子 a 変形 A 多項式で再帰的に計算、5_2：非特異性条件で一意化（非ファイバー結び目では再帰だけでは決まらない【明記 1167–1170】）。t 変形の 4_1 も super-A 多項式で計算 | **対称表現だけ**（x = q^r）。a あり。t は superpolynomial（FGS）から導入 | t 変形は superpolynomial（Poincaré 多項式）の解析接続で、index ではない。「t′ = t と同一視できれば BPS 状態の Poincaré 多項式」と記述（1500–1512） | BPS 状態の次元の生成関数と解釈【明記 1500–1505】。係数の正値性の主張は見当たらない（4_1 の f_k は符号混在、1330–1345） |
| **EGGKPS 2** *Branches, quivers, and ideals for knot complements*, 2110.13768 | 本文 `main.tex`：122（概要）、196、847–946（4_1, 5_2）、1031–、2118–2237 | 4_1, 5_2 で A 多項式の各枝の F_K、4_1 の quiver 形 | 対称表現（248 行）。a あり、t 変形版あり。sl_N の対称表現の R 行列（2118–） | 同上（superpolynomial 由来の t） | quiver 形 → DT 不変量の整数性 |
| **Kucharski** *Quivers for 3-manifolds*, 2005.13394 | 本文 `Quivers.tex`：425–440、509–540、1194–1300 | quiver は T(2,2p+1) だけ。4_1 は q^r = μ の代入が使えないと明記（429–434） | 対称表現。a あり、t 変形あり | なし | 結び目補空間の LMOV 型不変量 N_{r,i,j}（t 版 N_{r,i,j,k}）は DT 不変量なので整数【明記 533, 1294】。正値性の主張は無い |
| **KRSS** Kucharski–Reineke–Stošić–Sułkowski, *Knots-quivers correspondence*, 1707.04017 | 本文 `knots-quivers-long-arXiv.tex`：611、747、366 | 6 交点以下の全結び目（4_1, 5_2, 6_1, … を含む）とツイスト結び目の族で quiver を同定 | 対称表現。a・t あり（refined LMOV を quiver から） | 入力は対称彩色の superpolynomial（homological） | quiver の motivic DT 不変量は**正整数**（Efimov）【明記 611】。refined LMOV 不変量への言及で GKS15 と KN（Kameyama:2017ryw）を引用【明記 366, 747】 |
| **GKS** Garoufalidis–Kucharski–Sułkowski, *Knots, BPS states, and algebraic curves*, 1504.06327 | 概要のみ | ツイスト結び目と 10 交点までの結び目 | 対称表現（extremal）。refined BPS 縮退度を super-A 多項式から | 入力は superpolynomial | refined BPS 縮退度の整数性を確認【明記：概要】。正値性は【未確認】 |
| **FGS** Fuji–Gukov–Sułkowski, *Volume conjecture: refined and categorified*, 1203.2182／*Super-A-polynomial*, 1205.1515 | 概要のみ（1205.1515 のソースは取得済み、本調査では未読） | ツイスト結び目（4_1, 5_2 …） | 対称表現。a・t | 入力は彩色 HOMFLY ホモロジーの Poincaré 多項式（homological） | ホモロジー側なので本項の対象外。参考 |
| **DGG** Dimofte–Gaiotto–Gukov, *3-manifolds and 3d indices*, 1112.5179 | 本文 `indexbox.tex`：890、1072–1160、236–237 | 理想三角形分割から全ての双曲結び目補空間で定義。4_1 を計算 | SL(2)（電荷 (m,e)）。Young 図形の彩色・a・t なし。高階数は将来課題と明記（236–237） | 3d N=2 理論 T_M の超共形 index。KN の M2–M5 の index とは別物 | 整数係数。符号は混在（4_1 の I(1,1) = −q − q² + 2q³ + …【明記 1148】） |
| **Cherednik 2026** *Instanton slices and their superpolynomials*, 2607.25666（v3） | 本文 `instantons.tex`：100–118（概要）、315–350、374–378、640–668、700–712、805–884、905–911、2420–2459、3055–3083、3236–3250、4483–4523、6345–6505、6585–6615、784–795 | **K12n242**（λ = (3,2) など 4 つの図形）、**K12n725**（λ = (3,1,1)）、λ = (2,2) の双曲結び目、フック族。2 行の図形はトーラス結び目に当たる 2 系列を除いて全て双曲【明記 6352–6357】 | **minuscule の ω_m（列 [1^m]）だけ**【明記 378, 911】。a・t あり（ℋ^inst は (1+aq)…(1+aq^{ϱ−1}) で a が入る） | 𝔥^daha_λ は DAHA の coinvariant {W_λ[m] X_m}（Coxeter 結び目）。torus／iterated torus の DAHA と同じ枠組み（本研究の rCS と同族）。ℋ^inst は F_q 点の数え上げ（motivic） | 予想 (iii)：ℋ^inst の係数は全ての λ で非負【明記 4519–4523】。𝔥^daha と ℋ^inst の一致 ⇔ 𝔥^daha が非負 ⇔ ℋ^inst が superdual【明記 341–350, 4505–4523】。λ = (2,2) では 𝔥^daha が負【明記 3055–3071】。正の 𝔥_λ[m] は KhR（[1^m] 彩色）と一致と予想【明記 Conj. `conj:KhR`, 2844–2855】。**入力の段**の主張 |
| **Cherednik 2024** *Superpolynomials of algebraic links*, 2410.14703 | 本文 `superpolynoms.tex`：438–444、674–684、4936–4948、8051– | 代数的絡み目と非代数的ケーブルの例（Cab(7,2)T(4,3) など） | DAHA は任意の支配的重み。motivic 側は ω_m だけ | 同上 | 長方形・代数的結び目で正、非長方形は (1−t),(1−q) で割れば正と予想【明記 676–684】。非 minuscule では「motivic 側はこの一般性では知られていない」【明記 4946–4948】 |
| **Galashin–Lam** *Monotone links in DAHA and EHA*, 2307.16794 | 本文 `EHA.tex`：1727–1733（概要）、1824、2303–2320 | monotone link = Coxeter link（双曲を含む。Cherednik 910 行が「generally hyperbolic」と記す） | 無彩色。a・q・t（EHA の対称関数） | EHA（本研究の EHA と同じ代数、Morton–Samuelson 経由） | 凸な曲線で q,t-Schur 正と予想【明記 1732, 2304–2310】。t = q^{−1} で HOMFLY を回復【明記 1824】 |
| **Galashin–Lam** *Positroids, knots, and q,t-Catalan numbers*, 2012.09745 | 本文 `qtcat.tex`：368–370、425 | Richardson link／positroid link（双曲を含むかは本文で未確認） | 無彩色。KhR の**最高 a 次数**の部分だけ | positroid 多様体の混合 Hodge 多項式・F_q 点の数 | 幾何的に正（Betti 数）。一番上の a 次数だけ |
| **STZ** Shende–Treumann–Zaslow, *Legendrian knots and constructible sheaves*, 1402.0490 | 本文 `stz-inventiones-final.tex`：316–322、425–434、578–590 | 正のブレイドの閉包（双曲を含む） | 「colored by r」の HOMFLY ホモロジー（型は【未確認】）。F_q 点の数は HOMFLY の**最低 a 次数**の係数 | 数え上げ（ruling polynomial） | 点の数は最低 a 次数の項に限られる |
| **Trinh**, 2106.07444／**Trinh–Williams**, 2601.17293 | 概要のみ | 正のブレイドの閉包 | 無彩色。2601.17293 は任意の a 次数の HOMFLY-PT の式を与えると概要に明記 | braid Steinberg 多様体・Hecke trace | 【未確認】（本文未読） |
| **Ekholm–Shende** *Colored HOMFLYPT counts holomorphic curves*, 2101.00619（と *Skeins on branes* 1901.08027） | 本文 `cylinder.tex`：145–151、200–206、305 | **任意の結び目・絡み目**（双曲を含む）で Ooguri–Vafa の公式を証明 | **全てのヤング図形**（[2,1] を含む）。a あり（resolved conifold）。**t なし**（skein 値の数え上げ） | 正則曲線の skein 値の数え上げ。refined ではない | 彩色 HOMFLY-PT = 曲線の数え上げ【明記】。LMOV 不変量の整数性・正値性の主張は本文 grep では見当たらない【未確認】 |
| **DDP** Diaconescu–Donagi–Pantev, 1704.07412／**CDDP** 1311.3624 | 概要のみ | なし（トーラス絡み目） | refined、彩色の一般化を予想 | wild character variety のコホモロジー＝refined stable pair | トーラス絡み目の範囲 |
| **AENV** Aganagic–Ekholm–Ng–Vafa, 1304.5778 | 概要のみ | 任意の結び目（knot contact homology の augmentation 多様体） | 高次表現への一般化に言及（概要）。t なし | なし | 【未確認】 |
| **Kivinen–Trinh** *Hilb-vs-Quot*, 2310.19633 | 題名のみ（ソース取得済み、未読） | 平面曲線特異点（代数的） | — | — | 対象外 |

## 2. 構成ごとのメモ

### 2.1 3d–3d 対応・Ẑ・F_K
- **双曲結び目で計算があるもの**：sl_2 の F_K（4_1 は GM が数値的に、正のブレイド結び目・double twist などは Park 2020）。a 変形・(a,t) 変形は対称表現だけで、4_1・5_2（EGGKPS、2110.13768）。
- **高い階数**：Park 2019 の F_K^G はトーラス結び目だけ。Gruen–San Martín Suárez 2508.15171 が sl_3 の全既約表現の F_K を正のブレイド結び目で構成したが、計算例はトーラス結び目だけで、a（大きい N）と t は無い。
  sl_3 の [2,1] は x, y の 2 変数の級数の中の情報として含まれるが、F_K は彩色多項式の値そのものではなく（resurgence による解析接続）、x = q^λ を代入して [2,1] の多項式を取り出せるとは限らない（EGGKPS は 4_1 で q^r = x の代入が使えないと明記、Kucharski 429–434）。
- **t の精密化**：Seifert 多様体では U(1)_β から本来の refinement がある（GPV, GPPV）。結び目補空間（双曲を含む）では、t 変形は superpolynomial（DGR／FGS の Poincaré 多項式）を入力にして定義されている（EGGKPS 1252–1346）。**index としての t は双曲結び目補空間では定義されていない**【導出：GPPV 799 行の「Seifert 多様体に U(1)_β がある」と EGGKPS の t 変形の作り方から】。
- **正値性**：KN 872 行が引く「同じ現象」は GPV16 2808 行・GPPV17 1823–1906 行の、L(p,1) の refined homological block の正値性。これは S 変換で「分配関数」を「BPS 状態の数」に組み替えた後の量で、KN の G 展開（\hat f → \hat N）と同じ「数え上げの段に移すと正」という型の主張になっている【推測】。ただし Seifert 多様体に限る。

### 2.2 Cherednik の instanton slice（2607.25666）
- **構成**：F = F_q 上で A = F[[x,y]] の階数 n の部分加群 M ⊂ A^n で conductor 𝒞 を含むものの Quot スタック（instanton slice）の点を、t^{deg M}(1+aq)…(1+aq^{ϱ(M)−1}) で数える（概要 303–311 行）。𝒞 = I_λ（ヤング図形 λ の単項式イデアル）のとき、λ に対応する Coxeter 結び目 K_λ の「超多項式」になると予想。
- **双曲結び目**：λ = (3,2), (2,2,1), (4,1), (2,1,1,1) が K12n242、(3,1,1) が K12n725、(4,3,1) が Cab(13,2)T(3,2)（2455–2459 行）。2 行の図形は T(3k±1,3) に当たる 2 系列以外は全て双曲（6352–6357 行）。
  K12n242 の ℋ^inst と 𝔥^daha の式（6362–6418 行）、K12n725 の ℋ^inst（6600–6612 行。因数分解の形では −t⁴ があるが、展開すると係数は全て非負）。
- **彩色**：ここで λ は**結び目を決める図形**で、彩色ではない。彩色は加群の階数 n ↔ ω_n（列 [1^n]）で、minuscule に限る（378, 648–668, 911 行）。理由として、minuscule では非対称 Macdonald 多項式が E_{ω_r} = X_{ω_r} と単項式になり、DAHA 側が完全に組み合わせ的になることを挙げている（648–668, 3246 行）。
  2410.14703 4946–4948 行は、非 minuscule の重みでは DAHA 側が組み合わせ的に複雑になり、motivic 側はその一般性で知られていないと明記。**[2,1]（sl_N の ω_1+ω_{N−1}）への拡張の見通しは本文に書かれていない**。
- **a 依存性**：ℋ^inst の (1+aq)…(1+aq^{ϱ−1}) の因子で a が入る。𝔥^daha は A_{N−1} の a 安定化（a = −t^N）で得る（6450–6455 行）。
- **iterated torus の DAHA との関係**：代数的結び目では ℋ^inst と 𝔥^daha（＝CD の DAHA 超多項式）の一致は「Coincidence Conjecture」で、motivic な ℏ^mot を経由する（Conj. `conj:main` (i)、4488–4500 行）。Coxeter 結び目では 𝔥^daha_λ は DAHA の coinvariant {W_λ[m] X_m} で定義する。Galashin–Lam の EHA の量（m = 1）と一致すると予想（374–378, 2855 行）。
- **正値性**：ℋ^inst の係数の非負性は全ての λ で予想【明記 4519】。𝔥^daha は λ によって負になり（(2,2)）、非負のときに限り ℋ^inst と一致すると予想。正値な 𝔥_λ[m] は [1^m] 彩色の reduced KhR 多項式と一致すると予想（Conj. `conj:KhR`）。**いずれも入力の段（不変量そのもの）の主張で、BPS 段（LMOV／G 展開）の議論は無い**。

### 2.3 その他の数え上げ型の構成
- **点の数え上げ（STZ、Galashin–Lam、Trinh、Trinh–Williams）**：正のブレイドの閉包（10_139 などの双曲結び目を含む）で、HOMFLY-PT の一番低い（または高い）a 次数、Trinh–Williams 2601.17293 では任意の a 次数（概要）を、F_q 上の多様体の点の数で与える。無彩色で t は無い（混合 Hodge 構造で t を入れる話は Galashin–Lam のトーラス／positroid の範囲）。
- **Ekholm–Shende**：任意の結び目と任意のヤング図形で、彩色 HOMFLY-PT を resolved conifold の正則曲線の skein 値の数え上げで与えた（Ooguri–Vafa の公式の証明）。[2,1] を含み双曲結び目でも成り立つ唯一の数え上げの定理だが、t（refined）は無い。EGGKPS 1515–1520 行は、refined 版は「special Lagrangian（unknot conormal の分岐被覆など）」で定義されるだろうと推測しているだけ。
- **knots–quivers**：対称表現に限るが、双曲結び目（4_1, 5_2, 6_1 など）で LMOV 不変量を quiver の DT 不変量（正整数）で書く。refined（t）版の入力は superpolynomial。
- **wild character variety（DDP, CDDP）、ORS・Hilbert スキーム（Kivinen–Trinh）**：トーラス絡み目または代数的結び目に限る。

## 3. 双曲結び目の [2,1] に届きうる構成とその条件【推測】

| 候補 | 足りないもの | 届くための条件 |
|---|---|---|
| Cherednik の 𝔥^daha_λ を非 minuscule に拡張する | DAHA の定義自体は任意の支配的重みで書ける（2410.14703 の DAHA-Jones は任意の重み、640–668 行の説明）。Coxeter 結び目の {W_λ X_b} を b = ω_1+ω_{N−1} にして a 安定化すればよい形をしている。ただし、本文が Coxeter 結び目への拡張を minuscule に限っているのは、非対称 Macdonald 多項式 E_b を使う必要が生じるため（911 行）。 | (a) 非 minuscule の b で coinvariant の a 安定化が成り立つこと、(b) 結び目の不変性（同じ K12n242 を与える 4 つの λ で一致するか、6421–6505 行の [1] の確認を [2,1] で）、(c) t = q で彩色 HOMFLY-PT に戻ること。本研究の EHA の計算基盤（Morton–Samuelson／Galashin–Lam の monotone link の EHA 元）で λ ごとの EHA 元に [2,1] を作用させる方が近道かもしれない。いずれも計算が要る（W0 の範囲外） |
| Galashin–Lam の monotone link の EHA 元に [2,1] の彩色を入れる | 原典は無彩色。EHA 元から彩色版を作る規則（衛星／ケーブルの楕円 Hall 代数の作用）は Morton–Samuelson にある | 本研究の [2,1] の EHA 計算（代数的 iterated torus）と同じ枠組みに乗る。ただし monotone link の EHA 元が結び目の位相不変量を与えるかは、原典でも Coxeter link との対応を通してしか保証されていない |
| sl_3 の F_K（Gruen–San Martín Suárez） | a（大きい N）、t、[2,1] の多項式への特殊化 | 正のブレイド結び目（10_139, 10_152）に限る。sl_3 の随伴の彩色 Jones との関係を取り出せるか（resurgence の特殊化）は未解決。a・t が無いので本研究の判定（a・t の BPS 不変量）には直接使えない |
| Ekholm–Shende の数え上げを refine する | t（refined な曲線の数え上げ）の定義 | 原典にも EGGKPS にも定義は無い（推測のみ） |

## 4. 「BPS 不変量の段で正」と整合する主張

- **双曲結び目について、本研究と同じ型（入力は負になりうるが \hat N の段で正）の主張は見つからなかった。**
- 型が近いのは次の二つ。
  1. **GPV16・GPPV17**：分配関数を S 変換して「BPS 状態の数え上げ」に移した量（homological block）が、q と −t について正（Seifert 多様体）。KN 872 行が本研究の正値性と「同じ現象」として引用している。双曲結び目補空間では t（U(1)_β）が無い。
  2. **knots–quivers**：LMOV 不変量を quiver の DT 不変量（正整数）で書く。双曲結び目 4_1, 5_2 などを含むが対称表現だけ。正なのは quiver 変数の段で、KN の \hat N の段との対応（変数変換の符号 (−t)^{t_i} など）は確かめていない。
- **入力の段の正値性**（本研究の型とは違う）：Cherednik の ℋ^inst の非負性予想、Galashin–Lam の Schur 正値性予想。Cherednik の λ = (2,2) の例は、DAHA 型の量が双曲結び目の [1] で負になりうることを示している（ℋ^inst とは一致しない）。

## 5. 未確認の点

- Cherednik の λ = (2,2) の Coxeter 結び目の同定（本文に記載なし）。本調査で導いたブレイド (T₂T₁²T₂)²T₁T₂ が 10 交点の正の 3 本ブレイドであることから 10_139 または 10_152 の可能性があるが【推測】、HOMFLY-PT の照合はしていない。本研究は 10_139・10_152 を計算済みなので、照合すれば DAHA 型の [1] と homological な [1] の比較ができる。
- STZ の「colored by r」の彩色の型（対称／反対称／その他）。
- Trinh 2106.07444、Trinh–Williams 2601.17293、GKS 1504.06327、FGS 1203.2182、DDP 1704.07412、CDDP 1311.3624、AENV 1304.5778、Kivinen–Trinh 2310.19633 は概要（または題名）のみ。
- Ekholm–Shende 2101.00619 が LMOV 不変量の整数性を主張・証明しているか（本文 grep では見当たらない）。
- KRSS の refined LMOV 不変量（quiver の DT 不変量の変数変換）と KN の \hat N_{ρ;i,j,k} の符号の規約の対応。
- Cherednik 2607.25666 の ℋ^inst と本研究の rCS の変数の対応。代数的結び目では ℋ^inst(q,t,a;n) = 𝔥^daha(tq^n, 1/q, a; ω_n)（4494–4498 行の予想）と conventions.md の DAHA ↔ KN（a_DAHA = −a_KN）を合成すれば決まるはずだが、照合していない。
- 双曲結び目で F_K の t 変形を index（refined CS）として定義した文献の有無。arXiv の検索 API がレート制限で使えず、Web 検索（「refined / t-deformed F_K hyperbolic」「higher rank F_K figure-eight sl_3」）では 2005.13349, 2110.13768, 2007.13277（Chae、トーラス結び目の ADO の t 変形）, 2508.15171 以外は見つからなかった。
- 文献の DOI・巻・頁は全て未確認（GPV16 の JHEP 07 (2017) 071 は KN の文献リスト `refined-integrality.bbl:154-158` の記載で、原典では未確認）。

## 6. 取得したファイル（references/、git 管理外）

本調査で新たに取得（TeX ソース、`src.tar` と展開した `src/`）：1904.06057, 2005.13349, 2005.13394, 1602.05302, 1701.06567, 1909.13002, 2004.02087, 2607.25666, 2410.14703, 2307.16794, 2012.09745, 2101.00619, 1112.5179, 2110.13768, 2310.19633, 1402.0490, 2508.15171。1703.05408 はソースを再取得して展開（`src/refined-integrality.tex`）。
題名と著者は arXiv の API または概要ページのメタデータで確認した。references/README.md の表への追記は、既存の docs を編集しない指示に従い行っていない。
