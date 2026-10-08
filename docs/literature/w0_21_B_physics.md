# W0-B：物理・表現論側の構成の調査（evolution 法、knots–quivers、Aganagic らの構成）

作成：2026-10-08。別研究（research_plan_21.md）の W0 のうち物理・表現論側。計算はしていない。引用・転載はしていない（要約のみ）。

確認方法の記号：**A** = arXiv API のメタデータ（題名・著者・掲載誌）、**T** = TeX ソースを取得して本文を検索・該当節を読んだ（references/<id>/src/）、**W** = Web 検索の結果のみ。
掲載誌は arXiv の journal_ref に書かれたものだけを書く（DOI・頁は未確認）。

---

## 0. 結論（先に）

1. **t = q の [2,1] 彩色 HOMFLY-PT は、2 橋結び目・ツイスト結び目・3 本ブレイドの結び目でかなり計算済み**。
   - ツイスト結び目の全族 K_k（k ∈ ℤ）：Mironov–Morozov–Morozov 1408.3076 が evolution 法で閉じた式を与えた（6 つの固有値、6 つの係数 u_α(A,q) を明示）【T】。
   - 4_1 の [2,1]：AMMM 1211.6375（ケーブル法、9 本ブレイド）と Gu–Jockers 1407.5643（WZW の 6j 記号）の両方で計算され、一致【T】。
   - Gu–Jockers 1407.5643：8 交点以下の 2 橋結び目すべての [2,1]（4_1, 5_2, 6_1–6_3, 7_2–7_7, 8_1–8_4, 8_6–8_9, 8_11–8_14）【T】。
   - 3 本ブレイドの結び目（5_2, 6_2, 6_3, 7_3, 7_5, 8_2, 8_5, 8_7, 8_9, 8_10 ほか）：AMMM 1211.6375 付録、MMMS 1508.02870、Bai ほか 1801.09363（10 交点まで、データは knotebook.org。サイトの中身は未確認）【T/A】。
   - 反平行の double braid（2 パラメータ）の [2,1]：A. Morozov 1612.00422 の微分展開（予想、既知のデータに合う）【T】。
2. **refined（t ≠ q）の evolution 法で [2,1] を扱った文献は見つからなかった**。
   - refined の固有値の規則は DMMSS 1106.4305／MMM 1306.3197 にあるが、**多重度のない場合（対称表現）に限定**と本文に明記【T】。
   - 対称表現のツイスト結び目の refined 版（MMM 1306.3197）でも、係数は「自明な結び目・三葉結び目・**8 の字結び目**」の値から決めており、4_1 の superpolynomial は定義からではなく別の予想（IMMM, FGS など）から入れている【T】。
   - 4_1 の [2,1] superpolynomial は AMMM 1211.6375 の予想だけ（Z 因子の変形。条件 (A)–(F) では決まらない曖昧さが残り、正値性と自己双対点で選んだと本文に明記）【T】。
   - ホモロジー（Poincaré 多項式）の evolution には、t = −1 で消える**余分な固有値**、**部屋（chamber）ごとに違う evolution**、thick な領域での固有値の冪（"nimble"）が現れる（1802.09383, 1812.00858, 1904.10277、いずれも基本表現）【T】。KN の rCS（index）で同じことが起きるかは文献に無い。
3. **knots–quivers は対称表現だけ**（KRSS, EKL, Stošić 2402.03066, JKLNS など全て）。KRSS 1707.04017 は任意のヤング図形を「今後の課題」と明記【T】。
   任意の分割の彩色を扱うのは Ekholm–Shende の skein 値の曲線の数え上げ（2101.00619：Ooguri–Vafa の式の証明。t の精密化なし、計算例は一般論）と、Ekholm–Longhi–Nakamura 2407.09836（Hopf 絡み目の conormal のみ）【T/A】。
4. **Aganagic らの構成は [2,1] に届かない**。AV 1204.4709・AENV 1304.5778 は対称表現（古典極限・augmentation 多様体の段）、多行の表現は「平行な m 本のコピー」に帰着させる記述のみ【T】。
   Aganagic 2004.14518 / 2105.06039 / 2207.14104 はミニュスキュール表現に限定（sl_n では列 [1^k]）、階数固定で a 依存性なし、計算例は sl_2 の三葉結び目など【T】。
5. **W1 への含意**：固有値は t = q では文献で決まっている（6 個）。refined の固有値は文献の規則を [2,1] に外挿すれば候補が出る（§2.3、本調査の推測）。
   しかし係数 B_j（6 個）を決める「定義どおりに分かる点」はツイスト族では 2 点（自明な結び目と三葉結び目）しかなく、文献の t = q の計算でも 6〜7 個の結び目の値を**別の方法で計算して**入れている。refined ではそれに当たる値が無い。

---

## 1. 文献ごとの表

(i) 双曲結び目での定義と計算例、(ii) 彩色（[2,1] を含むか）と a 依存性、(iii) t の精密化と KN の rCS（index）との対応、(iv) 正値性の主張。

### 1.1 evolution 法・Racah 行列（t = q が中心）

| arXiv | 著者・題名（A で確認） | 確認 | (i) 双曲結び目 | (ii) 彩色・a | (iii) t・rCS との対応 | (iv) 正値性 |
|---|---|---|---|---|---|---|
| 1106.4305 | Dunin-Barkowski–Mironov–Morozov–Sleptsov–Smirnov, Superpolynomials for toric knots from evolution induced by cut-and-join operators（JHEP 03 (2013) 021） | A | トーラス結び目のみ | 任意（トーラス）、a あり | refined の固有値の規則（q^κ → q^{ν_{Q'}} t^{−ν_Q}）の出典 | — |
| 1306.3197 | Mironov–Morozov–Morozov, Evolution method and "differential hierarchy" of colored knot polynomials | A,T | ツイスト結び目・double braid・反平行 2 本の絡み目 | **対称表現のみ**（多重度のない場合に限ると §4 で明記）、a あり | refined の固有値：平行は (t/q)^i 倍、反平行は割る（対称表現）。ツイスト結び目の refined [r] の式。係数は自明な結び目・3_1・4_1 の値から決める | superpolynomial の正値性は前提として扱う |
| 1211.6375 | Anokhina–Mironov–Morozov–Morozov, Knot polynomials in the first non-symmetric representation | A,T | 4_1（HOMFLY は厳密、ケーブル法）、付録に 3 本ブレイドの結び目 10 個 | **[2,1]**、a あり | 4_1 の [2,1] superpolynomial を**予想**（Z 因子の変形）。条件 (A) t=q, (B) 正値性, (C) 鏡像, (D) A=t², q⁻² で [1] に簡約, (E) A=t/q, (F) 一次ずれ。下線の項の選び方は (F) 以外に効かず、α は正値性と自己双対点で選ぶ | 入力（superpolynomial）の段の正値性を条件に使う |
| 1307.2216 | Anokhina–An. Morozov, Cabling procedure for the colored HOMFLY polynomials（TMF 178 (2014)） | A,T(要旨) | |Q|·(ブレイドの本数) ≤ 12 の結び目 | 任意（ケーブル法）、a あり | t = q のみ | — |
| 1407.5643 | Gu–Jockers, A note on colored HOMFLY polynomials for hyperbolic knots from WZW models（CMP 338 (2015) 393） | A,T | 8 交点以下の 2 橋結び目すべて | **[2,1]**（多重度を含む 6j 記号を計算）、a あり | t = q のみ。3 本ブレイドの結び目で AMMM と一致と本文 | — |
| 1408.3076 | Mironov–Morozov–Morozov, On colored HOMFLY polynomials for twist knots（MPLA 29 (2014) 1450183） | A,T | **ツイスト結び目の全族**（閉じた式） | **[2,1]**、a あり | t = q のみ。固有値 6 個、係数 6 個を k = 0（自明）, 1（3_1）, −1（4_1）, 2（5_2）, −2（6_1）, 3（7_2）で決め、k = −3（8_1）で検証 | — |
| 1412.8432 | Mironov–Morozov–Sleptsov, Colored HOMFLY polynomials for the pretzel knots and links（JHEP 07 (2015) 069） | A,T(要旨) | 一般化 pretzel | 対称表現のみ。非対称への拡張は未解決と要旨 | t = q | — |
| 1504.00364 | Nawata–Ramadevi–Singh, Colored HOMFLY polynomials that distinguish mutant knots | A,T(要旨) | Kinoshita–Terasaka と Conway 結び目 | **[2,1]**、a あり | t = q | — |
| 1508.02870 | Mironov–Morozov–Morozov–Sleptsov, Colored knot polynomials. HOMFLY in representation [2,1]（IJMPA 30 (2015) 1550169） | A,T(要旨) | 3 本ブレイド＋double fat の 7 パラメータ族（10 交点以下の 8 割超と要旨） | **[2,1]**、a あり | t = q | — |
| 1601.04199 | Mironov–Morozov–Morozov–Ramadevi–Singh–Sleptsov, Tabulating knot polynomials for arborescent knots（J. Phys. A 50 (2017) 085201） | A,T(要旨) | arborescent 結び目の表（knotebook） | [21] の符号の処方に言及、a あり | t = q | — |
| 1612.00422 | A. Morozov, Factorization of differential expansion for non-rectangular representations | A,T | 反平行 double braid（ツイスト結び目を含む） | **[2,1], [3,1]**、a あり | t = q。[2,1] の微分展開（Z 因子 5 つ、うち「異常な」2 成分が 1 つの寄与）とツイスト族の F 関数を明示。Racah 行列の多重度部分は分離できないと明記 | F 関数は正係数に見えると付記（"amusing side remark"） |
| 1801.09363 | Bai–Jiang–Liang–Mironov–Morozov–Morozov–Sleptsov, Quantum Racah matrices up to level 3 and multicolored link invariants（J. Geom. Phys. 132 (2018) 155） | A,T(要旨) | 3 本ブレイドの結び目・絡み目、10 交点まで | |R| ≤ 3（[2,1] を含む）、a あり | t = q | — |
| 1903.00259 | A. Morozov, Extension of KNTZ trick to non-rectangular representations（PLB 793 (2019) 464） | A,T(要旨) | ツイスト・double braid | [r,1] | t = q | — |
| 2007.12532 | Bishler ほか, Distinguishing Mutant Knots（J. Geom. Phys. 159 (2021) 103928） | A | 変異結び目 | 非対称表現（詳細未確認） | t = q（未確認） | — |
| 2204.05977 | A. Morozov–Tselousov, Evolution properties of the knot's defect | A,T(要旨) | 一般 | 一般 | t = q（微分展開の defect の evolution） | — |
| 2609.09437 | Kononov–A. Morozov, Determinantal representation of colored HOMFLY for double braids | A | ツイスト・double braid | 長方形のみ。非長方形への拡張は未解決と要旨 | t = q | — |

### 1.2 refined／ホモロジーの evolution

| arXiv | 著者・題名 | 確認 | (i) | (ii) | (iii) | (iv) |
|---|---|---|---|---|---|---|
| 1201.3339 | Mironov–Morozov–Shakirov–Sleptsov（JHEP 2012 (2012) 70） | A | トーラス | 任意 | Macdonald と Hall–Littlewood 展開（トーラス） | — |
| 1503.01441 | Cherednik–Elliot, Refined composite invariants of torus knots via DAHA | A,T | トーラスのみ | **合成表現 [λ, μ]**（任意の組、a あり） | DAHA による refined の合成表現の不変量。合成表現の braiding 固有値（t = q）を明示 | 正値性は長方形で予想、と本文 |
| 1802.09383 | Anokhina–A. Morozov, Are Khovanov–Rozansky polynomials consistent with evolution in the space of knots?（JHEP 1804 (2018) 066） | A,T(要旨) | トーラス | 基本表現 | KR 多項式（Poincaré）でも evolution と簡約は保たれるが、n → −n の鏡像とは整合しないと要旨 | Poincaré 多項式 |
| 1812.00858 | Dunin-Barkowski–Popolitov–Popolitova, Evolution for Khovanov polynomials for figure-eight-like family of knots | A,T | 4_1 を含む 2 パラメータ族（T(2,·) とツイスト結び目を含む） | 基本表現、N = 2 | Khovanov では部屋が 4 つ＋孤立点 2 つ。(1+t) に比例し t = −1 で消える**余分な固有値**がある | Poincaré 多項式 |
| 1904.10277 | Anokhina–A. Morozov–Popolitov, Nimble evolution for pretzel Khovanov polynomials（EPJC 79 (2019) 867） | A,T(要旨) | 種数 1, 2 の pretzel | 基本表現、N = 2 | thin では固有値 1, q²T。thick では λ², …, λ^g が加わる。部屋の境界で切り替わる | Poincaré 多項式 |
| 2308.13095 | Anokhina–Lanina–Morozov, Towards tangle calculus for Khovanov polynomials | A | — | 基本表現 | 未読 | — |
| 1504.02620 | Arthamonov–Shakirov（a1_21_constructions.md による） | — | 種数 2 の pretzel・ツイスト | A₁ のみ | refined CS（index）。Poincaré 多項式と一致しないと本文（既存の調査） | — |

### 1.3 knots–quivers と skein 値の曲線の数え上げ

| arXiv | 著者・題名 | 確認 | (i) | (ii) | (iii) | (iv) |
|---|---|---|---|---|---|---|
| 1707.02991 | Kucharski–Reineke–Stošić–Sułkowski, BPS states, knots and quivers（PRD 96, 121902） | A | 6 交点以下など | **対称のみ** | 動機的 DT の生成関数 → refined LMOV（superpolynomial、Poincaré 側） | 整数性は自動 |
| 1707.04017 | 同, Knots-quivers correspondence（ATMP 23 (2019) 1849） | A,T | 6 交点以下の全結び目、無限族、thick な結び目の一部 | **対称のみ**。任意のヤング図形での LMOV の符号化は「望ましい課題」と本文 | quadruply-graded の refined LMOV の整数性も自動、と本文 | 整数性。正値性の主張は本文で未確認 |
| 1811.03110 | Ekholm–Kucharski–Longhi, Physics and geometry of knots-quivers correspondence（CMP 379 (2020) 361） | A,T | — | 「完全対称表現に限る」と本文 | 正則円板・skein | — |
| 1910.06193 | Ekholm–Kucharski–Longhi, Multi-cover skeins, quivers, and 3d N=2 dualities（JHEP 02 (2020) 018） | A,T | — | 対称のみと本文 | — | — |
| 2108.12645 | Ekholm–Kucharski–Longhi, Knot homologies and generalized quiver partition functions | A | — | 対称のみ（要旨） | HOMFLY-PT ホモロジーの Poincaré 多項式との関係を予想 | — |
| 2105.11806 | Jankowski–Kucharski–Larraguível–Noshchenko–Sułkowski, Permutohedra for knots and quivers（PRD 104, 086017） | A,T | — | 対称のみと本文 | — | — |
| 2212.04379 | 同著者, Quiver diagonalization and open BPS states | A,T(要旨) | — | 対称な quiver | LMOV 整数性の証明の一般化 | 整数性 |
| 2402.03066 | Stošić, Generalized knots-quivers correspondence | A,T | 9・10 交点の thick な結び目を含む多数 | **S^r のみ**（本文検索で非対称表現の記述なし） | colored differential が見える quiver 型の式 | — |
| 2505.05668 | Kucharski–Noshchenko, Knot-quiver correspondence: a brief review | A,T | — | 対称のみ | — | — |
| 2505.02059 | Sachdeva, A survey of knots and quivers | A | — | — | — | — |
| 2101.00619 | Ekholm–Shende, Colored HOMFLYPT counts holomorphic curves | A,T | 一般の結び目（計算例ではなく定理） | **任意の分割**（[2,1] を含む）、a あり | t の精密化なし。Ooguri–Vafa の式（彩色 HOMFLY = skein 値の曲線の数え上げ）の証明 | — |
| 2407.09836 | Ekholm–Longhi–Nakamura, The worldsheet skein D-module … Hopf link conormal | A | Hopf 絡み目のみ | 任意の分割（skein の段） | t なし | — |
| 2601.22922 | Chauhan–Ekholm–Longhi, Open strings on knot complements | A,T(検索) | ファイバー結び目（例はトーラス結び目） | skein 値 | t なし | — |

### 1.4 Aganagic らの構成

| arXiv | 著者・題名 | 確認 | (i) | (ii) | (iii) | (iv) |
|---|---|---|---|---|---|---|
| 1204.4709 | Aganagic–Vafa, Large N duality, mirror symmetry, and a Q-deformed A-polynomial for knots | A,T | 4_1 は対称表現の古典極限（Q 変形 A 多項式） | 対称のみ | refinement は可能と定性的に述べるだけ | — |
| 1304.5778 | Aganagic–Ekholm–Ng–Vafa, Topological strings, D-model, and knot contact homology（ATMP 18 (2014) 827） | A,T | augmentation 多様体（古典） | 多行の表現は m 本の平行コピーに帰着させると記述（古典の段） | t なし | — |
| 2004.14518 | Aganagic, Knot categorification from mirror symmetry, Part I（ATMP 28 (2024) 1151） | A,T | 一般論 | **ミニュスキュール表現に限定**、階数固定（a なし） | ホモロジー（Poincaré） | 定義上の次数付きベクトル空間 |
| 2105.06039 | Aganagic, Part II: Lagrangians | A,T | sl_2 の三葉結び目・自明な結び目などの計算例 | ミニュスキュールのみ | ホモロジー | 同上 |
| 2207.14104 | Aganagic, Homological knot invariants from mirror symmetry（ICM） | A,T | — | ミニュスキュール（任意の表現も原理的には現れうるが幾何ではミニュスキュールに取る、と本文） | ホモロジー | — |
| 1105.5117 | Aganagic–Shakirov（既存の調査） | — | トーラスのみ | 任意（トーラス） | refined CS ＝ 本研究の rCS の出発点 | — |

sl_n のミニュスキュール表現は列 [1^k] だけなので、sl_3 の随伴（= [2,1]）を含め [2,1] は対象外。

---

## 2. W1（refined evolution 法）への含意

### 2.1 t = q の固有値（文献で決まっている）
- ツイスト結び目のねじれの領域は**反平行**の 2 本なので、固有値は [2,1] ⊗ \overline{[2,1]} の既約成分（合成表現 [λ, μ]）で決まる（1408.3076, 1612.00422）【T】。
  research_plan_21.md W1 の「[2,1] ⊗ [2,1]」は平行の場合で、T(2, 2k+1) の族に当たる。ツイスト族では反平行を使う必要がある。
- 成分（1408.3076 の分解）：([2,1],[2,1])、([2],[2])、([2],[1,1])、([1,1],[2])、([1,1],[1,1])、([1],[1])×**2**（多重度 2）、∅。
- t = q の固有値（k を 1 つ増やすごと）：A⁶, q⁴A⁴, A⁴（([2],[1,1]) と ([1,1],[2]) が同じ値）, q⁻⁴A⁴, A², 1 の **6 個**（1612.00422 の Λ²）。
  多重度 2 の ([1],[1]) では R 行列は多重度空間に恒等的に作用するので、固有値の数は増えない（1408.3076 の ansatz は 6 項）。
- 平行（T(2,·) の族）では [2,1]⊗[2,1] の [3,2,1] が対称部分と反対称部分に 1 つずつ入るので、符号の違う固有値として分かれる（表現論の一般論。文献での確認はしていない）。

### 2.2 refined の固有値（文献には [2,1] の規則が無い）
- DMMSS／MMM の規則（q^{ν_{Q'}} t^{−ν_Q}、平行は (t/q)^i 倍、反平行は割る）は**多重度のない場合（対称表現）に限って**述べられている【T：1306.3197 §4】。[2,1] の合成表現への適用例は無い。
- Cherednik–Elliot 1503.01441 は DAHA で合成表現の refined 不変量を定義しているが、対象はトーラス結び目だけ【T】。ツイスト族の refined の固有値の出典にはならない（規則を作る材料にはなりうる）。
- ホモロジーの段（Poincaré 多項式）では、evolution が部屋ごとに変わり、余分な固有値（∝ (1+t)）や固有値の冪が現れる（1812.00858, 1904.10277）【T】。
  本研究の [2] の \hat f では p > 0 と p < 0 の違いは単項式 ε_p だけだった（REPORT_phase4.md §5.1）。index（rCS）が Poincaré 多項式のような部屋の構造をもたないかどうかは、文献では分からない【未解決】。

### 2.3 本調査での推測（仮説。文献には無い。W1 で検証すること）
REPORT_phase4.md §5.1 の [1], [2], [1,1] の入力の根（太字）1, X = 𝐚²𝐭², Y = 𝐚⁴𝐪⁴𝐭⁸, Ȳ = 𝐚⁴𝐪⁻⁴𝐭⁴ は、合成表現 [λ, μ]（|λ| = |μ| = k）について

  Λ²([λ, μ]) = (𝐚𝐭)^{2k} · (𝐪𝐭)^{2(n'(λ)+n'(μ))} · 𝐪^{−2(n(λ)+n(μ))}　（n(λ) = Σ(i−1)λ_i、n'(λ) = n(λᵀ)）

と書ける（DMMSS の q^{ν_{Q'}} t^{−ν_Q} を λ, μ ごとに当てはめた形。3 つの根に合わせて作ったので、ここまでは当てはめにすぎない）。これを [2,1] に外挿すると、候補は

  1, 𝐚²𝐭², 𝐚⁴𝐪⁴𝐭⁸（[2],[2]）, 𝐚⁴𝐪⁻⁴𝐭⁴（[1,1],[1,1]）, 𝐚⁴𝐭⁶（[2],[1,1] と [1,1],[2]）, 𝐚⁶𝐭¹⁰（[2,1],[2,1]）

の 6 個。t = q（𝐭 = −1、A² = 𝐚²、q = 𝐪）で 1, A², q⁴A⁴, q⁻⁴A⁴, A⁴, A⁶ に戻り、DGR の鏡像 (e_𝐚, e_𝐪, e_𝐭) ↦ (e_𝐚, −e_𝐪, e_𝐭 − e_𝐪) で集合として不変（Y ↔ Ȳ、他は不変）。
この形は λ ↔ μ で対称なので「異常な」2 成分は refined でも分かれないが、分かれる規則（gl_N の分割から作る規則など）を排除する根拠は無い【未確認】。
W1 の手順 1（T(2,2k+1) で根を確かめる）は平行の場合の確認なので、反平行の根はトーラス結び目では確かめられない。代わりに Cherednik–Elliot の合成表現の DAHA（トーラス結び目）で T 固有値の規則を確かめる手はある。

### 2.4 係数を決める情報
- ツイスト族 K_k で**定義どおりに分かる**のは k = 0（自明）と 3_1 の 1 点の計 2 点。t = q の文献でも残りは 4_1・5_2（ケーブル法）、6_1・7_2・8_1（Gu–Jockers）の値を**別の方法で計算して**入れている（1408.3076）【T】。
- refined の係数 B_j は 6 個（異常な成分が分かれれば 7 個）。2 点と、t = q での一致（全ての k で 1408.3076 の式に戻る）、鏡像、A = t², q⁻² での [1] への簡約（AMMM の条件 (D)）、Alexander（A = 1）の性質などは B_j に一次の条件を課すが、t についての自由度を消すかどうかは数えていない。
  第 3・4 期の結果（各結び目ごとに自由な因子 u が残る）から見ると、evolution は「結び目ごとの u_p」を「有限個の関数 c_j で u_p = Σ_j c_j Λ_j^p、ただし u_0 = u_{3_1} = 0」に縮めるだけで、c_j が Laurent 多項式であれば全ての p で Laurent 性は保たれる。したがって**Laurent 性と上の条件だけでは一意に決まらない見込みが高い**【推測】。全ての p での正値性（無限個の条件）でどこまで縮むかは未検討。
- 反平行の double braid（2 パラメータ、1606.06015・1612.00422）でも、定義どおりに分かるのは自明な結び目と 3_1 だけ（double twist 結び目は種数 1 以下なので T(2,5) 以上は含まない）。
  平行と反平行の領域をもつ族（1812.00858 の 4_1 型の族）は T(2,2k+1) を 1 本の直線として含むが、反平行の方向の情報は増えない。
- 1612.00422 の [2,1] の微分展開（因子化）を refined に持ち上げられれば強い構造条件になるが、その refined 版は 4_1 の AMMM の予想以外に無い。

---

## 3. 未確認の点
- knotebook.org のデータの中身（[2,1] の表の範囲、t の有無）。
- 2007.12532、2308.13095 の本文。
- refined knots–quivers の正値性を扱う論文（PRD 104, 106012 "Refined open topological strings revisited"）の arXiv 番号と著者（categorification_literature_map.md でも未確認）。
- §2.3 の規則は推測で、文献の裏付けは無い（対称表現の規則と 3 つの根に合うことだけを確かめた。数値計算はしていない）。
- 1408.3076 の k の符号と本研究の K_p の向きの対応（3_1 が k = +1 か −1 か）は conventions.md と照合していない。
