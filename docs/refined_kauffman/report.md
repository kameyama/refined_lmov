# 報告：refined Kauffman 多項式の large N duality（K0〜K2'、2026-10-14）

主張の分類：【導出】【計算（範囲）】【予想】【仮説】。記号：U = t^N（SO(2N) のランク N）、Q = U²/t（AS の Kähler 変数）、A = t^M（U(M)）。コードは src/rcsknots/kauffman_refined.py（閉じた式）、src/rcsknots/koornwinder.py（Koornwinder 多項式）、scripts/rk_k2_*.py、scripts/rk_k2p_hopf.py。

## 1. unrefined（K0）
Mariño 0904.1088 の構造を自前のコードで再現した【計算】（literature.md §1、tests/test_rk_unrefined.py）。原論文の誤植：(4.8) の \hat g_[2] の z⁹ の項、NRZ 1310.2240 (C.9) の表示式。

## 2. refined の入力（K2）
- AS の SO(2N) refined unknot dim_λ（D_N 型 Macdonald の主特殊化）は Koornwinder K_λ(x;q,t;1,−1,q^{1/2},−q^{1/2}) の主特殊化 x = (1, t, …, t^{N−1}) に等しい【計算：van Diejen–Sahi の評価式で 6 つの λ、q-差分作用素で求めた多項式でも N = 2, 3 で一致】。
- dim_λ は U の Laurent 多項式でなく、評価式の C⁺_λ(t₀t₁t₂t₃t^{2N−2}/q) = Π_{(i,j)∈λ}(1 − q^{λ_i+j−1} t^{2N−i−λ'_j}) の因子のうち分子で約分されないものを分母にもつ【計算：|λ| ≤ 4 の 8 例】。例：[2] は (1 − q t^{N−1})、[1,1] は (1 − q t^{2N−3})、[3] は (1 − q² t^{N−1})。t = q では約分される。計量 g_λ（AS (B.4)）は同じ因子を 2 乗でもつ。
- U(M) の composite [λ,μ] の Macdonald 次元（対称正規化）は A の有理式で、([1],[1]) は (1 − qA/t)、([2],[1]) は (1 − q²A/t)、([1],[1,1]) は (1 − qA/t²) を分母にもつ【計算：M = 8, 9 の直接計算と 13 例で一致】。

## 3. Mariño 型の構造は refined の unknot の 2 箱で閉じない（K2）【導出＋計算】
Z_G = Σ_λ dim_λ w_λ B_λ(v)、Z_R = Σ_{λ,μ} H_{(λ,μ)}(A) w_λ w_μ B_λ(v)B_μ(v)（w は q,t の重み、B は Schur / P_λ(v;q,t) / P_{λ^T}(−v;t,q)）とすると、F_G − c F_R の 2 次の部分の U の極は打ち消せない。composite ([1],[1]) の寄与は v の p₁² 方向だけにあるので、Kauffman 側の [2] と [1,1] が同じ極を比 f^{(1,1)}_{(1),(1)} でもつ必要があるが、D_N では極の位置が違い（U = t/q と U² = t³/q）、Koornwinder (1,−1,t^{1/2},−t^{1/2}) 型では [1,1] に極が無い。w = 1/g_λ でも (Uq² − t²), (Uq − t³) が残り、composite の 3 箱の極の位置に一致して次数が合わない。
積の形（refined Littlewood 恒等式）も w（5 通り）× B（3 通り）の全部で 3 次の p₁³, p₁p₂ の係数が 0 でない【計算：有理点 2 点】。t = q では成立する。

## 4. 双対 Cauchy（t-brane）型の定式化（K2'）
- Mimachi の恒等式 Σ_{λ⊆m^n}(−1)^{|λ|} K_{m^n−λ}(x;q,t;t̄) K_{λ'}(y;t,q;t̄) = Π_{i,j} x_i^{−1}(1 − x_i y_j^{±}) を q-差分作用素で作った Koornwinder 多項式で (n,m) = (1,1),(2,1),(1,2),(2,2),(2,3),(3,2) で確認した【計算】。
- これを refined orientifold の Ooguri–Vafa 核として読むと：ホロノミー U ∈ SO(2n) の側は AS の D_n（±1, ±q^{1/2}）、probe（y、m 個）の側は Koornwinder K_{λ'}(y; t, q; ±1, ±q^{1/2})、すなわち probe 自身の (q', t') = (t, q) に対する LRW 型 (±1, ±t'^{1/2}) である。**AS 型と LRW 型の境界パラメータは対立するものではなく、双対 Cauchy 核の両側に現れる**【導出】。
- 役割を入れ替えた読み（Mimachi は (x,n,q,t) ↔ (y,m,t,q)、λ ↔ λ' で対称）では、ホロノミー側が標準の表現 λ'（λ' ⊆ m^n）、probe 側が (det y)^n K_{n^m−λ}(y;t,q;t̄) = Rains の virtual Koornwinder \hat K_λ(y; t,q, T = t^n)：**probe 側の基底が Kähler 変数 U = t^n に依存する**。n → ∞ で (det y)^n K_{n^m−λ}(y) → P_λ(y;t,q) Π_i Π_r (t_r y_i;t)_∞/(y_i²;t)_∞ Π_{i<j}(q y_i y_j;t)_∞/(y_i y_j;t)_∞（Rains–Warnaar (2.6.9)）、すなわち A 型 Macdonald × (q,t)-Littlewood 型の普遍因子。
- unknot：Z(y) = ⟨Π x_i^{−1}(1 − x_i y_j^{±})⟩ = Π_{j}Π_{k=1}^{n} t^{−(n−k)}(1 − t^{n−k} y_j^{±}) = Π_j (y_j^{±};t)_∞/(U y_j^{±};t)_∞ × 単項式：AS (3.13) 型の第二量子化の積表示を自動的に満たす（½ も Σ_{d 奇数} も不要）【導出】。ただしこれは双対 Cauchy 恒等式の言い換えで、結び目の情報を含まない。
- Hopf 絡み目（probe 各 1 個）：AS の S_{αβ}/S_{00} = K_α(t^ρ)K_β(t^ρ q^α) から
  Z(y,z) = Σ_{k=0}^{n} (−1)^k K_{1^k}(t^ρ) K_{(n−k)}(y;t,q;t̄) · z^{−n}Π_i(1 − z x_i)(1 − z/x_i)|_{x = t^ρ q^{1^k}}
  を得た（直接の二重和と n = 2, 3 で一致【計算】）。y^n z^n Z(y,z) は n = 2, 3 で Q 上既約で、線形因子の積にならない。unrefined（t = q）でも既約なので、**有限 n での因数分解は判定基準にならない**【計算】。probe 側の K_{(n−k)}(y;t,q;t̄)（Askey–Wilson 多項式、次数 n−k）の低次の係数は Askey–Wilson の漸化式の係数の部分和で、t^n の有理式にならない。

## 5. 現時点の判断
- Mariño 型（bosonic、½、Σ_{d 奇数}、固定した小さい表現）の refined 化は、AS の D_N 理論の unknot で既に閉じない（§3）。BC 型には bosonic な Cauchy 恒等式が無いので、小さい表現で彩色した Kauffman 側の母関数が積になる理由も無い。
- 双対 Cauchy 型（§4）は数学的に整合し、第二量子化の積表示を unknot で満たすが、probe 側の基底が U に依存し（virtual Koornwinder）、Laurent 多項式性による整数性の言い方が使えない。整数性の判定には Z 自体の積表示（Z = Π(1 − q^• t^• U^• y^•)^{±1}）を無限積として扱う必要があり、Hopf 絡み目でそれを判定する方法は未確立。
- 残る道：(a) Askey–Wilson / Koornwinder の係数を T の関数として閉じた形にし（Rains の二項公式）、Hopf 絡み目の log Z を T の冪級数として整数性を見る；(b) Avetisyan–Mkrtchyan 2205.12832 の閉弦の構造（交差帽ごとの Kähler 変数のずれ、Klein 瓶の分母 (t^n + t^{−n})）を開弦の Hopf 絡み目に合わせる；(c) B_n / C_n（Sp）型や (±q^{1/2}, ±t^{1/2}) 型で同じ計算を繰り返し、分母の構造がより単純な型があるかを見る。いずれも 1 件 30 分以内の記号計算で始められる。
