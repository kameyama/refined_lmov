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
  を得た（直接の二重和と n = 2, 3 で一致【計算】）。y^n z^n Z(y,z) は n = 2, 3 で Q 上既約で、線形因子の積にならない。unrefined（t = q）でも既約なので、**有限 n での因数分解は判定基準にならない**【計算】。probe 側の K_{(n−k)}(y;t,q;t̄)（Askey–Wilson 多項式、次数 n−k）の係数は t^n の有理式になる（§4.1。当初「有理式にならない」と書いたのは漸化式の monic 化の誤りによる。2026-10-14 深夜に訂正）。

### 4.1 K2'(a)：U = t^N を形式変数とした Hopf 絡み目【計算】
- **probe 側の基底の閉じた形**：K_{(n)}(y; t,q; 1,−1,q^{1/2},−q^{1/2})（基底 t の Askey–Wilson 多項式、monic）の y^{n−2j} の係数は U_n = t^n の有理式
  c_{2j}(n) = t^j (q/t;t²)_j/(t²;t²)_j · (1 − U_n)(1 + U_n t^{−2j}) Π_{i=1}^{j−1}(1 − U_n² t^{−2i}) / Π_{i=1}^{j}(1 − q U_n² t^{−2i−1})
  【予想・計算】：有理数 (q,t) の 4 点で j ≤ 5、n ≤ 30 で漸化式と一致（scripts/rk_k2p_aw.py、tests/test_rk_askey_wilson.py）。U_n → 0 で (qy²;t²)_∞/(ty²;t²)_∞ の係数（Rains–Warnaar の極限）。分母 (1 − q t^{2n−2i−1}) は Askey–Wilson の abcd Q^{2n−•} 因子。
- **virtual 正規化**：Z~(y,z) := y^N z^N Z(y,z) は y, z の多項式で Z~(0,0) = 1。y^N K_{(N−k)}(y) = y^k Σ_j c_{2j}(N−k) y^{2j} は Rains の virtual Koornwinder \hat K_{(k)}(y;t,q,T=U)。y^a の次数 a には k ≤ a しか寄与しないので、Z~ mod (y^{A+1}, z^{B+1}) は N によらない有限和になり、K_{1^k}(t^ρ)（評価式）、c_{2j}、z 側の積 exp(−Σ_d p_d z^d/d)（p_d は等比級数の和）から **Q(q,t,U) の元として**計算できる（src/rcsknots/kauffman_hopf.py）。U = t^N（N = 3..8）での直接計算と a+b ≤ 8 で全係数一致【計算】。
- **unknot**：log Z~(y) = Σ_a L_{a0} y^a、Möbius 反転で f_{10} = (1 − U)(U + t)/(U(t − 1))、f_{a0} = 0 (a ≥ 2)。ν := (U²/t)^{1/2}（Mariño の ν、Q = U²/t）で f_{10} = −1 − (ν − ν^{−1})/(t^{1/2} − t^{−1/2})：**Mariño 型（(ν − ν^{−1})/z）の unknot と「−1」（交差帽、Mariño の g_{[1]} = 1 に対応）の和で、分母は (t^{1/2} − t^{−1/2}) のみ**【導出】（(y;t)_N (y/u;t)_N、u = t^{N−1} の log）。
- **Hopf 絡み目（a, b ≤ 4）**：f_{ab} := Σ_{d|(a,b)} (μ(d)/d) L_{a/d,b/d}(q^d,t^d,U^d) について
  - f_{11} = −f_{10} · z_q (ν~ − ν~^{−1})、z_q = q^{1/2} − q^{−1/2}、**ν~ = ν (q/t)^{1/2} = U q^{1/2}/t**。t = q では ν~ = ν で、Mariño の Hopf 絡み目 \hat g_{[1],[1]} = z(ν − ν^{−1})（0904.1088、Hopf 絡み目の例）と unknot 因子 f_{10} と符号を除いて一致する。refined では **Kähler 変数が q/t だけずれる**（Avetisyan–Mkrtchyan 型の量子シフトの開弦版）。
  - f_{1b} = f_{11} · χ_{b−1}(U/t)、χ_m(x) = x^m + x^{m−2} + … + x^{−m}（b ≤ 4）。
  - 全 (a,b)（a,b ≥ 1、≤ 4）で f_{ab} = (1 − q) × [q, t, U の整数係数 Laurent 多項式] / (1 − t)（(a,b) が共に偶数のときは /(1 − t²)）。**整数性は成り立つ**。t = q では f_{ab} は整数係数の Laurent 多項式（分母なし）。
  - (a,b) ≥ 2 は 1 変数の probe では Sp(2) 指標（Chebyshev）で混ざるので、Mariño の [2], [1,1] との比較には probe 2 個（y_1, y_2、2 変数 Koornwinder）が要る。
- 計算時間：A = B = 4 で約 1 分（sympy の有理関数体）。データ：data/refined_kauffman/hopf_f.pkl。

## 5. 現時点の判断
- Mariño 型（bosonic、½、Σ_{d 奇数}、固定した小さい表現）の refined 化は、AS の D_N 理論の unknot で既に閉じない（§3）。BC 型には bosonic な Cauchy 恒等式が無いので、小さい表現で彩色した Kauffman 側の母関数が積になる理由も無い。
- 双対 Cauchy 型（§4）は数学的に整合し、第二量子化の積表示を unknot で満たす。probe 側の基底は U に依存する（virtual Koornwinder）が、virtual 正規化 y^N Z で U を形式変数にでき（§4.1）、unknot と Hopf 絡み目（probe 各 1 個、a,b ≤ 4）で **f_{ab} ∈ (1−q)/(1−t) · Z[q^±,t^±,U^±]（(偶,偶) は /(1−t²)）の整数性が成り立つ**。unknot は Mariño 型 (ν−ν^{−1})/z + 交差帽の −1、Hopf 絡み目の f_{11} は Mariño の \hat g_{[1],[1]} の ν → ν(q/t)^{1/2} による refined 化になる。
- 残る未確認：(i) 2 箱以上の表現（probe 2 個）での整数性と Mariño の [2],[1,1] との対応、(ii) 結び目（trefoil）：AS の T 行列（framing 固有値）を双対 Cauchy 型に入れた計算、(iii) 正値性（refined BPS の抽出式は未定義）。
- 残る道（(a) は §4.1 で実施済み）：(a) Askey–Wilson / Koornwinder の係数を T の関数として閉じた形にし、Hopf 絡み目の log Z を T の冪級数として整数性を見る；(b) Avetisyan–Mkrtchyan 2205.12832 の閉弦の構造（交差帽ごとの Kähler 変数のずれ、Klein 瓶の分母 (t^n + t^{−n})）を開弦の Hopf 絡み目に合わせる；(c) B_n / C_n（Sp）型や (±q^{1/2}, ±t^{1/2}) 型で同じ計算を繰り返し、分母の構造がより単純な型があるかを見る。いずれも 1 件 30 分以内の記号計算で始められる。
