# framed knot の refined large N 双対：Mariño–Vafa と KN から自然に説明できるか（2026-10-09、ユーザーの指示）

問い：Mariño–Vafa（MV）による framed knot の large N 双対の説明と、KN による refined large N 双対の説明を組み合わせると、本研究の framed knot の refined LMOV（conventions.md §17）が自然に説明できるか。

印：【原典】原典で確認、【導出】、【計算】、【解釈】（物理的な読み）、【仮定】。
原典は、MV が hep-th/0108064 の TeX（framing.tex）、KN が references/1703.05408/src/refined-integrality.tex、AS が references/1105.5117/src、IKV が hep-th/0701156 の TeX（refinedvertex.tex）。

---

## 1. 二つの説明の要点

### 1.1 MV（非 refined）【原典】
1. **CS 側**：framing を p 単位変えると、W_R は exp(2πi p h_R) 倍になる。h_R は WZW の primary の共形次元（\naframing, \cweight, framing.tex:629-647）。
   U(N) では q^{pκ_R/2} λ^{pℓ/2} 倍になる。ここで κ_R = ℓ + Σ_i(l_i² − 2i l_i) = 2Σ content（\unframing:659）。
2. **弦側**：framing の整数は、Lagrangian L_K の open の平坦座標の取り方（AKV の整数の曖昧さ）にあたる。
   どの framing でも、別の IR の振る舞いをもつ幾何の domain wall の BPS 数として整数になるはずだ、と述べている（:807-816）。
3. **微妙な点**：平坦座標が û → û + iπp とずれる。そのため V → (−1)^p V となり、W_R は (−1)^{ℓp} 倍になる。MV はこの符号が整数性に不可欠だとしている（:833-845）。
   λ^{pℓ/2} は V に吸収され、f_R の全体の因子になる（:846-848）。結論は W^{(p)}_R = (−1)^{ℓp} q^{pκ_R/2} W_R（\framefin:854）。
4. 例の式 f^{(p)}_{[2]} = q^p W_{[2]} − ½(W_□² + (−1)^p W_□(q², λ²))（:858-866）では、Adams の写しで符号が 2 乗されない。符号は V の置換ではなく、W_R の係数として扱われている。

### 1.2 KN（refined、framing 0）【原典】
1. **CS 側**：refined CS は、M5（S³）と M5'（N_K）の間の M2 が作る 3d N=2 理論の refined index Tr(−1)^F q^{J_1−S_R} t^{S_R−J_2} として定義される。追加の U(1)_R が要り、これはトーラス結び目（Seifert）で存在する（:453-462）。
2. **弦側**：遷移の後、L_K に付く M2 の BPS 状態から、refined な Schwinger 計算で自由エネルギーの形が決まる（:708-716）。
   q-brane と \bar t-brane の 2 通りの式が \hat f を決める（(refined-GT1), (refined-GT2), :844-850）。\bar t-brane では色を転置し、(q,t) → (t^{−1}, q^{−1}) とする（:511）。
3. **正値性**：S_R が BPS 状態に自明に作用するという仮定から従う（:872-892）。
4. framing については、脚注で「U(1) 部分は framing と linking number の補正を与えるだけで、本論文では役割をもたない」としか述べていない（:462 の脚注）。

---

## 2. 組み合わせ：refined な framing の変化の因子

### 2.1 CS 側：refined の T 行列の固有値【原典 ＋ 導出】
- AS の refined T 行列は T_i = c_T q^{(λ_{R_i} + βρ)²/2}、t = q^β（1105.5117 の式 :1541-1543）。
- U(N) の λ（分割）について T_λ/T_∅ を計算する。SU(N) の部分 q^{−|λ|²/2N} を除き、a = t^N とおくと
  T_λ/T_∅ = (a^{1/2}(q/t)^{1/2})^{|λ|} · q^{n(λ')} t^{−n(λ)} 【導出。N = 3, 4, 5、|λ| ≤ 3 で指数を記号計算で確認】
  となる。つまり、framing を 1 単位変える因子は、refined の framing の固有値 q^{n(λ')} t^{−n(λ)}（(N2') の T_λ）に、箱 1 つあたりの単項式を掛けたものである。
  ただし、AS の本文の T_□ の表示（:2509）と向き・符号の規約が合うかは確かめていない（未確認）。
- 箱あたりの単項式 c^{|λ|} は、c = a^{1/2}, (q/t)^{1/2} のように Adams 作用と両立する（ψ_k(c) = c^k）。そのため \hat f_ρ を c^{|ρ|} 倍するだけで、β と 2J_r を一様にずらす（MV が λ^{pℓ/2} を V に吸収したのと同じ）。
  太字の変数では (q/t)^{1/2} = −𝐭 なので、全体の符号 (−1)^{|ρ|} も伴う。これは判定で「全体の符号は良い方」を選ぶことで吸収される。

### 2.2 弦側：refined topological vertex の framing の因子【原典】
- IKV（hep-th/0701156 の Appendix「Framing factors」、refinedvertex.tex:3305-3337）は、preferred direction に沿った framing の因子として次を与える。
  f_ν(t,q) = (−1)^{|ν|} t^{n(ν)} q^{−n(ν^t)}
  t = q では q^{−κ(ν)/2}。また、Nekrasov の分配関数への外部ゲージ場の結合（Tachikawa）からも、同じ t^{n(ν)} q^{−n(ν^t)} が出ると述べている。
- これは本研究の framing[λ] = (−(q/t)^{1/2})^{|λ|} q^{n(λ')} t^{−n(λ)} の逆数と、次の 2 点を除いて同じ形である。
  - 箱あたりの (q/t)^{1/2}（§2.1 により影響しない）。
  - 向き（x ↔ −x）と、(q,t) の対応の取り方。
  ただし IKV は、符号 (−1)^{|ν|} を「後の便宜のため」に入れており、導出はしていない。

### 2.3 符号 (−1)^{|λ|x}：MV の平坦座標のずれの refined 版【導出 ＋ 解釈】
- **太字の変数で見ると、符号は 𝐭 の冪になる**【導出】。framing[□] = −(q/t)^{1/2} = 𝐭（KN の (cov)）。一般の λ では
  framing[λ] = 𝐭^{|λ| + 2n(λ')} 𝐪^{2(n(λ') − n(λ))}
  となり、符号のない単項式である（[2] は 𝐪²𝐭⁴、[1,1] は 𝐪⁻²𝐭²、[2,1] は 𝐭⁵。conventions.md §17）。
- 非 refined の極限 t = q は 𝐭 = −1 にあたる。したがって **MV の符号 (−1)^{ℓp} は、refined では 𝐭^{ℓp} の t = q での影**と読める【解釈】。
  framing は、(彩色付き) HOMFLY-PT ホモロジーの次数を符号なしでずらす操作になる。ホモロジーとして自然な形である。
- 物理的には、framing の 1 単位が [1] の BPS 状態のスピン 2J_r を 1 ずらす（\hat N^{(x)}_{[1]} は 2J_r が x ずれるだけ、conventions.md §17）。
  非 refined の index では、それが (−1)^{2J} として符号に見える。MV の平坦座標のずれ û → û + iπp（V → −V）の refined 版は、V → 𝐭V = −(q/t)^{1/2}V と読める【解釈】。
- AS も、三葉結び目の前因子を「−6 単位の framing と √(q/t)」と書いている。そして √(q/t) は「L_K 上の brane の open の平坦座標 V のずれに関係するはずだ」と述べている（1105.5117 :2507-2513）【原典】。
  本研究の framing[□] = −(q/t)^{1/2} は、MV の −1 と AS の √(q/t) を合わせたもので、この予想を具体化したものにあたる【解釈】。
- Adams の写しでの符号の扱い（framing[□](q², t²) = −𝐭²。符号は 2 乗されない）も、MV の f^{(p)}_{[2]} の (−1)^p W_□(q², λ²) と同じである【原典で照合】。指数和の根 −1 は、ここから来る。

### 2.4 t = q での MV との一致【計算】
- framing[λ]^x の t = q の極限は (−1)^{ℓx} q^{xκ_R/2}（MV の κ_R）に一致する。λ は |λ| ≤ 3 の全て、x = −3, −1, 1, 2 で確認した。
  さらに、本研究の unknot（主特殊化）と MV の量子次元が、a = λ、q = q_MV で単項式倍を除いて一致する（q ↔ 1/q では一致しない）。したがって framing の向きも x = p で MV と同じである。
  テストは tests/test_framing.py の test_t_eq_q_limit_is_marino_vafa_framing と test_t_eq_q_unknot_matches_marino_vafa_variables。

---

## 3. KN の議論は framing を変えても成り立つか

| KN の議論の要素 | framing を変えたとき | 根拠 |
|---|---|---|
| refined index による CS の定義（U(1)_R が要る） | framing の変化は T 行列（モジュラー変換）の作用なので、Seifert の構造と U(1)_R を壊さないと考えられる | 【仮定】。AS は framing を T 行列の冪として、refined CS の中で扱っている（:2030 のコメント、:2509） |
| 遷移の後の Schwinger 計算の形（自由エネルギーの形） | 変わらない。framing は L_K の無限遠での境界条件（平坦座標）の選び方で、M2 の BPS 状態の数え方の構造は同じ | MV/AKV の議論そのまま【解釈】 |
| q-brane と \bar t-brane の 2 式が同じ \hat f を与えること | 保たれる。framing[λ](1/t, 1/q) = framing[λ^T](q, t) なので、\bar t-brane の規則（色の転置と (q,t) → (t^{−1}, q^{−1})）と両立する | 【導出】（tests/test_framing.py::test_framing_factor_transposition）。計算でも、全ての判定で f^{\bar t} = M(q^{−1})\hat f が成り立った（framing_check の検査）【計算】 |
| S_R が自明に作用 ⇒ 正値性 | 同じ仮定の下で、全ての x で \hat N^{(x)} ≥ 0 が予言される | 【仮定】。計算では、732 件の判定と、2 箱の 8 結び目の全ての x で成立 |
| 非トーラス結び目への拡張（ホモロジーを入力） | framing は結び目によらない単項式を掛けるだけなので、拡張の仕方は同じ | 太字の変数で符号のない単項式なので、ホモロジーの次数のずらしとして自然（§2.3） |

**結論（考察）**：framed knot の refined large N 双対は、次の 3 つの組み合わせとして自然に説明できる。
1. CS 側：AS の refined T 行列の固有値（本質は q^{n(λ')} t^{−n(λ)}）。
2. 弦側：IKV の refined vertex の framing の因子と同じ形。
3. MV の平坦座標のずれの refined 版：太字の 𝐭 で、t = q で MV の符号 (−1)^{ℓp} になる。

このとき、KN の抽出式（refined-GT1, GT2）とその物理的な根拠（Schwinger の計算、S_R の自明な作用）は、そのまま全ての framing に適用できる。
本研究の定義は、ユーザーの以前の計算の定義と同じもので、この説明と矛盾しない。

**自然には説明できない点・仮定の点**
1. 符号の起源：MV は非 refined で AKV の平坦座標から符号を導いたが、refined の平坦座標のずれ（V → −(q/t)^{1/2}V）は AS の示唆と本研究の読みにとどまり、導出していない。IKV も符号は便宜として入れている。
2. U(1)_R が framing の変化で保たれること（KN の正値性の議論の前提）は仮定である。
3. 正値性（整数性より強い主張）は MV の議論からは出ない。KN の S_R の仮定による。framing を変えても成り立つことは、計算で支持されている。
4. MV は、固定した種数で BPS 数が p の多項式になると観察している（unknot の Π(l + jp)、:1455）。本研究の「x の指数和」は、これを g_s で展開する前の形と読める【解釈、未確認】。
5. Cheng–Sułkowski（2104.00713）の toric brane の refined open BPS の framing との対応は未確認（文献調査に回した）。

---

## 4. 論文での位置づけ（outline.md §3 への反映）
- §3.1（定義）で、framing[λ] の 3 つの由来を述べる：AS の T 行列、IKV の framing の因子、MV の符号の refined 版（太字では 𝐭^{|λ|+2n(λ')}𝐪^{κ}）。
- §3.2〜§3.5 の「x によらない構造」は、上の §3 の表の【導出】の部分（転置との両立、Laurent 性、G 展開が閉じる条件）。
- 物理的な予言（全ての x で正）は、KN の S_R の仮定を framed に広げたものとして述べる。根拠は計算である。
