# 2018 年の試行錯誤（Dropbox「Mathematica/Refined Kauffman」）の整理とミスの調査

作成 2026-10-14。対象：Mathematica ノートブック 42 本（2018-07-10〜10-02、計約 1 MB のテキスト）と同梱の PDF 3 本（Cherednik 1111.6195、Cherednik–Elliot 1503.01441、「tensor product of SO(N).pdf」）。
ノートブックは scripts/nb2txt.py でテキスト化して読んだ。パッケージ `Macdonald\``、`rCS\``、`myfunc\`` の中身はノートブックに無いので、関数の挙動は出力から推定した（推定したものは【推定】と記す）。
読解は 4 つの補助エージェントに分担させ（batch2〜9）、疑わしい点のうち計算で確かめられるものは本セッションで再計算した（【計算】）。

## 1. 試行の全体像（時系列）

| 時期 | ノートブック | 内容 |
|---|---|---|
| 7/10〜7/11 | h_hat_fig8, g_hat fig8, refined kauffman check, unrefined large N Kauff | **unrefined** の Mariño の構造の確認：4_1 の \hat h_[2], \hat g_[2], \hat g_[1,1]（NRZ 1310.2240 付録 C の再現）、trefoil の refined Kauffman の候補（Cherednik の DAHA、AS (6.15)、GW (6.14)）の比較 |
| 7/20〜7/27 | refined Kaufman, unrefined Kaufman, refined Kauffman 2 | trefoil [1] の refined 候補と refined HOMFLY の差（c=1 部分）の genus 展開；SO(N) 指標と SU への分岐（so2su）；Macdonald を使った最初の LHS/RHS |
| 8/2〜8/9 | LR for Macdonald 1〜3, hat h 計算 (t-brane) 1〜2, transposition of composite, composite vs adjoint | Macdonald の LR 係数・skew Macdonald の実装と検査；composite 側（t-brane）の \hat h の計算 |
| 8/14〜8/15 | composite dimeniton, composite q-brane, hat h 計算 (q-brane), refined composite trefoil ([1],[1]), composite t-brane, check | q-brane 版と t-brane 版の composite unknot と \hat h；trefoil の composite ([1],[1]) |
| 8/18〜9/4 | composite R t-brane, framed hat g, difference R and comp, test, hat h (t-brane)R / Aganagic, composite R t-brane Aganagic, refined Kaufman t-brane, SO(even) to SU, hat g 計算 (t-brane), check LR, hat h (t-brane)3, …new skew Macdonald | t-brane での記号解（h_λ, g_λ、3 箱まで）、framing 付き unknot の \hat g, \hat h、SO(odd)/SO(even) 型の Kauffman unknot、skew Macdonald の正規化の変更 |
| 9/14〜10/2 | composite 4box, composite q,t-brane (framed), composite for adjoint, refined Kaufman t-brane 2 / after sum / (adjoint), composite (adjoint) q,t-brane | 4 箱の記号解；q-brane と t-brane の混合；**有限 N の Macdonald 次元から refined composite unknot の真の値を決定**（composite for adjoint）；それを入れた最終計算（分母 (t − qa) が残り整数性が破れて停止） |

## 2. 試行が置いていた仮説（最終段階の形）

記号：u_n = (a^{n/2} − a^{−n/2})/(t^{n/2} − t^{−n/2})（KN と同じ主特殊化。a は HOMFLY の a、t^{1/2} は Mariño の t）、P̃_λ(p) = P_λ(−p; t, q)（KN の t-brane の双対基底）、T_λ^{(p)} = [(−1)^{|λ|} q^{n(λ^T)} t^{−n(λ)} (q/t)^{|λ|/2}]^p（framing）、s_μ は Schur。

```
log[1 + Σ_λ K_λ P̃_{λ^T}] − ½ log[1 + Σ_{(λ,μ)≠(∅,∅)} H_{(λ,μ)} P̃_{λ^T} P̃_{μ^T}] = Σ_{d 奇数} Σ_μ (1/d) g_μ(a^d,q^d,t^d) s_μ(p_{d·})
log[1 + Σ_{(λ,μ)} H_{(λ,μ)} P̃_{λ^T} P̃_{μ^T}] = Σ_{d≥1} Σ_μ (1/d) h_μ(a^d,q^d,t^d) s_μ(p_{d·})
\hat g = M(q^{−1})^{−1} g、\hat h = M(q^{−1})^{−1} h
```
要請：\hat g, \hat h（に (t^{1/2} − t^{−1/2}) 等を掛けたもの）が正整数係数の Laurent 多項式。
入力（unknot）：
- Kauffman 側："odd" 型 K_λ = Σ_{ν 自己共役} (−1)^{(|ν| − rank ν)/2} P_{λ/ν}(q,t)|_{p_n → u_n}（SO(2n+1) の分岐則 o_λ = Σ (±) s_{λ/ν} の s → P 置換。unrefined では正しい）、"even" 型 K_λ = Σ_{ν ∈ (α+1|α)} (−1)^{|ν|/2} P_{λ/ν}（Littlewood の O(N) 公式の置換。+1 が無く Kauffman の unknot にならない）。
- composite 側：Koike 型 H_{(λ,μ)} = Σ_ν (−1)^{|ν|} P_{λ/ν}(q,t) P_{μ/ν^T}(q,t)|_{p_n → u_n}（冊により λ, μ の (q,t)↔(t,q) の割り当てが違う）、または有限 N の Macdonald 次元（composite のヤング図、a = t^N）。
- 結び目は trefoil（AS (6.15)、GW (6.14)、Cherednik）と 4_1（NRZ の対称表現の式）を unrefined で使っただけで、refined の結び目の \hat g は計算されていない。

要するに「Mariño の F^{SO} − ½F^{comp} = Σ_{d 奇数} … の構造を unrefined のまま置き、(i) 母関数の基底を KN の t-brane の双対基底に、(ii) 展開側は KN と同じ Schur に、(iii) Koike の composite 公式と Littlewood の SO 分岐則の s_{λ/ν} を P 型 skew Macdonald に置き換え、(iv) unknot は t 分母の主特殊化、framing は AS の T_λ に MV 符号を付けたもの」が仮説である。

## 3. 到達点

- うまくいった：記号的な h_λ, g_λ（|λ| ≤ 3；composite R t-brane, refined Kaufman t-brane）；framing 0 の unknot の composite 版（2 箱 −1、3 箱 0。ただし §4 (m) のとおり恒等式）；"odd" 型 Kauffman unknot の unrefined の整合；framed unknot の 2 箱で \hat g が正係数の Laurent 多項式になる例（framed hat g 計算：p = 2, 3, 5；refined Kaufman t-brane 2："odd" 2 箱 p = 4、"even" 3 箱 p = 3）。
- うまくいかなかった：refined で g_[1,1,1] ≠ 0（(q − t)/(1 − q²t) 型）；"even" 型の半整数・分母 (1 + q)；framed composite の分母 (q − 1)(t − 1) と負係数；真の composite（有限 N）を入れると 2 箱で既に分母 (aq − t) と ½（refined Kaufman t-brane (adjoint)、composite (adjoint) q,t-brane）；4 箱はどこでも完走していない。
- 破れた場所は正値性以前に **整数性（Laurent 多項式性）** で、破れる項はすべて (q − t) に比例するか t = q で unrefined の値に戻る（t = q への還元は全冊で整合）。

## 4. ミスと疑わしい点の一覧

計算で確認したものに【計算】、読解から判断したものに【読解】を付ける。

### 4.1 構造上のミス（結論に影響する）

(A) **refined Kauffman unknot の仮説が AS の SO(2N) refined CS と一致しない**【計算、scripts/rk_so_unknot_check.py】。"odd" 型 K_λ は [1] では AS の unknot（D_N Macdonald の主特殊化、a = t^{2N−1}）と一致するが、[2] で差 (q − t)(q − t^{N−1})/((qt − 1)(q t^{N−1} − 1))、[1,1] で差 −(q − t)(…)/(q t^{2N−3} − 1) が残る。試行は refined Kauffman の結び目の値に AS (6.15) を使おうとしていたので、unknot の正規化が理論と合っていない。文献から取った式ではなく、unrefined の分岐則の s → P 置換（仮説）であることがノートブックに明記されていない。

(B) **AS の refined unknot も composite unknot も a の Laurent 多項式でない**【計算、scripts/rk_refined_unknot_poles.py】。D_N の主特殊化は [2] で (1 − q t^{N−1})、[1,1] で (1 − q t^{2N−3}) を分母にもち、U(N) の composite ([1],[1]) の A 型 Macdonald 次元は (1 − q t^{N−1})（a = t^N で (t − qa)）を分母にもつ。試行の最終段階（composite for adjoint → refined Kaufman t-brane (adjoint)）は後者を発見して入力に使ったが、Kauffman 側は (A) の多項式の仮説のままだったので、両者の極が打ち消すはずもなく、\hat g に (aq − t) の分母が残って停止した。**極をもつ入力同士をどう組み合わせれば多項式になるか（組み合わさるべきか）が、この研究の中心の問いになる。**

(C) **Σ_{d 奇数} と ½ を unrefined のまま置いた**【読解＋AS §7 の導出】。AS の閉弦の結果では、oriented 部分の分母が (t^d − t^{−d}) で「(q, t²) の refined conifold」の形になり、½F_cov ＋ unoriented の形に書き直すと偶数 d の引き算に (1 + q^d)/(1 + t^d) の因子が付く。試行にはこの構造が無い。unknot の 1 箱までは影響しないが、2 箱以上の g_μ の定義式（½ の係数、d = 2 の項）に効く。

(D) **Koike 型の refined composite は恒等式で unrefined に帰着する**【補助エージェントの導出。本セッションで 2 箱まで記号計算で確認：log Z = 2u_1 p_1 + (u_2 p_2 − p_1²) + …】。λ 側 (q,t)、μ 側 (t,q) の割り当て（composite q,t-brane）では Macdonald の Cauchy 恒等式により F^{comp} = Σ_d (1/d)[2u_d p_d − p_d²] となり、「h_[2] = h_[1,1] = −1、3 箱以上 0」は検証ではなく恒等式。同じ (q,t) の割り当て（hat h 計算 (t-brane)3）でも unknot・framing 0 では −1, 0 が出るが、これも Π(1 − v_i v_j) の undeformed 因子と整合しているだけで refined の情報を含まない。

(E) **skew Macdonald の正規化と ν^T 側の (q,t)↔(t,q)**【読解】。パッケージの SkewMacdonald は P 型（b 因子込み）で、LR 係数から作る Q 型とは (q − t) に比例してずれる（check LR for macdonald）。Koike の公式 Σ_ν (−1)^{|ν|} s_{λ/ν} s_{μ/ν^T} の s → P 置換では、双対 Cauchy 恒等式 Σ_ν P_ν(x;q,t) P_{ν^T}(−y;t,q) = Π(1 − x_i y_j) に合わせて ν^T 側を (t,q) にするのが自然で、冊によって割り当てが違う（§4.2 (b)）。試行で出た (q − t) 比例の残留（g_[1,1,1]、difference R and comp、SO(even) to SU の rkauf3）はこの種の差。

(F) **"even" 型の Kauffman unknot は構造的に不整合**【読解＋導出】。+1 を含まないので g_[1] = 0 となり、unknot で F_K − ½F^{comp} = −Σ_{d 偶数} p_d/d となって奇数 d の仮説と矛盾する（2 箱の ±½ はこれ）。"odd" 型が Mariño の unknot（g_[1] = 1、\hat h_[2] = \hat h_[1,1] = −z^{−1}）に対応する。

### 4.2 実装上のミス・不整合（出力の信頼性に影響する）

(a) **`/.` の優先順位による brane 型の不整合**（composite q,t-brane framed、composite (adjoint) q,t-brane）：混合項 H_{(λ,μ)} の λ 因子が 2 回 swap されて (q,t) に戻り、単独項は (t,q) のまま。最終ノートブックに残っている。
(b) **dimcomp（Koike 型）の定義が 3 通り**：λ (q,t)・μ (t,q)、特殊化後に λ を swap、両方 (q,t)。どれが意図かは不明。
(c) **基底の変更**：9 月中旬までは ωP_{λ^T}(x;t,q)（IMacdonald の符号 (−1)^{n+1}）、最終段階は (−1)^{|λ|}P_{λ^T}(x;t,q) = P_{λ^T}(−x;t,q)（KN の t-brane 基底）。前者は Q_λ(x';q,t)（x' は p_n → (1−q^n)/(1−t^n) p_n）に等しく KN のどちらの brane でもない。両者の \hat g は転置と符号で結ばれ、直接比較できない。
(d) **M 行列は unrefined の M(q^{−1})**：KN の f^{\bar t} = M(q^{−1}) \hat f の構造には合っている。
(e) **framing の規約の混在**：arm/leg 版（箱あたり (q/t)^{1/2} 込み）と db 版、符号 (−1)^{|λ|} を R 側・q,t 部分を dimcomp 側に分けてグローバル記号 p で評価する冊（hat h 計算 (t-brane)R）。composite は (T_λ T_μ)^p で λ, μ 同符号。評価順序に依存する出力（unframed のはずの test1 に √q が混入）。
(f) **a の同定**：SO 側（"odd" で a = t^{N−1} 相当）と U(N) 側（a = t^N）で同じ a を使う。unknot では相殺して見えないが、結び目では効く。本来は Mariño の ν（共通の Kähler 変数）に両者を揃える必要がある。
(g) **古い出力の混在（stale）**：composite R t-brane 冒頭、composite R t-brane Aganagic の 3 箱 Solve、hat h (t-brane) Aganagic 全体、hat g 計算 (t-brane) の "odd" Framed 節（$Aborted なのに出力がある）、hat h (t-brane)3 と new skew Macdonald の Framed 節（同一出力）、composite q,t-brane の h1/test1、composite for adjoint の comp[{1},{1,1}]（comp[{2},{1,1}] の式と一致）。ノートブックの出力だけで「一致した」と判断してはいけない箇所。
(h) **検査の抜け**：3 箱の R ↔ H 関係（p_1³, p_1p_2, p_3 の係数）未検査（test）；dimso の虚数単位（framed hat g）；composite R t-brane Aganagic の Solve が解なし；composite 4box の連立が 1 次元退化（h[{4},1] が自由、p_2² の方程式が欠落）；composite for adjoint の comp[{1},{1,1}], comp[{2},{1,1}] は hi = 1 の確認なし。
(i) **genus 展開の変数**：genusExpansion は単一の g = Q − 1/Q で展開（KN の (q^{1/2} − q^{−1/2})^g (t^{−1/2} − t^{1/2})^g と違う）。hat g 計算 (t-brane) の p = 1 では内部変数が漏れて失敗。
(j) **4 箱はどこでも完走していない**。

### 4.3 8 月前半（LR for Macdonald 1〜3、hat h 計算 (t-brane) 1〜2、transposition of composite、composite vs adjoint）の要点【読解】

- **Macdonald の LR 係数 f^λ_{μν}(q,t)** は P_μ P_ν = Σ f^λ_{μν} P_λ を 3 箱まで Solve で求めたもの（例：f^{(1,1)}_{(1),(1)} = (1−q)(1+t)/(1−qt)、f^{(2,1)}_{(2),(1)} = (1−q²)(1−qt²)/((1−qt)(1−q²t))）。独立な検査は {1}×{1} だけ。
- **composite の最初の実装 nakami2[R,S,U,V,W] = (−1)^{|U|} f^R_{UV} f^S_{U^T W} P_V(t^ρ) P_W(t^ρ)** は Koike の公式の s_{R/U} を Σ_V f^R_{UV} P_V と実装したもので、Macdonald では Q_{R/U} = Σ f P_V... ではなく P_{R/U} = Σ_V (b_U b_V/b_R) f^R_{UV} P_V なので、P と Q の正規化が混ざっている（2 箱では b が相殺し影響なし、3 箱で効く）。
- **同一ノートブック内で composite ([1],[1]) の重みが 2 通り**（(q/t) b_1² と b_2）、2 箱は b 付き・3 箱は b 無しのコピー、など規約が定まっていない。plethystic log は Z を Macdonald 基底で展開しつつ log の係数を Schur 基底で取っている（h_11 = f^P_11 + (1 − c) f^P_2）。
- **M 行列**：Mq2 = M^{MV}(t^{−1})/(t^{−1/2} − t^{1/2})、Mt2 = M^{MV}(q)/(q^{1/2} − q^{−1/2})（Mariño–Vafa の M を q または 1/t で評価して全体因子を落としたもの）。(q,t) 変形した M は試していない。
- **検査はすべて t = q**（Mariño の構造の再現：2 箱 √q/(q−1)、3 箱 0）。refined で「0 になる」「整数になる」例は無い（hat h 計算 (t-brane)2 の 2 箱 = −1 は、composite の LR 係数を Schur に落としたハイブリッドで、定数なので情報が無い）。
- **hat h 計算 (t-brane) の 3 箱はバグ**（R[…] で書いた h3 に H[…] 用の置換を適用し、出力に R が残る）。
- **transposition of composite**：composite の転置対称性 (λ,μ;q,t) → (λ^T,μ^T;1/t,1/q) を unknot と trefoil で検査したが、比較した組の多くは転置対になっていない（(3,1)↔(2,1,1) だけが単項式 q^{20} を返し、整合的）。dimcomp[{1},{1}] は t のみに依存するので対称性を満たさない。
- **composite vs adjoint（8/9）が破綻点を最初に特定**：Koike 型 dimcomp[{1},{1}] = P_1(t^ρ)² − 1 は t = q では SU(N) の随伴の量子次元と一致するが、refined では Macdonald の随伴 P_{(2,1^{N−2})}(t^ρ) = P_1 P_{1^{N−1}} − f^{(1^N)}_{(1),(1^{N−1})}、f = (1−q)(1−a)/((1−t)(1−qa/t))（a = t^N）と一致しない。Koike の交代和の係数 (−1)^{|U|} が refined では a 依存の (q,t) 係数に置き換わる、というのが正しい読みで、「LR 係数を Macdonald にするだけ」の持ち上げでは composite の refined 不変量は得られない。9/30 の composite for adjoint はこれを有限 N の計算で確定させた（§4.1 (B)）。
- trefoil の composite ([1],[1])（unrefined）の照合では Stevan と Mariño が一致し、Ramadevi の式とは q の奇数冪の符号だけ違う（q^{1/2} → −q^{1/2} の規約差か転記ミス、未解決）。

### 4.4 8 月中旬（composite dimeniton、composite q-brane、check、hat h 計算 (q-brane)、refined composite trefoil ([1],[1])、composite t-brane）の要点【読解＋補助エージェントの再計算】

- **q-brane 版**：Z_H = Σ_{λ,μ} H_{(λ,μ)} m_λ m_μ P_λ(x;q,t) P_μ(x;q,t)、m_λ = (q/t)^{|λ|/2}/⟨P_λ,P_λ⟩（KN の (refined-GT1) の測度を λ, μ に別々に掛ける）。log の展開は Schur 基底、d は全部。ただし「refined g」の節では R_λ = Σ f^λ_{μν} H_{(μ,ν)} を LR で合成してから m_λ を掛けており、**composite の母関数の定義が同じノートブック内で 2 通り**（t = q でだけ一致）。
- **composite 次元の候補**：(i) AS 流の Cauchy 核（dimAganagic、i ≤ 5 で打ち切り）は t = q でも一致せず不成立、(ii) Koike 型 P-skew（dimkanno = 後の dimcomp）、(iii)「積 − 跡」P_2 P_1 − P_1。(ii) と (iii) は t ≠ q で一致しない。
- **composite unknot の h（q-brane）**：h_1 = 2(a^{1/2} − a^{−1/2})/(q^{1/2} − q^{−1/2})、h_2 = h_11 = −m_1²、3 箱 0【補助エージェントが独立実装で 3 点で再現】。すなわち Z_H(unknot) = exp(Σ_d (1/d)[2u_q(d) p_d − m_1(d)² p_d²])（Koike の母関数の持ち上げ、4 箱以上未確認）。\hat h_[2] = \hat h_[1,1] = (t^{1/2} − t^{−1/2})/(q^{1/2} − q^{−1/2}) で t = q で 1 だが Laurent 多項式でない。
- **trefoil の composite ([1],[1])**：Rosso–Jones の Adams 作用を p_n → p_{kn} の**素朴なプレシズム（(q,t) を変えない）**で持ち上げた（docs/framed_knot/conventions.md §6 で誤りと確定している規約）。t = q では unrefined の CompositeHomfly と一致するが、t ≠ q では分母に (1 + t), (qt − 1)⁴ が残る。refined の値を独立な情報源（Cherednik–Elliot の composite 超多項式）と比べた箇所はない。
- **Kauffman unknot（rKauff）の実装**：rkauf_λ = Σ_{μ∈{∅,[1]}} Σ_ν f^λ_{μν} u_ν（LR 係数を使用）。s_λ(z',1) = Σ_μ s_{λ/μ}(z') の持ち上げに必要なのは P-skew（b 因子付き）で、[2] で係数 1 と (1+q)(1−t)/(1−qt) の違いが出る（9 月の check kauf では SkewMacdonald（P-skew）に変わっている）。符号 (−1)^{(|μ|−ℓ(μ))/2} は μ = [2] で破綻。[1] の dimSO との比較は `/.` の優先順位で未判定。
- **refined g（unknot）**：g_1 = m_1 = (t^{1/2} − t^{−1/2})/(q^{1/2} − q^{−1/2})、g_2 は (q − t) に比例して非零（分母 2, (1 + q), t）。補助エージェントの変種の試験では、K_λ を Littlewood の母関数 Σ_λ o_λ s_λ = Π_{i≤j}(1 − x_i x_j) Σ_ν s_ν(z',1) s_ν(x) の各因子を s → m_λ P_λ で持ち上げたもの（P-skew と同値）にし、composite の測度を λ, μ に別々に掛けると、unknot の F_K − ½F_H が 3 箱まで厳密に Σ_{d 奇数} (1/d) m_1(d) p_d になる【1 点＋t = q で確認、2 箱は恒等式 b_1²/b_2 + b_11/b_1² = 2 で導出】。これは Cauchy 型の恒等式で閉じているだけで（§4.1 (D) と同種）、refined の理論的入力（AS の D_N 次元）とは §4.1 (A) のとおり 2 箱で一致しない。
- **t-brane 版（composite t-brane）**：基底 P_{λ^T}(−x;t,q)（測度なし、KN の (refined-GT2)）。h の式を導いただけ。
- 評価順の問題：dimcomp の二重定義、未定義の IMt3 に出力がある（3 箱が 0 なのに非零の列ベクトル。古い出力）、trefoilcomp11 の参照。

### 4.5 unrefined の部分（7 月）は正しい【計算】
h_hat_fig8・g_hat fig8 の h_[2], \hat h_[2], \hat g_[2], \hat g_[1,1]（4_1）は NRZ 付録 C と一致し、本セッションの再計算（scripts/rk_unrefined_check.py）とも一致する。M_{RS}^{−1} の形（(t h_2 + t^{−1} h_11)/(t² − t^{−2})）、[1,1] を [2] から ν → −ν、t → t^{−1} で作る対称性も正しい（trefoil の Mariño (4.8) で検証）。

## 5. 判断

- 試行は「unrefined の恒等式（Koike、Littlewood、Σ_{d 奇数}、½、M(q^{−1})）の s → P 置換」という一つの仮説の下で整合的に進められており、計算の不整合（§4.2）は結論を変えるほどではない。**本質的な失敗の原因は §4.1 (A)〜(C)**：refined の unknot（Kauffman 側・composite 側とも）が a の Laurent 多項式でないという理論側の事実を、Kauffman 側の仮説（多項式）が反映していなかったこと、そして refined の orientifold では ½ と Σ_{d 奇数} の構造自体が変わる（AS §7）ことを取り込んでいなかったこと。
- 試行の「うまくいった」部分（unknot で 2 箱 −1・3 箱 0、F_K − ½F_H が Σ_{d 奇数} m_1(d) p_d/d に閉じる）は、いずれも Cauchy 型の恒等式の s → P 置換で自動的に成り立つもので、refined の理論的入力を含まない。結び目で検証された refined の量は無い（trefoil の composite は素朴なプレシズムで、独立な照合なし）。
- 同じ方向で進めるなら、出発点は「unknot（AS の D_N Macdonald 次元と U(N) の composite Macdonald 次元、どちらも閉じた式）で refined の orientifold の開弦の構造を決める」ことである（research_plan.md K2）。2018 年の仮説のうち残す価値があるのは、KN の q-brane / t-brane の基底と測度の扱い、composite の母関数の 2 つの定義（別々の測度か LR 合成後か）、Littlewood・Koike の母関数の因子ごとの持ち上げ、という「候補の集合」であり、それらを AS の次元で篩にかける。

## 6. 手書きノート（GoodNotes、2018-10-13「Refined Kauffman for unknot」3 頁、2018-10-26「Refined Kauffman」8 頁）【読解、2026-10-14】

場所：`~/Library/CloudStorage/Dropbox/GoodNotes Backup/研究/`。poppler が無いので scripts/pdf2png.py（CoreGraphics を ctypes で呼ぶ）で PNG にして読んだ。後者は前者を含む。

- p.1：Mariño 0904.1088 の構造の要約と、AS 1202.4456 の (4.15) \hat h（composite）・(4.16) \hat g（unoriented）・(7.9) 閉弦の自由エネルギーの書き写し。refined の \hat g_R の形の推測：Σ (N̂^{c=1}(q^{1/2}−q^{−1/2})^g (t^{1/2}−t^{−1/2})^g a^Q (q/t)^{?} + N̂^{c=2} …(q^{1/2}−q^{−1/2})^g (t^{1/2}−t^{-1/2})^{g+1}…)「(q or t)」。
- p.2：composite の母関数の基底を Schur から g_λ M_λ(x;q,t)（q-brane の測度）に置き換える方針。Stevan の dim^comp(R,S) = Σ_U (−1)^{|U|} N^R_{UV} N^S_{U^T W} dim V dim W、Kanno の s_{[λ,μ]}、Stevan の SO(2N+1) 指標 ch_λ = Σ_{μ ⊂ λ, μ = μ^T} (−1)^{(|μ|−ℓ(μ))/2} s_{λ/μ} を書き、「s を Macdonald に」と注記（§2 の仮説 (iii) の出所）。
- p.3：refined large N duality（t-brane）：Z_comp = Σ rCS_{(λ,μ)} M_{λ^T}(−x;t,q) M_{μ^T}(−x;t,q) = exp(Σ_d (1/d) h_λ(a^d,q^d,t^d) s_λ(x^d))、Z_kauff = Σ rCS^{SO}_λ M_{λ^T}(−x;t,q)、F_kauff − ½F_comp = Σ_{d 奇数}…（§2 の式そのもの）。「refined quantum dim」として rCS_{(λ,μ)}(○) := Σ_ν (−1)^{|ν|} M_{λ/ν}(x;q,t) M_{μ/ν^T}(x;q,t)、rCS^{SO}_λ(○) := Σ_μ (−1)^{(|μ|−ℓ(μ))/2} M_{λ/μ}(x;q,t) と**定義**している。「unknot 以外も計算したいが composite / Kauffman の Rosso–Jones の公式が必要。refine すると γ 因子が必要で難しい。modular S 行列があれば計算できるかも」「refined Kauffman は AS, GW, Cherednik が計算しているが彼らの結果はダメ。composite も計算しているがダメ」。
- p.4：Macdonald の LR 係数と R_λ = Σ N^{λ^T}_{μ^T ν^T}(t,q) H_{(μ,ν)}。
- p.5：unrefined の unknot（標準 framing）：Cauchy・双対 Cauchy・skew の母関数の公式から、composite の母関数 exp(2Σ p_d^* p_d/d + Σ (−1)^d p_d²/d) と Kauffman の母関数 Σ_{τ=τ^T} (−1)^{(|τ|−rank τ)/2} exp(Σ p_d^* p_d/d) s_τ(x) を導き、F_g − ½F_R = Σ_{d 奇数} p_d/d を「Mathematica で確認」。
- p.6：refined の unknot の公式（Awata–Kanno 付録 B の Cauchy 恒等式）。composite の母関数を exp(2Σ p_d^* p_d/d − Σ p_d²/d) の refined 版として、Σ_τ (−1)^{|τ|} P_τ(x;q,t) P_{τ^T}(x;t,q) を挿入して skew Macdonald で書き直し、dim^comp_{λμ} = Σ_τ (−1)^{|τ|} P_{λ/τ^T}(p^*;q,t) P_{μ/τ}(p^*;t,q) を得る。「i と −x の位置は動かすことができる」。
- p.7：**仮説の出所**。「もし refined の \hat h_λ(○), \hat g_λ(○) が unrefined と同じ（h_□ = 2p_1^*、h_{□□} = h_{□,□} = −1、g_□ = 1、他は 0）なら」F^{ref}_comp = 2Σ p_d^*p_d/d − Σ p_d²/d、F_g = Σ p_d/d + Σ p_d^* p_d/d − ½Σ p_d²/d となり、Z_g = exp(Σ p_d^* p_d/d) Σ_{τ=τ^T} (−1)^{…} s_τ(x) から dim^{SO}_λ := Σ_{τ=τ^T} (−1)^{(|τ|−rank)/2} P_{λ/τ}(p^*;q,t) を**定義**し、残る問いを log(Σ_{τ=τ^T} (−1)^{…} P_τ(x;q,t)) が Σ_{d 奇数} p_d/d の refined 版になるか、としている。
  → §4.1 (A)・(D) の読みを裏付ける：Kauffman 側の refined unknot は理論（AS）から取ったものではなく、「unknot の reformulated 不変量が unrefined と同じ」という要請から逆算した定義で、「2 箱 −1・3 箱 0」はその要請の言い換えである。
- p.8：composite vs adjoint：unrefined の composite 図形 [λ,μ]（a = q^N）と、GNS³ (2.13) の M_λ(p_n = (a^{n/2}−a^{−n/2})/(t^{n/2}−t^{−n/2})) = Π_s (t^{l'} − a q^{a'})/(1 − q^a t^{l+1}) を書き、「M_{[λ,μ]}(a = t^N) を composite 図形で定義すると分母が変わりそう」「さらに分母に (1 − a^• q^• t^•) のような因子が入るため変」と注記。
  → §4.1 (B) の極（composite の Macdonald 次元の N 依存の分母）に 2018 年 10 月の時点で気づいていた記録。K2 で確定させた内容（極の位置が Kauffman 側と合わない）はここから先。
