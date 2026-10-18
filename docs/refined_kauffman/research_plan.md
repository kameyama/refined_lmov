# 研究計画：refined Kauffman 多項式の large N duality（第 3 の研究、2026-10-14 作成）

規則は README.md。文献は literature.md、2018 年の試行の整理は past_attempt.md。

## 0. 問い
- Q1 unrefined の Kauffman の large N duality（Mariño 0904.1088）の構造は何か。自分のコードで再現できるか。【K0、済】
- Q2 2018 年の試行はどこで、なぜ破れたか。【K1、済（past_attempt.md）】
- Q3 refined では unknot の時点で何が起きるか。AS の SO(2N) refined unknot（D_N Macdonald 次元）と U(N) の composite の Macdonald 次元（どちらも閉じた式）から、refined の orientifold の開弦の構造（½ と Σ_{d 奇数} の refined 版、基底、極の打ち消し）を決められるか。【K2】
- Q4 その構造の下で、結び目（trefoil, T(2,5), 4_1）の refined Kauffman 不変量（AS の SO(2N) refined CS、または NRZ/GW の Kauffman ホモロジー）と refined composite 不変量（Cherednik–Elliot の composite DAHA 超多項式、または U(N) refined CS の composite 彩色）から抽出した \hat g_ρ, \hat h_ρ は整数係数の Laurent 多項式か。正値性は何に課せるか。【K3〜K5】

## 1. 段階

### K0：unrefined の再現【済、2026-10-14】
- Mariño の構造（g_R, h_R, M_{RS}、予想）を literature.md §1 に整理。scripts/rk_unrefined_check.py で trefoil の \hat g_[2], \hat g_[1,1]（Mariño (4.8)）と 4_1 の \hat h, \hat g（NRZ 付録 C）を全係数で再現。誤植 2 点（Mariño (4.8) の z^9、NRZ (C.9)）を特定。
- 残り：tests/ に移す（test_rk_unrefined.py）。

### K1：2018 年の試行の監査【済（草稿）、2026-10-14】
- past_attempt.md。構造上のミス (A)〜(F)、実装上の不整合 (a)〜(j)。

### K2：refined unknot から構造を決める（本計算に入る前の停止点）
入力（すべて閉じた式、記号計算）：
- Kauffman 側：AS の SO(2N) refined unknot dim_{q,t}^{D_N}(λ)（scripts/rk_refined_unknot_poles.py の eval_DN）。|λ| ≤ 3（必要なら 4）。
- composite 側：U(N) の composite 表現 [λ,μ]_N の A 型 Macdonald 次元（(q,t)-hook の式、N を記号として扱えないので N = 4..8 で計算して a = t^N の有理式に直す）。
- 変数の同定：両者を共通の Kähler 変数で表す（unrefined では Mariño の ν。refined では AS の Q = t^{2N−1} と U(N) の a = t^{N_U} の関係を unknot の [1] で決める）。
検討する仮説（優先順）：
1. 基底：Kauffman 側・composite 側とも KN の t-brane 基底 P_{λ^T}(−x;t,q)（2018 年の最終形）か、composite の片側を q-brane にするか。
2. 展開：Σ_d と Σ_{d 奇数} の refined 版。AS §7 の閉弦の構造（oriented 部分が (q, t²) の形、偶数 d の引き算に (1 + q^d)/(1 + t^d)）を開弦に持ち上げた形を第一候補にする。
3. 極の打ち消し：D_N 次元の分母 (1 − q t^{N−1}) 等と composite の分母 (1 − q t^{N_U−1}) 等が、どの同定と組み合わせで打ち消すかを 2 箱で調べる。
判定：unknot の \hat g, \hat h が（unrefined の \hat g_[1] = 1、\hat h_[2] = \hat h_[1,1] = −z^{−1}、他 0 を refined 化した形で）Laurent 多項式になる組み合わせがあるか。
停止点：K2 の結果（できた／できない、どの仮説が残るか）を報告して確認を求める。

### K2'：K2 の否定的な結論を受けた再定式化（2026-10-14 追加、deep research の報告に基づく）
- 判定基準を第二量子化に移す：AS (3.13) に倣い、Kauffman 側の母関数（unknot）が Π(1 − q^• t^• Q^•)^{整数} の形（向き付け部分は t の偶数冪、固定面の部分は t の昇冪なし）に書けるかを判定する。½ と Σ_{d 奇数} は使わない。
- 基底：A 型 Macdonald の代わりに lifted Koornwinder \tilde K_λ(x;q,t,T;t̄)（T = Kähler 変数）と Rains の lifted Cauchy 恒等式を使う。D_N の分母が C⁺ 因子（T の極）であること【計算】を利用して、分母を基底と計量に吸収できるかを調べる。
- Kähler 変数のずれ：クロスキャップの数ごとに Q を量子シフトする（Avetisyan–Mkrtchyan 2205.12832 の閉弦の三重化の開弦版）【仮説】。
- 境界パラメータ：(±1, ±q^{1/2})（AS）と (±1, ±t^{1/2})（LRW）、(±q^{1/2}, ±t^{1/2})（LRW Conjecture 1.4）を unknot の 2〜3 箱で比べる。
- 停止点：unknot の 2 箱で積表示が得られるか否かが決まったとき。
- K2'(a)（2026-10-14 深夜、実施済み）：probe 側の Askey–Wilson 係数を U = t^N の閉じた形にし、virtual 正規化 y^N z^N Z で Hopf 絡み目の log Z を Q(q,t,U) の上で計算、Möbius 反転した f_{ab} の分母と整数性を見る。
- K2'(b)（2026-10-14、実施済み）：成分 1 の probe を無限個にしても S 行列の対称性で 1 変数の \hat K_{(k)}(z) と冪和に縮約できる（多変数 Koornwinder 不要）。|R| ≤ 4、b ≤ 2 で Schur 基底の整数性を確認（report.md §4.2）。Mariño の \hat g_{[2],[1]} との係数ごとの対応は未実施。
- K3'（次、trefoil）：AS 1202.4456 (6.15), (6.17) の T(2,2m+1) の [1] を入力に、f_1 = −W_{[1]} の構造（refined の h/g 分解：½h^{ref}_[1] = A 型 refined CS の [1] に A = U²/t を入れたもの、g^{ref}_[1] = f_1 − ½h^{ref}_[1] が Laurent 多項式か＝refined Rudolph 性）を調べる。framing は refined Rosso–Jones（p_2 の Koornwinder 展開）で AS の式と合わせる。

### K3：結び目の refined 入力の整備
- refined Kauffman：AS の SO(2N) refined CS の S, T 行列（D_N Macdonald 多項式、小さい N で数値的に）から trefoil の [1], [2], [1,1]。NRZ の対称表現の Kauffman ホモロジー（trefoil, 4_1）との関係（AS は「トーラス結び目では変数変換で結べない」）を記録。
- refined composite：Cherednik–Elliot の composite DAHA 超多項式（trefoil の ([1],[1]) は本文に明示）。unreduced 化に必要な次元は K2 の結果に従う。
- refined HOMFLY：既存（src/rcsknots）。

### K4：結び目での検証
- trefoil（2 箱）で \hat g, \hat h の整数性。破れるなら、どの入力の組み合わせでも破れるかを分類する。
- 正値性は [1] から（CLAUDE.md の規則どおり、G 展開の前後）。

### K5：報告
- report.md に結果と導出、novelty を literature.md §3 に。

### K4'：refined のクロスキャップ展開（c = 1, 2 の refined 版）と正整数性（2026-10-14 追加、ユーザーの指示）
- 定義の場所：unknot と Hopf 絡み目（双対 Cauchy 型の f_{R,S}、Schur 基底、仮説なしで得られる）。
- 分解：f_{R,S} = ½ f^{cov}_{R,S} + g_{R,S}。f^{cov} は被覆（U(N) refined CS の composite 表現、A = ν² = U²/t）の連結自由エネルギーで、
  単一変数の probe では Z^{cov}(y,z) = Σ_λ C_λ(y) P_λ(t^ρ) Π_i (1 − z x_i^{(λ)})(1 − z/x_i^{(λ)})（λ は composite (1^a; \overline{1^b})、C_λ は Pieri 係数）。
  [1],[1] は解析的：½f^{cov}_{11} = (ν−ν^{−1})(ν~−ν~^{−1}) z_q/z_t（ν~ = ν(q/t)^{1/2} は q^{1/2}t^{(N−1)/2} として自然に現れる）、g_{11} = z_q(ν~ − ν~^{−1})。
- クロスキャップ展開：g = Σ N^{c=1}_{g,β,J} G^g w^{2J} ν^β + X · Σ N^{c=2}_{g,β,J} G^g w^{2J} ν^β、G = (q^{1/2}−q^{−1/2})(t^{1/2}−t^{−1/2})、w = (q/t)^{1/2}、
  奇数因子 X の候補は z_q = q^{1/2}−q^{−1/2}、z_t = t^{1/2}−t^{−1/2}、w − w^{−1}。どの X で割り切れるかを unknot・Hopf・trefoil [1] で決める。
- 判定：N^{c} の整数性と符号（正値性）、および DGR 変数（𝐪 = √t、𝐭 = −√(q/t)、𝐚 = ν）での正整数性。
- 順序（効果の高い順）：(1) unknot と Hopf [1],[1]（定義と X の決定）、(2) trefoil [1]（結び目での正値性の最初の判定。K3' の単項式の仮説に依存）、(3) Hopf の (a),(b)（a+b ≤ 4、被覆の composite の計算が必要）、(4) 結び目の 1^k（SL(2,Z) の実装が必要）。

## 2. 停止点
- K2 の後（refined unknot の構造）。
- research_plan の各段階で 1 件 30 分を超える計算が必要になったとき（CLAUDE.md の準備の規則）。

## 3. 状態
| 段階 | 状態 |
|---|---|
| K0 | 済（2026-10-14） |
| K1 | 済・草稿（2026-10-14） |
| K2 | 済・否定的（2026-10-14、progress.md） |
| K2' | 済・条件付き（2026-10-14、report.md §4–5） |
| K2'(a) | 済・肯定的（2026-10-14、report.md §4.1）：Hopf 絡み目（probe 各 1 個、a,b ≤ 4）で f_{ab} の整数性。次：probe 2 個、trefoil |
| K2'(b) | 済・肯定的（2026-10-14、report.md §4.2）：probe 無限個、Schur 基底で |R| ≤ 4、b ≤ 2 の整数性。Macdonald 基底では破れる |
| K4' | 一部済（2026-10-14/17、report.md §7）：被覆の実装、Hopf の refined Rudolph 性（a,b ≤ 4）、(1,b) の閉じた形（ν → ν~）。(a) 偶数×偶数の (1+t)：全 d（第二量子化）の定義 g = ½E/(1+t) + C で説明され、2(1+t)g は整数 Laurent、(4,4) に (1+t²) なし、unknot の −½ は消える（2026-10-17）。E と C の分離：セクター別の被覆と t = q の条件で (2,2),(2,4) は一意（積の形、E_24 = E_22 χ_1）、(4,4) は不定（2026-10-18）。Ĩ = 入れ替え（Adams）を |R| ≤ 4、b ≤ 4 で検証、ε_{(1,1),(1)} は一意。両成分 2 箱の E は Adams でも決まらない。結び目の 1^k は未 |
| K3'（trefoil [1]） | 済（2026-10-14/17/18、report.md §6）：g^{ref}_[1] の整数性は単項式 T_[1]^{−n}(q/t)^{δ/2} のもとで成立。T(2,2m+1)（m ≤ 5）で、この単項式は Laurent 性と t = q の Kauffman 型の条件だけから決まり、AS の SU(N) の因子と同じ形（2026-10-18）。(6.15) は N = 2, 3 で GN（A 型）により厳密に検証。1^k（k ≤ N）を N = 2, 3 で計算（2026-10-18）：N 依存の極は打ち消し、f_a は整数 Laurent（framing は極の上でしか決まらない）。2 箱（probe 2 個、Schur |R| ≤ 2）を N = 2, 3 で計算（2026-10-18/19）：framing は「極 + t = q で 1」で決まり f は整数 Laurent。被覆 U(3), U(5) とあわせた g は 2g 整数。ĝ = M^{-1}(q^{-1})g の Laurent 性は m_(1) = t/q を要求し、K3' の (t/q)^{1/2} と両立しない（未解決）。3 箱（|R| ≤ 3、2026-10-19）：SO の f̂ = M^{-1}(q^{-1})f は m_(1) = t/q のときだけ整数 Laurent、framing は全部 (t/q)^{|μ|/2}、g（|R| = 3）は整数、被覆は混ざった方向で hatted が閉じない。一様な規格化（data × (t/q)^{δ|μ|/2}、2026-10-19）で、全方向の hatted が閉じ、K3' の (q/t)^{1/2} は被覆側の規格化のずれと特定（m ≤ 3）。判定式も一様な規格化で述べ直した（SO(2) の点で 1、𝐭 について偶、c = 2 の符号は m によらない、2026-10-19）。一般の N は未（D_N の SL(2,Z)/DAHA が必要） |
| B（クロスキャップ BPS 数の定義） | 一部済（2026-10-18、report.md §8）：定義の提案と、T(2,2m+1) の [1] での符号の規則 |
| K3〜K5 | 未着手 |
