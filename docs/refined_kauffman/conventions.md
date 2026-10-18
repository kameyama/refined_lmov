# 規約（第 3 の研究：refined Kauffman）

作成 2026-10-14。unrefined は Mariño 0904.1088 の規約を基準にする（literature.md §1.1）。refined の規約は、計算に入る前（研究計画 K2）に確定する。

## unrefined（Mariño）
| 項目 | 規約 |
|---|---|
| 変数 | t, ν。z = t − t^{−1}。U_q(sl_N)：t = q^{1/2}, ν = t^N。U_q(so_N)：t = q^{1/2}, ν = q^{N−1} = t^{2N−2} |
| unknot | H = (ν − ν^{−1})/z、G = 1 + (ν − ν^{−1})/z |
| 彩色の unknot | H_{(R,S)} = dim_q (R,S)、G_R = dim_q^{SO(N)} R |
| framing | H_R, G_R → (−1)^{f|R|} t^{f κ_R} ×（標準 framing から f 単位） |
| composite | (R,S) = Σ_{U,V,W} (−1)^{|U|} N^R_{UV} N^S_{U^T W} (V ⊗ \bar W) |
| 母関数 | Z_R = Σ R_R s_R(v)、Z_G = Σ G_R s_R(v)、R_R = Σ N^R_{R₁R₂} H_{(R₁,R₂)} |
| reformulated | F_R = Σ_d Σ_R (1/d) h_R(t^d,ν^d) s_R(v^d)；F_G − ½F_R = Σ_{d 奇数} Σ_R (1/d) g_R(t^d,ν^d) s_R(v^d) |
| 基底変換 | \hat h = M^{−1} h、\hat g = M^{−1} g、M_{RS} = Σ_μ z_μ^{−1} χ_R(C_μ) χ_S(C_μ) Π_i (t^{μ_i} − t^{−μ_i})/(t − t^{−1}) |
| 整数性 | \hat h_R = z^{−1} Σ N^{c=0} z^{2g} ν^Q、\hat g_R = Σ (N^{c=1} z^{2g} + N^{c=2} z^{2g+1}) ν^Q |

KN（本研究）の変数 (a, q, t) との対応：KN の unrefined は 𝐭 = −1 で a ↔ ν、q ↔ t（Mariño の t は KN の q に当たる）。本研究の docs/conventions.md を参照。

## refined（未確定）
- 候補となる入力：AS 1202.4456 の SO(2N) refined CS（変数 q, t、N）、GW/NRZ の Kauffman 超多項式（λ, q, t）、Cherednik の B/C 型 hyperpolynomial、Cherednik–Elliot の composite 超多項式。
- 注意（AS §7）：refined の orientifold では「偶数 d の引き算」の分母が (q^{d/2} − q^{−d/2})(t^{d/2} + t^{−d/2}) になり、q と t が混ざる。Σ_{d 奇数} を素朴に Macdonald 多項式に持ち上げる仮説は、この構造を見落とす可能性がある（2018 年の試行の検証項目）。

## refined（双対 Cauchy 型、K2'、2026-10-14）
- ホロノミー側：AS の D_N Macdonald = Koornwinder K_λ(x; q,t; 1,−1,q^{1/2},−q^{1/2})、unknot は x = t^ρ（x_i = t^{N−i}）での評価。U = t^N。
- probe 側（probe m 個、変数 y_1..y_m）：K_{λ'}(y; t,q; 1,−1,q^{1/2},−q^{1/2})。核は Mimachi の双対 Cauchy 恒等式 Π_{i,j} x_i^{−1}(1 − x_i y_j)(1 − x_i/y_j)。
- virtual 正規化：Z~ := (Π_j y_j)^N Z。probe 側の基底は (Π y_j)^N K_{N^m−λ}(y) = Rains の virtual Koornwinder \hat K_λ(y; t,q, T = U)。Z~ は y の多項式で、U を形式変数として log を展開する。
- 自由エネルギー：log Z~ = Σ L_{a…} y^a…、f = Möbius 反転（f_{ab} = Σ_{d|(a,b)} (μ(d)/d) L_{a/d,b/d}(q^d,t^d,U^d)）。
- Kähler 変数：ν = (U²/t)^{1/2}（Mariño の ν、Q = U²/t）。Hopf 絡み目では ν~ = ν (q/t)^{1/2} = U q^{1/2}/t が現れる。z_q = q^{1/2} − q^{−1/2}。

## trefoil（K3'、2026-10-14、仮説を含む）
- AS 1202.4456 (6.15)/(6.17) は unknot で正規化した [1] の値。t = q で Mariño の F_{3_1} の鏡像 × ν^{−6}（ν = 1/𝐚）【計算】。
- 標準 framing の unnormalized W_SO = (6.15) × r^{−5} s^{11} U^{−6} × unknot_SO（r = q^{1/2}, s = t^{1/2}）【仮説：t = 1 の極の打ち消しと最短性から】。T_[1] = q^{1/2} t^{N−1} に対し r^{−5}s^{11}U^{−6} = T_[1]^{−6}(q/t)^{1/2}。
- oriented 部分：½h^{ref}_[1] = 𝒫_DGR の鏡像(𝐚_U,𝐪,𝐭) × (ν−ν^{−1})/(𝐪−𝐪^{−1})、𝐚_U² = A (t/q)^{1/2}、A = ν² = U²/t（unknot から）。

## クロスキャップ展開（K4'、2026-10-14）
- 被覆：U(N_A) refined CS、A = ν² = U²/t = t^{N_A}（N_A = 2N_D − 1）、核 Π_i(1 − y x_i)(1 − y/x_i)（fermionic、composite）。oriented 部分 = ½F^{cov}。
- クロスキャップ部分（2026-10-17 に全 d の定義へ変更、report.md §7）：g_{ab} := f_{ab} − ½ f^{cov}_{ab}（f, f^{cov} は全 d の Möbius 反転）。第二量子化（AS §3）の構造 g = ½E/(1+t) + C（E：Ĩ で固定される oriented 状態、(偶,偶) のみ；C：凍結した unoriented 状態）をもち、2(1+t)g は整数 Laurent。
  Mariño 型の Σ_{d 奇数} の定義 g^{odd}_{ab} = L_{ab} − ½L^{cov}_{ab} − Σ_{d 奇数>1}(1/d) g^{odd}_{a/d,b/d}(x^d) とは g_{ab} = g^{odd}_{ab} − ½ g^{odd}_{a/2,b/2}(x²) で結ばれ、t = q では g^{odd} = Mariño の ĝ。refined では g^{odd} の (偶,偶) に (1+t^{d/2}) の分母が残るので、整数性の判定には全 d の g（と f 自身）を使う。
  理由：refined の軌道体では z₂ → −z₂ が probe（t-brane）に沿い、固定される oriented 状態のモード射影が ½/(1+t) を生む（閉弦では AS §7 の (t^{d/2}+t^{−d/2}) の項）。E と C の分離は未確定（E(t=−1) だけが分離に依らない）。
  E の決め方（2026-10-18）：E = P − 2Q（P = (1−t)f^{cov}_{固定セクター}、Q は Ĩ 固有値 −1 の状態の部分和）、t = −1 で E = 2(1+t)g、t = q で C = ĝ（Mariño）。(2,2),(2,4) で一意、(4,4) は未確定（report.md §7）。
- 変数：z_q = q^{1/2}−q^{−1/2}、z_t = t^{1/2}−t^{−1/2}、G = z_q z_t、w = (q/t)^{1/2}、ν~ = ν w。DGR：𝐪 = t^{1/2}、𝐭 = −w、𝐚 = ν、ν~ = −𝐚𝐭。

## refined のクロスキャップ BPS 数（2026-10-18、report.md §8、【仮説】）
- 奇数のセクターの [1]：g^{ref}_[1] の v = (qt)^{1/2} の整数冪の部分を c = 1、半整数冪の部分を c = 2 とする。c = 1 は Σ_β ν^β Σ N^{c=1}_{g,β,J} G^g w^{2J}（G = z_q z_t、w = (q/t)^{1/2}、𝐭 = −w）。
- 正値性は G 展開の後で、ν^β ごとの符号のそろい方として判定する。t = q で消えるブロックと c = 2 は普遍的な因子 (𝐭² − 1)（c = 2 は (𝐚² − 𝐪²)(𝐭² − 1)）を除いてから判定する。
