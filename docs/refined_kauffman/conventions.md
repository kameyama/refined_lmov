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
