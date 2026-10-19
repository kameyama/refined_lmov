# 文献調査：Kauffman 多項式の large N duality（unrefined と refined）

作成 2026-10-14。書誌は arXiv API（export.arxiv.org）で確認した。TeX ソースは references/<id>/src/ にある（git の外、`just fetch-refs` で再取得）。

## 1. unrefined：Mariño の予想（0904.1088）

Mariño, "String theory and the Kauffman polynomial", Commun. Math. Phys. 298 (2010) 613–643, arXiv:0904.1088（TeX で確認）。

### 1.1 規約（Mariño §2）
- 変数 t, ν、z = t − t^{−1}。HOMFLY の skein：X − Y = (t − t^{−1}) I、カール = ν。self-writhe で正規化（標準の HOMFLY に ν^{2 lk} を掛けたもの）。
- HOMFLY 不変量 H(L) = P_L · H(unknot)、H(unknot) = (ν − ν^{−1})/(t − t^{−1})。量子群 U_q(sl_N) との対応は t = q^{1/2}, ν = t^N。
- Kauffman（Dubrovnik 版）の skein：X − Y = (t − t^{−1})(I − ∞)、カール = ν。G(L) = F_L · G(unknot)、G(unknot) = 1 + (ν − ν^{−1})/(t − t^{−1})。U_q(so_N) との対応は t = q^{1/2}, ν = q^{N−1}（**ν = t^{2N−2}**）。
- 彩色：HOMFLY は annulus の skein の基底 P_{(R,S)}（Hadji–Morton）で composite 表現 (R,S) の彩色 H_{(R,S)}；Kauffman は基底 y_R（Beliakova–Blanchet）で G_R。unknot の値はそれぞれ dim_q (R,S)、dim_q^{SO(N)} R。
- composite 表現：(R,S) = Σ_{U,V,W} (−1)^{|U|} N^R_{UV} N^S_{U^T W} (V ⊗ \bar W)（Koike の式。Mariño (2.19)）。dim_q (R,S) も同じ式で V, W の量子次元の積に展開する。C_{(R,S)} = C_R + C_S。
- framing の変更：H_R → (−1)^{f|R|} t^{f κ_R} H_R（Kauffman も同じ規則）。κ_R = Σ l_i(l_i − 2i + 1)。
- 例：trefoil の P = 2ν² − ν⁴ − z²ν²、F = 2ν² − ν⁴ + z(−ν³ + ν⁵) + z²(ν² − ν⁴)。Hopf 絡み目の F = z^{−1}{ν − ν^{−1} + z + z²(ν − ν^{−1})}。

### 1.2 構造の式（Mariño §3）
- composite 不変量：R_R(K) = Σ_{R₁,R₂} N^R_{R₁R₂} H_{(R₁,R₂)}(K)。向きに依らない。例：R_[1] = 2H_[1]、R_[2] = 2H_[2] + H_{([1],[1])}、R_[1,1] = 2H_[1,1] + H_{([1],[1])}。
- 母関数：Z_R(v) = Σ_R R_R s_R(v)、Z_G(v) = Σ_R G_R s_R(v)、F = log Z。
- reformulated 不変量：
  - F_R(v) = Σ_{d≥1} Σ_R (1/d) h_R(t^d, ν^d) s_R(v^d)
  - F_G(v) − (1/2) F_R(v) = Σ_{d 奇数} Σ_R (1/d) g_R(t^d, ν^d) s_R(v^d)
- 2 箱までの g_R（Mariño (3.14)）：
  - g_[1] = G_[1] − H_[1]
  - g_[2] = G_[2] − ½G_[1]² − H_[2] + H_[1]² − ½H_{([1],[1])}
  - g_[1,1] = G_[1,1] − ½G_[1]² − H_[1,1] + H_[1]² − ½H_{([1],[1])}
  - 3 箱は Mariño (3.15)（g_[1](t³,ν³) の項を含む）。
- 2 箱までの h_R（NRZ 1310.2240 付録 C；Mariño の定義から導ける）：h_[1] = 2H_[1]、h_[2] = 2H_[2] + H_{([1],[1])} − 2H_[1]² − H_[1](t², ν²)、h_[1,1] = 2H_[1,1] + H_{([1],[1])} − 2H_[1]² + H_[1](t², ν²)。
- 基底変換：M_{RS} = Σ_μ (1/z_μ) χ_R(C_μ) χ_S(C_μ) Π_i (t^{μ_i} − t^{−μ_i})/(t − t^{−1})（|R| = |S| のとき。LMOV と同じ）。\hat h_R = Σ_S M^{−1}_{RS} h_S、\hat g_R = Σ_S M^{−1}_{RS} g_S。2 箱では M = [[t, −t^{−1}], [−t^{−1}, t]]（順に [2], [1,1]）、M^{−1} = (t² − t^{−2})^{−1} [[t, t^{−1}], [t^{−1}, t]]。
- **予想（Mariño Conjecture 3.3、結び目）**：
  - \hat h_R ∈ z^{−1} Z[z², ν^{±1}]：\hat h_R = z^{−1} Σ_{g,Q} N^{c=0}_{R;g,Q} z^{2g} ν^Q
  - \hat g_R ∈ Z[z, ν^{±1}]：\hat g_R = Σ_{g,Q} (N^{c=1}_{R;g,Q} z^{2g} + N^{c=2}_{R;g,Q} z^{2g+1}) ν^Q
  - N^{c} は整数（正とは限らない）。絡み目は z の冪が L−2, L−1 だけずれる。framed でも成り立つとされる。
- 帰結：\hat g_[1] = g_[1] に z^{−1} の項が無い ⇔ k₀ = p₀（Lickorish の定理）。g_[2] の整数性 ⇔ G² ≡ H_{([1],[1])} mod 2（Rudolph の定理）。g_[4] の整数性 ⇒ G_[2]² ≡ H_{([2],[2])} mod 2（Morton–Ryder）。
- 証拠（Mariño §4）：unknot（\hat h_[1] = 2(ν − ν^{−1})z^{−1}、\hat h_[2] = \hat h_[1,1] = −z^{−1}、\hat g_[1] = 1、他は 4 箱まで 0）、trefoil の \hat g_[2], \hat g_[1,1]（(4.8) に明示。c=2 の項が非零）、Hopf 絡み目。トーラス結び目の H_{(R,S)} は knot operator の式 (4.2) で composite にも使える。
- 文字列論の解釈（Mariño §5）：Z_{SO/Sp} = ½ Z_cov + Z_unor。covering の寄与は、conormal bundle の変形 N_{K,±ε} の 2 組の brane から R と \bar S の Wilson line が出て、composite 表現 (R,S) の HOMFLY になる（テンソル積 R ⊗ \bar S ではなく composite であることは unknot の topological vertex 計算（BFM hep-th/0411227 (3.10)）で確認）。unoriented の寄与は X/I の crosscap c = 1, 2 の曲線を数える。d の和が奇数だけなのは unoriented の多重被覆の一般的な性質（SV, BFM）。SO と Sp の違いは c = 1 の項の全体の符号。

### 1.3 関連文献（unrefined）
| arXiv | 著者・題名 | 役割 |
|---|---|---|
| hep-th/0405083 | Bouchard–Florea–Mariño, Counting higher genus curves with crosscaps in Calabi–Yau orientifolds, JHEP 0412:035 | 閉弦の orientifold の整数性（c=1, 2、d 奇数） |
| hep-th/0411227 | Bouchard–Florea–Mariño, Topological open string amplitudes on orientifolds, JHEP 0502:002 | 開弦の orientifold。unknot と Hopf の topological vertex 計算、T(2,2m+1) の Kauffman の式。c=1 を parity で分離（c=2 は未決定）。M_{RS} の規約が Mariño と異なる（N_{R;g,Q} → (−1)^{|R|−1} N_{R^T;g,Q}） |
| hep-th/0505008 | Borhade–Ramadevi, SO(N) reformulated link invariants from topological strings, Nucl. Phys. B727 (2005) 471 | SO(N) の reformulated 不変量 |
| 1003.5282 | Paul–Borhade–Ramadevi, Composite invariants and unoriented topological string amplitudes | framed トーラス結び目・絡み目での Mariño の予想の検証 |
| 1003.2861 | Stevan, Chern–Simons invariants of torus links, Ann. Henri Poincaré 11 (2010) 1201 | トーラス絡み目の colored HOMFLY/Kauffman の knot operator の式と Mariño の予想の検証 |
| 1007.1656 | Lin Chen–Qingtao Chen, Orthogonal quantum group invariants of links | orthogonal LMOV 予想の定式化、列（column）の場合の証明、トーラス結び目の colored Kauffman の明示式 |
| 1310.2240 | Nawata–Ramadevi–Zodinmawia, Colored Kauffman homology and super-A-polynomials, JHEP 1401 (2014) 126 | 付録 C：4_1 の \hat h_[2], \hat h_[1,1], \hat g_[2], \hat g_[1,1] と BPS 整数の表（非トーラス結び目での最初の検証）。対称表現の colored Kauffman（trefoil, 4_1）の式 |
| 1702.06316 | Mironov–Morozov–Morozov–Sleptsov ほか（Web 検索の結果。原典未確認） | arborescent 結び目（8 交点まで）での Mariño の予想の検証 |

Web 検索（2026-10-14）の範囲では、Mariño の予想の完全な証明は見つからなかった（Chen–Chen が部分的に証明。未確認：2017 年の論文も「未証明」と記す）。

## 2. refined：既知の構成要素

### 2.1 Aganagic–Schaeffer（1202.4456）
"Orientifolds and the refined topological string", JHEP 09 (2012) 084（TeX で確認）。
- **閉弦の orientifold の refined の自由エネルギー**（AS §3 (3.7)）：unrefined では F(X/I) = ½F(X) + F_unor、F_unor は c=1（1/(q^{d/2}−q^{−d/2})^{1−2g}）と c=2（(q^{d/2}−q^{−d/2})^{2g}）の項、d は奇数か偶数だけ。
- **強い整数性**（AS §4）：Z = e^F が M 理論の指数であるためには、½ の因子が打ち消す必要があり、N^g_β ≡ N^{g/2,c=1}_{β/2}（g 偶数）、N^g_β ≡ N^{(g−1)/2,c=2}_{β/2}（g 奇数）mod 2。開弦（結び目）では \hat N_{([1],[1]),g,β} ≡ \hat N^{c=1}_{[1],g/2,β/2}（g 偶数）、≡ \hat N^{c=2}_{[1],(g−1)/2,β/2}（g 奇数）mod 2。unknot・trefoil・T(2,5) で確認済み（AS）。
- **SO(2N) refined CS**（AS §5）：D 型 Macdonald 多項式で S, T 行列を定義。S_{ij} = S_{00} M_{λ_i}(t^{−ρ}) M_{λ_j}(t^{−ρ} q^{−λ_i})、T_{ii} = g_i q^{(λ,λ)/2} t^{(λ,ρ)} ×定数、計量 g_λ（D_N の箱の積の式は AS 付録 B）。(q,t) LR 係数で knot operator を作り、トーラス結び目 Z(T(n,m)) = Σ K_{0k} N^k_{ij} (K^{−1})^j_l S^l_0。
- **refined Kauffman 不変量**（AS §6）：unknot S_{·V}/S_{00} = (t^{(2N−1)/2} − t^{−(2N−1)/2})/(t^{1/2} − t^{−1/2}) + 1。Hopf：q t^{N−1} + 1 + q^{−1}t^{−(N−1)} + (t^{(2N−3)/2} − t^{−(2N−3)/2})/(t^{1/2} − t^{−1/2})。変数変換 𝐪 = √t、𝐭 = −√(q/t)、𝐚 = t^{(2N−1)/2} で Hopf は Kauffman ホモロジー（GW）と一致。trefoil（AS (6.15)）：𝐭⁶𝐪²𝐚⁴ + 𝐭⁴𝐪^{−2}𝐚⁴ + 𝐭⁴𝐪𝐚³ − 𝐭⁴𝐪^{−1}𝐚³ − 𝐭⁴𝐪²𝐚² + 𝐭²𝐚² − 𝐭²𝐪^{−2}𝐚² − 𝐭²𝐪𝐚 + 𝐭²𝐪^{−1}𝐚。**一般のトーラス結び目では GW の Kauffman ホモロジーと変数変換で結べない**（AS の結論：SO(2N) refined CS は新しい不変量）。T(2,2m+1) の一般式 (6.17)、T(3,4) の SO(4), SO(6)。
- **large N**（AS §7）：Z = S_{00} の自由エネルギーは
  F = Σ_d (1/d) Q^d q^{−d/2} / ((q^{d/2} − q^{−d/2})(t^d − t^{−d})) + Σ_d (1/d) Q^{d/2} t^{d/2} q^{−d/2} / (q^{d/2} − q^{−d/2})、Q = t^{2N−1}。
  同じものを ½F_cov ＋ unoriented の形に書くと、第 3 項 −Σ_d (1/2d) Q^d t^{d/2} q^{−d/2} / ((q^{d/2} − q^{−d/2})(t^{d/2} + t^{−d/2})) が現れる（t = q で Σ_{d 奇数} の「偶数 d の引き算」に当たる項。**refined では分母が (q^{d/2} − q^{−d/2})(t^{d/2} + t^{−d/2}) で q と t が混ざる**）。orientifold は SU(2)_L × SU(2)_R を U(1)×U(1) に破るので、unoriented の BPS は spin の完全な多重項をなさなくてよい。
- 第 1 量子化（自由エネルギーの g, c による展開）の refined 版は AS でも未定義。開弦の refined の整数性は「Z の M 理論の指数としての整数性」だけが要請されている。

### 2.2 Gukov–Walcher（hep-th/0512298）
"Matrix factorizations and Kauffman homology"（TeX で確認。雑誌未確認）。
- 変数 λ, q, t。\bar F(λ,q) − \bar P(λ,q) = Σ N^{c=1}_{g,Q}(q − q^{−1})^{2g} λ^Q + Σ N^{c=2}_{g,Q}(q − q^{−1})^{2g+1} λ^Q（[1] の Mariño の g_[1] と同じ。トーラス結び目では N^{c=2}_{[1]} = 0）。
- Kauffman 超多項式：\bar F(λ,q) = (q − q^{−1})^{−1} Σ (−1)^r λ^Q q^s D_{Q,s,r}、\bar𝓕(λ,q,t) = Σ λ^Q q^s t^r D_{Q,s,r}。予想 1：HSO_N = \bar𝓕(λ = q^{N−1})、HSp_N(\bar K) = \bar𝓕(λ = q^{−N−1})。
- trefoil の unreduced 超多項式 (5.22)（符号あり）、reduced (6.14)（正）、T(2,2k+1) の一般式 (5.23)、4_1 の (5.25)。SO(4) = Khovanov²、Sp(2) = Khovanov(q²) で決める。
- 微分 d_N（次数 (−1, N−1, −1)）、universal な d_→, d_← で HOMFLY ホモロジーに落ちる（§6）。

### 2.3 Cherednik（1111.6195）、Cherednik–Elliot（1503.01441）
- Cherednik, "Jones polynomials of torus knots via DAHA"（arXiv:1111.6195、IMRN に掲載と Web 検索の結果で推定、未確認）：DAHA による超多項式。B, C 型の hyperpolynomial（Kauffman の refined 化、追加の変数 1 つ）。
- Cherednik–Elliot, "Refined composite invariants of torus knots via DAHA"（arXiv:1503.01441）：composite DAHA 超多項式（ヤング図の対）。composite HOMFLY-PT（annulus の skein）の refined 化。（TeX の読解は未了：本文の見出しが grep に掛からない。要確認。）

### 2.4 Nawata–Ramadevi–Zodinmawia（1310.2240）
- 対称表現の colored Kauffman ホモロジー（四重次数）。trefoil と 4_1 の [r] 彩色の refined Kauffman の式（"triple Kauffman"：3 変数 (a, q, t) の [r] 彩色）。2018 年のノートブックの tripleKauffmanTrefoil, TripleKauffmanFig8 はこれ（未確認、照合予定）。
- 付録 C（unrefined）：4_1 の \hat h, \hat g と表。

## 3. 既知の結果の範囲と新規性（暫定）
- 既知：Mariño の予想（定式化と検証）、AS の SO(2N) refined CS と refined の large N（閉弦の自由エネルギーのみ）、GW の Kauffman 超多項式、NRZ の対称表現の refined Kauffman、Cherednik–Elliot の composite 超多項式。
- Web 検索（2026-10-14、refined Kauffman large N / orientifold refined topological string knot / colored Kauffman homology annulus composite）の範囲では、**結び目の refined Kauffman 不変量に対する Mariño 型（g_R, h_R）の refined 版の定式化と検証は見つからなかった**。AS 自身が「refined では第 1 量子化の定義が無い」と書いている。
- 本研究の新規性の候補：refined の g_ρ（unoriented 部分）と h_ρ（composite 部分）の定義、整数性と（何らかの）正値性の検証。

## 4. U(N)/SU(N) 以外の refined CS の計算の範囲（2026-10-19 の調査）
確認の度合い：【TeX】は references/ の原典の TeX で確認、【抄録】は arXiv の抄録（Web 検索の結果）で確認、【未確認】はそれ以外。
| 文献 | 群 | 計算されているもの | 確認 |
|---|---|---|---|
| Aganagic–Schaeffer 1202.4456 | 単純レース型（ADE）で定義。計算は SO(2N) だけ | unknot、Hopf（V, V）、trefoil の V（一般の N、(6.15)）、trefoil のスピノル S（SO(4), SO(6), SO(8)）、T(2,2m+1) の V の一般式、T(3,4) の V（SO(4), SO(6)）、大きな N の閉弦の自由エネルギー。非単純レース型は脚注で「(2,0) 理論の外部自己同型のひねりで得られるはず」とするだけ | 【TeX】 |
| Cherednik 1111.6195 | 任意の被約ルート系と C^∨C_1 | DAHA-Jones 多項式（トーラス結び目）。§Kauffman：C–D（部分的に B）型の hyper-polynomial（4 変数、予想と例、「最小の重みでも正値性はほとんど無い」）。§Exceptional：例外型（G_2 など）は t = q の値だけ、最も簡単なトーラス結び目、最小・準最小の重み | 【TeX】 |
| Cherednik–Danilenko 1509.08351 | B, C, D | 「B, C, D 型の DAHA-hyperpolynomial の存在と対称性は未解決」、D 型の t = q が Kauffman 多項式に一致することは示された | 【TeX】 |
| Cherednik–Elliot 1503.01441 | A（composite）、例外型 | composite DAHA 超多項式。最後の節で例外型の Deligne–Gross 系列の DAHA-hyperpolynomial を実験的に構成 | 【抄録】 |
| Elliot–Gukov 1505.01635 | 例外型 | 例外型の結び目ホモロジーの微分の構造（トーラス結び目）、E_6 の hyperpolynomial | 【抄録】 |
| Avetisyan–Mkrtchyan 2205.12832（v2） | A, B, C, D | S³ の分配関数（閉弦）の多重サイン関数表示と GV 型の展開。D, C に c = 1（奇数の n）と「Klein 瓶」（c = 2、全 n、分子 q/t − t/q、y = 1 で 0）、B は c = 1 だけで分母に (t^{n/2} + t^{−n/2})。Kähler 変数の「三重化」τ₀, τ₁, τ₂。結び目は無い。本研究で Klein 瓶の項が全 n の凍結した状態＋不変な状態（AgSc §6 の形）の並べ替えであることを確認した（theory_orientifold_index.md §8） | 【TeX】 |
| Avetisyan–Mkrtchyan 2107.08679、2302.14319 | 任意（非単純レース型を含む） | S³ の分配関数（閉弦）。結び目は無い | 【抄録】（TeX は取得済み、未読） |
| Piazzalunga–Uranga 1405.6019（v2） | — | real topological string の M 理論の解釈。Ĩ で不変な状態だけが奇数の巻き付きで寄与し、real BPS 不変量は偶奇で重みを付けた差 n⁺ − n⁻。被覆の取っ手一つの捻ったトレースは −2 sinh（z 一つ）。refined は脚注で触れるだけ | 【TeX】 |
| Kim–Wei 2201.12264（v2） | D 型の箙（5d） | ON 面（O5 の S 双対）のある 5-brane web の refined 位相的頂点。向き付けられない面は数えない（本研究の照合の対象外） | 【TeX、序論のみ】 |
| Bishler–Mironov 2504.13831 | 単純レース型 | Macdonald 次元（q^ρ での値）の Vogel 普遍性 | 【抄録】 |
| NRZ 1310.2240、Gukov–Walcher hep-th/0512298 | SO/Sp | Kauffman ホモロジー（refined CS とは別物と明記） | 【TeX】 |
- 非トーラス結び目、高い彩色の一般の N の refined SO/Sp の結び目不変量、Sp の refined CS の結び目不変量は、今回の検索の範囲では見つからなかった。
- 2605.12469（"A note on universality in refined Chern-Simons theory"）は題名だけを見た【未確認】。
