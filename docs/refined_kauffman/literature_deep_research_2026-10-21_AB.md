# refined SO(2N) Chern–Simons の large N 双対におけるクロスキャップ部分：問題 A・B の文献調査報告

結論から書く。問題 A（向き付けられない曲線の refined 多重項）と問題 B（Ĩ 不変状態の数 E を決める原理）のどちらにも、直接答える文献は見つからなかった。いちばん近いのは数学側の orientifold Donaldson–Thomas 理論（Young、Franzen–Young、Bu）である。そのうち Young の arXiv:1603.05401 は、「orientifold の不変量は通常の DT 不変量からは計算できず、ℤ₂ 同変な精密化が要る」と原典の本文で明言している。これは、E が指数から 2(1+t) を法としてしか決まらないというこの研究の結果を、数学の側から言い換えたものと読める（解釈）。

## TL;DR

- **A（refined の向き付けられない多重項）：見つからなかった。** refined な real GV 不変量、refined real vertex、向き付けられない面を数える refined O5/ON 頂点は、いずれも確認できなかった。現存するのは unrefined の M 理論的解釈（Piazzalunga–Uranga の後続 1508.06644）、数学の局所 real GV 公式（Georgieva–Ionel 1812.02505、Klein TQFT）、motivic（すなわち refined）な orientifold DT 不変量（Young 1408.4888・1603.05401、Franzen–Young 1605.06596、Bu 2503.20667）である。A3 については、Koornwinder の duality（q, t を固定し境界パラメータを変える）は σ ではない。σ に対応しうるのは、Rains の lifted Koornwinder 多項式に対する Macdonald 対合の作用である（抜粋で存在だけ確認）。
- **B（E の決定）：見つからなかった。** 境界をもつ曲線への Piazzalunga–Uranga 型 t⁺ − t⁻ の拡張も、「差 z_q が和 z_q⁺ になる」という境界の捻ったトレースも見当たらない。最も近いのは Young 1603.05401 である。そこでは、CoHA 上の自由加群でない CoHM（cohomological Hall module）の分解 Σ_e [H_Q(e)]·Ω^σ_e = [M_Q] と、loop quiver では Ω_Q が ℤ₂ 同変 DT を決めるが一般には決めないという主張が、原典の本文で確認できる。Hopf では E が一意に決まり結び目では決まらないという観察と、構造が合う（解釈）。
- **検証できる予言：** (A) σ は太字変数では 𝐪 ↦ −(𝐪𝐭)^{−1}、𝐭 ↦ 𝐭、𝐚 ↦ 𝐚 である。A 型の類推からは、σ が対応づけるのは同じ R どうしではなく R と R^T だと予想される。したがって σ(ĝ^{c=1}_{[2]}) と ĝ^{c=1}_{[1,1]} を比べることを勧める。また普遍因子のうち (𝐚² − 𝐪²) は U = ±t（N = 1、SO(2)）で消えるので、refined だけの部分は D_1 の自明性から来ていると読める。(B) Hopf の方法（固定セクターの指数、t = −1、t = q の条件）は、quiver が一頂点・loop 型に当たる unknot の多箱セクターでも E を一意に決めるはずである。失敗が始まるのは、矢でつながった多頂点 quiver に当たる結び目からだと予言できる。

## 調査の範囲と確認度の表記

確認度は次の 3 段階で書く。
- 「原典本文で確認」：arXiv の abs ページ、または本文 HTML/PDF を取得して読んだ箇所。
- 「本文抜粋（検索経由）」：検索結果に表示された原典 PDF/HTML の本文の断片だけを見た箇所。
- 「抄録のみ」／「二次情報」：抄録だけ、または他論文の参考文献リスト・講演スライド・紹介サイトによる情報。

書誌のうち著者や掲載誌を今回見ていないものは「未確認」と書く。本調査で使えた検索は 18 回である。加えて、サブエージェント 1 件に Kauffman skein 値の曲線計数と orientifold quiver の追跡を依頼した。

---

## 問題 A：凍結した状態（向き付けられない面）の refined の多重項

### A1. 向き付けられない曲線に巻き付いた M2 の refined 多重項、捻った Lefschetz トレース、refined real GV 不変量

**結論：見つからなかった。** 反正則対合のもとでの sl₂（Lefschetz）作用の捻ったトレースを使って refined real GV 不変量を定義した文献は、物理側・数学側のどちらにも確認できなかった。近い文献を以下に挙げる。

**(1) Hayashi–Piazzalunga–Uranga, arXiv:1508.06644, "Towards a gauge theory interpretation of the real topological string"（Phys. Rev. D 93, 066001, DOI 10.1103/PhysRevD.93.066001）**
- (b) 主張：非コンパクトなトーリック CY 上の real topological string を、O4 面と、その上に固定された 1 枚の D4 をもつ IIA の orientifold として扱う。M\[1\] 理論への持ち上げから real GV 不変量（巻き付いた M2 の BPS 状態）を多くの幾何で計算する。flop 転移で real topological string の振幅が決まったパターンの符号を拾い、それが real BPS 不変量の連続性と整合することを示す。orientifold は 8 超電荷のゲージ理論の新種の surface operator とみなせる。\[1\]
- 当てはめ（解釈）：Piazzalunga–Uranga（1405.6019）の t⁺ − t⁻ を、計算できる例に展開した直接の後続である。ただし unrefined で、SU(2)_L × SU(2)_R の破れ方や取っ手あたりの因子の refined 版は抄録には出てこない。この研究にとっての価値は二つある。一つは「flop を越えたときの符号の規則」が、Ĩ の固有値 ε の符号の決め方の手本になりうること。もう一つは、O4 と固定された D4 の系が、refined にする余地のある gauge theory（Nekrasov 型）の言葉で書かれていることである。
- (c) 確認度：抄録のみ（arXiv abs で確認）。

**(2) Narain–Piazzalunga–Tanzini, arXiv:1612.07544, "Real topological string amplitudes"（掲載誌：JHEP と二次情報にあるが、巻号は未確認）**
- (b) 主張：type I（あるいは orientifold つきの type II）の物理的超弦の相関関数のうち、real topological string の振幅を計算するものを論じる。世界面の Euler 標数 χ を固定した real 振幅 G_χ の正則微分を、orientifold 不変部分に制限した N = 2 Weyl 多重項の g′ = −χ + 1 乗の結合として同定する。orientifold\[2\] の作用は X_D への ℤ₂ 反転を伴う。\[3\]
- 当てはめ（解釈）：unrefined の z^{2g−2+c} という冪の数え方を、Weyl 多重項の冪 g′ = −χ + 1 として物理的に導いている。refined 化（Ω 背景で ε₊ ≠ 0）の扱いは抜粋には見えない。refined 版の取っ手あたりの因子を超弦の振幅から導く場合、この論文が出発点になる。
- (c) 確認度：本文抜粋（検索経由）。

**(3) Georgieva–Ionel, arXiv:1812.02505, "A Klein TQFT: the local Real Gromov–Witten theory of curves"（Adv. Math. と二次情報にあるが、巻号は未確認）**、および同著者の arXiv:2005.05928, "Splitting formulas for the local real Gromov–Witten invariants"（J. Symplectic Geom. と二次情報）
- (b) 主張：実曲線上の局所 3-fold L ⊕ c*L → Σ の real GW 理論が、向き付けられない面の圏の拡張の上の 2 次元 Klein TQFT になることを示し、表現論のデータによる閉じた公式を与える（Bryan–Pandharipande の拡張）。帰結として、連結な\[4\] real GW 不変量を整数の不変量で表す局所 real Gopakumar–Vafa 公式（Theorem 1.2／§10 の Theorem 10.1）を得る。クロスキャップ（Möbius\[5\] の帯）の像は特別な元 U で、その係数 U_ρ は構造定数の正の平方根である。向きを反転する tube も生成元に入る。\[5\]
- 当てはめ（解釈）：数学側で「real GV 不変量」を定義・証明した最も確かな文献だが、unrefined で閉弦だけを扱う。この研究が問う「Jacobian のコホモロジーへの反正則対合の作用」は使われておらず、TQFT の切り貼りで済ませている。それでも、クロスキャップの係数が「構造定数の平方根」という符号の決まった量で与えられる点は示唆的である。refined の c=2 部分が「−(𝐚²−𝐪²)(𝐭²−1) × 正の多項式」になるという観察を説明する、数学側の候補になりうる。
- (c) 確認度：抄録と本文抜粋（検索経由）。

**(4) Young, arXiv:1408.4888, "Self-dual quiver moduli and orientifold Donaldson–Thomas invariants"（Commun. Number Theory Phys. 9(3) (2015) 437–475。掲載誌は所属機関のページによる二次情報）**
- (b) 主張：orientifold での BPS 状態の数え上げを動機に、反変対合をもつ quiver の自己双対表現のモジュライを調べる。Hall module の方法で、半安定な自己双対表現のモジュライスタックの有限体上の点の数を計算する。安定性を変えたときの壁越え公式が、物理文献で提案された orientifold DT/BPS 不変量の原始的な壁越えを再現する。有限型では、量子トーラスの表現に作用する量子二重対数の恒等式として理解できる。\[6\]
- (c) 確認度：抄録と本文抜粋（検索経由）。

**(5) Young, arXiv:1603.05401, "Representations of cohomological Hall algebras and Donaldson–Thomas theory with classical structure groups"（Commun. Math. Phys. 380 (2020)。頁は 272–322 と二次情報にある）**
- (b) 主張（序論を原典本文で確認）：自己双対表現のスタックのコホモロジー M_Q = ⊕_e H^•(M_e^σ){ℰ(e)/2} を CoHA H_Q 上の加群（CoHM）として構成する。motivic な orientifold DT 不変量 Ω^σ_{Q,e} を、その最小生成部分対象の類として定義し、σ 対称な quiver について整数性を証明する（Theorem A = Theorem 3.4）。M_Q\[7\] は H_Q 上の自由加群ではない。Conjecture A（Conjecture 3.7）は、部分代数 H_Q(e) ⊂ H_Q による分解 ⊕_e H_Q(e) ⊠ W^prim_{Q,e} ≅ M_Q を予想し、Grothendieck 群では Σ_e [H_Q(e)]·Ω^σ_{Q,e} = [M_Q] となる。さらに本文は「一般に Ω_Q を知っても [H_Q(e)] は決まらず、Ω^σ_Q は Ω_Q から直接計算できない。代わりに Ω_Q の ℤ₂ 同変な精密化が要る」と述べる。CoHM は「orientifold 弦理論の BPS 状態の空間の数学的モデル」と位置づけられている。\[7\]\[8\]
- 当てはめ（解釈）：これが A1 で求めている「refined real GV 不変量」に最も近い数学的対象である。motivic（Hodge の重みの次数つき）なので、unrefined ではなく refined な情報（スピンに当たる次数）をもつ。ただし対象は quiver（と潜在関数）に限られ、曲線の Jacobian の言葉には翻訳されていない。knots–quivers 対応（後述）と組み合わせれば、結び目の orientifold セクターをこの枠組みで扱う道がありうるが、そうした文献は見つからなかった（サブエージェントの調査でも同じ）。
- (c) 確認度：原典本文で確認（序論と主定理の記述）。

**(6) Franzen–Young, arXiv:1605.06596, "Cohomological orientifold Donaldson–Thomas invariants as Chow groups"（Selecta Math. (N.S.) 24 (2018) 2035–2061）**
- (b) 主張：σ 対称な quiver について、cohomological orientifold DT 不変量が σ 安定な自己双対表現のモジュライ空間の有理 Chow 群と同型であることを証明する。応用として、古典 Lie 環の安定な m 組のモジュライの Chow Betti 数を数値的に計算できることを示す。半安定 Hall module の cohomological 壁越え公式も与える。改訂版で、ODT 不変量に壁越えが起きない結果（Thm. 4.3）が加わった。\[9\]
- 当てはめ（解釈）：「refined の orientifold 不変量は固定点モジュライの正の Betti 数（Chow 群）である」という幾何学的な正値性の定理である。この研究の「[1] の c=2 部分は (𝐚²−𝐪²)(𝐭²−1) × 正の多項式」といった正値性の観察に、数学的な裏付けの型を与える。
- (c) 確認度：抄録のみ（arXiv abs で確認）。

**(7) Bu, arXiv:2503.20667, "Orthosymplectic Donaldson–Thomas theory"（前身は arXiv:2302.00038, "Enumerative invariants in self-dual categories. I. Motivic invariants"。著者本人のページで「superseded」とされている）**
- (b) 主張：線形圏の直交・シンプレクティックな対象を数える DT 不変量を構成する（構造群 GL(n) を O(n)、Sp(2n) に置き換える）。motivic Hall 代数とその orthosymplectic 版（motivic Hall module）で定義し、壁越え公式を証明する。例として、CY3 上の直交・シンプレクティックな perfect complex を数える DT 不変量と、潜在関数つき自己双対 quiver の自己双対表現を数える DT 不変量を定義する。\[10\]
- 当てはめ（解釈）：Young の quiver 版を、CY3 の層の圏（resolved conifold の orientifold を含みうる）に広げる一般論である。ただし、対合が「反正則な I と時空の z₂ → −z₂ の組」（Aganagic–Schaeffer §3）に当たるかどうかは、抄録からは判断できない。Bu の自己双対構造は双対関手（層の双対）であり、反正則な対合とは異なる型の orientifold である可能性が高い（解釈、未確認）。
- (c) 確認度：抄録のみ。

**(8) 既読文献 Aganagic–Schaeffer 1202.4456 のうち、依頼書に書かれていない箇所（本文抜粋、検索経由）**
本文には次の三つの記述が見える。
- 「I で固定される点こそ、固定された実の部分に巻き付いた BPS 状態である」\[11\]
- 「これは "real"（向き付けられない）GV 不変量と向き付けられた GV 不変量が 2 の倍数だけ異なることを意味する」\[11\]
- orientifold の固定面が、実の部分 RP² 上の実直線束に巻き付く例、および Klein 瓶の項（種数 0、クロスキャップ 2）の分子の因子に関する記述\[11\]

解釈：前者の「2 を法とした合同」は、問題 B で E が 2(1+t) を法としてしか決まらないことの unrefined の原型である。refined では、法が 2 から 2(1+t) に「1 + t 倍」されていることになる。

### A2. 2012 年以降の orientifold の refined 開弦 BPS 展開（第一量子化）

**結論：見つからなかった。** refined の開弦について第一量子化の BPS 展開を提案した文献は確認できなかった。閉弦・行列模型・5-brane web の側に近い文献がある。

**(1) Krefl–Walcher, arXiv:1207.1438, "ABCD of beta ensembles and topological strings"（JHEP 11 (2012) 111）**
- (b) 主張：B_N、C_N、D_N の固有値測度をもつ β アンサンブルと refined topological string の関係を調べる。古典固有値アンサンブルとその位相的弦の対応物が、変形・特殊化・量子シフト・離散商で互いに関係することを示し、対応する局所幾何上の refined orientifold との関係を議論する。\[12\]
- 当てはめ（解釈）：refined の閉弦の orientifold を、行列模型（β アンサンブル）側から扱った 2012 年の文献である。開弦の BPS 展開ではないが、SO(2N) refined CS の行列模型表示（Aganagic–Schaeffer の SO(2N) refined CS の S³ 上の測度は D 型 Macdonald 測度）と直結する。検索結果には、Krefl–Walcher 系列の論文（どの論文かは特定できず）に「orbifold 半径 β = 2 の理論は β = 1 の理論の orientifold と同一視できる場合がある」という一文も見えた。refined\[13\] の変数で orientifold を「β の特殊値」として捉える見方であり、帰属は未確認である。
- (c) 確認度：抄録のみ（掲載誌のページ）。

**(2) Krefl, arXiv:1001.5031, "Wall Crossing Phenomenology of Orientifolds"** — 他論文（arXiv:1506.03907）の参考文献リストで番号と題名を見ただけで、内容は未確認（二次情報）。Young\[14\] 1408.4888 が「物理文献で提案された orientifold DT/BPS 不変量の原始的な壁越え」と呼ぶものの候補である（解釈）。

**(3) 5-brane web の O 面頂点（unrefined または ON のみ）**
- Kim–Yagi, arXiv:1709.01928, "Topological vertex formalism with O5-plane"（Phys. Rev. D 97, 026011）：(p,q) 5-brane web に O5 面がある場合の新しい頂点規則を提案し、5d Sp(1)（N_f = 0, 1, 8）で検証する。本文は、O5\[15\] で鏡映された部分図では頂点の Young 図の順序が逆になることを述べ、real topological string の頂点（Walcher、Krefl–Walcher）を「別の定式化」として挙げる。unrefined。確認度は抄録と本文抜粋。\[16\]\[17\]
- 既読の arXiv:2201.12264（ON 面の refined 頂点）のうち、依頼書に書かれていない箇所：本文（抜粋、検索経由）は、O5 面による鏡映の像では 𝔱 と 𝔮 の割り当てが入れ替わり、そのため refined の頂点で扱うのが難しく、以後は ON 面に限ると明言している。解釈：O5\[18\] 型の鏡映で q ↔ t が入れ替わるという観察は、この研究の σ（q ↦ t^{−1}, t ↦ q^{−1}）と同じ型の対合が refined の orientifold で自然に現れることを示す。Aganagic–Schaeffer の Ĩ は z₂ → −z₂ で一方の平面だけを反転するので、O5 の鏡映（両平面の入れ替えを伴う）とは別物だが、σ 非対称の起源を考えるうえで比べる価値がある。
- arXiv:2412.19655, "O-vertex, O7⁺-plane, and Topological Vertex"（JHEP04(2025)182）と arXiv:2510.15662, "Topological Vertex for Symmetric matter" は、いずれも Sung-Soo Kim, Xiaobin Li, Futoshi Yagi, Rui-Dong Zhu による。前者は O5 面をもつ 5d SO(N) の O-vertex を扱い、抄録によれば unrefined の分配関数を対象とする。後者は O7⁺ 面に対して「Z₂-vertex」と「FD-vertex」を導入し、その Appendix B は「real topological vertex との比較」で、O5/O7 面つきの 5-brane web の頂点は real topological vertex とは別物だと明言している（本文抜粋）。向き付けられない面の数え上げではない。

**(4) unrefined の閉弦の後続（参考）**
- arXiv:2309.12046, Mariño–Schwick（Univ. Genève）, "Non-Perturbative Real Topological Strings"（SIGMA 22 (2026) 039, doi:10.3842/SIGMA.2026.039）：real topological string の resurgence。抄録によれば、円板を数える整数不変量が resurgent 構造の Stokes 定数として現れ（"the integer invariants counting disks appear as Stokes constants in the resurgent structure"）、local P² でその証拠を示す。本文抜粋では、追加の Stokes 定数は D-brane と orientifold 面がある場合の BPS 状態の数え上げに関係する、と述べる。
- arXiv:2602.04539, "Crosscap resurgence and Riemann–Hilbert problems for orientifolded conifolds"（著者未確認、v3 は 2026-09-26）：SO/Sp CS の large N 双対の自由エネルギーを、resolved conifold の部分と奇の射影を受けたクロスキャップの部分に分ける。有限のクロスキャップ射線に台をもつ階数 5 の有理 BPS 構造を構成し、正規化された double sine 関数による解を与える（本文抜粋）。解釈：double\[19\] sine 関数は二つの周期をもつので refined 化と相性がよいが、抜粋からは refined を扱っているかどうか判断できない。
- Avetisyan–Mkrtchyan 2205.12832（既読）の参考文献には、"Mkrtchyan, Refined Chern-Simons / topological string for classical gauge groups: conifold singularity, in preparation" とある（本文抜粋）。この後続が公開されたかどうかは確認できなかった。\[20\]

### A3. σ: q ↦ t^{−1}, t ↦ q^{−1} で対称でない量を展開する Koornwinder／D 型 Macdonald の枠組み

**結論：σ 非対称な量を自然に展開する基底は見つからなかった。** 関係する構造は三つある。

**(1) Koornwinder の duality は σ ではない。**
- Sahi, arXiv:q-alg/9710032, "Nonsymmetric Koornwinder polynomials and duality"（Annals of Math. 1999 と二次情報）：6 パラメータの DAHA を定義し、非対称 Koornwinder 多項式を導入して、Macdonald–Koornwinder の duality 予想を証明する（原典の冒頭を抜粋で確認）。\[21\]
- Philippe Di Francesco と Rinat Kedem の arXiv:2112.09798, "Macdonald duality and the proof of the Quantum Q-system conjecture"（抄録によれば、証明は Macdonald 多項式と Koornwinder 多項式の duality に強く依拠する）の式 (3.11)–(3.12)（本文抜粋）：duality は DAHA の反対合から来る。q, t を固定したまま、境界パラメータを a* = (abcd/q)^{1/2}、b* = −(q ab/(cd))^{1/2}、c* = (q ac/(bd))^{1/2}、d* = −(q ad/(bc))^{1/2} に替え、P_λ^{(a,b,c,d)}(q^μ t^{ρ*})/P_λ(t^{ρ*}) = P_μ^{(a*,b*,c*,d*)}(q^λ t^ρ)/P_μ(t^ρ) とする。
- 解釈：これは境界パラメータの対合であり、(q, t) を入れ替える σ とは別物である。SO(2N) に当たる特殊化でこの duality がどう働くかは、この研究の「super-duality が無い」こととは独立の問題である。

**(2) Macdonald 対合（(q,t) ↔ (t,q) と λ ↔ λ′）は、lifted Koornwinder 多項式に作用する。**
- Rains, arXiv:math/0112035, "BC_n-symmetric polynomials"（Transformation Groups 10 (2005) 63–132）：Okounkov の BC_n 対称補間多項式と Koornwinder 多項式についての基本文献で、§6–7 で lifted interpolation functions と lifted Koornwinder 対称関数 K̃_λ(q, t, T; t₀, t₁, t₂, t₃) を導入する（他論文 1405.7035 による要約、二次情報）。T = t^n で通常の Koornwinder に戻る。\[22\]\[23\]
- Rains・Vazirani（PDF 冒頭の抜粋で著者表記を一部確認）, arXiv:math/0606204, "Quadratic transformations of Macdonald and Koornwinder polynomials"：本文抜粋に「lifted Koornwinder polynomials への Macdonald 対合の作用は、これを双対化して…」という一文がある。lifted\[24\] Koornwinder に Macdonald 対合が作用する公式があることまでは確認したが、具体形は未確認。
- 関連：arXiv:1506.02755, "Bounded Littlewood identities"（著者未確認）は lifted Koornwinder と virtual Koornwinder K̂_λ(q, t, Q; t₀..t₃) を復習し、Q = 0 で virtual Koornwinder を Macdonald 多項式で表す（[108, Corollary 7.21] の引用）。Chul-hee\[25\] Lee, Eric M. Rains, S. Ole Warnaar の arXiv:2007.03174, "An elliptic hypergeometric function approach to branching rules"（SIGMA 16 (2020) 142, doi:10.3842/SIGMA.2020.142。抄録によれば、古典的な分岐則の多くについて Macdonald 型の変形を証明する）は、Mimachi の Cauchy 恒等式（Mimachi の論文の Theorem 2.1）と lifted Koornwinder を使い、普遍指標 sp_λ, o_λ, so_λ を lifted Koornwinder の特殊化として得る。t = q, t = (q^{1/2}, −q^{1/2}, q, −q) で C_n 測度になる例も挙げている。
- 解釈：large N の変数 T（この研究の U = t^N に当たる）を独立変数として保ったまま Macdonald 対合を作用させられる唯一の既知の枠組みが、lifted Koornwinder である。Macdonald 多項式では P_λ(x; q^{−1}, t^{−1}) = P_λ(x; q, t) なので、σ は「(q,t) ↔ (t,q)」と実質同じである。したがって A 型の super-duality（σ と転置 λ ↔ λ′ の組）の BC 型の類似を探すなら、ω_{q,t} の lifted Koornwinder への作用（パラメータ t₀..t₃ と T の変換を伴うはず）が候補になる。ただしこれは「右辺の基底の取り替え」と同じ方向に見える。依頼書で失敗とされた基底（Hall–Littlewood、Macdonald P/Q、modified Macdonald）には lifted Koornwinder（dual パラメータつき）は含まれていないが、重複する可能性は指摘しておく。

**(3) σ を太字変数で書くと「𝐪 ↦ −(𝐪𝐭)^{−1}, 𝐭, 𝐚 固定」になる（計算による確認、解釈）。**
𝐪 = t^{1/2} ↦ q^{−1/2} であり、𝐪𝐭 = −q^{1/2} だから σ(𝐪) = −(𝐪𝐭)^{−1}。また σ(𝐭) = −(t^{−1}/q^{−1})^{1/2} = 𝐭 で、ν は固定されている。A 型の色付き超多項式では、転置 R ↔ R^T と 𝐪 の反転を組み合わせた「mirror 対称性」が知られている（Gukov–Stošić, arXiv:1112.0030 "Homological algebra of knots and BPS states" が典拠とされる。番号は他論文の参考文献による二次情報で、式は未確認）。unrefined の SO/Sp では SO(±N) ↔ Sp(∓N) の対応がある（Astorino, arXiv:1005.3861, "Kauffman knot invariant from SO(N) or Sp(N) Chern-Simons theory and the Potts model"、Phys. Rev. D 81, 125026、抄録で確認。摂動の 2 次まで）。解釈：SO\[26\]\[27\] 側の σ が同じ R の中で閉じず、R と R^T（さらに N の符号反転、すなわち Sp 側）を対応づけると考えると、「[1] では σ 対称、多箱では σ 非対称」というこの研究の観察は自然に説明できる。[1]^T = [1] だからである。

### A4. 2013 年以降の refined SO(2N) CS と Kauffman ホモロジーの関係

**結論：refined SO(2N) CS と Kauffman ホモロジーを直接結ぶ新しい文献は見つからなかった。**

- Nawata–Ramadevi–Zodinmawia, arXiv:1310.2240, "Colored Kauffman homology and super-A-polynomials"（本文抜粋で確認）：序論に「SO(N)/Sp(N) では refined CS 不変量は Kauffman ホモロジーの Poincaré 多項式と異なる [25, 26]」とある。彩色 Kauffman ホモロジーは refined exponential growth を使って決めている。\[28\]
- Aganagic–Schaeffer の該当箇所（他文献の引用で見た本文、二次情報）：「SO(2N) では、Hopf 絡み目については（SU(N) と同様に）単純な変数変換が存在して Kauffman ホモロジーにつながるが、一般のトーラス結び目にはそのような変数変換が存在しない」。解釈：Hopf\[29\] 絡み目だけ変数変換が効くという事情は、問題 B で Hopf だけ E が一意に決まった事情と並行しているように見える。どちらも、Hopf の BPS 状態が「打ち消しのない」単純な構造をもつことの表れかもしれない。
- Wang–Yang–Zhang–Nawata の arXiv:2012.12008, "On knots, complements, and 6j-symbols"（Ann. Henri Poincaré, DOI 10.1007/s00023-021-01033-4。サブエージェントが本文抜粋で確認）：薄い結び目について、[r²] 彩色 HOMFLY-PT ホモロジーから [r] 彩色 Kauffman ホモロジーへの次数の付け替え規則を与え、抄録によればこの規則は同型 (so_6,[r]) ≅ (sl_4,[r²]) に由来する。ほかに（t 次数つきの）SO 型 super-A-多項式と SO(N) の 6j 記号を与える。refined SO CS との比較はなく、quiver 公式もない。
- Petrou–Hikami（OIST）の arXiv:2412.04933, "The HOMFLY-PT polynomial and HZ factorisation"（本文抜粋）：HZ 分解可能な結び目では、HOMFLY-PT と Kauffman の関係が向き付けられない弦の 2 クロスキャップ BPS 不変量の消滅 N̂^{c=2} = 0 と同値である（抄録："this relation is equivalent to the vanishing of the two-crosscap BPS invariants in topological string theory"）。解釈：c = 2 部分が消える結び目族の判定基準として、この研究の c=2 部分の正値性（−(𝐚²−𝐪²)(𝐭²−1) × 正）を refined で比べる材料になる。
- Mariño 0904.1088（既読）の §6（結論、本文抜粋）のうち依頼書にない箇所：N^{c=0,1,2} は BPS 状態に付随するコホモロジーの Euler 標数であり、それらが彩色 Kauffman 不変量の圏化を与えると予想するのが自然だと述べる。さらに、絡み目や高次の表現では Kauffman 不変量と composite 表現の HOMFLY 不変量の両方が関わるはずだと書いている。\[30\]\[31\]

---

## 問題 B：Ĩ で不変な状態の数 E の決定

### B1. 向き付けられた M2（被覆の環状の弦、商では境界を一つもつ Möbius の帯）の BPS 状態への Ĩ の固有値をモジュライの幾何から決める文献

**結論：見つからなかった。** Piazzalunga–Uranga の t⁺ − t⁻ を、境界をもつ曲線（開弦）に広げた文献は確認できなかった。近い文献と状況は次のとおり。

- **Young 1603.05401（原典本文で確認）**：先述の「Ω^σ は Ω_Q から直接計算できず、ℤ₂ 同変な精密化が要る」という主張に加えて、本文（抜粋）に二つの記述がある。「通常の場合と違い、例によっては\[8\] orientifold DT 不変量がすでに無限個ある」こと、「loop quiver は Ω_Q が ℤ₂ 同変 DT 不変量を決めるという特別な性質をもつ」ことである。\[8\]
  - 当てはめ（解釈）：この研究の三つの観察と一対一に対応する。(i) E は指数（Ω_Q に当たる）から 2(1+t) を法としてしか決まらない。(ii) Hopf では決まったが結び目では決まらない。(iii) T(2,5) では指数の上限の中に解がなく、打ち消しで見えない状態が要る（orientifold 側の不変量が通常の不変量の台に収まらない）。Young の枠組みで決め手になるのは、自己双対表現のモジュライ M_e^σ 上の ℤ₂ 同変コホモロジー、すなわち「Ĩ の固定点のモジュライの幾何」である。これは B1 の「モジュライ空間上の対合」に当たる。
- **Franzen–Young 1605.06596**：σ 安定な自己双対表現のモジュライの Chow 群として、正の Betti 数が決まる（抄録）。解釈：E\[9\] の「正しい代表元」を選ぶ原理として、「E の各係数が固定点モジュライの Chow Betti 数（ε ごとの符号つき）になる」という正値性の条件が考えられる。
- **Ekholm–Shende の skein 値の曲線計数（サブエージェントの調査）**：Ekholm の 2021 年 3 月の講演スライド（Stanford NCSGS、Shende・Kucharski・Longhi・Georgieva・Ng との共同研究の報告）に、「一般の曲線には HOMFLY skein を、Lagrangian を固定する対合で不変な曲線には Kauffman skein を使う」「real な曲線は通常の GW 理論と同様に数えられる」とある（二次情報：講演スライドの抜粋）。しかし、Kauffman\[32\] skein 値の曲線計数の論文、またはそれによる Mariño 0904.1088 の整数性予想の証明は、arXiv では見つからなかった。サブエージェントは Shende の arXiv 著者リストに該当する題名がないことを確認したが、Ekholm と Georgieva の著者リストはレート制限で取得できなかった。したがって「存在しない」ではなく「見つからなかった」である。既存の skein 値の論文はすべて HOMFLYPT である。
  - arXiv:1901.08027 Ekholm–Shende "Skeins on Branes"：Theorem 1.2、証明は §8。Ooguri–Vafa の予言の証明で、v6 は 2026-04-24。\[33\]\[34\]
  - arXiv:2406.00890 Ekholm–Shende "Counting bare curves"：v3 で §10 が加わり、多重被覆を含む skein 値の計数を扱う（Corollary 10.1）。\[35\]\[36\]\[37\]
  - 題名・抜粋だけ見たもの：arXiv:2012.15366、2101.00619、2212.05835、2412.15454（Ekholm–Longhi–Shende）、2510.19041（Ekholm–Longhi–Park–Shende）、2407.09836、2512.07762（Hu–Shende）。
- **Hayashi–Piazzalunga–Uranga 1508.06644**：flop で振幅の符号が規則的に変わる（抄録）。解釈：固有値\[1\] ε の大域的な符号の選び方（向きのデータ）が壁越えで変わる可能性があることを示している。

### B2. Ĩ が二つの境界 L と I(L) を入れ替えるときの捻ったトレースで「差 z_q が和 z_q⁺ になる」こと

**結論：見つからなかった。** 近い文献と、この研究の式の読み方を挙げる。

- 最も近い文献：Georgieva–Ionel 1812.02505 の Klein TQFT には、「向きを反転する tube」と「クロスキャップ元 U」が生成元として入っている（抄録・抜粋）。境界を入れ替える対合を\[5\] TQFT の射として扱う枠組みはあるが、refined の境界のゼロモードの捻ったトレースは扱っていない。
- この研究の式の読み方（解釈、計算は自明な恒等式）：z_q · z_q⁺ = (q^{1/2} − q^{−1/2})(q^{1/2} + q^{−1/2}) = q − q^{−1} = z_{q²} である。したがって z_q⁺ = z_{q²}/z_q となり、E_22 = −z_q (z_q⁺)²(ν̃ − ν̃^{−1})(t q^{−1/2}ν̃ + q^{1/2}ν̃^{−1}) は E_22 = −(z_{q²})² z_q^{−1}(ν̃ − ν̃^{−1})(t q^{−1/2}ν̃ + q^{1/2}ν̃^{−1}) と書ける。入れ替えで不変な二つの境界のゼロモードの組の捻ったトレースは、一般に Tr_{V⊗V}(swap · g) = Tr_V(g²) という Adams 型の規則に従う。差が和に替わる現象は、「境界因子 z_q の 2 個組が q ↦ q² の境界因子 z_{q²} 1 個に置き換わる」と読める。ただし Adams 型の E = F·(ψ²[g_[1]] + 1) は依頼書で既に失敗とされているので、ここでは境界因子だけに限った読み方として記すにとどめる。

### B3. 第二量子化の積表示の M（向き付けられた、Ĩ 不変）と M̃（凍結）を独立に決める原理

**結論：開弦の refined 積表示で M と M̃ を独立に決める原理を論じた文献は見つからなかった。** 原理の候補となる文献は次のとおり。

- Young 1603.05401 の Conjecture A（Conjecture 3.7）：H_Q(e) ⊠ W^prim_{Q,e} → M_Q の分解。\[8\]解釈：積表示の二つの因子の型（Ĩ で入れ替わる対の塔 1/(1−t) と、z₂ = 0 に凍結した状態）は、CoHA が作用して生成する部分（H_Q(e)、対の塔）と原始生成元 W^prim（凍結状態と不変状態）に対応する。M_Q が自由加群でないこと、すなわち H_Q(e) が e によって変わることが、E → E + 2(1+t)X、C → C − X という不定性の代数的な起源にあたる。この不定性は、生成元の取り方（どの状態を「生成元」とみなし、どれを「対の作用の像」とみなすか）が指数だけでは決まらないことに対応する。
- Franzen–Young 1605.06596：Chow 群による同定と、壁越えのない結果（Thm. 4.3）。保護された量として「固定点モジュライの Chow Betti 数」を使う原理を与える。
- Bu 2503.20667：orthosymplectic DT の壁越え公式。refined の壁越えで M と M̃ を分ける可能性を与える一般論である（抄録のみ）。\[10\]
- Krefl 1001.5031（内容未確認、二次情報）。

### B4. composite（annulus の skein）の refined 不変量を brane–image の入れ替えで分解し、固有空間ごとの不変量を定義した文献。環状 HOMFLY ホモロジーでのトーラス結び目の composite 彩色の計算例

**結論：どちらも見つからなかった。**

- Queffelec–Rose–Sartori, arXiv:1802.04131, "Annular evaluation and link homology"：他論文による要約（二次情報）では、annulus の skein を環状 Khovanov–Rozansky ホモロジーで圏化し、三重次数の KhR ホモロジーが underived horizontal trace を経由することを示したとされる。トーラス結び目の\[38\] composite（adjoint）彩色の具体的な計算例は、今回の検索では見えなかった。
- Queffelec–Rose, arXiv:1506.08188, "Sutured annular Khovanov–Rozansky homology"（Trans. Amer. Math. Soc. 370 (2018) 1285–1319 と二次情報）：環状 sl(n) ホモロジーの導入（抄録）。\[39\]\[40\]
- arXiv:2512.20597, "Doubly-graded exponential growth of colored torus knot homology"（著者未確認）：彩色トーラス結び目ホモロジーの計算ソフトウェアを参照しているが、参考文献の抜粋から見る限り対称彩色が中心で、composite は確認できなかった。
- 参考文献リストで題名だけ見たもの（番号未確認）：Gorsky ほか "Evaluations of annular Khovanov–Rozansky homology"、Huerfano–Khovanov "A category for the adjoint representation"（sl₂ の adjoint）。\[41\]\[42\]
- 解釈：Mariño 0904.1088 の予想では、向き付けられない不変量に composite 表現の HOMFLY、すなわち annulus の HOMFLY skein 全体が入る（サブエージェントも本文抜粋で確認）。Cherednik–Elliot\[29\]\[43\] の言うとおり composite の DAHA 超多項式には正値性がないので、正値性は「brane–image 入れ替えの固有空間ごと」に初めて現れると考えるのが自然である。しかし、その分解を定義した文献は見つからなかった。

### 補足：knots–quivers 対応の Kauffman 版

サブエージェントの調査では、knots–quivers 対応はすべて HOMFLY-PT/超多項式に限られており、彩色 Kauffman 多項式や SO 型の BPS 不変量の quiver 公式（refined・unrefined とも）は見つからなかった。Young の orientifold DT（1408.4888、1603.05401）を彩色 Kauffman 多項式と結びつけた文献も見つからなかった。確認された主な文献は次のとおり。
- arXiv:1707.04017 Kucharski–Reineke–Stošić–Sułkowski "Knots-quivers correspondence"（abs を確認）\[44\]\[45\]
- arXiv:1811.03110 Ekholm–Kucharski–Longhi
- arXiv:1910.06193、2105.11806、2405.10161、2408.01832、2505.05668、2505.02059、2412.10885

---

## 次に読むべき文献と、検証できる予言

### 問題 A

**次に読むべき文献（3 件）**
1. Rains, arXiv:math/0112035, "BC_n-symmetric polynomials"（§6–7 の lifted Koornwinder）と、Rains–Vazirani, arXiv:math/0606204 の、Macdonald 対合が lifted Koornwinder に作用する公式。σ（＝(q,t) ↔ (t,q)）が large N の変数 T とともにどう作用するかを特定し、Mimachi の dual Cauchy 核の右辺で何が対応するかを見る。
2. Young, arXiv:1603.05401（§3.3–3.4 と §4.2 の loop quiver）。motivic（refined）な orientifold 不変量の具体例と、Hodge の重み（スピン）の付き方を、この研究の c=1・c=2 の多重項と比べる。
3. Georgieva–Ionel, arXiv:1812.02505（§7–10、Klein TQFT とクロスキャップ元 U_ρ）。unrefined の局所 real GV 公式のクロスキャップ係数の符号と正値性を、この研究の c=2 部分の「−(𝐚²−𝐪²)(𝐭²−1) × 正」と比べる。

**検証できる予言（いずれも解釈に基づく）**
- (A-i) σ と転置の組：太字変数での σ は 𝐪 ↦ −(𝐪𝐭)^{−1}、𝐭 ↦ 𝐭、𝐚 ↦ 𝐚 である。A 型の類推からは、
  σ(ĝ^{c=1}_{[2]})(𝐚, 𝐪, 𝐭) = ± 𝐚^{k}𝐪^{l}𝐭^{m} · ĝ^{c=1}_{[1,1]}(ε𝐚^{±1}, 𝐪, 𝐭)（ε = ±1）
  のように、σ が R と R^T を対応づける関係が予想される。3 箱なら [3] ↔ [1,1,1]、[2,1] ↔ [2,1] である。[2,1] は自己転置なので、3 箱の中でも [2,1] だけは σ 対称になるはずだ、というのが鋭い予言になる。これが成り立てば、「多箱で σ 対称でない」のは同じ R の中で σ を見ていたためと説明できる。
- (A-ii) 普遍因子の (𝐚² − 𝐪²) の起源：𝐚² − 𝐪² = ν² − t = t^{−1}(U − t)(U + t) なので、この因子は U = t（N = 1、すなわち SO(2) ≅ U(1)）で消える。D_1 には根がなく、Macdonald 理論は t によらない。したがって
  f_R(q, t, U = t) と g_R(q, t, U = t) は、ν の置き換えを除けば t = q の値と一致する
  と予言できる。これは「refined だけの部分が (𝐚² − 𝐪²) で割り切れる」ことの説明になる。U = −t でも消えることの意味（SO(2) の別の分岐、あるいは N → −N 型の対応）も確かめる価値がある。
- (A-iii) Franzen–Young の正値性の型：c=2 部分で (𝐚²−𝐪²)(𝐭²−1) を割った後の正の多項式は、何らかの固定点モジュライの Poincaré（Chow）多項式であるはずである。したがって、𝐪 について回文的（Poincaré 双対性）であることが期待される。T(2,2m+1) の [1] で、𝐪 ↦ 𝐪^{−1}（適当な単項式の正規化を含む）の回文性を確かめることを勧める。

### 問題 B

**次に読むべき文献（3 件）**
1. Young, arXiv:1603.05401（Conjecture 3.7 と、§4.2 の loop quiver の場合に Ω_Q が ℤ₂ 同変 DT を決める箇所）。「E が指数から決まる場合」の代数的な特徴づけを得る。
2. Franzen–Young, arXiv:1605.06596（Thm. 4.3 と、Chow 群との同型）。E の代表元を選ぶ幾何学的な正値性の原理を得る。
3. Bu, arXiv:2503.20667（自己双対 quiver と潜在関数の例、壁越え公式）。結び目の quiver（KRSS）に自己双対構造を入れたときの不変量の計算法を得る。

**検証できる予言（いずれも解釈に基づく）**
- (B-i) loop quiver の場合の一意性：Young によれば、loop quiver では通常の DT が ℤ₂ 同変 DT を決める。knots–quivers 対応では unknot は一頂点（loop 型）の quiver に当たる。したがって、Hopf で成功した方法（固定セクターの指数を状態の和とみなす仮定と、t = −1・t = q の条件）は、unknot の [2]・[3]（p₂ 方向のセクター）でも E を一意に決めるはずである。失敗が初めて起きるのは、矢でつながった多頂点の quiver に当たる結び目（trefoil は KRSS で 3 頂点）からだと予言する。
- (B-ii) 境界の因子：入れ替えで不変な境界対ごとに z_q² が z_q · z_{q²}/z_q · … に置き換わる、という読み方に従えば、結び目の p₂ 方向のセクターでは、
  E_{[2]} + 2(1 + t)X ≡ 0 (mod z_q⁺)
  を満たす整数係数の X が存在するはずである（Hopf の (2,2) セクターでは (z_q⁺)² で割り切れることが確認済み）。これは X の選び方を絞る条件として、整数計画の前段で使える。
- (B-iii) 不定性の代数的な意味：E → E + 2(1+t)X、C → C − X の不定性は、Young の CoHM が自由加群でないこと（H_Q(e) が e に依存すること）に対応する。したがって、固定点モジュライの Chow Betti 数が非負であるという条件、すなわち E の各係数を ε ごとに分けて非負にする条件を課すと、X は一意に決まるか、少なくとも有限個に絞られるはずである。T(2,5) で解が指数の上限の外に出るという観察は、Young の「orientifold DT 不変量は無限個になりうる」という記述と整合する。

## 注意点

- 今回の検索は 18 回で予算が尽きたため、Ekholm と Georgieva の arXiv 著者リスト、Rains の lifted Koornwinder に対する Macdonald 対合の具体式、Queffelec–Rose–Sartori の本文は確認できていない。「見つからなかった」は「存在しない」を意味しない。
- Bu の orthosymplectic DT の自己双対構造（層の双対）が、Aganagic–Schaeffer の Ĩ（反正則な I と z₂ → −z₂ の組）に当たるかどうかは未確認である。Young の quiver の反変対合と Ĩ の対応も、この研究の側で確かめる必要がある。
- 予言 (A-i)〜(B-iii) はすべて類推に基づく解釈で、文献に書かれた主張ではない。

## 出典

1. [\[1508.06644\] Towards a gauge theory interpretation of the real topological string](https://arxiv.org/abs/1508.06644)
2. [Real topological string amplitudes](https://www.researchgate.net/publication/386964743_Real_topological_string_amplitudes)
3. [Real topological string amplitudes K.S. Naraina N. Piazzalungab,c A. Tanzinic](https://arxiv.org/pdf/1612.07544)
4. [\[1812.02505\] A Klein TQFT: the local Real Gromov-Witten theory of curves](https://arxiv.org/abs/1812.02505)
5. [arXiv:1812.02505v1 \[math.SG\] 6 Dec 2018](https://arxiv.org/pdf/1812.02505)
6. [Motivic Donaldson-Thomas invariants and Kac conjecture](https://www.researchgate.net/publication/50353625_Motivic_Donaldson-Thomas_invariants_and_Kac_conjecture)
7. [Representations of cohomological Hall algebras and Donaldson-Thomas theory with classical structure groups](https://arxiv.org/pdf/1603.05401)
8. [Representations of cohomological Hall algebras and Donaldson-Thomas theory with classical structure groups](https://arxiv.org/html/1603.05401)
9. [\[1605.06596\] Cohomological orientifold Donaldson-Thomas invariants as Chow groups](https://arxiv.org/abs/1605.06596)
10. [\[2503.20667\] Orthosymplectic Donaldson-Thomas theory](https://arxiv.org/abs/2503.20667)
11. [Orientifolds and the Refined Topological String](https://arxiv.org/pdf/1202.4456)
12. [ABCD of beta ensembles and topological strings](<https://link.springer.com/article/10.1007/JHEP11(2012)111>)
13. [ABCD of Beta Ensembles and Topological Strings](https://www.researchgate.net/publication/228333665_ABCD_of_Beta_Ensembles_and_Topological_Strings)
14. [Exact Chern-Simons / Topological String duality](https://arxiv.org/pdf/1506.03907)
15. [Topological vertex formalism with O5-plane](https://link.aps.org/doi/10.1103/PhysRevD.97.026011)
16. [Topological vertex formalism with O5-plane](https://arxiv.org/html/1709.01928)
17. [Topological vertex formalism with O5-plane](https://doi.org/10.1103/PhysRevD.97.026011)
18. [Refined topological vertex with ON-planes](https://arxiv.org/pdf/2201.12264)
19. [CROSSCAP RESURGENCE AND RIEMANN–HILBERT PROBLEMS FOR ORIENTIFOLDED CONIFOLDS](https://arxiv.org/pdf/2602.04539)
20. [On refined Chern-Simons / topological string duality for classical gauge groups](https://arxiv.org/pdf/2205.12832)
21. [arxiv.org](https://arxiv.org/pdf/q-alg/9710032v1)
22. [Refined Cauchy/Littlewood identities and six-vertex model partition functions: II. Proofs and new conjectures](https://arxiv.org/pdf/1405.7035)
23. [BC\_n-symmetric polynomials - NASA/ADS](https://ui.adsabs.harvard.edu/abs/2001math.....12035R/abstract)
24. [arXiv:math/0606204v1 \[math.RT\] 9 Jun 2006](https://arxiv.org/pdf/math/0606204v1)
25. [Bounded Littlewood identities](https://arxiv.org/pdf/1506.02755)
26. [\[1005.3861\] Kauffman Knot Invariant from SO(N) or Sp(N) Chern-Simons theory and the Potts Model](https://arxiv.org/abs/1005.3861)
27. [Kauffman knot invariant from or Chern-Simons theory and the Potts model](https://link.aps.org/doi/10.1103/PhysRevD.81.125026)
28. [Colored Kauffman Homology and Super-A-polynomials](https://arxiv.org/pdf/1310.2240)
29. [Matrix Factorizations and Kauffman Homology](https://www.researchgate.net/publication/2065294_Matrix_Factorizations_and_Kauffman_Homology)
30. [String theory and the Kauffman polynomial](https://arxiv.org/pdf/0904.1088)
31. [Archive ouverte UNIGE https://archive-ouverte.unige.ch Article scientifique](https://access.archive-ouverte.unige.ch/access/metadata/56824fb7-148e-4acf-b7d6-6d9135d1117f/download)
32. [Skein module curve counts and recursion Tobias Ekholm](http://math.stanford.edu/~ionel/ncsgs/NCSGSMar2021_TobiasEkholm.pdf)
33. [\[1901.08027\] Skeins on Branes](https://arxiv.org/abs/1901.08027)
34. [Skeins on Branes](https://arxiv.org/pdf/1901.08027)
35. [Counting bare curves](https://arxiv.org/html/2406.00890v1)
36. [Counting bare curves](https://arxiv.org/pdf/2406.00890)
37. [Vivek Shende's articles on arXiv](https://arxiv.org/a/shende_v_1.html)
38. [Sutured annular Khovanov-Rozansky homology](https://www.researchgate.net/publication/279311591_Sutured_annular_Khovanov-Rozansky_homology)
39. [Link splitting deformation of colored Khovanov--Rozansky homology](https://www.researchgate.net/publication/353354195_Link_splitting_deformation_of_colored_Khovanov--Rozansky_homology)
40. [\[1506.08188\] Sutured annular Khovanov-Rozansky homology](https://arxiv.org/abs/1506.08188)
41. [Computing colored Khovanov homology](https://arxiv.org/pdf/2505.03916)
42. [Annular Khovanov homology and knotted Schur-Weyl representations](https://arxiv.org/pdf/1505.04386)
43. [(PDF) String Theory and the Kauffman Polynomial](https://www.researchgate.net/publication/24169632_String_Theory_and_the_Kauffman_Polynomial)
44. [Knots-quivers correspondence](https://arxiv.org/html/1707.04017v1)
45. [\[1707.04017\] Knots-quivers correspondence](https://arxiv.org/abs/1707.04017)
