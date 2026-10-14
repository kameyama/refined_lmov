# 進捗記録（第 3 の研究：refined Kauffman の large N duality）

## 現在の状況（再開用、2026-10-14）
- K0（unrefined の構造と再現）と K1（2018 年の試行の監査）は済。結果は literature.md、conventions.md、past_attempt.md、research_plan.md。
- 次は K2（refined unknot から構造を決める）。停止点：K2 の後。ユーザーの確認待ち。
- ブランチ：claude/refined-kauffman-large-n-00d58e（worktree）。push と PR は未（外部送信の許可が必要）。
- テスト：tests/test_rk_unrefined.py（約 12 秒）。

## 2026-10-14
- ディレクトリ docs/refined_kauffman/ を作成（README.md、progress.md）。
- 過去のノートブック 42 本（2018-07-10〜10-02）を scripts/nb2txt.py でテキスト化して読解を開始。同梱の PDF は Cherednik 1111.6195（DAHA-Jones、B/C 型の hyperpolynomial）と Cherednik–Elliot 1503.01441（refined composite invariants of torus knots via DAHA）（arXiv API で題名を確認）。
- 文献の取得：justfile の papers に 0904.1088（Mariño）、hep-th/0405083・hep-th/0411227（Bouchard–Florea–Mariño）、hep-th/0505008（Borhade–Ramadevi）、1003.5282（Paul–Borhade–Ramadevi）、1003.2861（Stevan）、hep-th/0512298（Gukov–Walcher）、1202.4456（Aganagic–Schaeffer）、1111.6195、1503.01441、1310.2240（Nawata–Ramadevi–Zodinmawia）、1007.1656（Chen–Chen）、hep-th/0108064（Mariño–Vafa）を追加して `just fetch-refs` を実行（題名・内容は原典で確認してから literature.md に書く）。
- **unrefined の構造の整理**：literature.md §1（Mariño の g_R, h_R, M_{RS}, 予想、証拠、弦の解釈）、conventions.md（Mariño の規約表）。
- **refined の既知の構成要素**：literature.md §2（AS の SO(2N) refined CS と large N の自由エネルギー、GW の Kauffman 超多項式、Cherednik–Elliot の composite DAHA 超多項式、NRZ の対称表現の refined Kauffman と 4_1 の unrefined の表）。
- **2018 年の refined SO unknot の検算**【計算】（scripts/rk_so_unknot_check.py）：ノートブック check kauf の refined SO(N) unknot（nakami[λ,μ] = (−1)^{(|μ|−rank μ)} P_{λ/μ} を p_n → (a^{n/2}−a^{−n/2})/(t^{n/2}−t^{−n/2}) で評価し、μ を自己共役な分割で足したもの）を、AS の SO(2N) refined CS の unknot（D_N 型 Macdonald 多項式の主特殊化 M_λ(t^ρ)、a = t^{2N−1}）と比べた。
  - AS (5.5) の t = q^β の積と、一般 t の評価式 t^{−⟨λ,ρ⟩} Π_{α>0} Π_{i=0}^{⟨λ,α⟩−1} (1 − q^i t^{⟨ρ,α⟩+1})/(1 − q^i t^{⟨ρ,α⟩}) は N = 3, 4、β = 2, 3、λ = [1], [2], [1,1] で一致。
  - [1]：ノートブックと AS は一致（N = 5, 6 で記号的に 0）。
  - [2]：差 = (q − t)(q − t^{N−1})/((qt − 1)(q t^{N−1} − 1))。[1,1]：差 = −(q − t)(1 + t + … + t^{2N−3} の因子)/(q t^{2N−3} − 1)。いずれも t = q で消える。
  - したがって、2018 年の試行の「refined SO(N) unknot」は、AS の SO(2N) refined CS の unknot と 2 箱で異なる。試行の入力（AS の trefoil (6.15) など）が AS の理論なら、unknot の正規化が整合していなかったことになる（ノートブック全体の読解で確認する）。
- **unrefined の再現**【計算】（scripts/rk_unrefined_check.py、sympy で厳密）：
  - trefoil：colored HOMFLY は src/rcsknots/homfly.py の Rosso–Jones（KN 変数 a = ν², Q = t、sign = −1）、colored Kauffman [2] は NRZ 1310.2240 (6.11) の t = −1、[1,1] は対称性 F_{[1,1]}(ν,t) = F_{[2]}(−ν, t^{−1})、composite ([1],[1]) は Mariño (4.4)。\hat g_[1,1] は Mariño (4.8) と全係数一致。\hat g_[2] は z^9 の項だけ不一致：Mariño の (4.8) の z^9 の項「−11ν^11 + 78ν^9 − 79ν^7 + 12ν^5」は z^8 の項の複写（奇数 z の項は偶数 ν 冪のはず）で、計算値は z^9(11ν^12 − 78ν^10 + 78ν^8 − 11ν^6)。**(4.8) の z^9 は原論文の誤植と判断**（他の全項は一致）。
  - Mariño の本文の P(3_1) = 2ν² − ν⁴ − z²ν² は量子群の標準値 2ν² − ν⁴ + z²ν²（= ν²(t² + t^{−2}) − ν⁴）と z² の符号が違う。(4.8) の再現には標準値を使った（Mariño の F(3_1) は NRZ の r = 1 と一致）。本文の P の符号は誤植の可能性（未確認）。
  - 4_1：\hat g_[2], \hat g_[1,1] は NRZ (C.10) と全係数一致。\hat h_[2], \hat h_[1,1] の z 展開から読んだ N^{c=0} の表は NRZ の表と完全一致（NRZ (C.9) の表示式は計算値と合わず、表示の誤植と判断）。
  - 2018 年 7 月のノートブック（h_hat_fig8, g_hat fig8）はこの NRZ の計算をなぞったもので、数値も一致する。
- **refined unknot の極**【計算】（scripts/rk_refined_unknot_poles.py）：AS の SO(2N) refined unknot（D_N Macdonald の主特殊化）は [2] で分母 (1 − q t^{N−1})、[1,1] で (1 − q t^{2N−3})、[3] で (1 − q² t^{2N−3}) などの **N に依存する分母**をもち、a = t^{2N−1} の Laurent 多項式にならない（t = q では約分される）。U(N) の composite ([1],[1]) = (2,1^{N−2}) の A 型 Macdonald 次元も分母 (1 − q t^{N−1})（a = t^N で (t − qa)）をもつ。2018 年 9 月 30 日の composite for adjoint ノートブックが見つけた「真の refined composite unknot の分母 (t − qa)」はこれと一致【補助エージェントの読解＋本計算で確認】。
- **AS §7 の閉弦の構造から読める示唆**【導出】：AS の自由エネルギーの oriented 部分 Σ_d (1/d) Q^d q^{−d/2}/((q^{d/2} − q^{−d/2})(t^d − t^{−d})) は、分母が (t^d − t^{−d}) = ((t²)^{d/2} − (t²)^{−d/2}) なので「(q, t²) の refined conifold」の形をしている（½ は現れない）。½F_cov ＋ unoriented の形に書き直すと、偶数 d の引き算に (1 + q^d)/(1 + t^d) の因子が付く（単純な d → 2d の置換ではない）。2018 年の試行は Σ_{d 奇数} と ½ を unrefined のまま置いたので、この構造を含んでいない。

## 2026-10-14（K2、途中）
- **閉じた式**【計算】（src/rcsknots/kauffman_refined.py）：AS の SO(2N) refined unknot dim_λ と計量 g_λ（D_N、AS (B.4)）を U = t^N の有理式として実装（dim は N = 5, 6 の積と一致を確認）。U(M) の composite [λ,μ] の Macdonald 次元（対称正規化）を、長方形の部分の積を列の対で伸縮させて A = t^M の有理式にし、M = 8, 9 の直接計算と 13 例で一致。
  - dim_[2] の分母 (Uq − t)、dim_[1,1] の分母 (U²q − t³)；g_[2] は (Uq − t)²、g_[1,1] は (U²q − t³)² を分母にもつ。composite ([1],[1]) の分母 (Aq − t)、([2],[1]) は (Aq² − t)、([1],[1,1]) は (Aq − t²)。
  - 訂正（2026-10-14 夕）：[3] の分母は (Uq² − t)（q²t^{N−1} = 1）であって (1 − q²t^{2N−3}) ではない。[2,1] は (U²q² − t³)、[1,1,1] は (U²q − t⁵)、[4] は (Uq² − t)(Uq³ − t)、[2,2] は (Uq − t²)(U²q² − t³)(U²q³ − t³)。
- **極の打ち消しの不可能性**【導出】：Z_G = Σ_λ dim_λ w_λ B_λ(v) で w_λ が U に依らなければ、2 次の [2] 成分の極 (Uq − t) と [1,1] 成分の極 (U²q − t³) は、composite ([1],[1]) の極 (Aq − t) 一つ（A は一つの値）では同時に打ち消せない。w_λ = 1/g_λ でも、(Uq² − t²), (Uq − t³) が composite の 2 箱の項には現れない。したがって「D_N の次元 × U に依らない重み」と「U(M) の composite 次元」の組み合わせでは、unknot の 2 箱で Laurent 多項式にならない。
- **積の形の検査**【計算】（scripts/rk_k2_littlewood.py、有理点 2 点）：w ∈ {1, b^A_λ, (q/t)^{|λ|/2} b^A_λ, 1/g_λ, b^A_λ/g_λ}、基底 ∈ {Schur, P_λ(v;q,t), P_{λ^T}(−v;t,q)} の 15 通りすべてで、log Z_G の 3 次の p_1³, p_1p_2 の係数が 0 でない（積の形でない）。t = q では w = 1・Schur で 0（Littlewood の恒等式）になることを確認（検査の妥当性）。
- **文献**：Rains–Warnaar 1506.02755（Memoirs AMS 270 (2021) 1317）と Lee–Rains–Warnaar 2007.03174（SIGMA 16 (2020) 142）を取得。
  - (q,t)-Littlewood 恒等式（Macdonald の 4 つと Kawanaka の 1 つ）、有界版（長方形の D_n, B_n, C_n Macdonald 多項式を A 型 Macdonald で展開）。
  - 辞書：P^{(D_n,D_n)}_λ(x;q,t) = P^{(B_n,B_n)}_λ(x;q,t,1)（ℓ(λ) < n）= Koornwinder K_λ(x;q,t; −1, −q^{1/2}, 1, q^{1/2})。AS の D_N Macdonald はこれ。
  - LRW の GL(2n) ⊃ O(2n) 型の (q,t)-分岐則（長方形の λ）は K_μ(x;q,t; 1, −1, t^{1/2}, −t^{1/2}) を使う：**AS の D_N（±1, ±q^{1/2}）と境界パラメータの q ↔ t が違う**（t = q で一致）。
  - 普遍指標の (q,t)-版は lifted Koornwinder 多項式 \tilde K_λ(x;q,t,T;t̄)（Rains 2005）。T = t^n で K_λ に戻る。T に極をもつ（Eq (2.?)：t_0t_1t_2t_3 = q^{1−â(s)} t^{l̂(s)+1} T^{−2}）。Cauchy 恒等式は双対型（Σ (−1)^{|λ|} \hat K_λ(x;t,q,T) \tilde K_{λ'}(y;q,t,T) = Π(1 − x_iy_j)）。
- 次：Koornwinder 多項式の主特殊化（評価式）を両パラメータ（±1,±q^{1/2} と ±1,±t^{1/2}）で求め、AS の D_N の式を再現することを確認した上で、極の構造を composite 側と比べる。
- **Koornwinder の主特殊化**【計算】（kauffman_refined.koornwinder_eval、van Diejen–Sahi の式）：パラメータ (1, −1, q^{1/2}, −q^{1/2})・評価点 (1, t, …, t^{n−1}) で AS の D_N の unknot と 6 つの λ で一致（辞書 P^{(D_n,D_n)} = K(±1, ±q^{1/2}) と実装の両方を確認）。LRW 型 (1, −1, t^{1/2}, −t^{1/2}) では [1^k] は極を持たず、[2] は qQ = 1、[3] は q³Q = 1、[2,1] は qQ = t に極（Q = U²/t）。
- **2 箱の極の打ち消し**【計算、scripts/rk_k2_poles.py】：composite ([1],[1]) の寄与は v の p_1² 方向だけにあるので、打ち消しには Kauffman 側の [2] と [1,1] が同じ極を残数の比 f^{(1,1)}_{(1),(1)} でもつ必要がある。D_N では極の位置が違い（U = t/q と U² = t³/q）、LRW 型では [1,1] に極が無い。重み 1/g_λ（D_N の計量）と A = U/t では (Uq − t²) が共有されるが、(Uq² − t²), (Uq − t³) が残り、これらは composite の 3 箱の極の位置に一致する（次数が合わない）。**結論：AS の D_N 理論でも LRW 型でも、U(M) の composite の Macdonald 次元を v の積の基底で組み合わせる限り、unknot の 2 箱で Laurent 多項式にならない**【導出＋計算】。
- K2 の停止点に到達（研究計画 §2）。

## 2026-10-14（夕）deep research の報告の読解
- ユーザーが deep research で得た文献調査の報告を docs/refined_kauffman/literature_deep_research_2026-10-14.md に保存（外部の要約、未検証）。要点：refined かつ開弦の orientifold の整数性を定式化した文献は無い（AS 自身が第一量子化は未定義と明言）。U(N) では KN（1703.05408）と Cheng–Sułkowski（2104.00713）が分母を Macdonald 基底・計量・refined holonomy に吸収する方式を確立。Avetisyan–Mkrtchyan（2205.12832、閉弦）は交差帽の数ごとに Kähler 変数の量子シフトが違う「三重化」を提案し、Klein 瓶の項の分母 (t^n + t^{−n}) を与える。Cherednik–Danilenko 1509.08351 は B/C/D の DAHA-hyperpolynomial の存在と対称性を未解決としている。Koornwinder 型（C^∨C_n）の refined CS は文献に無い。
- **D_N の unknot の N 依存の分母は Koornwinder の評価式の C⁺ 因子（lifted Koornwinder の T の極）に由来する**【計算】：C⁺_λ(t₀t₁t₂t₃ t^{2n−2}/q) = Π_{(i,j)∈λ}(1 − q^{λ_i+j−1} t^{2n−i−λ'_j})（t₀t₁t₂t₃ = q）の因子のうち分子 (t^n, t₀t_r t^{n−1}; q,t)_λ で約分されないものが、[1]〜[4], [2,1], [1,1,1], [2,2] の分母と一致する（上の訂正を含む）。報告の「Q3 の推測」は [3] の分母の誤記（本セッションの以前の記録に由来）を除いて正しい。
- 報告の推奨のうち、本研究で採用する順：(1) 整数性の判定基準を AS (3.13) 型の第二量子化（積表示）に移す、(2) Kauffman 側の母関数の基底を lifted Koornwinder（T を Kähler 変数）にして、Rains の lifted Cauchy 恒等式 Σ_λ (−1)^{|λ|} \hat K_λ(x;t,q,T) \tilde K_{λ'}(y;q,t,T) = Π(1 − x_i y_j) の上で組み直す、(3) 交差帽ごとの Kähler 変数のずれ（Avetisyan–Mkrtchyan）を開弦に持ち込む。いずれも【仮説】として研究計画 K2' に追加する。

## 2026-10-14（夜）K2'
- Koornwinder 多項式を q-差分作用素の固有関数として実装（src/rcsknots/koornwinder.py、n = 2, 3）。van Diejen–Sahi の評価式と一致、Mimachi の双対 Cauchy 恒等式を (n,m) ≤ (3,2) で確認【計算】。tests/test_rk_koornwinder.py。
- 双対 Cauchy 型の refined orientifold OV 核の定式化：AS 型（ホロノミー側）と LRW 型（probe 側）のパラメータが両側に現れること、unknot が第二量子化の積表示を自動的に満たすこと、Hopf 絡み目の明示式（scripts/rk_k2p_hopf.py、n = 2, 3 で直接の和と一致）。y^n z^n Z は n = 2, 3 で既約（t = q でも既約）。report.md §4–5 に整理。
- 停止点：K2' の判定（unknot の 2 箱で整合する定式化）は「双対 Cauchy 型なら積表示として整合、Laurent 多項式性の意味では未定義」。次の方向は report.md §5。

## 2026-10-14（深夜）K2'(a)：U を形式変数とした Hopf 絡み目
- 訂正：Askey–Wilson 漸化式（src/rcsknots/askey_wilson.py）の monic 化が誤っていた（KS (3.1.4) の漸化式は p̃_n = p_n/(a^{−n}(ab,ac,ad;Q)_n) に対するもの）。修正後、固有関数として解いた Koornwinder と n ≤ 4 で一致。昨夜の report.md §4 の「係数は t^n の有理式にならない」はこの誤りによるもので、撤回（§4 に訂正を記した）。
- probe 側の上位係数 c_{2j}(n) は U_n = t^n の有理式（有理補間、次数 (2j,2j)）。閉じた形（report.md §4.1）を 4 つの有理点で j ≤ 5、n ≤ 30 で確認【計算】。記号的な証明は未（【予想】扱い）。
- Hopf 絡み目（probe 各 1 個）を virtual 正規化 Z~ = y^N z^N Z で Q(q,t,U) の上で計算（src/rcsknots/kauffman_hopf.py、scripts/rk_k2p_hopf_U.py、a,b ≤ 4、約 1 分）。U = t^N（N = 3..8）の直接計算と一致【計算】。
- 結果【計算】：f_{10} = −1 − (ν−ν^{−1})/(t^{1/2}−t^{−1/2})、f_{a0} = 0 (a ≥ 2)；f_{11} = −f_{10} z_q(ν~−ν~^{−1})、ν~ = U q^{1/2}/t（t = q で Mariño の \hat g_{[1],[1]} = z(ν−ν^{−1}) と unknot 因子・符号を除いて一致）；f_{1b} = f_{11} χ_{b−1}(U/t)；全 f_{ab} は (1−q)/(1−t)（(偶,偶) は /(1−t²)）× 整数係数 Laurent 多項式。t = q では分母なし。
- テスト：tests/test_rk_askey_wilson.py（3 件）、tests/test_rk_hopf.py（2 件、A = B = 2、N = 3, 4 との照合と f_{10}, f_{11}, f_{12} の式）。
- 停止点：K2'(a) の結果（Hopf 絡み目で整数性が成り立つ）。次の候補：(i) probe 2 個（2 変数 Koornwinder、Mimachi の (n,m) = (N,2)）で 2 箱の表現を分離し Mariño の [2],[1,1] と比較、(ii) trefoil（AS の T 行列）、(iii) 閉じた形 c_{2j}(n) の証明（Rains の二項公式）。
- ユーザーの指示（2026-10-14）：(i) → (ii) の順で進める。

## 2026-10-14（深夜 2）K2'(b)：成分 1 に probe 無限個
- 訂正：K2'(a) の報告で「t = q の 1 変数 probe は Sp(2) 指標で混ざる」と書いたのは誤り。SO(2N) では Λ^k V（k < N）が既約なので、双対 Cauchy から K_{(m)}(y;q,q;t̄) = y^m + y^{−m}（D_1 指標）になり、virtual 正規化では y^k。Askey–Wilson の直接計算（K_{(2)} = y² + y^{−2}）と閉じた形（t = q で c_{2j} = 0）とも整合。report.md §4.1 を訂正。
- S 行列の対称性で成分 1 の和を縮約すると、probe の数によらず 1 変数の \hat K_{(k)}(z) と冪和だけで Z~(P;z) が書ける（多変数 Koornwinder 不要）。|R| ≤ 4、b ≤ 2 で f_{R,b} を Schur 基底と Macdonald 基底で計算（約 20 秒）。
- 結果【計算】：Schur 基底では全 f_{R,b} が (1−q)/(1−t)（b=2 かつ |R| 偶数は /(1−t²)）× 整数係数 Laurent 多項式。Macdonald 基底 P_R(y;t,q) では (1−qt), (1−q²t), (1−qt²) などの分母が出て整数性が破れる。f_{(a),b} = 1 変数の f_{ab}、f_{R,0} = 0（|R| ≥ 2）。report.md §4.2。
- テスト：tests/test_rk_hopf.py::test_multi_probe_consistency（A = B = 2、約 3 秒）。
- 次：(ii) trefoil。AS 1202.4456 §6 の (6.15)（T(2,3)、[1]、unknot で正規化）と T(2,2m+1) の一般式を使う。t = q で AS の値は Mariño の F_{3_1}(z,ν) の鏡像（ν → 1/ν、z → −z）× ν^{−6}（framing 6 単位）に一致【計算】。

- 未解決：f_{ab} の (偶,偶) の分母 (1−t²) の意味（refined LMOV の分母の構造との対応）。正値性の判定基準（refined BPS の抽出式）は未定義。

## 2026-10-14（深夜 3）K3'：trefoil の [1]
- AS 1202.4456 §6 の T(2,2m+1) の一般式を実装（src/rcsknots/kauffman_trefoil.py）。m = 1 で (6.15) に一致。t = q で Mariño の F_{3_1} の鏡像 × ν^{−6}【計算】。
- 古典的 Rosso–Jones（D_N 指標）は t = q で AS の値に係数 1 で一致するが、素朴な refined Rosso–Jones（Koornwinder 展開、T_μ = q^{(μ,μ)/2}t^{(μ,ρ)}）は AS と単項式倍でも一致しない【計算、n = 2, 3】。p_2 の Koornwinder 展開：K_{(1,1)} の係数は A 型と同じ (1+q)(1−t)/(1−qt)、K_{(1)} の係数 0、定数項は N 依存（n = 2 で −a11²）。
- refined の h/g 分解（仮説）：g^{ref}_[1] = W_SO − ½h^{ref}（½h^{ref} は DGR の超多項式の鏡像に 𝐚_U² = ν²(t/q)^{1/2} を入れたもの × (ν−ν^{−1})/(𝐪−𝐪^{−1})）。t = 1 の極の打ち消しから (6.15) に掛ける単項式は r^{−5} s^{d} U^{−6}（d 奇数）に限られ、d = 11（T_[1]^{−6}(q/t)^{1/2}）で最短の 11 項。整数係数。t = q で (ν²−ν⁴)(𝐪²+1+𝐪^{−2}) + ν⁶（z の偶数冪のみ、N^{c=2}_{[1]} = 0 と整合）。report.md §6。
- テスト：tests/test_rk_trefoil.py（3 件、約 1 秒）。
- 停止点：trefoil [1] の整数性は仮説の単項式のもとで成立。未解決：余分な (q/t)^{1/2} の由来、refined Rosso–Jones が D 型で AS と合わない原因、1^k（k ≥ 2）の refined SO 不変量の計算法（SL(2,Z)、有限レベル）。

## 2026-10-14（深夜 4）K4'：refined の交差帽展開
- ユーザーの指示：refined の交差帽展開（c = 1, 2）を unknot・Hopf で定義し、正値性と DGR 変数での正整数性を見る項目を計画に追加（research_plan.md K4'）。効果の高い順に (1) unknot・Hopf [1],[1]、(2) trefoil [1]、(3) Hopf (a),(b)、(4) 結び目の 1^k。
- 被覆（composite 表現の U(N) refined Hopf）を実装（kauffman_cover.py）。Pieri 係数の形を数値で特定、記号版（A = ν²）は有限 N と一致。unknot・[1],[1] の解析式と一致【計算】。
- Hopf の g_{ab} = f_{ab} − ½f^{cov}_{ab} は a,b ≤ 4 で ν の Laurent 多項式（有理補間）【計算】。
- 記号計算：(1,b) は g_{1b} = z_q(ν~−ν~^{−1})χ_{b−1}(ν~/q^{1/2})、t は ν~ 経由のみ（refined 交差帽 = unrefined with ν → ν~）。(2,3),(3,3) などは Laurent・整数だが t 依存。(2,2),(2,4),(4,4) は 2(1+t), 4(1+t)(1+t²) の分母が残る（AS §7 の (t^{d/2}+t^{−d/2}) の開弦版）。閉弦の (1+q^d)/(1+t^d) 規則では解消しない。report.md §7。
- 交差帽展開（G = z_q z_t の冪 + 奇数因子 X）を unknot・Hopf [1],[1]・trefoil [1] に適用（scripts/rk_k4_crosscap.py）：unknot は c=1 の −1、Hopf [1],[1] は X = z_q の c=2（±1）、trefoil の奇数部分は G 環に入らない（仮説依存）。
- テスト：tests/test_rk_cover.py（5 件、約 6 秒）。
- 未解決：refined の偶数 d の引き算の正しい形（(2,2) の留数を記録）、trefoil の単項式、結び目の 1^k。
