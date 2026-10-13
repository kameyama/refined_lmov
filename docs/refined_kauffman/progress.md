# 進捗記録（第 3 の研究：refined Kauffman の large N duality）

## 現在の状況（再開用、2026-10-14）
- 開始。ユーザーの指示に従い、(1) unrefined の Kauffman の large N duality の文献調査、(2) 2018 年の試行錯誤（Dropbox の Mathematica ノートブック 42 本）のミスの調査、の順に進める。
- ブランチ：claude/refined-kauffman-large-n-00d58e（worktree）。

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
