# 主張の格付け表（2026-10-09、論文執筆の準備：docs/paper/outline.md §3 の 1）

論文（outline.md の章立て）に載せる主張を、すべて次のどれかに割り当て、根拠と前提を対応させる。

**種別**
- **定理 / 命題 / 補題（Thm / Prop / Lem）**：【導出】の主張に限る。前提に予想式や仮説を含むときは、命題の文の中に明記する。
- **計算結果（CR）**：【計算】の主張。範囲を必ず書く。
- **予想（Conj）**：論文として提出する予想。
- **仮説（Hyp）**：文献から取れず、検証で裏付けた規則。
- **既知（K）**：再現・引用だけで、成果として書かない。
- **注意（Rem）**：観察・解釈。

列の意味
- 前提：仮説 H_sat（refined 衛星公式）、予想式（KNTZ 1902.02275 の refined の式。GNSSS 1512.07883 の予想式と一致。NRZS 1209.1409、Stošić 2402.03066 ＋ KRSS 1707.04017）、仮定 S_R（KN の U(1)_R の自明な作用）など。
- 根拠：導出のある文書、データ・ログ、機能テスト（`results` は保存済みの大きな結果を判定し直すテスト）。
- 新規性：novelty_report.md の判断（「新規」は「新規と思われる」の意味）。

---

## 1. 枠組み：framed knot への拡張（outline §3）

| ID | 主張 | 種別 | 前提 | 根拠 | 新規性 |
|---|---|---|---|---|---|
| F1 | 定義：rCS^{(x)}_λ = framing[λ]^x rCS_λ、framing[λ] = (−(q/t)^{1/2})^{\|λ\|}q^{n(λ')}t^{−n(λ)}。これから \hat f^{(x)}_ρ、\hat𝐍^{(x)}、\hat N^{(x)} を定める | 定義 | — | conventions.md §17、scripts/framing_check.py、tests/test_framing.py::test_framing_factor_matches_notebook_definition, test_framed_is_monomial_rescaling | KN の \hat f への拡張は新規（§10）。因子の形は既知（F3） |
| F2 | t = q で framing[λ]^x = (−1)^{ℓx}q^{xκ_λ/2}。Mariño–Vafa の framing に x = p で一致する（変数の向きも含む） | 命題（自明な導出） | — | framing_large_n.md §2.4、test_t_eq_q_limit_is_marino_vafa_framing, test_t_eq_q_unknot_matches_marino_vafa_variables | 照合（新規性は主張しない） |
| F3 | framing[λ] は Taki の refined vertex の framing の因子 f_λ(t,q) の逆数にちょうど一致する。IKV の形は (t/q)^{\|ν\|/2} が欠けている | 注意（導出による照合） | — | framing_large_n.md §2.2（Taki 0710.1776 の式 (modification)、原典で確認） | 因子は既知（Taki、Cheng–Sułkowski）。IKV の誤りは Taki 自身が修正済み。成果として書かない |
| F4 | AS の refined T 行列の固有値の比は、SU(N) の部分を除いて (a^{1/2}(q/t)^{1/2})^{\|λ\|}·q^{n(λ')}t^{−n(λ)} | 補題 | — | framing_large_n.md §2.1（記号計算 N = 3..5、\|λ\| ≤ 3） | AS の式からの直接の計算。新規性は主張しない |
| F5 | 箱あたりの Adams 作用と両立する単項式（a^{1/2}, (q/t)^{1/2}）は、\hat f_ρ を c^{\|ρ\|} 倍するだけで、β・2J_r の一様なずれと全体の符号しか変えない | 補題 | — | framing_large_n.md §2.1、§2.5（三葉結び目の 3 箱 x = ±1, ±3 で、(q/t)^{1/2} の有無の判定が項数まで一致） | 標準的（MV の λ^{pℓ/2} の吸収と同じ）。証明は短く書く |
| F6 | 太字の変数では framing[λ] = 𝐭^{\|λ\|+2n(λ')}𝐪^{κ_λ} で、符号のない単項式になる。framing は（彩色付き）HOMFLY-PT ホモロジーの次数のずらしとして働く | 命題（導出）＋ 注意（解釈） | — | conventions.md §17、framing_large_n.md §2.3 | 新規（観察として） |
| F7 | 𝐭^{−nx}\hat f^{(x)}_ρ は x の指数和。2 箱の根は {𝐪²𝐭², 𝐪⁻², 1, −1}、3 箱は 7 根。x の偶奇で分けると根は係数 1 の単項式 | 命題 | — | conventions.md §17。test_two_box_framing_is_exponential_sum（三葉結び目、当てはめに使わない x で予測が一致）。3 箱は自明な結び目の x = 7..10, −1, −2 で確認 | 新規 |
| F8 | framing[λ](1/t, 1/q) = framing[λ^T](q, t)。転置対称性、σ 不変性（G 展開が閉じる条件）、q-brane と \bar t-brane の整合性は x によらない | 命題 | — | test_framing_factor_transposition。計算でも全ての判定で f^{\bar t} = M(q^{−1})\hat f（framing_check の検査） | 新規 |
| F9 | 2 箱の公式の分母の各既約因子の零点で、残る 2 入力の framing の因子は一致する。したがって Laurent 性は x によらない（framing 0 で Laurent ⇔ 全ての x で Laurent） | 命題 | — | test_two_box_denominator_loci_framing_factors_agree | 新規 |
| F10 | MV の符号 (−1)^{\|λ\|x} を外すと、x が奇数のとき、一般の入力で 2 箱の \hat f が Laurent でなくなる（零点 q = −1, t = −1）。符号は KN の形に必要 | 命題 ＋ 計算結果（三葉結び目 x = 1 で分母 (1+q)(1+t)） | — | test_marino_vafa_sign_is_needed_for_laurent、framing_large_n.md §2.5 | 新規（MV の「符号は整数性に不可欠」の refined 版） |
| F11 | G 展開の補題：G 展開が閉じ（各 𝐚 行が σ 不変で 𝐪 の指数が偶数）、展開の前が非負係数なら、\hat N も非負整数。逆は成り立たない（K5c） | 補題 | — | REPORT_knot_table.md §2、REPORT_phase4.md §5.1。test_lemma_chebyshev_expansion_is_positive、test_K5c_before_and_after_G_differ | 新規と思われる（構成要素は初等的） |
| F12 | 命題 1：入力が非負係数で DGR の対称性をもてば、\hat f_{[1]} の G 展開は整数で閉じ、前後とも非負。全ての x で成り立つ（framing[□]^x = 𝐭^x なので 2J_r が x ずれるだけ） | 命題 | DGR の対称性は Gorsky–Hogancamp–Mellit 2103.01212 の定理（ホモロジーを入力とする場合） | REPORT_knot_table.md §2。test_positive_symmetric_polynomials_give_positive_bps、test_fundamental_framing_is_t_shift | 新規（構成要素は既知、§6） |
| F13 | 系：HOMFLY-PT ホモロジーの Poincaré 多項式を入力とする限り、任意の結び目・任意の x で [1] は正 | 系 | F12 と同じ | 同上 | 新規 |
| F14 | [1] の 2J_r の偶奇が (g, β) ごとにそろう ⇔ 𝐚 の次数が連続する範囲ごとに e_𝐭 + e_𝐚/2 の偶奇が一定（parity knot の条件） | 命題 ＋ 計算結果（801 個中 9 個） | 𝐚 次数の台に隙間が無いこと（801 個すべてで成立）。CG の定義と同値（§7 の 4） | REPORT_knot_table.md §2、REPORT_phase4.md §7。test_knot_table_all_positive_and_parity_criterion | 条件自体は既知（Chandler–Gorsky 2209.13058）。BPS との同値性と、11 交点以下の一覧は新規 |
| F15 | 命題 2：\hat f_{[2]} が Laurent なら、非負の入力では dim H^{[2]} = (dim H^{[1]})² が必要（[1,1] も同様）。この条件は x によらない | 命題 | — | REPORT_knot_table.md §3.3、F9 | 導出は新規。9_42 の次元 401（GGS）と、9_42 が KN の形に乗らないこと（KN）は既知 |
| F16 | 全ての x の判定法：x の偶奇・向きごとにニュートン形を作り、係数の一様な符号と、起点と 0 の間の直接判定で、全ての x を有限個の判定に帰着する | 方法（補題としての正しさは導出） | — | scripts/framing_newton.py、test_trefoil_two_box_all_framings_newton | 方法。V4 のツイスト結び目の方法の転用 |
| F17 | (N2') の彩色依存部分はちょうど F（= pq）単位の framing：M_λ = (M_□(−(q/t)^{1/2})^{−F})^{\|λ\|}framing[λ]^F。基準の framing を偶数ずらすのは x のずれに吸収され、奇数ずらすと F10 により Laurent 性が壊れる | 命題 | — | framing_large_n.md §2.6、conventions.md §17 | 新規（正規化の意味づけ） |
| F18 | framed knot の refined large N 双対の物理的な説明：MV（平坦座標と framing）＋ KN（refined index、Schwinger、S_R）＋ AS の T 行列 ＋ Taki の framing の因子。抽出式（refined-GT1, GT2）は全ての x でそのまま適用できる | 注意（解釈） | 仮定：framing を変えても U(1)_R が保たれる。符号の幾何的な導出（V → −(q/t)^{1/2}V）はしていない | framing_large_n.md §3 | 新規（議論として） |
| F19 | **予想 A（framed refined LMOV の正値性）**：U(1)_R をもつ結び目（トーラス・代数的な iterated torus）と、ホモロジーを入力とする結び目で、全ての x、全ての ρ について、全体の符号を選べば \hat𝐍^{(x)}_ρ と \hat N^{(x)}_ρ が非負整数 | 予想 | 仮定 S_R（KN）＋ F18 | 証拠：§2〜§4 の計算結果（R1〜R14） | 新規の予想（KN の予想の framed への拡張） |

## 2. 入力（outline §4）

| ID | 主張 | 種別 | 前提 | 根拠 | 新規性 |
|---|---|---|---|---|---|
| I1 | KN のノートブックの 27 件（トーラス結び目 \|λ\| ≤ 3）を EHA（GN の構成）から全係数で再現し、BPS の表 13 個も再現した（規約の固定） | 既知の再現 | — | test_gn_reproduces_notebook、test_saved_tables_match_paper、test_kn_tables_reproduced | 既知。成果として書かない |
| I2 | ケーブルの refined 不変量：DAHA の定義を companion の段で分解した計算法が、CD14 の明示例 10 件（第 1 期の 8 件に K7, K8 を追加）と全係数一致する。素朴な Rosso–Jones の t 変形は誤り | 計算結果（方法の検証） | — | REPORT.md §2、conventions.md §7。test_matches_daha_example | 素朴な形の失敗は既知。正しい形は定義から従うので、新規性は弱い（novelty_report §5 (N4)）。「計算法として明示し検証した」と書く |
| I3 | ケーブルの t = q の極限が、ブレイドから計算した彩色 HOMFLY-PT と一致する | 計算結果 | — | test_cable_formula_matches_braid、test_t_equals_q_is_homfly | 照合 |
| I4 | 正規化 (N1)(N2')：[1] は t = q で HOMFLY-PT に一致し最低 𝐭 次数 0。彩色付きは [1] の単項式と framing の固有値の積 | 仮説（ケーブル）。仮説の中身は「t = q で決まらない彩色ごとの (q/t)^{c_λ} が framing の固有値に従う」の一点（§7 の 2） | — | conventions.md §5。トーラス結び目では確定（27 件）。ケーブルでは K5c で (N2) が破れ (N2') で回復（転置・Laurent 性）。test_normalization_N2prime | 新規（規則の提案） |
| I5 | M_□ = (−q/t)^{b−1}(a(q/t)^{1/2})^δ が保存済みの 22 結び目で成立する | 計算結果 | — | conventions.md §13。test_box_monomial_formula | 新規 |
| I6 | I5 の導出：H_sat とトーラス結び目の事実 (I)〜(IV) から | 定理（条件つき） | **仮説 H_sat**。(IV) の辞書は、Mellit の Thm (dyck paths)（定理）と、標準の正規化の Cor (daha) の形 𝒫 = −(a(qt)^{−1/2})^δ(1/(1−t))(Q·1)[a/q−1]/(a/q−1) に依存する。後者は Thm (daha) から d_± の表現の元として導出でき（2026-10-09、原典で確認）、標準の正規化での形は計算（12 例、総当たりで唯一）。d_± の元と BGLX の Q の対応の照合が残る | conventions.md §13.1。test_framing_operator_shifts_slope、test_bottom_coefficient_is_rational_catalan | 新規（条件つき）。GN の Catalan 予想が Mellit の定理の系であることは既知と考えられるので、成果にしない |
| I7 | 結び目表 801 個の [1] の入力（yui-kr の HOMFLY-PT ホモロジー）が、DGR の対称性・ブレイドの HOMFLY-PT・トーラス結び目の rCS と整合する | 計算結果（入力の検証） | — | test_dgr_symmetry、test_euler_characteristic_is_braid_homfly、test_torus_knots_match_rcs | 照合 |
| I8 | 2 箱の入力：KNTZ（1902.02275）の refined の式（ツイスト・ダブルツイスト）。GNSSS の予想式と全ての (m, n) で一致（各象限の 3×3 の格子で 144 件一致 ＋ 指数和の根の数）。t = q では Fang–Zhou（2610.00303）の定理。t_c 次数の変換は GGS の 𝐐 → 𝐪^R で、三葉結び目の EHA と較正 | 既知（入力の予想式） | **予想式** | conventions.md §15。scripts/kntz_vs_gnsss.py、logs/kntz_vs_gnsss.log、test_kntz_equals_gnsss_on_grids、test_trefoil_matches_eha、test_twist_knots_match_phase2_inputs | 入力は既知の予想式。KNTZ と GNSSS の両方を引用 |
| I9 | thick な結び目の t_c 次数：Stošić の quiver ＋ KRSS の (quad-C) の組み合わせが、8_19 で EHA の [2]〜[1,1,1] と全項一致する | 計算結果（方法の較正） | **予想式** | open_problems.md B4、logs/b4_819_quad.log | 部分的に新しい（予想式の検証例） |
| I10 | 厳密な計算：素数 1〜3 個の補間と CRT、補間に使っていない素数での検証、G 展開の高速化 | 方法 | — | tests/test_speedup.py（7 件） | 方法（付録） |

## 3. 検証 I：[1]（outline §5）

| ID | 主張 | 種別 | 前提 | 根拠 | 新規性 |
|---|---|---|---|---|---|
| R1 | 11 交点以下の全素結び目 801 個で、\hat f_{[1]} の G 展開が閉じ、前後とも非負（命題 1 の確認） | 計算結果 | 入力は yui-kr | data/knot_table_summary.md、logs/knot_table_t1.log、test_knot_table_all_positive_and_parity_criterion | 新規（表として） |

## 4. 検証 II：2 箱（outline §6）

| ID | 主張 | 種別 | 前提 | 根拠 | 新規性 |
|---|---|---|---|---|---|
| R2 | ツイスト結び目 K_p：全ての p ≠ 0 で、2 箱は G 展開の前後とも正（x = 0） | 定理（導出 ＋ 有限個の厳密な判定） | **予想式（KNTZ の refined の式＝GNSSS）**、t_c 次数 δ' = −S | REPORT_phase4.md §5.1。c_2 を使わない LP 証明書あり（data/v4/v4_ck5_cert_*）。test_v4_twist.py（7 件）、logs/v4_closed.log, v4_newton_*.log | 新規（\hat f の段階） |
| R3 | ツイスト結び目の 2 箱は、全ての p で sl_2 指標の正係数の和（単峰） | 命題（導出 ＋ 有限個の判定） | 予想式 | REPORT_phase4.md §5.1（c_0..c_2 は計算） | 新規。注意として載せる |
| R4 | ダブルツイスト K_{m,n}（GNSSS の pretzel K_{s,t}）：m, n が同符号なら全て正。異符号は min(\|m\|, n) ≤ 20 で証明、\|m\| ≤ 30 かつ n ≤ 40 の範囲で直接確認（\|m\|, n > 20 の 200 点。閉じた形での判定で、閉じた形は 54 点で直接計算と一致）。残りは 1 変数の 3 系列の非正値性に帰着（u ≤ 25 で成立） | 定理（同符号）＋ 定理（異符号の一部）＋ 計算結果 | **予想式（KNTZ の refined の式＝GNSSS）** | REPORT_phase4.md §5.2、REPORT_phase5.md §2。data/v4/ の証明書。test_v4_pretzel.py | 新規。**残りが未証明であることを定理の文に明記**（d_1 の LP は打ち切り、対角線は未着手） |
| R5 | thick な非トーラス結び目 10_139, 10_152 の 2 箱は G 展開の前後とも正（次元 441 = 21², 841 = 29²） | 計算結果 | **予想式（Stošić ＋ KRSS）** | RESULTS.md §3、logs/b4_s2_thick.log、data/b4/ | 新規（knots–quivers の正値性からは従わない、novelty_report §9） |
| R6 | 正値性は自動的でない：5_2 で t_c 次数の別の割り当て 12 通りは全て負の項をもつ | 計算結果 | 予想式 | REPORT_knot_table.md §4、logs/t3_nonunique.log | 新規。注意として載せる |
| R7 | 2 箱は、全ての整数 x で G 展開の前後とも正。対象は自明な結び目、T(2,3), T(2,5), T(2,7), T(3,4), 4_1, K1, K3 | 定理（導出 ＋ 有限個の厳密な判定） | 入力：EHA/DAHA、4_1 は KN のデータ（ホモロジー由来） | logs/framing_newton.log、F16、test_trefoil_two_box_all_framings_newton | 新規 |
| R8 | ツイスト結び目の (p, x)：p と x が同符号の向きでは全て正。逆符号は \|p\| ≤ 6、\|x\| ≤ 14 の直接計算で全て正（696 件） | 定理（同符号）＋ 計算結果（逆符号） | 予想式 | logs/framing_twist.log, framing_twist_direct.log | 新規。**逆符号は証明書なし**と明記 |
| R9 | framing を振った直接の判定：計 1821 件（(x, ρ) の組）で G 展開の前後とも負の項 0（9 件は \hat f = 0）。3 箱まで：自明な結び目と T(2,3) は \|x\| ≤ 6、T(2,5) と 4_1 は ≤ 4、T(2,7) と T(3,4) は ≤ 3、T(3,5), T(2,11), T(3,7), T(4,5) と K6a, K7, K8, K9 は ≤ 10、K3 は ≤ 2、K1 は ≤ 2 と ±6。K10 は 2 箱まで ≤ 10、3 箱は ≤ 2。ツイスト結び目 12 個は 2 箱まで ≤ 4 | 計算結果 | 入力は R7・R8 と同じ | logs/framing_batch.log（732 件）、data/framed/ と data/framed_summary.md（1083 件 ＋ K1 の 6 件、main の PR #6）、docs/paper/generated/framing_table.tex（scripts/paper_framing_table.py）、test_framing_big_matches_framing_check。以前の計算（未発表）の出力 9 件と全係数一致 | 新規 |

## 5. 検証 III：3 箱（outline §7）

| ID | 主張 | 種別 | 前提 | 根拠 | 新規性 |
|---|---|---|---|---|---|
| R10 | トーラス結び目 T(3,4), T(3,5), T(2,11), T(3,7), T(4,5) の \|ρ\| ≤ 3 は G 展開の前後とも正（x = 0）。T(2,3), T(2,5) は KN の表の再現 | 計算結果 | — | RESULTS.md §1、logs/u4_torus.log, u4_T3_4.log、data/lmov/ | KN の予想の確認（新しい表だが、主張の新規性は弱い） |
| R11 | 代数的な iterated torus 結び目 7 個（K1, K3, K6a, K7, K8, K9, K10）の \|ρ\| ≤ 3 の全 6 ρ で、G 展開の前後とも正（x = 0）。K7 は 3 本ケーブル、K10 は 2 段の反復（3 箱は約 183 万項、係数の最大 1.1×10¹⁹） | 計算結果 | 入力は DAHA の定義どおり。正規化は I4（仮説） | RESULTS.md §1、data/cables/, data/lmov/, data/bps/、logs/u4_cables*.log, v0_K10_*.log。test_u4_three_box_positivity（results）、test_u4_bps_saved、test_positive_before_and_after_G | 新規（K7, K10 の \|λ\| = 3 の値自体も文献に無い、novelty_report §10） |
| R12 | 入力の rCS_{[2,1]} 自体は負の係数をもつのに、BPS の段で正になる | 注意（計算での観察） | — | REPORT_phase3.md §3 | 観察 |
| R13 | 3 箱を x ≠ 0 で判定：R9 の範囲で全て正。最大は K10 の 3 箱 x = 2（G 展開の前に約 189 万項）。全ての x の証明書は未取得（6 根のニュートン形では見つからない） | 計算結果 | 4_1 の [2,1] は KN の候補 | R9 と同じ | 新規 |
| R14 | **予想 B**：代数的な iterated torus 結び目は、全ての ρ について KN の正値性を満たす（予想 A の x = 0、結び目のクラスを明示した版） | 予想 | KN の U(1)_R の議論は直接は当てはまらない（補空間が Seifert でない）。動機は CD14 の正値性予想（代数的な結び目・長方形の DAHA の値は正）と Connection Conjecture (ii)（無彩色では KhR と一致）。[1] は CD14 (ii) ＋ 命題 1 から従う。[2,1] は CD14 の範囲外（§7 の 5） | R10, R11 | 新規の予想（KN の予想のクラスの拡張） |
| R15 | 双曲結び目の [2,1]：構造的な条件では決まらない。解の差は「結び目によらない因子 × 自由な多項式」 | 命題（計算 ＋ 導出） | — | REPORT_phase3.md §2、REPORT_phase4.md §3 | 曖昧さ自体は既知（AMMM 1211.6375）。解の空間の構造の明示は新規 |
| R16 | 弱い判定：正値性を満たす [2,1] は 4_1（KN の候補）と 6_1（813 項の整数解）に存在し、5_2 は有理数解のみ（整数は未決）。σ = −1 は 3 結び目とも解なし | 計算結果 | 入力の [1]〜[3] は予想式（GNSSS） | REPORT_phase4.md §6、logs/v5_weak.log | 新規。**予想の検証ではなく「矛盾しない」ことの確認**と明記 |
| R17 | 族として：t = q ではダブルツイスト結び目の [2,1] の Morozov の式がケーブルの値と一致する。refined の長方形は KNTZ の式を実装して照合した。refined の [2,1] は、族の仮定を置いても一意に決まらない | 計算結果 ＋ 命題 | 族の仮定（Morozov の構造、KNTZ の正則な行） | REPORT_overview.md §6、tests/test_x0_double_twist.py, test_x1_refined.py, test_x2_family.py | t = q は既知の式の検証。refined の否定的な結果は新規。要点だけ書く |

## 6. 付録に回す主張

| ID | 主張 | 種別 | 根拠 | 扱い |
|---|---|---|---|---|
| A1 | 非代数的なケーブル（DAHA の値）では正値性が破れ、破れ方は 3 つの型に分かれる。DAHA の値はスペクトル系列の必要条件を満たさず、ホモロジーではない | 計算結果 ＋ 命題 | REPORT.md §5.2, §5.3、test_negativity_cause.py | 付録 B。CLAUDE.md・open_problems D5 により深追いしない |
| A2 | 非代数的なケーブルの [1] で、DGR 型の superpolynomial の候補からは前後とも正の refined LMOV が得られる | 計算結果 | REPORT.md §5.3、test_constructed_superpolynomial_constraints | 付録 B（候補であって定義ではない） |
| A3 | 文献の誤植と読み方：GS 1112.0030 の 9_42 の S²、1612.00422 の分母、1903.00259 の正規化、KNTZ の P*、KN 付録 B の f^t、BGLX (5.26) の内積 | 注意 | test_gs_942_s2_misprint、test_ft111_matches_notebook など | 付録 C。著者への連絡は別途判断 |
| A4 | f^q, f^{\bar t} の係数 | 計算結果（補足） | REPORT.md §4 | 付録 D（refined LMOV ではない） |
| A5 | 2J_r の偶奇の一致は、rCS の指数の偶奇の一様性から自動的に従う（ケーブル） | 命題 | REPORT.md §5.5、test_uniform_parity | 付録または脚注 |

---

## 7. 書く前に解決・明記が必要な点（2026-10-09 に検証済み）

1. **I6 の Mellit の正規化**【解決（影響の範囲を確定）】
   - 原典の確認（conventions.md §13.1 (IV)）：Cor (daha) の文面は、d_± の表現の元と読めば Thm (daha) と整合する【導出】。標準の正規化では [a/q−1]・符号 −1 の形だけが成り立つ【計算、総当たりで唯一】。
   - 影響の範囲：(IV) が決めるのは M_□ の (q/t) の冪 c だけである。(q/t)^{c|λ|} は Adams 作用と両立するので、\hat f_ρ を (q/t)^{c|ρ|} 倍し、2J_r を一様にずらすだけ（F5、conventions.md §5）。
     **したがって、(IV) の未了の部分（d_± の元と BGLX の Q の対応）は、正値性・整数性の結果に影響しない。** 影響するのは閉じた式 I5 の (q/t) の指数の導出（I6）だけで、I5 自体は 22 結び目で計算により成立している。
   - 書き方：I6 を「H_sat の下で、(I)〜(III) は導出、(IV) は Mellit の Thm (dyck paths) と、標準の正規化での Cor (daha) の形（計算で確認）の下で導出」とする。
2. **I4 の (N2') の仮説の範囲**【解決（範囲を訂正）】
   - 前の記述「仮説として残るのは M_□」は不正確だった。M_□ は [1] の t = q（(N1)、HOMFLY-PT との一致）と (N2) で決まる。(N2) の自由度 (q/t)^c は正値性に影響しない。
   - 彩色付きの単項式 ±a^α q^β t^γ のうち、t = q で決まるのは符号・α・β + γ である（彩色付きの (N1) を検査として確認：test_t_equals_q_is_homfly）。**t = q で決まらないのは、彩色ごとの (q/t)^{c_λ} の自由度だけ**。
   - (N2') はこれを「c_λ は |λ| に比例する部分（正値性に無関係）と、refined の framing の固有値 T_λ^F（F = pq）の分だけ」と定める。**仮説の中身はこの一点**。
   - 根拠は 3 つ：(i) H_sat からの導出（conventions.md §13.1 (A)）。(ii) K5c で、[2] を (q/t)^u 倍する 1 パラメータ族のうち、Laurent 性と転置対称性を満たすのは (N2') の u だけ（conventions.md §5）。(iii) 全ケーブルで (N2') の形が保存値と一致（test_normalization_N2prime、9 結び目）。
   - F17 により、この仮説は「彩色によって変わる部分は Seifert framing の F 単位の framing である」と言い換えられる。
3. **R4・R8 の未証明部分**【検証済み、1 件訂正】
   - R8：ログと一致する。p と x が同符号の 2 象限は、テンソル積のニュートン形と、起点と 0 の間の判定で全て正。逆符号は |p| ≤ 6、|x| ≤ 14 の直接判定 696 件で負の項 0（logs/framing_twist_direct.log）。
     ログの「符号 0（一様でない）」は、x を固定したときの p のニュートン形の係数が一様でないこと（証明法が使えない）を表す。負の項ではない。
   - R4：**範囲の記述に誤りがあった**。|s|, t > 20 で直接判定したのは 21 ≤ |s| ≤ 30・21 ≤ t ≤ 40 の 200 点だけである（logs/v4_pretzel_np2_check_a.log。「264 点」は、その中の 21..28 の正方形 64 点を重ねて数えていた）。
     正しくは「min(|s|, t) ≤ 20 で証明、|s| ≤ 30 かつ t ≤ 40 で確認」。RESULTS.md、REPORT_overview.md、REPORT_phase5.md、outline.md、R4 を訂正した。
     この直接判定は閉じた形を使っている。閉じた形は、作るのに使っていない 54 点で直接計算と全係数一致している。
4. **F14 の parity knot**【解決：同値】
   - Chandler–Gorsky（2209.13058、main.tex:334-336）の定義は「HHH(K) 全体で t/2 の偶奇が一定（Nakagane–Sano の次数）」で、同値な言い換えは「t_DGR − a/2 の偶奇が一定」。大域的な条件である。
   - 本研究の条件は「𝐚 の次数が連続する範囲ごとに e_𝐭 + e_𝐚/2 の偶奇が一定」。𝐚 の次数は偶数なので ±a/2 は偶奇が同じ。**𝐚 次数の台に隙間が無ければ、二つの条件は同じ**【導出】。
   - 801 結び目すべてで台に隙間が無く、該当する結び目の集合は一致した（9 個：3_1, 5_1, 7_1, 9_1, 11a_367, 8_19, 10_124, 10_125, 10_140）。CG の元の定義（NS の t/2）で数えても同じ【計算】。CG の本文には一覧が無い。
   - 書き方：「[1] の 2J_r の偶奇がそろう ⇔ parity knot（𝐚 次数の台に隙間が無い結び目で。11 交点以下の全素結び目がこれに当たる）」。
5. **F18・F19・R14 の仮定**【一部解決、一部は仮定として明記】
   - CS 側（解決）：AS（1105.5117 :646-648, :1236-1239）では、U(1)_R は M が半自由な U(1) 作用をもつ（Seifert）ことから存在する。framing の変化は refined CS の T 行列（K と K T^f の違い。AS :2030 のコメントアウトされた注に「framing の曖昧さ」と明記）として、同じ S, T の演算子の形式の中で扱われる。
     したがって、トーラス結び目では framing を変えても CS 側の U(1)_R は保たれる【原典 ＋ 解釈】。
   - 弦側（仮定のまま）：遷移の後の L_K の境界条件を変えたときに、BPS 状態への U(1)_R の作用が自明なままであることは導出できない。予想 A の前提として明記する。符号の幾何的な導出も無い（必要性は F10 で導出済み）。
   - **新たに見つかった点（R14）**：代数的な iterated torus 結び目の補空間は一般に Seifert ではない（graph manifold）。そのため、KN の U(1)_R による正値性の議論は x = 0 でも直接は当てはまらない。
     予想 B の動機は KN の物理的な議論ではなく、DAHA の値（定義どおりの refined 不変量）とホモロジーの関係として述べる。
   - **CD14 の原典での確認（2026-10-09、references/1408.4348/src/cable.tex）**：
     - 要旨（:55-70）：DAHA superpolynomial は「係数が非負の場合、reduced stable Khovanov–Rozansky 多項式と一致すると思われる」。代数的な結び目では、ORS の予想を一般の dominant weight に広げたものと一致すると期待する。
     - 正値性の予想（:1072-1076）：「一般の正値性予想は、代数的な結び目と**任意の長方形のヤング図形**について、DAHA superpolynomial の全ての係数が正であると主張する」。この正値性をもつ結び目を pseudo-algebraic と呼ぶ（:1076-1077、:2892-2894。代数的な結び目とその鏡像を、長方形について全て含むと予想）。
     - Connection Conjecture（Conj. CONCONJ、:1488-1550）：(i) t ↦ q, a ↦ −a で HOMFLY-PT に一致（任意の iterated knot と weight）。(ii) **無彩色**の pseudo-algebraic 結び目では、DAHA superpolynomial ＝ tilde 正規化した reduced stable KhR。(iii) ORS の Conjecture 2 と合わせて、平面曲線の特異点の Jacobian factor の Betti 数。
       (ii)(iii) を一般の weight で述べない理由は、彩色付きの KhR と P_alg が無いことだとしている（:236-241, :1086-1089）。標準の変数の対応（qtareli、:1494-1499）は conventions.md §1.2 と同じ。
     - 本研究との関係：
       - [1]：CD14 の (ii) の下で、代数的な iterated torus 結び目の DAHA の [1] は KhR の Poincaré 多項式になる。命題 1（F12）の入力の条件（非負係数と DGR の対称性）がホモロジーとして満たされ、[1] の正値性が従う。つまり**[1] については、予想 B は CD14 (ii) ＋ 命題 1 から従う**。
       - **(N1)(N2) と tilde 正規化の対応（2026-10-09 に確認）**：
         - CD14 の tilde 正規化は、KhR を最小の a_st の冪と q_st^•t_st^• で割り、a_st = 0 の部分が定数項 1 の非負係数の多項式になるようにしたもの（:1525-1531）。標準の変数（qtareli）は太字の変数と同じ。
         - 【計算】係数が非負の保存値 19 個（トーラス 9、ケーブル 10）で、本研究の太字の P は **𝐚^{2δ}𝐪^{−2δ} × (tilde 正規化)** になっている。最低 𝐚 の行の最低の項がちょうど 𝐚^{2δ}𝐪^{−2δ}（係数 1）。
           CD14 の明示例との比も、正のもの（K1, K3, K5b, K7, K8）で全て 𝐚^{2δ}𝐪^{−2δ}（𝐭 のずれなし、符号 +）。どの P も DGR の対称性をもち、𝐪 の指数は偶数。
           scripts/check_cd14_tilde.py、logs/check_cd14_tilde.log、tests/test_normalization_formula.py::test_n1n2_is_cd14_tilde_times_corner。
         - 【導出】CD14 (ii) を仮定すると、本研究の P = m·P_KhR（m は単項式、P_KhR は DGR の正規化の reduced KhR の Poincaré 多項式）。(N1) により両者は 𝐭 = −1 で同じ reduced HOMFLY-PT になるので、m = ±𝐭^k で ±(−1)^k = 1。両者とも係数が非負なので m = 𝐭^{2l}。
           𝐭^{2l} は DGR の対称性 σ で不変なので、P_KhR の対称性（Gorsky–Hogancamp–Mellit）がそのまま P に移る。よって命題 1 の条件が満たされ、[1] の正値性が従う。正規化のずれは 2J_r の一様なずれ（𝐭^{2l}）以外に無い。
         - 計算した結び目については CD14 (ii) を使わずに済む：保存値の P が非負で対称であることを直接確かめてあるので、命題 1 が直接適用できる。CD14 (ii) が要るのは「全ての代数的な iterated torus 結び目で [1] が正」という一般の主張だけ。
       - 長方形（[2], [1,1], [3], [1,1,1]）：CD14 は、代数的な結び目の長方形の DAHA の値が正であることを予想し、ホモロジーとの一致を期待している（予想としては述べていない）。予想 B の長方形の部分は、この期待を KN の BPS の段に移したものと位置づけられる。
       - [2,1]：CD14 は、長方形でない weight の正値性を主張していない。本研究の観察 R12（rCS_{[2,1]} は負の係数をもつが BPS の段で正）は、CD14 の範囲の外での新しい現象である。
       - 非代数的なケーブル（付録 A1）：CD14 は、正値性の無い非代数的な結び目（non-pseudo-algebraic、:1638, :3068）があることを認めている。本研究で DAHA の値の [1] に負の係数が出た型 II（REPORT.md §5.2）は、これに当たる。
6. **R11 の正規化の依存**【解決】
   - 2 により、x = 0 の判定が依存するのは「彩色ごとの (q/t)^{c_λ} が framing の固有値に従う」という一点の仮説だけである。符号と t = q の部分は HOMFLY-PT から決まる。F17 により、この仮説の下で正値性の主張は「F 単位の framing を基準にした列」の 1 点についてのものと書ける。
   - K6a〜K10 の x ≠ 0 は未計算（表の空欄。このセッションでは計算しない）。
7. **Zhu の用語**【確認済み】
   - Zhu 2105.02037 の Remark（HOMFLYPT-skein-theory3.tex:673-679）で、「refined LMOV conjecture」の refinement は KN の refined とは意味が違う、と Zhu 自身が明記している。本文で区別を書く。
8. **KN の書誌**【解決：食い違いは無い】
   - Crossref と INSPIRE で確認した：J. Knot Theory Ramifications **32** (2023) no. 06, 2041001、DOI 10.1142/S0218216520410011、オンライン公開 2020-07-08、印刷 2023-05。
   - DOI の「20」はオンライン公開の年に当たり、矛盾ではない。novelty_report.md §5 の記載（32 (2023) no.06, 2041001）は正しい。

## 8. 種別ごとの数（論文の骨格の確認）

- 定理（導出 ＋ 有限個の判定を含む）：R2, R4（同符号・異符号の一部）, R7, R8（同符号）。いずれも入力の前提つき。
- 命題・補題（導出）：F2, F4〜F12, F14, F15, F17, R3, R15。
- 計算結果：I2, I3, I5, I7, I9, R1, R5, R6, R9〜R13, R16, R17。
- 予想：F19（予想 A）、R14（予想 B）。
- 仮説：I4（(N2') のケーブルへの適用）、H_sat（I6 の前提）。
- 既知（成果として書かない）：I1, I8, F3（因子の形）, R10 の T(2,3), T(2,5)、R17 の t = q。
