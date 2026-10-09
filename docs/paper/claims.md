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
- 前提：仮説 H_sat（refined 衛星公式）、予想式（GNSSS 1512.07883、NRZS 1209.1409、Stošić 2402.03066 ＋ KRSS 1707.04017）、仮定 S_R（KN の U(1)_R の自明な作用）など。
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
| F14 | [1] の 2J_r の偶奇が (g, β) ごとにそろう ⇔ 𝐚 の次数が連続する範囲ごとに e_𝐭 + e_𝐚/2 の偶奇が一定（parity knot の条件） | 命題 ＋ 計算結果（801 個中 9 個） | Chandler–Gorsky の parity knot との規約の照合は未確認 | REPORT_knot_table.md §2、REPORT_phase4.md §7。test_knot_table_all_positive_and_parity_criterion | 条件自体は既知（Chandler–Gorsky 2209.13058）。BPS との同値性は新規 |
| F15 | 命題 2：\hat f_{[2]} が Laurent なら、非負の入力では dim H^{[2]} = (dim H^{[1]})² が必要（[1,1] も同様）。この条件は x によらない | 命題 | — | REPORT_knot_table.md §3.3、F9 | 導出は新規。9_42 の次元 401（GGS）と、9_42 が KN の形に乗らないこと（KN）は既知 |
| F16 | 全ての x の判定法：x の偶奇・向きごとにニュートン形を作り、係数の一様な符号と、起点と 0 の間の直接判定で、全ての x を有限個の判定に帰着する | 方法（補題としての正しさは導出） | — | scripts/framing_newton.py、test_trefoil_two_box_all_framings_newton | 方法。V4 のツイスト結び目の方法の転用 |
| F17 | (N2') の彩色依存部分はちょうど F（= pq）単位の framing：M_λ = (M_□(−(q/t)^{1/2})^{−F})^{\|λ\|}framing[λ]^F。基準の framing を偶数ずらすのは x のずれに吸収され、奇数ずらすと F10 により Laurent 性が壊れる | 命題 | — | framing_large_n.md §2.6、conventions.md §17 | 新規（正規化の意味づけ） |
| F18 | framed knot の refined large N 双対の物理的な説明：MV（平坦座標と framing）＋ KN（refined index、Schwinger、S_R）＋ AS の T 行列 ＋ Taki の framing の因子。抽出式（refined-GT1, GT2）は全ての x でそのまま適用できる | 注意（解釈） | 仮定：framing を変えても U(1)_R が保たれる。符号の幾何的な導出（V → −(q/t)^{1/2}V）はしていない | framing_large_n.md §3 | 新規（議論として） |
| F19 | **予想 A（framed refined LMOV の正値性）**：U(1)_R をもつ結び目（トーラス・代数的な iterated torus）と、ホモロジーを入力とする結び目で、全ての x、全ての ρ について、全体の符号を選べば \hat𝐍^{(x)}_ρ と \hat N^{(x)}_ρ が非負整数 | 予想 | 仮定 S_R（KN）＋ F18 | 証拠：§2〜§4 の計算結果（R1〜R14） | 新規の予想（KN の予想の framed への拡張） |

## 2. 入力（outline §4）

| ID | 主張 | 種別 | 前提 | 根拠 | 新規性 |
|---|---|---|---|---|---|
| I1 | KN のノートブックの 27 件（トーラス結び目 \|λ\| ≤ 3）を EHA（GN の構成）から全係数で再現し、BPS の表 13 個も再現した（規約の固定） | 既知の再現 | — | test_gn_reproduces_notebook、test_saved_tables_match_paper、test_kn_tables_reproduced | 既知。成果として書かない |
| I2 | ケーブルの refined 不変量：DAHA の定義を companion の段で分解した計算法が、CD14 の明示例 8 件と全係数一致する。素朴な Rosso–Jones の t 変形は誤り | 計算結果（方法の検証） | — | REPORT.md §2、conventions.md §7。test_matches_daha_example | 素朴な形の失敗は既知。正しい形は定義から従うので、新規性は弱い（novelty_report §5 (N4)）。「計算法として明示し検証した」と書く |
| I3 | ケーブルの t = q の極限が、ブレイドから計算した彩色 HOMFLY-PT と一致する | 計算結果 | — | test_cable_formula_matches_braid、test_t_equals_q_is_homfly | 照合 |
| I4 | 正規化 (N1)(N2')：[1] は t = q で HOMFLY-PT に一致し最低 𝐭 次数 0。彩色付きは [1] の単項式と framing の固有値の積 | 仮説（ケーブル） | — | conventions.md §5。トーラス結び目では確定（27 件）。ケーブルでは K5c で (N2) が破れ (N2') で回復（転置・Laurent 性）。test_normalization_N2prime | 新規（規則の提案） |
| I5 | M_□ = (−q/t)^{b−1}(a(q/t)^{1/2})^δ が保存済みの 22 結び目で成立する | 計算結果 | — | conventions.md §13。test_box_monomial_formula | 新規 |
| I6 | I5 の導出：H_sat とトーラス結び目の事実 (I)〜(IV) から | 定理（条件つき） | **仮説 H_sat**。(IV) の辞書は Mellit の Cor (daha) の正規化に依存し、その原典での確定は**未了**（符号と a → a/q） | conventions.md §13.1。test_framing_operator_shifts_slope、test_bottom_coefficient_is_rational_catalan | 新規（条件つき）。GN の Catalan 予想が Mellit の定理の系であることは既知と考えられるので、成果にしない |
| I7 | 結び目表 801 個の [1] の入力（yui-kr の HOMFLY-PT ホモロジー）が、DGR の対称性・ブレイドの HOMFLY-PT・トーラス結び目の rCS と整合する | 計算結果（入力の検証） | — | test_dgr_symmetry、test_euler_characteristic_is_braid_homfly、test_torus_knots_match_rcs | 照合 |
| I8 | 2 箱の入力：GNSSS の予想式（ツイスト・pretzel）。t_c 次数の変換は、トーラス結び目と既知の彩色 HOMFLY-PT で較正した | 既知（入力の予想式） | **予想式** | conventions.md §12。test_tc_conversion_reproduces_torus_knots、test_cable_inversion_reproduces_known_colored_homfly | 入力は既知の予想式 |
| I9 | thick な結び目の t_c 次数：Stošić の quiver ＋ KRSS の (quad-C) の組み合わせが、8_19 で EHA の [2]〜[1,1,1] と全項一致する | 計算結果（方法の較正） | **予想式** | open_problems.md B4、logs/b4_819_quad.log | 部分的に新しい（予想式の検証例） |
| I10 | 厳密な計算：素数 1〜3 個の補間と CRT、補間に使っていない素数での検証、G 展開の高速化 | 方法 | — | tests/test_speedup.py（7 件） | 方法（付録） |

## 3. 検証 I：[1]（outline §5）

| ID | 主張 | 種別 | 前提 | 根拠 | 新規性 |
|---|---|---|---|---|---|
| R1 | 11 交点以下の全素結び目 801 個で、\hat f_{[1]} の G 展開が閉じ、前後とも非負（命題 1 の確認） | 計算結果 | 入力は yui-kr | data/knot_table_summary.md、logs/knot_table_t1.log、test_knot_table_all_positive_and_parity_criterion | 新規（表として） |

## 4. 検証 II：2 箱（outline §6）

| ID | 主張 | 種別 | 前提 | 根拠 | 新規性 |
|---|---|---|---|---|---|
| R2 | ツイスト結び目 K_p：全ての p ≠ 0 で、2 箱は G 展開の前後とも正（x = 0） | 定理（導出 ＋ 有限個の厳密な判定） | **予想式（GNSSS）**、t_c 次数 δ' = −S | REPORT_phase4.md §5.1。c_2 を使わない LP 証明書あり（data/v4/v4_ck5_cert_*）。test_v4_twist.py（7 件）、logs/v4_closed.log, v4_newton_*.log | 新規（\hat f の段階） |
| R3 | ツイスト結び目の 2 箱は、全ての p で sl_2 指標の正係数の和（単峰） | 命題（導出 ＋ 有限個の判定） | 予想式 | REPORT_phase4.md §5.1（c_0..c_2 は計算） | 新規。注意として載せる |
| R4 | pretzel K_{s,t}：s, t が同符号なら全て正。異符号は min(\|s\|, t) ≤ 20 で証明、\|s\|, t ≤ 40 の範囲で直接確認。残りは 1 変数の 3 系列の非正値性に帰着（u ≤ 25 で成立） | 定理（同符号）＋ 定理（異符号の一部）＋ 計算結果 | **予想式（GNSSS）** | REPORT_phase4.md §5.2、REPORT_phase5.md §2。data/v4/ の証明書。test_v4_pretzel.py | 新規。**残りが未証明であることを定理の文に明記**（d_1 の LP は打ち切り、対角線は未着手） |
| R5 | thick な非トーラス結び目 10_139, 10_152 の 2 箱は G 展開の前後とも正（次元 441 = 21², 841 = 29²） | 計算結果 | **予想式（Stošić ＋ KRSS）** | RESULTS.md §3、logs/b4_s2_thick.log、data/b4/ | 新規（knots–quivers の正値性からは従わない、novelty_report §9） |
| R6 | 正値性は自動的でない：5_2 で t_c 次数の別の割り当て 12 通りは全て負の項をもつ | 計算結果 | 予想式 | REPORT_knot_table.md §4、logs/t3_nonunique.log | 新規。注意として載せる |
| R7 | 2 箱は、全ての整数 x で G 展開の前後とも正。対象は自明な結び目、T(2,3), T(2,5), T(2,7), T(3,4), 4_1, K1, K3 | 定理（導出 ＋ 有限個の厳密な判定） | 入力：EHA/DAHA、4_1 は KN のデータ（ホモロジー由来） | logs/framing_newton.log、F16、test_trefoil_two_box_all_framings_newton | 新規 |
| R8 | ツイスト結び目の (p, x)：p と x が同符号の向きでは全て正。逆符号は \|p\| ≤ 6、\|x\| ≤ 14 の直接計算で全て正（696 件） | 定理（同符号）＋ 計算結果（逆符号） | 予想式 | logs/framing_twist.log, framing_twist_direct.log | 新規。**逆符号は証明書なし**と明記 |
| R9 | framing を振った直接の判定：自明な結び目・トーラス結び目 4 個・4_1 は 3 箱まで、ツイスト結び目 12 個は 2 箱まで、K1, K3 は 3 箱まで。\|x\| ≤ 2〜6 の 732 件で負の項 0（9 件は \hat f = 0） | 計算結果 | 入力は R7・R8 と同じ | logs/framing_batch.log、test_unknot_and_trefoil_three_boxes_positive。以前の計算（未発表）の出力 9 件と全係数一致（logs/framing_notebook_check.log、test_trefoil_framing_2_matches_notebook） | 新規 |

## 5. 検証 III：3 箱（outline §7）

| ID | 主張 | 種別 | 前提 | 根拠 | 新規性 |
|---|---|---|---|---|---|
| R10 | トーラス結び目 T(3,4), T(3,5), T(2,11), T(3,7), T(4,5) の \|ρ\| ≤ 3 は G 展開の前後とも正（x = 0）。T(2,3), T(2,5) は KN の表の再現 | 計算結果 | — | RESULTS.md §1、logs/u4_torus.log, u4_T3_4.log、data/lmov/ | KN の予想の確認（新しい表だが、主張の新規性は弱い） |
| R11 | 代数的な iterated torus 結び目 7 個（K1, K3, K6a, K7, K8, K9, K10）の \|ρ\| ≤ 3 の全 6 ρ で、G 展開の前後とも正（x = 0）。K7 は 3 本ケーブル、K10 は 2 段の反復（3 箱は約 183 万項、係数の最大 1.1×10¹⁹） | 計算結果 | 入力は DAHA の定義どおり。正規化は I4（仮説） | RESULTS.md §1、data/cables/, data/lmov/, data/bps/、logs/u4_cables*.log, v0_K10_*.log。test_u4_three_box_positivity（results）、test_u4_bps_saved、test_positive_before_and_after_G | 新規（K7, K10 の \|λ\| = 3 の値自体も文献に無い、novelty_report §10） |
| R12 | 入力の rCS_{[2,1]} 自体は負の係数をもつのに、BPS の段で正になる | 注意（計算での観察） | — | REPORT_phase3.md §3 | 観察 |
| R13 | 3 箱を x ≠ 0 で判定：自明な結び目 \|x\| ≤ 6、T(2,3) \|x\| ≤ 6、T(2,5) と 4_1 \|x\| ≤ 4、T(2,7), T(3,4) \|x\| ≤ 3、K1, K3 \|x\| ≤ 2 で全て正。全ての x の証明書は未取得（6 根のニュートン形では見つからない） | 計算結果 | 4_1 の [2,1] は KN の候補 | logs/framing_batch.log、logs/framing_newton_first.log | 新規 |
| R14 | **予想 B**：代数的な iterated torus 結び目は、全ての ρ について KN の正値性を満たす（予想 A の x = 0、結び目のクラスを明示した版） | 予想 | 仮定 S_R | R10, R11 | 新規の予想（KN の予想のクラスの拡張） |
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

## 7. 書く前に解決・明記が必要な点（格付けで見つかったもの）

1. **I6 の Mellit の正規化**：原典での確定が未了（outline §3 の 6）。確定できなければ、I6 は「H_sat と (I)〜(III) の下で導出、(IV) は辞書の仮定の下で」と書き分ける。
2. **I4 の (N2') は、ケーブルでは仮説**：F17 で framing として意味づけたので、「彩色依存部分は framing、仮説として残るのは M_□（I5・I6）」と書ける。仮説の範囲を狭めて書く。
3. **R4・R8 の未証明部分**：pretzel の異符号の残りと、ツイスト結び目の逆符号の領域。定理の文で範囲を限定する。
4. **F14**：Chandler–Gorsky の規約との照合が未確認。「同じ形の条件」と書くにとどめるか、照合してから同値と書く。
5. **F18・F19 の仮定**：U(1)_R が保たれることと、符号の幾何的な導出が無いこと。予想 A の動機の節で明記する。
6. **R11 の正規化の依存**：x = 0 の判定は (N2') に依存するが、F17 により、正値性の主張は「F 単位の framing を基準にした列」の中の 1 点についてのものだと書ける。K6a〜K10 の x ≠ 0 は未計算（表の空欄、このセッションでは計算しない）。
7. **新規性で注意する用語**：Zhu の "refined LMOV for framed links" は非 refined の概念。本文で区別を明記する。
8. **KN の書誌**：INSPIRE（JKTR 32 (2023) 2041001）と DOI（2020 年の巻を示唆）が食い違う。引用の前に確定する。

## 8. 種別ごとの数（論文の骨格の確認）

- 定理（導出 ＋ 有限個の判定を含む）：R2, R4（同符号・異符号の一部）, R7, R8（同符号）。いずれも入力の前提つき。
- 命題・補題（導出）：F2, F4〜F12, F14, F15, F17, R3, R15。
- 計算結果：I2, I3, I5, I7, I9, R1, R5, R6, R9〜R13, R16, R17。
- 予想：F19（予想 A）、R14（予想 B）。
- 仮説：I4（(N2') のケーブルへの適用）、H_sat（I6 の前提）。
- 既知（成果として書かない）：I1, I8, F3（因子の形）, R10 の T(2,3), T(2,5)、R17 の t = q。
