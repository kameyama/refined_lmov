# 主要な結果（2026-10-07 時点）

本研究（KN arXiv:1703.05408 の refined large-N 双対性の、トーラス結び目以外への拡張）の主要な結果の一覧。
詳しい議論は各期の REPORT（第 5 期のまとめは docs/REPORT_phase5.md）、規約は docs/conventions.md、未解決の問題は docs/open_problems.md。

**結果の信頼性の根拠**：結果はコードで計算し、計算時に判定した値をデータとログに保存している。
結果が正しいと判断する根拠は、コードと機能テスト（tests/、既定の実行）が信頼できること。
機能テストは、既知の値の再現（1703.05408 のノートブックと表、DAHA の明示例、ブレイドから計算した HOMFLY-PT）と、各アルゴリズムの小さな入力での照合からなる。
保存済みの大きな結果を判定し直すテストには `results` の印を付け、内容や理論に問題が出たときだけ実行する（`pytest -m results`）。

印：【導出】解析的に導出、【計算】計算機で確認（範囲を明記）、【予想】、【仮説】。判定は refined LMOV 不変量の G 展開の前（太字の係数 \hat{𝐍}）と後（\hat N）の両方。
「正」は、全体の符号を一つ選べば全ての項が非負整数であること。

---

## 1. 3 箱（|ρ| ≤ 3）の正値性：定義どおりに計算できる結び目

[2,1] まで定義どおり（EHA/DAHA）に refined CS 不変量が得られる結び目では、|ρ| ≤ 3 の全 6 つの ρ で G 展開の前後とも正【計算】。

| 結び目 | 種類 | 最大の項数（3 箱、前） | データ | ログ |
|---|---|---|---|---|
| T(2,3), T(2,5) | トーラス（KN の refined CS と BPS の表を全項で再現） | — | data/torus/, data/bps/, data/lmov/ | logs/stage6_run*.log |
| T(3,4), T(3,5), T(2,11), T(3,7), T(4,5) | トーラス（3 箱は本研究で計算） | 約 2 万（T(4,5)） | data/torus/, data/lmov/ | logs/u4_torus.log, logs/u4_T3_4.log |
| K1 = T(2,3) の (2,13)-ケーブル | 代数的 | — | data/cables/, data/bps/ | logs/stage6_run*.log |
| K3 = T(2,3) の (2,15)-ケーブル | 代数的 | — | data/cables/, data/bps/ | 同上 |
| K6a = T(2,5) の (2,21)-ケーブル | 代数的 | — | data/cables/, data/bps/ | 同上 |
| K9 = T(2,7) の (2,31)-ケーブル | 代数的 | 約 23 万 | data/cables/K9_*, data/lmov/K9_* | logs/u4_cables*.log |
| K8 = T(3,4) の (2,25)-ケーブル | 代数的 | 約 26 万 | data/cables/K8_*, data/lmov/K8_* | logs/u4_cables3.log |
| K7 = T(2,3) の (3,19)-ケーブル | 代数的（3 本ケーブル） | 約 35 万 | data/cables/K7_*, data/lmov/K7_* | logs/u4_cables3.log |
| K10 = K1 の (2,53)-ケーブル | 代数的（2 段の反復） | 約 183 万 | data/cables/K10_*, data/lmov/K10_*（3 箱は K10_3.json.xz） | logs/v0_K10_*.log |

- 代数的な iterated torus 結び目への KN の正値性予想の拡張を支持する【予想】。K10 も |ρ| ≤ 3 の全 6 つの ρ で G 展開の前後とも非負【計算】（2026-10-08、REPORT_phase4.md §1.1）。
- 4_1（KN が [2,1] の候補を構成、付録 C）の 3 箱までの表は全項で再現した（既知の結果の再現で、本研究の成果ではない。data/bps/4_1_*）。
- 非代数的なケーブル（DAHA の値を入力）では [1] またはそれ以降で破れる。破れ方は 3 つの型に分かれる（docs/REPORT.md §5.2）。
  DAHA の値は非代数的なケーブルではホモロジーではないので、これ以上は扱わない（CLAUDE.md、open_problems.md D5）。

## 2. [1]：ホモロジーを入力とする任意の結び目で正

- **命題 1【導出】**：入力が非負係数で DGR の対称性をもてば、\hat f_{[1]} の G 展開は整数で閉じ、前後とも非負。DGR の対称性は Gorsky–Hogancamp–Mellit の定理なので、ホモロジーを入力とする限り任意の結び目で [1] は正。
- 11 交点以下の全素結び目 801 個（yui-kr）で確認【計算】（data/knot_table_summary.md、logs/knot_table_t1.log）。
- [1] の 2J_r の偶奇がそろう条件【導出】：𝐚 の次数が連続する範囲ごとに e_𝐭 + e_𝐚/2 の偶奇が一定（Chandler–Gorsky の parity knot と同じ形）。801 個中 9 個。

## 3. 2 箱：必要条件と、thin な族での証明

- **命題 2【導出】**：\hat f_{[2]} が Laurent 多項式なら、非負の入力では dim H^{[2]} = (dim H^{[1]})² が必要。9_42（401 ≠ 81）が KN の形に乗らない原因（docs/REPORT_knot_table.md §3）。
- **ツイスト結び目 K_p**：全ての p ≠ 0 で 2 箱は G 展開の前後とも正【導出（入力は GNSSS の予想式）】（docs/REPORT_phase4.md §5.1、`scripts/v4_closed.py`, `scripts/v4_newton.py`）。
- **pretzel K_{s,t}**：s, t が同符号なら全て正【導出（同上）】。異符号では min(|s|, t) ≤ 20 で正【導出 ＋ 有限個の判定】。|s|, t > 20 は 1 変数の 3 系列の非正値性に帰着し（u ≤ 25 で成立。d_1 と対角線の LP 証明書は打ち切り、2026-10-08）、|s|, t ≤ 40 の範囲で直接確認【計算】（REPORT_phase4.md §5.2、REPORT_phase5.md §2、`scripts/v4_pretzel*.py`、証明書 `data/v4/`）。
- 証明の要点：\hat f は撚りの数の指数多項式（根 1, 𝐚²𝐭², 𝐚⁴𝐭⁴, 𝐚⁴𝐪⁴𝐭⁸, 𝐚⁴𝐪⁻⁴𝐭⁴）。差分商（ニュートン形）の係数が一様な符号をもち、G 展開の後は「展開が閉じる ＋ 展開前に非負 ⇒ \hat N ≥ 0」から従う。
- **thick な非トーラス結び目 10_139, 10_152**（2026-10-07）：2 箱の \hat f_{[2]}, \hat f_{[1,1]} は G 展開の前後とも正【計算、予想式に基づく】。
  - 入力：Stošić 2402.03066 の quiver 型の式を展開し、KRSS 1707.04017 の予想式 (quad-C) で S²（t_c 次数）を作った。この方法は 8_19 で EHA の正解（[2] から [1,1,1] まで）と全項一致することを確かめてある（open_problems.md B4）。
  - 確認：次元は 441 = 21², 841 = 29²（命題 2 の必要条件を満たす）。t = q の極限は、(2,s)-ケーブルから逆に解いた [2], [1,1] 彩色 HOMFLY-PT と一致。展開の順序によらない。
  - 意味：命題 2 の必要条件を満たす thick な非トーラス結び目で、2 箱の正値性が成り立つ例。文献の範囲では見当たらず、knots–quivers の DT 不変量の正値性からも従わない（docs/novelty_report.md §9）。
  - 方法の較正のうち、8_19 の t_c 次数と refined CS の一致そのものは既知（GNSSS 1512.07883）。
  - `scripts/b4_s2_thick.py`、`logs/b4_s2_thick.log`、データ `data/b4/`。
- 2 箱の正値性は自動的ではない：5_2 で t_c 次数の別の割り当てを 12 通り試すと全て負の項が出る【計算】。

## 4. ケーブルの正規化の閉じた式

- 正規化 (N1)(N2') の単項式は M_□ = (−q/t)^{b−1}(a(q/t)^{1/2})^{δ}（b はブレイドの本数、δ は種数型の量）。保存済みの全 22 結び目で成立【計算】。
- refined 衛星公式（仮説 H_sat）と、トーラス結び目についての定理（Mellit の rational shuffle と torus knot homology、Franks–Williams）から導出【導出、H_sat と演算子の正規化の対応（計算で確認）の下で】（docs/conventions.md §13, §13.1）。

## 5. [2,1]（本研究では打ち切り、別研究へ）

- 双曲結び目の [2,1] は、構造的な条件では決まらない（解の差は結び目によらない因子 × 自由な多項式）【計算＋導出】。
- 弱い判定：正値性を満たす [2,1] は 4_1 と 6_1 に存在、5_2 は有理数解のみ【計算】。
- 定義に基づく構成は文献にもソフトにも無い。別研究の計画：docs/research_plan_21.md。
- 別研究 X0（2026-10-08）：t = q では、ダブルツイスト結び目の [2,1] の Morozov の式（1612.00422、1903.00259 の 𝓑^{[2,1]}）が、ツイスト結び目 |p| ≤ 10 と K_{±2,2}, K_{±3,2} の 3 本ケーブルの逆解きと全係数一致【計算】。
  長方形は KNTZ の式と一致（[2,2] は 3_1, 4_1, 5_2, 6_1）。データ data/x0/colored_t_eq_q.txt、ログ logs/x0_check_*.log、テスト tests/test_x0_double_twist.py。
- 別研究 X2（2026-10-08、停止点 1）：[2,1] の族の仮定（Morozov の構造、正則な行は KNTZ の refined の 𝐓𝐰）の下で、EGP・U(2) とは矛盾しないが、t = q・EHA の三葉結び目・Laurent 性・転置・ミラーを加えても一意に決まらない【計算＋導出】。
- 別研究 X1（2026-10-08）：KNTZ の refined の長方形の式を実装し、著者の補助ファイル・GNSSS・三葉結び目の EHA（[2,2] まで）と一致【計算】。tests/test_x1_refined.py、conventions.md §15。

---

## 計算基盤（結果の再現に使うもの）

| 部分 | コード | 主な機能テスト |
|---|---|---|
| Macdonald 多項式 | src/rcsknots/macdonald.py | tests/test_stage1_macdonald.py |
| 楕円 Hall 代数・GN の構成 | src/rcsknots/ehall.py, gn.py | tests/test_stage2_torus.py（ノートブックの 27 件を全項で再現） |
| 厳密な補間（mod p） | src/rcsknots/exact.py, interp.py | tests/test_speedup.py |
| refined LMOV と G 展開 | src/rcsknots/bps.py, extract.py | tests/test_stage3_bps.py（KN の表）, tests/test_speedup.py（G 展開） |
| ケーブル | scripts/stage5_cables.py | tests/test_stage5_cables.py（DAHA の明示例）, tests/test_stage4_cable_homfly.py（ブレイドの HOMFLY-PT） |
| 正規化 | src/rcsknots/normalize.py, scripts/renormalize_cables.py | tests/test_normalization_formula.py |
| 結び目表の入力 | src/rcsknots/knot_table.py | tests/test_knot_table.py |
| 2 箱の入力（GNSSS の式） | scripts/t3_pretzel.py | tests/test_v4_twist.py, tests/test_v4_pretzel.py |

訂正の履歴：pretzel の入力スクリプトの浮動小数点のバグ（2026-10-07 修正、REPORT_phase4.md §8）。修正前のログは logs/*_float_bug.log。
