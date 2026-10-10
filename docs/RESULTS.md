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

- G 展開の後の \hat N_{ρ,g,β,J_r} は data/bps/{K}_{ρ}.json（\hat f も同梱）。T(3,4)〜T(4,5), K7, K8, K9 の分は 2026-10-09 に保存済みの \hat f から `scripts/export_bps.py` で書き出した（K7〜K9 の 3 箱は .json.xz のみ。計算時の判定 48 件と項数・負の数が一致、tests/test_phase3.py::test_u4_bps_saved）。K10 も 2026-10-09 に書き出した（同じ方法、全 6 ρ で負の項 0、計算時の判定と項数が一致、196 秒）。K10 の 2・3 箱は .json.xz のみ（3 箱は各 87〜90 MB）。
- 代数的な iterated torus 結び目への KN の正値性予想の拡張を支持する【予想】。K10 も |ρ| ≤ 3 の全 6 つの ρ で G 展開の前後とも非負【計算】（2026-10-08、REPORT_phase4.md §1.1）。
- 4_1（KN が [2,1] の候補を構成、付録 C）の 3 箱までの表は全項で再現した（既知の結果の再現で、本研究の成果ではない。data/bps/4_1_*）。
- 非代数的なケーブル（DAHA の値を入力）では [1] またはそれ以降で破れる。破れ方は 3 つの型に分かれる（docs/REPORT.md §5.2）。
  DAHA の値は非代数的なケーブルではホモロジーではないので、これ以上は扱わない（CLAUDE.md、open_problems.md D5）。

## 1.2 トーラス結び目の 4 箱以上と 4_1 の 4 箱（2026-10-10）

- トーラス結び目 8 個で |ρ| = 4〜9（結び目ごとの範囲は data/highbox_summary.md）、計 340 件の (結び目, ρ) で G 展開の前後とも負の項 0【計算】。
  T(2,3) は 9 箱（30 個の ρ）、T(2,5), T(2,7), T(2,9) は 8 箱、T(3,4), T(3,5) は 6 箱、T(3,7), T(4,5) は 5 箱まで。KN の表（T(2,3) の 3 箱まで、T(2,5) の [1]）を超える範囲での KN の予想の確認。
  計算法は conventions.md §18（bpsgen.py、保存済みの 3 箱を全係数で再現）。ログ logs/highbox_torus.log, logs/highbox_batch.log。
- 4_1 の 4 箱：[3,1] の refined の値が未知のため弱い判定。σ = −1 は G 展開の前の正値性で解なし、σ = +1 は前後の正値性を満たす有理数解が存在、整数解は 50 分の MILP で未決【計算】。
  KN の [2,1] の候補、U(2) の仮説、未知数の窓に依存する（progress.md 2026-10-10、scripts/f41_box4.py）。

## 2. [1]：ホモロジーを入力とする任意の結び目で正

- **命題 1【導出】**：入力が非負係数で DGR の対称性をもてば、\hat f_{[1]} の G 展開は整数で閉じ、前後とも非負。DGR の対称性は Gorsky–Hogancamp–Mellit の定理なので、ホモロジーを入力とする限り任意の結び目で [1] は正。
- 11 交点以下の全素結び目 801 個（yui-kr）で確認【計算】（data/knot_table_summary.md、logs/knot_table_t1.log）。
- [1] の 2J_r の偶奇がそろう条件【導出】：𝐚 の次数が連続する範囲ごとに e_𝐭 + e_𝐚/2 の偶奇が一定（Chandler–Gorsky の parity knot と同じ形）。801 個中 9 個。

## 3. 2 箱：必要条件と、thin な族での証明

- **命題 2【導出】**：\hat f_{[2]} が Laurent 多項式なら、非負の入力では dim H^{[2]} = (dim H^{[1]})² が必要。9_42（401 ≠ 81）が KN の形に乗らない原因（docs/REPORT_knot_table.md §3）。
- **ツイスト結び目 K_p**：全ての p ≠ 0 で 2 箱は G 展開の前後とも正【導出（入力は GNSSS の予想式）】（docs/REPORT_phase4.md §5.1、`scripts/v4_closed.py`, `scripts/v4_newton.py`）。
- **pretzel K_{s,t}**：s, t が同符号なら全て正【導出（同上）】。異符号では min(|s|, t) ≤ 20 で正【導出 ＋ 有限個の判定】。|s|, t > 20 は 1 変数の 3 系列の非正値性に帰着し（u ≤ 25 で成立。d_1 と対角線の LP 証明書は打ち切り、2026-10-08）、|s| ≤ 30 かつ t ≤ 40 の範囲で直接確認【計算】（|s|, t > 20 の 200 点。2026-10-09 に範囲の記述を訂正）（REPORT_phase4.md §5.2、REPORT_phase5.md §2、`scripts/v4_pretzel*.py`、証明書 `data/v4/`）。
- 証明の要点：\hat f は撚りの数の指数多項式（根 1, 𝐚²𝐭², 𝐚⁴𝐭⁴, 𝐚⁴𝐪⁴𝐭⁸, 𝐚⁴𝐪⁻⁴𝐭⁴）。差分商（ニュートン形）の係数が一様な符号をもち、G 展開の後は「展開が閉じる ＋ 展開前に非負 ⇒ \hat N ≥ 0」から従う。
- **thick な非トーラス結び目 10_139, 10_152**（2026-10-07）：2 箱の \hat f_{[2]}, \hat f_{[1,1]} は G 展開の前後とも正【計算、予想式に基づく】。
  - 入力：Stošić 2402.03066 の quiver 型の式を展開し、KRSS 1707.04017 の予想式 (quad-C) で S²（t_c 次数）を作った。この方法は 8_19 で EHA の正解（[2] から [1,1,1] まで）と全項一致することを確かめてある（open_problems.md B4）。
  - 確認：次元は 441 = 21², 841 = 29²（命題 2 の必要条件を満たす）。t = q の極限は、(2,s)-ケーブルから逆に解いた [2], [1,1] 彩色 HOMFLY-PT と一致。展開の順序によらない。
  - 意味：命題 2 の必要条件を満たす thick な非トーラス結び目で、2 箱の正値性が成り立つ例。文献の範囲では見当たらず、knots–quivers の DT 不変量の正値性からも従わない（docs/novelty_report.md §9）。
  - 方法の較正のうち、8_19 の t_c 次数と refined CS の一致そのものは既知（GNSSS 1512.07883）。
  - `scripts/b4_s2_thick.py`、`logs/b4_s2_thick.log`、データ `data/b4/`。
- 2 箱の正値性は自動的ではない：5_2 で t_c 次数の別の割り当てを 12 通り試すと全て負の項が出る【計算】。

## 3.5 framed knot（framing を変えた refined 不変量）でも正（2026-10-09）

- framing x の rCS_λ = framing[λ]^x rCS_λ（framing[λ] = (−(q/t)^{1/2})^{|λ|} q^{n(λ')} t^{−n(λ)}、conventions.md §17）で、自明な結び目・トーラス結び目 4 つ・4_1 は 3 箱まで、ツイスト結び目 12 個は 2 箱まで、代数的なケーブル K1, K3 は 3 箱まで、framing を |x| ≤ 2〜6 で振った 732 件で、G 展開の前後とも負の項 0【計算】（9 件は \hat f = 0）。
- 基準：ユーザーの以前の計算（未発表）の出力 9 件と全係数一致。[1] は 𝐭^x 倍だけで自明【導出】。
- データ・コード：logs/framing_batch.log、logs/framing_notebook_check.log、scripts/framing_check.py、tests/test_framing.py。
- 拡張（2026-10-09）：T(3,5), T(2,11), T(3,7), T(4,5) は 3 箱まで x ∈ [−10, 10]、代数的なケーブル K6a, K7, K8, K9 は 3 箱まで x ∈ [−10, 10]、K10 は 2 箱まで x ∈ [−10, 10] と 3 箱の x ∈ {−2, −1, 1, 2}（x = 0 は framing 0 の既存の判定）。
  計 1083 件（x × ρ）で G 展開の前後とも負の項 0【計算】（2026-10-09 訂正：以前は 1167 件と書いていたが、data/framed/ の判定を数えると 1083 件）。ほかに K1 の 3 箱 x = ±6 の 6 件も負の項 0。前の 732 件と合わせて 1821 件。最大の \hat f は K10 の 3 箱 x = 2 の約 189 万項（G 展開の前）。表は data/framed_summary.md。
  データ・コード：data/framed/（(結び目, x, n) ごとの判定）、scripts/framing_big.py（\hat f だけを補間）、scripts/framing_big_batch.sh、scripts/framing_big_K10_n3.sh、logs/framing_big_*.log。
- framing を記号として扱う（2026-10-09）：𝐭^{−nx}\hat f^{(x)} は x の指数和で、命題 1・命題 2・Laurent 性・G 展開の閉じ方は framing によらない【導出】。
  2 箱は自明な結び目、T(2,3), T(2,5), T(2,7), T(3,4), 4_1, K1, K3 で**全ての整数の framing** で G 展開の前後とも正【導出 ＋ 有限個の厳密な判定】（logs/framing_newton.log）。
  ツイスト結び目の 2 箱は p と framing が同じ符号の向きで全て正（逆向きは |p| ≤ 6、|x| ≤ 14 の直接計算のみ）。3 箱の全ての framing の証明書は未取得。

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
| 任意の箱数の \hat f | src/rcsknots/bpsgen.py, scripts/highbox_torus.py | tests/test_bpsgen.py, tests/test_f41_box4.py |

訂正の履歴：pretzel の入力スクリプトの浮動小数点のバグ（2026-10-07 修正、REPORT_phase4.md §8）。修正前のログは logs/*_float_bug.log。
