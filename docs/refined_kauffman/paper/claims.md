# 主張の格付け表（2026-10-20、論文執筆の準備：outline.md §3 の 1）

論文に載せる主張を種別に割り当て、根拠と前提を対応させる。種別は docs/paper/claims.md と同じ：
- 定理・命題・補題（Thm / Prop / Lem）：【導出】に限る。
- 計算結果（CR）：【計算】。範囲を書く。
- 予想（Conj）、仮説（Hyp）、既知（K、成果として書かない）、注意（Rem）。

テストは tests/test_rk_*.py。`results` は保存済みの大きな結果を判定し直すテスト（既定では実行しない）。

---

## 1. 枠組み（§2）

| ID | 主張 | 種別 | 前提 | 根拠 | 新規性 |
|---|---|---|---|---|---|
| S1 | AgSc の SO(2N) の refined unknot dim_λ は Koornwinder K_λ(x;q,t;1,−1,q^{1/2},−q^{1/2}) の主特殊化 x = t^ρ に等しい | CR（6 つの λ、N = 2, 3） | — | test_rk_refined_unknot::test_koornwinder_evaluation_is_DN, test_dim_DN_matches_generic_evaluation | D_N の Macdonald = Koornwinder の特殊化は既知（van Diejen など、原典は未確認）。照合として書く |
| S2 | Mariño 型（bosonic 核、½、Σ_{d 奇数}）の素朴な refined 化は unknot の 2 箱で U の極が打ち消さない | CR | — | report.md §3、test_rk_refined_unknot::test_poles | 新規（否定的な結果として短く） |
| S3 | Mimachi の双対 Cauchy 恒等式の両側に AgSc 型（±1, ±q^{1/2}）と、(q,t) を入れ替えた LRW 型の境界パラメータが現れる | Prop（恒等式の読み替え） | Mimachi の恒等式（文献、原典の番号は未確認） | report.md §4、test_rk_koornwinder::test_mimachi | 読み替えは新規と思われる |
| S4 | Askey–Wilson 係数 c_{2j}(n) の閉じた形 | Conj（有理点 4 つ、j ≤ 5、n ≤ 30） | — | test_rk_askey_wilson::test_top_coefficients_closed_form | 既知かどうか未調査（Askey–Wilson の明示式から従う可能性） |
| S5 | 一様な規格化 W_μ = (t/q)^{δ\|μ\|/2} R_μ P_μ(t^ρ) は行列式のずらしと荷電共役で不変（δ は [1] の superpolynomial の δ 次数） | CR（trefoil 16 組、T(2,5)・T(2,7) 各 13 組、二重ツイスト 7 結び目） | data の規格化（KN の (N1)+(N2)） | test_rk_trefoil::test_trefoil_uniform_normalization, test_torus_family_uniform_normalization, test_double_twist_offset_and_41_fixed_N | 新規 |
| S6 | 一様な規格化では、固定 N の framing の単項式（極の打ち消し＋t = q で 1）はすべて 1 になる。DAHA-Jones（鉛直 framing）とも framing 固有値の冪を除いて一致する | CR（N = 2, 3、\|R\| ≤ 3、T(2,3), T(2,5), T(2,7)） | — | test_rk_trefoil::test_T25_fixed_N_uniform_hatted、test_rk_daha::test_dn_daha_jones_matches_data_and_as | 新規 |
| S7 | 一様な規格化が framing 0 の値であること | Hyp | — | S5, S6 と DAHA の一致 | — |
| S8 | reformulation は Mariño の M(q^{−1})。M(t), M(t^{−1}) では Laurent にならず、M(q) は R ↔ R^T | CR | — | report.md §8、test_rk_cover::test_hopf_hatted_crosscap_invariants | M = KN の M(x) の x = q^{−1} の同定は新規と思われる |

## 2. orientifold の refined 指数（§3）

| ID | 主張 | 種別 | 前提 | 根拠 | 新規性 |
|---|---|---|---|---|---|
| O1 | 商の一粒子状態を Ĩ の対・不変（ε = ±1、モードの射影）・凍結に分けると g = f − ½f^{cov} = ½E/(1+t) + C、E = Σ_不変 ε m、C = Σ_凍結 m | Prop | AgSc §3 の M 理論の定義（積表示 (3.13)）、t-brane（M5 が z₂ 平面を張る） | theory_orientifold_index.md §2 | 新規（開弦の全 d の構造の導出） |
| O2 | E は Ĩ で自分に移るセクター（(偶,偶)、p₂ 型の方向）にだけ現れる。E_(2) = −E_(1,1) | 系 | O1 | theory §2、Hopf と結び目の計算（O4, D3） | 新規 |
| O3 | E と C の分離は指数からは決まらない（凍結した状態 m ≡ 不変な状態 m と mt の組）。指数から決まるのは E(t = −1) だけ | Prop | O1 | theory §3 | 新規 |
| O4 | 検証：Hopf（a, b ≤ 4）で全 d の g_{ab} の分母は (偶,偶) で 2(1+t)・単項式、他は単項式。Σ_{d 奇数} の 4(1+t)(1+t²) は全 d で消える | CR | — | test_rk_cover::test_hopf_crosscap_all_d_even_sector | 新規 |
| O5 | 閉弦：D_N の積表示の −log Z と AgSc §6 の 3 項は、N の 1 次式を除いて一致（N = 2..7） | CR | — | test_rk_closed_orientifold::test_agsc_large_n_matches_product_formula | AgSc §6 の式は既知（K）。積表示との照合は検算 |
| O6 | AM の D/C 型の「奇数 n の c = 1 ＋ Klein 瓶（c = 2）」は、全 n の凍結した状態＋不変な状態（ε = +1、z₂ の偶数モード）に恒等的に等しい。Klein 瓶の項は −½ψ²C と ½E/(1+t) の差 | Prop（有理式の恒等式） | AM (3.4)–(3.7) | test_rk_closed_orientifold::test_am_klein_bottle_is_frozen_plus_invariant、theory §8 | 新規と思われる（AM と AgSc の整理の関係） |
| O7 | 正値性は E を含まない方向（p₁ⁿ、奇数のセクター）で E の不定性によらず判定できる | 系 | O1, O2 | theory §4 | 新規 |
| O8 | AgSc §4 の「差 × 和」は ψ²z で、取っ手一つの因子ではない（PU の捻ったトレースは同じ変数で z 一つ） | Rem（文献の読み） | — | theory §5 の訂正 | — |

## 3. unknot と Hopf 絡み目（§4）

| ID | 主張 | 種別 | 前提 | 根拠 | 新規性 |
|---|---|---|---|---|---|
| H1 | unknot：f_{10} = −1 − (ν − ν^{−1})/(t^{1/2} − t^{−1/2})、f_{a0} = 0（a ≥ 2）。Mariño 型の unknot とクロスキャップ −1 の和 | Prop | 双対 Cauchy 核 | report.md §4.1、test_rk_hopf::test_unknot_and_hopf_fundamental | 新規（refined SO の開弦の unknot） |
| H2 | Hopf：f_{ab} ∈ (1 − q)/(1 − t)·Z[q^±, t^±, U^±]（(偶,偶) は /(1 − t²)）、a, b ≤ 4。Schur 基底の f_{R,b}（\|R\| ≤ 4、b ≤ 2）でも同じ | CR | — | test_rk_hopf::test_symbolic_matches_finite_N, test_multi_probe_consistency | 新規 |
| H3 | f_{11} = −f_{10} z_q(ν~ − ν~^{−1})、ν~ = ν(q/t)^{1/2}、f_{1b} = f_{11} χ_{b−1}(U/t) | CR（b ≤ 4） | — | 同上 | 新規 |
| H4 | g_{1b} = z_q(ν~ − ν~^{−1}) χ_{b−1}(ν~/q^{1/2})。t は ν~ だけを通して入る | CR（b ≤ 4） | — | test_rk_cover::test_hopf_crosscap_fundamental_sector | 新規 |
| H5 | ĝ_{R,[1]}（\|R\| ≤ 4）：hook (n−d,1^d) で (−1)^d z_q(ν~ − ν~^{−1})ν~^{2d−n+1}、それ以外で 0。t = q で Mariño の ĝ_{[2],[1]}, ĝ_{[1,1],[1]} に R ↔ R^T と符号で一致 | CR | — | test_rk_cover::test_hopf_hatted_crosscap_invariants | 新規 |
| H6 | E_22, E_24 は (i) t = −1 と (ii) t = q の条件で一意に決まり、積の形。(4,4) は一意でない | CR | 固定セクターの部分和の仮定 | test_rk_cover::test_E22_unique_from_t_minus1_and_t_q | 新規 |
| H7 | Ĩ は「固定セクターの入れ替え＋Adams」として働く：E_{R,S} = Σ c^R_{2;ρ} c^S_{2;σ} ε_{ρσ}（\|R\| ≤ 4、\|S\| ≤ 2 の 10 セクター、4×4 で t = −1 の 25 組） | CR | — | test_rk_cover::test_multi_probe_cover_and_adams_b2、results の 2 件 | 新規 |
| H8 | Ĉ_{R,S} の refined だけの部分は (ν~² − q)(q − t) で割り切れる（10/10）。冪ごとの符号は 8/10 でそろい、転置の対 (3,1)/(1,1), (2,1,1)/(2) で破れる | CR | — | test_rk_cover::test_saved_42_sectors_adams_and_hatted_factor（results） | 新規 |

## 4. トーラス結び目の [1]（§5）

| ID | 主張 | 種別 | 前提 | 根拠 | 新規性 |
|---|---|---|---|---|---|
| T1 | h/g 分解（½h^{ref} = 𝒫(𝐚_U,𝐪,𝐭)(ν − ν^{−1})/(𝐪 − 𝐪^{−1})） | Hyp | DGR の superpolynomial | report.md §6 | — |
| T2 | g^{ref} が Laurent になる単項式は T_[1]^{−2(2m+1)}(q/t)^{δ/2} に一意（d は Kauffman 多項式の形で決まる）。(q/t)^{δ/2} は被覆の規格化のずれ | CR（m = 1..5） | T1 | test_rk_trefoil::test_torus_family_framing_monomial, test_torus_family_uniform_normalization | 新規。AS 1105.5117 §7 の SU(N) の √(q/t) と同じ形（既知の観察との照合） |
| T3 | g^{true}_[1]：(N) SO(2) の点で 1、(P) 𝐭 について偶、(c=2) −(𝐚²−𝐪²)(𝐭²−1)𝐚^{−(4m+1)}𝐪^{−(2m−1)}𝐭^{−4m}Π_m、(c=1) G 展開のブロックが 2 種 | CR（m = 1..5） | T1, S7 | test_rk_trefoil::test_crosscap_criterion_uniform_normalization, test_refined_crosscap_c1_blocks | 新規 |
| T4 | Π_m の閉じた形（係数 +1 の m(m+1)/2 項） | Conj（m ≤ 3 から推測、m = 4, 5 で確認） | — | test_rk_trefoil::test_torus_family_refined_only_part | 新規 |
| T5 | AgSc (6.15)（N = 2, 3）は Gorsky–Negut の SU(2)²・SU(4) の値と一致 | CR（有理点 2 つで厳密） | — | test_rk_trefoil::test_as_615_matches_gn_at_N2_N3 | 照合 |
| T6 | 判定式（正値性の形）：(c=1) G 展開の後でブロックごとに同符号か (𝐭²−1)×同符号、(c=2) 普遍因子 × 正 | Conj | T1, S7 | T3 | 新規 |

## 5. 一般の N の多箱（§6）

| ID | 主張 | 種別 | 前提 | 根拠 | 新規性 |
|---|---|---|---|---|---|
| D1 | D_N・GL_N の DAHA の多項式表現での DAHA-Jones の直接計算（一般の t）。DAHA の関係式、Macdonald 作用素の固有関数を検証 | CR（N = 2, 3, 4） | Cherednik 1111.6195 の定義 | test_rk_daha::test_dn_relations_and_macdonald_eigen | 定義は既知（K）。D 型の一般の t の計算は Cherednik にも例がある（範囲は要確認） |
| D2 | D_N：JD = T_λ^{−2n} × 鏡像（固定 N の W_λ、一様な規格化）。N = 4, 5 の [1] は AgSc (6.15) | CR | — | test_rk_daha::test_dn_daha_jones_matches_data_and_as | 新規（照合） |
| D3 | GL の単項式の規則 M_b = q^{−n(b,b)} t^{−2n(b,ρ)+(m+1)\|b\|(N−1)}（n = 2m+1） | CR（T(2,3), T(2,5)、N = 2..5） | — | test_rk_daha::test_gl_daha_jones_monomial_rule_and_adjoint | 新規 |
| D4 | 一般の N の補間：trefoil の (1),(2),(1,1),(3),(2,1),(1,1,1)、T(2,5) の (1),(2),(1,1)、被覆の随伴と composite (2;1),(11;1)。補間に使わない N で照合 | CR | 補間の窓（s の次数の傾き） | test_rk_daha::test_interpolation_vector_reproduces_as、results の test_general_N_three_box_tq_and_T25_universal_factor | 新規（refined SO(2N) の多箱の一般の N） |
| D5 | 大きな N の 2 箱：g_(2), g_(1,1) の分母は 2·単項式·(1+t) だけ。g_(2)+g_(1,1) と (1+t)(g_(2)−g_(1,1)) は整数 Laurent。U = t², t³ で固定 N と一致 | CR（T(2,3), T(2,5)） | — | test_rk_daha::test_large_N_two_box_crosscap、results の crosscap_structure | 新規 |
| D6 | 大きな N の 3 箱（trefoil）：g_R は整数 Laurent で (1+t) の極なし（E = 0） | CR | — | results の test_large_N_three_box_crosscap | 新規 |
| D7 | t = q で Mariño との辞書：W は Mariño の標準 framing の鏡像に比 1、ĝ は奇数 d に直して −鏡像(ĝ^Mar)。0904.1088 (4.8) の ĝ_[2] の z⁹ の誤植を訂正 | CR | — | test_rk_daha::test_large_N_two_box_marino_dictionary、test_rk_unrefined::test_trefoil_marino_48 | 誤植の指摘は新規 |
| D8 | KN の q 枠と t 枠の hatted は 2 箱で一致 | CR | — | test_rk_daha::test_large_N_two_box_kn_frames | 新規 |
| D9 | 普遍因子：ĝ − ĝ\|_{t=q}(ν → ν w^k) が (q − t)(ν² − t) で割り切れるのは k = +1 だけ（[1]、2 箱、3 箱、T(2,5) の 2 箱） | CR（k = −3..3） | — | test_rk_closed_orientifold::test_shift_direction_unique_for_trefoil_one_box、logs/refined_kauffman/b_shift_direction.log、results の universal_factor | 新規 |
| D10 | 多箱の c = 1 部分は σ 対称でなく G 展開できない。SO の R_λ に super-duality に当たる対称性が無い | CR | — | report.md §9、scripts/rk_b_2box_structure.py | 新規（否定的） |

## 6. 正値性（§7）

| ID | 主張 | 種別 | 前提 | 根拠 | 新規性 |
|---|---|---|---|---|---|
| P1 | [1]（T(2,2m+1)、m ≤ 5）は判定式 T6 を満たす | CR | T1, S7 | T3 | 新規 |
| P2 | 多箱の普遍因子で割った商の符号は大部分のブロックで混在（2 箱の和 5/10、T(2,5) 6/18、3 箱 2/15, 3/14, 2/13） | CR | — | progress.md 続き 28〜33 | 新規（否定的） |
| P3 | 基底の取り替え・二次量子化の指数・被覆の brane/image の分割でも符号はそろわない。被覆の brane–brane は KN の正値性を満たし、brane–image は (𝐭² − 1) を含んで混在 | CR | — | progress.md 続き 27〜31 | 新規 |
| P4 | 固定セクターの整数計画（状態の部分和、捻ったトレース）で E を決める試みは、T(2,5) で解なし（指数の上限では少なくとも 2 組の隠れた対が要る）。結び目の E の境界因子は和 z_q^+ | CR | 捻りの規則（Hyp） | test_rk_twisted_multiplet（results）、theory §9 | 新規 |
| P5 | 被覆の随伴 R_adj はどの偶奇の符号規則でも正にならない（Cherednik–Elliot の「composite の DAHA 超多項式は正値性をもたない」と一致） | CR ＋ K | — | scripts/rk_b_adj_positivity.py、1503.01441 序論 | 照合 |
| P6 | 多箱の正値性は、凍結した状態の多重項の基底と E・C の分離が決まるまで判定できない | Rem（現状の結論） | O3, P2〜P5 | — | — |

## 7. 既知（成果として書かない）

- Mariño 0904.1088 の構造と trefoil の数値、NRZ 1310.2240 の 4_1 の表（K0 の再現：test_rk_unrefined）。
- AgSc 1202.4456 の T(2,2m+1) の [1] の一般式、Hopf、大きな N の閉弦（§6）。
- AM 2205.12832 の閉弦の GV 型展開。
- Cherednik の DAHA-Jones の定義。Cherednik–Danilenko の D 型の t = q。
- KN 1703.05408 の M(x)、q 枠と t 枠。
