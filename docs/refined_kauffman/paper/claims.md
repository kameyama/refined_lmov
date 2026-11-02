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
| O2 | E は Ĩ で自分に移るセクター（(偶,偶)）にだけ現れる（導出）。p₂ 型の方向に入り E_(2) = −E_(1,1) となること（計算） | 注意＋CR | O1 | theory §2、草稿の注意 rem:Esector、Hopf と結び目の計算（H7, D5） | 新規 |
| O3 | E → E + 2(1+t)X, C → C − X で g は不変。g から決まるのは E の 2(1+t) を法とした類だけ（E(−1) は決まる）。状態の言葉では「凍結 m ＋ 対 {m}, {mt}」と「ε = +1 の不変 m, m, mt, mt」が f と f^{cov} の両方で区別できない | Prop | O1 | 草稿の命題 prop:EC（2026-10-20 の照合で訂正。theory §3 の「m と mt の組」は誤り） | 新規 |
| O4 | 検証：Hopf（a, b ≤ 4）で全 d の g_{ab} の分母は (偶,偶) で 2(1+t)・単項式、他は単項式。Σ_{d 奇数} の 4(1+t)(1+t²) は全 d で消える | CR | — | test_rk_cover::test_hopf_crosscap_all_d_even_sector | 新規 |
| O5 | 閉弦：D_N の積表示の −log Z と AgSc §6 の 3 項は、N の 1 次式を除いて一致（N = 2..7） | CR | — | logs/refined_kauffman/b_closed_orientifold.log（N = 2..7）、test_rk_closed_orientifold::test_agsc_large_n_matches_product_formula（N = 2..5） | AgSc §6 の式は既知（K）。積表示との照合は検算 |
| O6 | AM の D/C 型の「奇数 n の c = 1 ＋ Klein 瓶（c = 2）」は、全 n の凍結した状態＋不変な状態（ε = +1、z₂ の偶数モード）に恒等的に等しい。Klein 瓶の項は −½ψ²C と ½E/(1+t) の差 | Prop（有理式の恒等式） | AM (3.4)–(3.7) | test_rk_closed_orientifold::test_am_klein_bottle_is_frozen_plus_invariant、theory §8 | 新規と思われる（AM と AgSc の整理の関係） |
| O7 | 正値性は E を含まない方向（p₁ⁿ、奇数のセクター）で E の不定性によらず判定できる | 系 | O1, O2 | theory §4 | 新規 |
| O8 | AgSc §4 の「差 × 和」は ψ²z で、取っ手一つの因子ではない（PU の捻ったトレースは同じ変数で z 一つ） | Rem（文献の読み） | — | theory §5 の訂正 | — |
| O9 | SO(2N) 版の trace 表示：g = ½Tr_{H₁^{cov}}(Ĩ·) + C。E は Witten の ζ 挿入のトレースの型、C は商の理論にしかない状態 | Prop | O1 | theory §11.2、草稿 prop:trace。Witten 1101.3216 §3.4・§5.1.4、AgSc 付録 A（原典で確認） | 新規（読み替え） |
| O10 | 境界ゼロモードへの Ĩ の作用：E の境界因子 z_q^{n−1}(z_q^+)^n、P は z_q^{2n−1}。系：E は p₂ 型の方向だけ | Prop ＋ CR（Hopf (2,2),(2,4)、T(2,3),(2,5)、framed unknot x=2,4） | KN の境界ゼロモードの描像、Ĩ が L と I(L) の境界を入れ替えること | theory §11.3、test_rk_boundary_mode（2 件）、logs/b_boundary_mode_rule.log、b_framed_unknot.log | 新規 |
| O11 | refined の強い整数性：全 d の f が整数 ⟺ E ≡ P_inv (mod 2) | Prop | O1 | theory §11.7 (c) | 新規 |
| O12 | σ は SU(2)_ℓ（Lefschetz）の Weyl 反射。凍結した状態には Lefschetz が無いので C は σ 対称でなくてよい。[1] の σ 対称性は未説明 | Rem ＋ Conj | — | theory §11.4 | 新規（観察） |

## 3. unknot と Hopf 絡み目（§4）

| ID | 主張 | 種別 | 前提 | 根拠 | 新規性 |
|---|---|---|---|---|---|
| H1 | unknot：f_{10} = −1 − (ν − ν^{−1})/(t^{1/2} − t^{−1/2})、f_{a0} = 0（a ≥ 2）。Mariño 型の unknot とクロスキャップ −1 の和 | Prop | 双対 Cauchy 核 | report.md §4.1、test_rk_hopf::test_unknot_and_hopf_fundamental | 新規（refined SO の開弦の unknot） |
| H2 | Hopf：f_{ab} ∈ (1 − q)/(1 − t)·Z[q^±, t^±, U^±]（(偶,偶) は /(1 − t²)）、a, b ≤ 4。Schur 基底の f_{R,b}（\|R\| ≤ 4、b ≤ 2）でも同じ | CR | — | test_rk_hopf::test_symbolic_matches_finite_N, test_multi_probe_consistency | 新規 |
| H3 | f_{11} = −f_{10} z_q(ν~ − ν~^{−1})、ν~ = ν(q/t)^{1/2}、f_{1b} = f_{11} χ_{b−1}(U/t) | CR（b ≤ 4） | — | 同上 | 新規 |
| H4 | g_{1b} = z_q(ν~ − ν~^{−1}) χ_{b−1}(ν~/q^{1/2})。t は ν~ だけを通して入る | CR（b ≤ 4） | — | test_rk_cover::test_hopf_crosscap_fundamental_sector | 新規 |
| H5 | ĝ_{R,[1]}（\|R\| ≤ 4）：hook (n−d,1^d) で (−1)^d z_q(ν~ − ν~^{−1})ν~^{2d−n+1}、それ以外で 0。t = q で Mariño の ĝ_{[2],[1]}, ĝ_{[1,1],[1]} に R ↔ R^T と符号で一致 | CR | — | test_rk_cover::test_hopf_hatted_crosscap_invariants | 新規 |
| H6 | E_22, E_24 は (i) t = −1 と (ii) t = q の条件で一意に決まり、積の形。(4,4) は一意でない | CR | 固定セクターの部分和の仮定 | test_rk_cover::test_E22_unique_from_t_minus1_and_t_q | 新規 |
| H7 | Ĩ は「固定セクターの入れ替え＋Adams」として働く：E_{R,S} = Σ c^R_{2;ρ} c^S_{2;σ} ε_{ρσ}（\|R\| ≤ 4、\|S\| ≤ 2 の 10 セクター、4×4 で t = −1 の 25 組のうち 4 組から解き残り 21 組で一致 組） | CR | — | test_rk_cover::test_multi_probe_cover_and_adams_b2、results の 2 件 | 新規 |
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
| D3 | GL の単項式の規則 M_b = q^{−n(b,b)} t^{−2n(b,ρ)+(m+1)\|b\|(N−1)}（n = 2m+1） | CR（m = 1：[1],[2],[1,1] で読み [3],[2,1],[1,1,1] と随伴で予言。m = 2：N = 2..4 の 9 例と随伴） | — | logs/refined_kauffman/c_gl_mono_m2.log、test_rk_daha::test_gl_daha_jones_monomial_rule_and_adjoint | 新規 |
| D4 | 一般の N の補間：trefoil の (1),(2),(1,1),(3),(2,1),(1,1,1)、T(2,5) の (1),(2),(1,1)、被覆の随伴と composite (2;1),(11;1)。補間に使わない N で照合 | CR | 補間の窓（s の次数の傾き） | test_rk_daha::test_interpolation_vector_reproduces_as、results の test_general_N_three_box_tq_and_T25_universal_factor | 新規（refined SO(2N) の多箱の一般の N） |
| D5 | 大きな N の 2 箱：g_(2), g_(1,1) の分母は 2·単項式·(1+t) だけ。g_(2)+g_(1,1) と (1+t)(g_(2)−g_(1,1)) は整数 Laurent。U = t², t³ で固定 N と一致 | CR（T(2,3), T(2,5)） | — | test_rk_daha::test_large_N_two_box_crosscap、results の test_large_N_two_box_crosscap_structure | 新規 |
| D6 | 大きな N の 3 箱（trefoil）：g_R は整数 Laurent で (1+t) の極なし（E = 0） | CR | — | results の test_large_N_three_box_crosscap | 新規 |
| D7 | t = q で Mariño との辞書：W は Mariño の標準 framing の鏡像に比 1、ĝ は奇数 d に直して −鏡像(ĝ^Mar)。0904.1088 (4.8) の ĝ_[2] の z⁹ の誤植を訂正 | CR | — | test_rk_daha::test_large_N_two_box_marino_dictionary、test_rk_unrefined::test_trefoil_marino_48 | 誤植の指摘は新規 |
| D8 | KN の q 枠と t 枠の hatted は 2 箱で一致 | CR | — | test_rk_daha::test_large_N_two_box_kn_frames | 新規 |
| D9 | 普遍因子：ĝ − ĝ\|_{t=q}(ν → ν w^k) が (q − t)(ν² − t) で割り切れるのは k = +1 だけ（[1]、2 箱、3 箱、T(2,5) の 2 箱） | CR（k = −3..3） | — | test_rk_closed_orientifold::test_shift_direction_unique_for_trefoil_one_box、logs/refined_kauffman/b_shift_direction.log、results の test_general_N_three_box_tq_and_T25_universal_factor | 新規 |
| D10 | 多箱の c = 1 部分は σ 対称でなく G 展開できない。SO の R_λ に super-duality に当たる対称性が無い | CR | — | report.md §9、scripts/rk_b_2box_structure.py | 新規（否定的） |
| D11 | C_N の DAHA の実装：関係式、u = 1 で D_N と一致、Cherednik の ℋ^C_{3,2}(ω₁) と一致。ℋ^C_{5,2}(ω₁) は新規 | CR | Cherednik の C_n の規約 | test_rk_daha_cn、logs/d_cn_check.log、d_hyper_fit.log | 実装と T(5,2) は新規 |
| D12 | refined Sp(2N) = C_N(u = t²)（AM の定義）の [1]：g^{Sp}(𝐚) = −g^{SO}(−𝐚)（T(2,3), T(2,5)）。c = 1 の符号反転が refined でも成立 | CR | SO の規約の較正、Sp の単項式の決め方 | test_rk_daha_cn::test_sp_one_box_is_minus_so_with_a_flipped、logs/d_sp_one_box.log、theory §11.12 | 新規（refined Sp の結び目不変量の最初の値） |
| D13 | refined Sp(2N) の trefoil の 2 箱（C_N の P_{2ε₁}, P_{ε₁+ε₂}、評価公式、一般の N）：整数性は SO と同じ全 d の構造。E^{Sp}(t=−1) = E^{SO}(t=−1)。hatted の p₁² 成分で Sp(𝐚) = +SO(−𝐚)、p₂ 成分と g_R そのものには単純な符号規則が無い | CR | C_N の実装、T_λ^{−2n} の規格化（[1] で較正済み） | test_rk_sp_two_box、logs/d_sp_two_box_compare.log、theory §11.14 | 新規 |
| D14 | SO と Sp の 2 箱の差：f^{Sp}_(2)(𝐚) − f^{SO}_(2)(−𝐚) = ψ²[δ]、(1,1) は −ψ²[δ]、δ = even_𝐚(g^{SO}_1) + 1 = (𝐚²−𝐪²)(𝐚²𝐭²−1)(𝐚²𝐪²𝐭²−1)/(𝐚⁶𝐪²𝐭⁴)。W の段階では W_λ^{Sp}(𝐚) = (−1)^{|λ|}W_λ^{SO}(−𝐚)、冪和の係数は p₂ の定数項だけ −2 違う（→ 定数 1）、Adams 項 ½ψ²[f₁] が規則と可換でない（→ ψ²[even g₁]） | CR | trefoil の m=1 | test_rk_sp_so_adams、theory §11.15 | 新規 |

## 6. 正値性（§7）

| ID | 主張 | 種別 | 前提 | 根拠 | 新規性 |
|---|---|---|---|---|---|
| P1 | [1]（T(2,2m+1)、m ≤ 5）は判定式 T6 を満たす | CR | T1, S7 | T3 | 新規 |
| P2 | 多箱の普遍因子で割った商の符号は大部分のブロックで混在（2 箱の和 5/10、T(2,5) 6/18、3 箱 2/15, 3/14, 2/13）。[1] の商は T(2,3)（3/3）、T(2,5)（5/5）でそろう | CR | — | progress.md 続き 28〜33 | 新規（否定的） |
| P3 | 基底の取り替え・二次量子化の指数・被覆の brane/image の分割でも符号はそろわない。被覆の brane–brane は KN の正値性を満たし、brane–image は (𝐭² − 1) を含んで混在 | CR | — | progress.md 続き 27〜31 | 新規 |
| P4 | 固定セクターの整数計画で E を決める試み：状態の部分和は trefoil で解なし（T(2,5) は未実行）。捻ったトレースは T(2,5) で解なし（指数の上限では少なくとも 2 組の隠れた対が要る）。結び目の E の境界因子は和 z_q^+ | CR | 捻りの規則（Hyp） | test_rk_twisted_multiplet（results）、theory §9 | 新規 |
| P5 | 被覆の随伴 R_adj はどの偶奇の符号規則でも正にならない（Cherednik–Elliot の「composite の DAHA 超多項式は正値性をもたない」と一致） | CR ＋ K | — | scripts/rk_b_adj_positivity.py、1503.01441 序論 | 照合 |
| P6 | 多箱の正値性は、凍結した状態の多重項の基底と E・C の分離が決まるまで判定できない | Rem（現状の結論） | O3, P2〜P5 | — | — |
| P7 | framed unknot：(i) と「t = q で C が整数」だけで x = 1, 2, 3 の E が一意（E(x) = (−1)^x(ν²+1)(q+1)Σ(−1)^k(ν²q)^{x−1−k}t^k）。x = 4 は (ii) の形 A で一意 | CR | — | logs/b_framed_unknot.log、test_rk_framed_unknot（x = 1 全数、x = 2 は results）、theory §11.8・§11.9 | 新規 |
| P8 | Ĩ の符号は電荷だけの関数でない：unknot・Hopf (2,2) で成り立つ (−1)^{(j+l)/2} の規則は Hopf (2,4) と結び目で (i)(ii) を満たさない。E の決定には t 電荷の J₂/S_R への分解が要る | CR ＋ Rem | — | logs/b_sign_rule.log、theory §11.8 | 新規（否定的） |
| P9 | [1] の c = 1 の ν^β ブロックの G の次数は brane–image の ν^{2β} ブロックの半分（全 β、T(2,3), T(2,5)）。[1] の凍結した状態は brane–image の環状の弦の商で、商の取っ手が G を担う | CR ＋ Conj | — | logs/b_one_box_vs_brane_image.log、theory §11.9 | 新規 |
| P10 | 多箱の c=1 の σ 非対称性の分解：(i) ν の偶奇で定義した c=1 には t=q で既に c=2（z 奇）が混ざる（p₂ 方向・(2,1)、閉じた形）、hatted では c=1 ↔ ν^β、β ≡ |R|−1 (mod 2)（[1] 偶、2 箱 奇、3 箱 偶）。冪和では p₂ を含む方向だけ逆。3 箱の ĝ_λ は U 一般で整数 Laurent。(ii) 正しい c=1 の σ 奇部分は 𝐪²(𝐭²−1)(V⁴−1) で割り切れ t=q で消える。σ 偶・奇の分解は半整数で、trefoil では σ 奇側の G 展開が全ブロックでそろう（交互）が T(2,5) では崩れる | CR（解釈は予想） | trefoil・T(2,5) の 2 箱、trefoil の 3 箱 | theory §11.16、b_multibox_sigma*.log | 新規 |
| P11 | ORS のセル分解（Γ=⟨2,3⟩）の実装が T(2,3) の DGR の unreduced superpolynomial を再現（l ≤ 8）。trefoil の固定セクターの指数 P = (1−t)x = −(𝐚²−1)(𝐚²−𝐪²)(𝐪²𝐭²−1)(𝐚²𝐪²𝐭²−1)N₈/(𝐚¹²𝐪⁴𝐭¹⁰)（framed unknot と同じ境界因子）。P は被覆の [1] の ψ² ではない | CR | T(2,3) | theory §11.18、test_rk_planB_cells | 再現＋新規 |
| P12 | Krefl–Walcher の「β=2 の理論 = β=1 の orientifold」の結び目版：unrefined の Kauffman [1] は refined U(N′) の [1] の β = 1/2, 1, 2 への特殊化と、Wilson loop でも自由エネルギーの形でも一致しない（a = λ^k q^j、単項式倍の範囲） | CR（否定的） | trefoil、T(2,5) | theory §11.19、e_krefl_walcher.log | 新規（否定的） |

## 6′. E の決定：n = 1 の原理と n = 2（§7 の続き、2026-10-28〜31 の続き 75〜98）

| ID | 主張 | 種別 | 前提 | 根拠 | 新規性 |
|---|---|---|---|---|---|
| P17 | 符号の規則（下の系列 k ≤ 2m で (−1)^q、上で (−1)^t）は T(2,3)・T(2,5) では (i)+(ii) で E を一意に決めるが、T(2,7) では (i) だけでも解がない。破れるのは上の系列の最初のブロック k = 2m+1 | CR（否定的） | — | theory §12.10〜§12.11、続き 88 | 新規（否定的） |
| P18 | n = 1 の E は状態の空間の同定から条件なしに作れる：N のブロックは下の系列 Sym^k V（V = 1 + t + qt、k < 2m）、移り目 Sym^{2m}V − S ⊗ S（S は簡約 superpolynomial）、上の系列 −t^{K−2j}[Sym^j V + (q − t)R_j]。Ĩ は V の対合 diag(1, 1, −1)（上は双対）、S ⊗ S の入れ替え（トレース ψ²S）。T(2,3)・T(2,5) の一意な E を全単項式で再現し、T(2,7)〜T(2,15) の E は予言として (i)・(ii) を満たす | CR（m ≤ 7）＋ Conj | — | theory §12.23・§12.29〜§12.31、test_principle_reproduces_unique_E、twisted_E_n1_m{3,...,7}.pkl | 新規 |
| P19 | 補正の形：R_{m+j′} = t^{m−1−j′} F_{j′}、F_{j′} は m によらず（j′ ≤ 2 で m = 4, 5）、母関数 Σ F_{j′} z^{j′} = [1 + (2t + t² + qt + 2qt²)z + qt³z²]/[(1 − z)(1 − qtz)(1 − t²z)(1 − q²t²z)]（F₀〜F₅ を再現、F₄ = 165 状態は T(2,13)、F₅ = 286 状態は T(2,15) で予言が的中、状態の数 C(2j′+3, 3)）。訂正：続き 89 の R_{m+1} = t^{m−2}(S′² + t²)（S′ = T(2,2m−3) の S）は T(2,9) で破れ、S′ → V | CR（j′ ≤ 5）＋ Conj | — | theory §12.29〜§12.31、correction_F、test_correction_generating_function | 新規 |
| P20 | 補正の Ĩ：E の補正は (−1)^{m+j′+1} t^{m−1−j′} ψ²(X_{j′})、X_{j′} = Sym^{j′}V + t·Sym^{j′−1}V（母関数 (1 + tz)/((1 − z)(1 − tz)(1 − qtz))、(j′+1)² 状態）、F_{j′} の残りは入れ替わる対。j′ = 2 は T(2,9) で (i) の行と台により 32 通りの変種から一意、j′ = 3 は T(2,11)、j′ = 4 は T(2,13)、j′ = 5 は T(2,15) で (i)+(ii) と整合 | CR（j′ ≤ 5）＋ Conj | P18, P19 | theory §12.29 追記・§12.31、diag_X、test_diag_X | 新規 |
| P21 | (q − t) の補正の対の符号は条件 (ii)（t = q）に見えず、(i) だけでは決まらない（n = 1 の j′ = 2 で 2/15、n = 2 の λ = (2) で 3/17。深い段の行を加えても同じ）。決めるには状態の空間の構造（P20）が要る | Prop ＋ CR | O1 | theory §12.27〜§12.28 | 新規 |
| P22 | n = 2（trefoil）：N_λ = Wilson loop なしの背景（閉じた形 −r⁷s³(q − t)ν⁶/[(A − t)(Aq − 1)(Aq − t)(Aq − t²)(Aq² − t)(A²q² − t²)]、真空 (q − t) と 6 セルの Fock 空間、塔なし）＋ 結び目の部分（λ = (2) の ν^{−10} から −(q − t)(1 + qt) ⊗ (1 + qt)）。上端〜ν^{−12} の Ĩ は「背景のセルに対角、平方に全因子の入れ替え」で (i)+(ii) と整合し λ = (1,1) は一意。λ = (2) の補正の Ĩ は未決定（n = 1 の型「対角は単項式 × ψ² の像」はどの平方の分け方でも載らない） | CR ＋ 未解決 | — | theory §12.26〜§12.28、test_rk_n2_rterms | 新規 |
| P23 | n = 1 の N と E の閉じた式：ブロックは普遍的なデータ V、F_{j′}（母関数）、X_{j′} と、結び目のデータ S_m・m だけで書ける（n_k = Sym^k V、n_{2m} = Sym^{2m}V − S_m²、n_{K−j} = −t^{K−2j}[Sym^j V + (q − t)t^{m−1−j′}F_{j′}]、E は対応する Ĩ のトレース）。T(2,3)〜T(2,15) で N は全係数で成立し E は (i)・(ii) を満たす | Conj ＋ CR（m ≤ 7） | P18〜P20 | theory §12.31 末尾、principle(m) の assert | 新規 |
| P24 | n = 1 の母関数の plethystic log：X_{j′} の母関数は有限（V ⊕ t の Fock 空間に関係 t²z²、PL = (1 + 2t + qt)z − t²z²）。F_{j′} の母関数と N・E の上の系列の母関数は無限で、上の系列の生成元は V z と (q − t) t^{m−1−j′} c_{j′} z^{m+j′}（c₀ = 1、c₁ = V² − t(1 + qt)、c₂ から負の係数）。上の系列は正の BPS 生成元の Fock 空間としては読めない | CR | P19, P20, P23 | theory §12.31 末尾、rk_f_plog.py、test_rk_plog | 新規 |
| P25 | n = 1 の原理は 2 本撚りに固有：T(3,4)（γ = τ₊τ₋³）では K = 9 = rs − 3、下の系列 Sym^k V（k ≤ 5）、移り目 Sym⁶V − S̃²（S̃ = (1 + qt)(1 + t + t² + qt² + q²t²)、10 状態、簡約 superpolynomial（11 状態）ではない）、上の系列は −t^{K−2j}[Sym^j V + (q − t)R_j] の形にならない（j = 0 で −t⁸(q + t)）。下の系列を固定した整数計画は (i)+(ii) で解あり（対 25、一意でない単項式 29/111） | CR | P18〜P20, P23 | theory §12.32、progress 続き 102、rk_f_twisted_tensor.py（blocks_raw・n1_upper(k, kt)）、test_rk_t34_blocks・test_rk_torus_tag | 新規 |
| P26 | n = 1 の N の母関数：x = ν^{−2}、[1] の鏡像の行 R(y) = Σρ_j y^j、随伴の行 A(y) = Σα_j y^j とすると Σ n_k x^k = [1 − x^{2δ}R(qx)²]/((1−x)(1−tx)(1−qtx)) + x^{2δ−1}G(qx)/(1−tx)、A = R² + (q−1)(t−1)G（導出：N の式の x 展開。Laurent ⟺ G が Laurent）。移り目の平方の因子は S̃ = ρ₀ − ρ₁（[1] の最上 2 行）で、g₀ = 0・g₁ = (2ρ₀ρ₁ − ρ₁²)/q と同値。T(2,3)〜T(2,15)・T(3,4)・T(3,5)・T(3,7)・T(3,8)・T(3,10)・T(4,5)・T(4,7) の 14 結び目で成立し、T(3,5)（τ₊τ₋τ₊τ₋）・T(3,7)・T(3,8)・T(3,10)・T(4,5)・T(4,7) は予言として的中（随伴と [1] の被覆データだけで N が決まる）。2 本撚りの S = 簡約 superpolynomial は行が 2 つだからで、T(3,4) の S̃ = S − t³ は 3 行目の切り捨て。最深ブロック n_K = −t^{K−s+2}h_{s−2}(q,t)（s = 2, 3, 4） | Prop（母関数）＋ CR（14 結び目）＋ Conj（g₀, g₁、n_K） | — | theory §12.33〜§12.34、rk_f_rows_gen.py、daha_dn.GAMMA_WORDS（τ の一般の語）、test_rk_rows_gen | 新規 |
| P27 | 3 本撚りの上の系列：U_j = −n_{K−j}/t^{K−2j} = (1 + q/t)Sym^jV + (q − t)C_j が j ≤ δ − 2 で成り立ち（C_j は結び目に依る）、j = δ − 1 で t = q の完全平方 −([⌊r/2⌋]_{t²}[r′]_t)²/t^{δ−1}（r′ は r 以下の最大の奇数。T(3,4)・T(3,5)・T(3,7)・T(3,8)・T(3,10)）。4 本撚り（T(4,5)・T(4,7)）では前因子 h₂/t²、平方の位置 j = r − s + 1、因子 ρ₀(T(s−1, r−s+1))·S(T(2,r′))|_{t=q}。refined の平方の因子は未同定 | CR ＋ 観察 | P26 | theory §12.34、f_rows_gen_upper*.log | 新規 |
| P28 | 移り目の Ĩ の規則 E_{2δ} = Sym^{2δ}V（セルの符号）− ψ²S̃ は 3 本撚りにも延びる：T(3,4)・T(3,5) で下の系列と移り目を固定し上を自由にした整数計画が (i)+(ii) に解をもつ（上の一意な単項式 22/31、24/78。最深 2 ブロックは一意、E_K = −t⁸(q − t)、+t¹¹(q − t)）。E の (R, G) の読み：E = F^ε(1 − x^{2δ}ψ²R) + x^{2δ−1}E_G、2 本撚りでは G の行の Ĩ は j ≥ δ+2 で q の規則（Ĝ_j = g_j(−q,t)(−q)^j）、j ≤ δ+1 で対の補正 D_{δ−i} = t^{a(m,i)}D̃_i（D̃_i は i と m の偶奇だけで決まる）。3 本撚りでは最上の行だけ q の規則（符号は結び目に依る） | CR ＋ Conj | P18〜P20, P26 | theory §12.34〜§12.35、rk_f_rows_E.py、test_rows_E_transition_rule_t25 | 新規 |
| P29 | composite の大きな N の因子化：R_{[α;β]} = R_α R_β + (q − 1)(t − 1) G_{αβ}（最上行は積、差は全行で (q−1)(t−1) で割れる）が trefoil の [2;1]・[1,1;1]・[2;2]・[2;1,1]・[1,1;1,1] で成り立ち、随伴では 14 結び目。G の行は α, β の次元型の因子 (1+q)・(1+t) を含み、最深行は (q−1)(t−1) × 次元型の因子。読み：(q−1)(t−1) は brane–image の環状の弦の refined の測度で、N の整数性はその連結相関の LMOV 型の整数性と同値。g₁ の式の一般化（ρ₀ρ₀ − S̃S̃）は他の composite で不成立 | CR（trefoil の 6 例）＋ Conj（一般）＋ 否定的な結果 | P26 | theory §12.36、rk_f_comp_factor.py、test_rk_comp_factor | 新規（論文には未反映） |
| P30 | n = 2（trefoil）の行の母関数（§12.40）での G の役割：λ = (1,1) の ν^{−12} の結び目の部分は K₁₂ = −(s³/r)[Y² + (q−t)(t/q⁴)a₂²]、Y = X + Z = (1+q)(1+t)(1+qt)/q³ で、G の項（＋結び目の Adams）がちょうど −X² を −Y² に補完し、(q−t) の補正は記号ブロックの a₂² の項そのもの（n = 1 の −ρ₀² + q g₁ = −S̃² の n = 2 版。§12.27 の形の導出）。ν^{−14} の最上行×次の行の交差項は両 λ で 2X₀X₁（行は級数として入る）。λ = (2) の ν^{−12} は行だけでは平方と補正に分かれない（c 単項式・[5]_q の refinement の探索で解なし）。ν^{−14} の最上行の 2 次形式は階数 3 | 導出（(1,1)、交差項）＋ CR（否定的、(2)） | P26, P29 | theory §12.41、rk_f_n2_grows.py、test_rk_n2_grows | 新規（論文には未反映） |
| P31 | n = 2 の N_λ の λ 依存性：2 箱の行と G の重みの付け替え σ（b → q b、c → c/t、G_{αβ} → w_αw_β G）で N_(2)∘σ − N_(1,1) = (q−t) × 随伴セクター（R₁(ν)²κ_R + G_adj(ν)κ_G + 背景）が ν^{−20} まで成立。σ 共変な部分は両 λ の ν^{−12} で完全平方で、Y_(2) = [(1+q)(1+t)(1+qt)² + qt(q−t)]/q⁴（§12.27 の 18 通りの曖昧さの解消）、補正は (q−t)[(t/q⁴)a₂² + δ]、δ = [tV′²a₂²/q³ − 2a₂S̃ + (V′−t)gadj6]/q³。κ_G の係数は階段型 T_n = Σ_{0≤j≤i≤n} q^i t^j × U_n | 導出（記号ブロック、ν^{−20} まで）＋ 観察（κ_G の形）＋ 予想（読み） | P30 | theory §12.42、rk_f_n2_grows.py part F、test_rk_n2_grows | 新規（論文には未反映） |
| P32 | E への適用：Y_(2) で λ = (2) の ν^{−12} は K₁₂ = −r^{−17}s³[Ŷ² + (q−t)R]、R は 15 状態（積型の Y の最小 17 より小さい）。平方の部分の Ĩ は +ψ²Ŷ（テンソル平方の入れ替え）だけが可能、補正の対の最小は 11（Y₂ では 20）。補正を行の部分に分けると (B) ψ²・(C) 入れ替え（0）は共通、(A) は ψ²/±恒等、G_adj の部分 (D) は (i)+(ii) で決まらない | CR（整数計画）＋ 否定的（(D)） | P31, P17〜P20 | theory §12.43、rk_f_n2_E_Y2.py、test_rk_n2_E_Y2 | 新規（論文には未反映） |
| P33 | G_adj の行の Ĩ：n = 1 の trefoil で Ĝ₁ = −ψ²(ρ₁)（g₁ = (2ρ₀ρ₁ − ρ₁²)/q の交差項が入れ替えで消える）、Ĝ₂ = −t²(q−1)(t+1)。n = 2 の λ = (2) の ν^{−12} には「平方・行の積すべてに入れ替え、随伴の塔の核 V′−t に F^ε 型の符号 (1+q)(1−t)」の一様な E が存在し (i)+(ii) を満たす（核の符号 16 通り中 1 通り）。この E は §12.27 の族（補正の単項式ごとに |e| ≤ |R|）の外にある | CR（整数計画）＋ Conj（入れ替えの原理） | P32, P18〜P20 | theory §12.44、rk_f_n2_EG_adj.py、test_rk_n2_EG_adj | 新規（論文には未反映） |
| P34 | λ = (1,1) の ν^{−12} の E の余分な (q−t)qt² は背景の E（B₁₂(1,1) = −r^{−13}s³(q−t)qt²、E_B = −B）：E₁₂(1,1) = r^{−13}s³[ψ²Ŷ + (q−t)tψ²(1+qt)] − B₁₂ で、結び目の部分は平方も補正も入れ替え。λ = (2) の ν^{−12} も全部入れ替えで符号まで一意（σ_Y = σ_A = σ_B = +1、核 (1+q)(1−t)）。ブロックの符号は σ = (−1)^{k/2+1}（n = 1 の移り目、(2) の ν^{−10}、両 λ の ν^{−12}） | 導出（(1,1)）＋ CR（(2)）＋ 観察（符号） | P33, P17〜P20 | theory §12.45、test_rk_n2_E11_bg | 新規（論文には未反映） |
| P35 | ν^{−14} の E は (i)+(ii) からは 2〜3 単項式しか決まらない（下のブロックを固定、113 未知数、対の最小 48 → 936）。brane/image を区別した記号ブロックで、ν^{−12} の平方 −X² は文字どおり −(s³/r)X^b X^i（brane の X ⊗ image の X、両 λ で 9 単項式が一致）。ν^{−14} の対角の単項式の係数は brane × image の核の積ではない | CR（否定的）＋ 導出（X^b ⊗ X^i） | P30, P34 | theory §12.46、rk_f_n2_E14.py、rk_f_n2_rows_lab.py、test_rk_n2_rows_lab | 新規（論文には未反映） |
| P36 | 核（次元・composite の次元・冪和の係数・φ）まで brane/image に分けた記号ブロックで、G を含まない対角の単項式に ψ² を掛けたトレースは λ によらず Laurent でもなく E ではない（否定的）。E の λ 依存性は Schur 基底でしか見えず、入れ替えのトレースは（brane の λ の対象 Ω_λ）の ψ² で、ν^{−12} の Ŷ_λ がその Ω_λ。原理 E_λ = ψ²(Ω_λ)、N_λ の積の部分 = Ω_λ² は予想 | CR（否定的）＋ 導出（理由）＋ Conj | P34, P35 | theory §12.47、rk_f_n2_rows_ker.py | 新規（論文には未反映） |
| P37 | n = 2 の混合キュムラントに純粋な brane×image の積は無い（共有の核も記号にした記号ブロックで、純粋な積の双線形形式は全ブロックで階数 0）。ν^{−12} の平方 Ŷ_λ² は連結の核（C_{αβ}・c₁・κ）を係数にもつ行の積が担う。σ 共変な部分は ν^{−14} で Ŷ_λ で割れず、級数の平方根 Ω_λ(ν) は得られない | CR（否定的）＋ 観察 | P35, P36 | theory §12.48、rk_f_n2_rows_ker.py --shared、rk_f_n2_omega.py | 新規（論文には未反映） |
| D15 | 2 箱の R_λ(U) と随伴 R_adj(A^{1/2}) は、r の次数ごとに U の冪の範囲を絞る補間（flint の約分）で N ≤ 8〜11 から一意に決まり、使わない N（≤ 20〜34）の全部で照合、Newton 補間の保存値と全係数一致。T(2,9)〜T(2,13) の 2 箱のデータ | CR | 帯の仮定（rk_c_bounded_interp） | rk_c_bounded_save.py、logs/c_bounded_save_*.log、続き 94〜97 | 方法 |

## 7. 既知（成果として書かない）

- Mariño 0904.1088 の構造と trefoil の数値、NRZ 1310.2240 の 4_1 の表（K0 の再現：test_rk_unrefined）。
- AgSc 1202.4456 の T(2,2m+1) の [1] の一般式、Hopf、大きな N の閉弦（§6）。
- AM 2205.12832 の閉弦の GV 型展開。
- Cherednik の DAHA-Jones の定義。Cherednik–Danilenko の D 型の t = q。
- KN 1703.05408 の M(x)、q 枠と t 枠。
| P13 | 凍結した状態 ↔ 被覆の Ĩ 不変な曲線（β ↔ 2β）。framed unknot（𝖿 = 1..4）の 2 箱の p₁² 方向の凍結した状態は結び目の [1] と同じ判定式を満たす（c=1 は σ 対称で G 展開の符号が同符号か (𝐭²−1)×同符号、c=2 は −(𝐚²−𝐪²)(𝐭²−1)×正）。Hopf では両端の β で P_{2β} = 単項式 × C_β²。結び目の多箱の σ 非対称部分は Lefschetz の重み m ごとに (𝐭²−1) × 整数多項式で、trefoil の 2 箱では全ブロックが基底 {G^g w^j} ∪ {(𝐭²−1) q_ℓ^m w^j} で同符号（β で交互）、3 箱では非対称部分だけ同符号、T(2,5) の中間ブロックは混在。非対称部分は t の有限の塔を含む | CR（解釈は Conj） | trefoil・T(2,5) の 2 箱、trefoil の 3 箱、framed unknot、Hopf | theory §12、f_frozen_basis.log、f_hopf_sqrt.log、test_rk_frozen_basis | 新規 |
| P14 | 有限個の z₂ モードをもつ Ĩ 不変な状態は、指数では E（J₂ の塔：½ε[L]_{−t}）か凍結した状態と区別できない形（S_R の塔の対）に吸収され、新しい型にならず、p₁² 方向には寄与しない。p₁² 方向の凍結した状態の σ 非対称部分は (1 − q/t) × t の指標 t^k[L]_t の和（trefoil 全列、T(2,5) 26/31 列）で、最低のブロックは −(1 − q/t) w^{−(8m+1)} t([2n]_t + [2n−1]_t)。T(2,5) の混在列は符号の変化 1 回 | Der ＋ CR | trefoil・T(2,5) の 2 箱 | theory §12.5、f_frozen_towers.log、test_rk_frozen_basis | 新規 |
| P15 | trefoil の被覆の brane 2・image 2 の固定セクター P（大きな N、composite [2;2], [2;11], [11;11] の一般の N から）は (ν²−1)(ν²−t)(ν²q−1) z_q³ で割り切れ、4P/z_q³ の ν のブロックはすべて σ 対称で G 展開できる（符号は中間で混在）。凍結した状態の σ 非対称な t の塔は被覆の側に無く、§12.5 の予想は素朴な形では支持されない。Hopf の両端の平方根の構造は 2 箱に持ち上がらない。G の次数の 2 倍の関係は 6 ブロック中 5。最低のブロックは −s^{20}r^{−41}Σ_{i,j<12} q^{max}t^{min} | CR（否定的を含む） | trefoil | theory §12.6、f_cover_fixed22.log、test_rk_daha_fast_lp | 新規 |
| P16 | 被覆の固定セクター（brane n・image n、n = 1, 2）への Ĩ の作用の分け方（E = P − 2Q 型の整数計画、t = q の条件 E = −(1+q)ψ²[F]）：n = 1 は (i)+(ii) で一意（§9 の再現）、n = 2 は (ii) だけでは決まらない（最低のブロックで 78 中 74 が自由）。refined の ψ² 型の分け方（t の塔 → t² の刻み）は、t = q で ψ² に戻る 9 通りの写像と単項式のずれのどれでも n = 1 で 7 中 2〜3、n = 2 で 12 中 0〜1 のブロックでしか条件を満たさない（否定的） | CR（否定的） | trefoil | theory §12.7、f_twisted22.log、test_rk_frozen_basis | 新規 |
