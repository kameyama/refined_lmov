"""境界ゼロモードへの Ĩ の作用の規則（theory_orientifold_index.md §11.3）：
Ĩ 不変な曲線の境界の対が n 組なら、E の境界因子は z_q^{n−1}(z_q^+)^n、挿入なしの P は z_q^{2n−1}。
Hopf の (2,2), (2,4) は n = 2（保存済みの E と固定セクターの f^{cov} で確認）。"""
import sys
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import rk_b_boundary_mode_rule as B


def test_hopf_E_has_one_zq_and_two_zq_plus():
    assert B.hopf() == {(2, 2): (1, 2), (2, 4): (1, 2)}


def test_hopf_P_has_three_zq_and_no_zq_plus():
    assert B.hopf_P() == {(2, 2): (3, 0), (2, 4): (3, 0)}
