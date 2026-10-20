"""C_N の DAHA（src/rcsknots/daha_cn.py）と refined Sp(2N) の [1]（案 D、progress.md 続き 50）。
(1) C_3 の関係式、(2) u = 1 で D_N の DAHA-Jones に一致、(3) Cherednik 1111.6195 (super-c-2-3) の ℋ^C_{3,2}(ω_1) に一致、
(4) trefoil：g^{Sp}_[1](𝐚, 𝐪, 𝐭) = −g^{SO}_[1](−𝐚, 𝐪, 𝐭)（c = 1 の符号の反転が refined でも成り立つ）。"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp


def test_cn_relations_and_u1_is_dn():
    import rk_d_cn_check as K
    assert K.relations(3)
    for m in (1, 2):
        assert sp.cancel(K.vs_DN(2, m) - 1) == 0


def test_cn_matches_cherednik_hyperpolynomial():
    import rk_d_cn_check as K
    res = K.vs_cherednik(3)
    assert res["そのまま"][0]


def test_sp_one_box_is_minus_so_with_a_flipped():
    import rk_d_hyper_fit as C
    import rk_d_sp_one_box as SPO
    r_, s_ = C.r_, C.s_
    H = SPO.load_H(1)
    gso, nu = C.g_known(1)
    Wc = C.cover_W1(1, nu)
    assert sp.cancel(SPO.W_group(1, H, "SO", nu) * (s_ / nu)**2 - Wc - gso) == 0
    gsp = sp.cancel(SPO.W_group(1, H, "Sp", nu) * s_**2 / nu**2 - Wc)
    A_, Q_, T_ = sp.symbols("A Q T")
    bold = lambda f: sp.expand(sp.cancel(f.subs({r_: -T_ * Q_, nu: A_}).subs(s_, Q_)))
    assert sp.expand(bold(gsp) + bold(gso).subs(A_, -A_)) == 0
