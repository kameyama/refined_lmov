"""多箱の c=1 の σ 非対称性（theory §11.16、trefoil の 2 箱、largeN_2box_m1.pkl）：
(1) t=q で p₁² 方向の ν 偶は q↔1/q 対称、(2) p₂ 方向の ν 偶の反対称部分は閉じた形 −2zδ(q¹²−1)/((q−1)ν⁶q^{11/2})、
(3) refined で hatted の ν 奇（c=1）の σ 非対称部分は (q−t)(qt−1) で割り切れる。"""
import sys
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rk_b_multibox_sigma import load, sigma, c1, r_, s_, nu


def test_two_box_sigma_structure():
    g, _ = load(1)
    g2, g11 = g["(2,)"], g["(1, 1)"]
    c2 = lambda e: sp.cancel(e - c1(e))
    x = sp.cancel(c1(g2 + g11).subs(s_, r_))
    assert sp.cancel(x - x.subs(r_, 1 / r_)) == 0
    x = sp.cancel(c1(g2 - g11).subs(s_, r_))
    anti = sp.cancel(x - x.subs(r_, 1 / r_))
    q = r_**2
    delta = (nu**2 - q) * (nu**2 * q - 1) * (nu**2 * q - q) / (nu**6 * q**2)
    closed = -2 * (r_ - 1 / r_) * delta * (q**12 - 1) / ((q - 1) * nu**6 * q**sp.Rational(11, 2))
    assert sp.cancel(anti - closed) == 0
    X = c2((g2 + g11) / (r_ - 1 / r_))
    a = sp.cancel(X - sigma(X))
    assert a != 0
    n_, _ = sp.fraction(a)
    P = sp.Poly(sp.expand(n_), r_, s_, nu)
    assert sp.rem(P, sp.Poly(r_**2 - s_**2, r_, s_, nu)).is_zero
    assert sp.rem(P, sp.Poly(r_**2 * s_**2 - 1, r_, s_, nu)).is_zero


def test_three_box_hatted_parity():
    """trefoil の 3 箱（Schur に直した hatted）：整数 Laurent、t=q で ν 偶が q↔1/q 対称・ν 奇が反対称（theory §11.17）。"""
    from rk_b_3box_hatted_norm import load3
    from rcsknots.bps import M_matrix, partitions
    u_ = sp.Symbol("u", positive=True)
    P, S = load3()
    parts = partitions(3)
    hat = M_matrix(3, 1 / r_).inv() * sp.Matrix([S[l] for l in parts])
    for h in hat:
        h = sp.cancel(h)
        n_, d_ = sp.fraction(h)
        assert sp.Poly(d_, r_, s_, u_).is_monomial and sp.Poly(d_, r_, s_, u_).coeffs()[0] == 1
        x = sp.cancel(h.subs(u_, nu * s_).subs(s_, r_))
        ev, od = c1(x), sp.cancel(x - c1(x))
        assert sp.cancel(ev - ev.subs(r_, 1 / r_)) == 0 and sp.cancel(od + od.subs(r_, 1 / r_)) == 0
