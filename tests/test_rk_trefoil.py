"""refined Kauffman K3'：trefoil [1]（AS 1202.4456 §6）と refined の h/g 分解。"""
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rcsknots.kauffman_trefoil import (Z_V_T2, Z615, r, s, U, marino_F31, g_ref, is_laurent, W_SO, h_half_ref)  # noqa: E402


def test_general_formula_reproduces_615():
    assert sp.simplify(Z_V_T2(1) - Z615) == 0


def test_unrefined_is_mirror_of_marino_with_framing():
    """t = q：AS (6.15) = F_{3_1}(ν → 1/ν, z → −z) × ν^{−6}、ν = 1/𝐚。"""
    Zu = Z615.subs(r, s)
    FM = marino_F31(s / U, 1 / s - s)
    assert sp.cancel(Zu / FM - U**6 / s**6) == 0


def test_refined_rudolph_with_identified_monomial():
    """g^{ref}_[1] = W_SO − ½h^{ref} は単項式 r^{−5} s^{11} U^{−6} で Laurent 多項式（11 項、整数係数）、
    T_[1]^{−6} = r^{−6} s^{12} U^{−6} では Laurent にならない。"""
    g = g_ref()
    assert is_laurent(g)
    num, den = sp.fraction(g)
    P = sp.Poly(num, r, s, U)
    assert len(P.terms()) == 11 and all(c.is_integer for c in P.coeffs())
    assert not is_laurent(g_ref(r**-6 * s**12 * U**-6))
    # unrefined でも Laurent
    assert is_laurent(g.subs(r, s))


def test_as_615_matches_gn_at_N2_N3():
    """SO(4) = SU(2)×SU(2)、SO(6) = SU(4)：AS (6.15) は、A 型の Gorsky–Negut（楕円 Hall 代数）で計算した
    SU(2) の [1] の 2 乗（N = 2）、SU(4) の [1,1]（N = 3）の正規化した trefoil に、単項式 q^{−1}、q^{−1}t^{−4} 倍で一致する。"""
    from fractions import Fraction as F
    from rcsknots.gn import GN
    for rr, ss in ((F(3, 7), F(5, 11)), (F(2, 3), F(7, 13))):
        qq, tt = rr**2, ss**2
        R = lambda x: sp.Rational(x.numerator, x.denominator)
        G2 = GN(qq, tt, tt**2, lambda x: x)
        v2 = (G2.torus((1,), 2, 3) / G2.torus((1,), 1, 0))**2
        assert R(v2) == Z615.subs({r: R(rr), s: R(ss), U: R(ss)**4}) * R(rr)**-2
        G4 = GN(qq, tt, tt**4, lambda x: x)
        v4 = G4.torus((1, 1), 2, 3) / G4.torus((1, 1), 1, 0)
        assert R(v4) == Z615.subs({r: R(rr), s: R(ss), U: R(ss)**6}) * R(rr)**-2 * R(ss)**-8


def test_trefoil_fixed_N_pole_cancellation():
    """(b)：N を固定した refined SO(2N) の trefoil（1^k、k ≤ N）。Z~ の N 依存の極が打ち消すには、reduced 不変量が極の上で単項式（倍）でなければならない。
    SO(4) = SU(2)²：R_[2] は q t = 1（r = ±1/s）で s^{−4}；m_2 = s^4 で f_2 は整数 Laurent、m_2 = 1 では Laurent でない。
    SO(6) = SU(4)：R_[2,1,1] は q t³ = 1（r = ±s^{−3}）で s^{−16}。R は data/refined_kauffman/torus の鏡像（AS の向き）。"""
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "scripts"))
    import rk_k5_framing_scan as Fs
    from rcsknots.kauffman_knot import free_energy, r as r_, s as s_
    for sign in (1, -1):
        assert sp.expand(Fs.Rdata("2", 2).subs(r_, sign / s_) - s_**-4) == 0
        assert sp.expand(Fs.Rdata("211", 4).subs(r_, sign * s_**-3) - s_**-16) == 0
    parts = Fs.setup(2)
    mono = lambda e: sp.Poly(sp.fraction(sp.cancel(e))[1], r_, s_).is_monomial
    for m2, expect in ((s_**4, True), (sp.Integer(1), False)):
        W = [sp.Integer(1), parts[1][0], sp.cancel(m2 * parts[2][0])]
        f = free_energy(2, W, 2)
        assert mono(f[2]) == expect
