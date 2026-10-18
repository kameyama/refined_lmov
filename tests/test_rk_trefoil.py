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


def test_torus_family_framing_monomial():
    """A：T(2,5) の [1]。g^{ref} を Laurent にし、t = q で Kauffman 型（r → −1/r 不変）にする単項式は
    r^{−8} s^{18} U^{−10} = T_[1]^{−10} (q/t)^{δ/2}（δ = 2、T_[1] = q^{1/2}t^{N−1} = r U s^{−2}）。AS 1105.5117 の SU(N) の因子と同じ形。
    t = q の g^{ref} は z の偶数冪だけ（N^{c=2} = 0）。"""
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "scripts"))
    import rk_k3_torus_family as Tf
    assert Tf.kauffman_q_power(2) == 0
    Tm = r * U / s**2
    assert sp.simplify(r**-8 * s**18 * U**-10 - Tm**-10 * (r / s)**2) == 0
    g = Tf.g_family(2, 18)
    assert is_laurent(g)
    # d の偶奇を変えると Laurent でない
    assert not is_laurent(sp.cancel(Tf.Z_V_T2(2) * r**-8 * s**19 * U**-10 * Tf.unknot_SO - Tf.h_half(2)))
    Q, Tb, A = sp.symbols("Q T A")
    gq = sp.expand(sp.cancel(g.subs({r: -Tb * Q, U: A * Q}).subs(s, Q)).subs(Tb, -1))
    assert sp.expand(gq - gq.subs(Q, -1 / Q)) == 0


def test_torus_family_refined_only_part():
    """A：g^{ref}_[1](T(2,2m+1)) の 𝐚 の奇数冪（t = q で消える refined の部分）は
    (−1)^{m+1}(𝐚²−𝐪²)(𝐭²−1) 𝐚^{−(4m+1)}𝐪^{−(2m−1)}𝐭^{−3m} Π_m、Π_m = Σ_{j<m} 𝐚^{2j}𝐪^{2(m−1−j)}𝐭^{2j} Σ_{i≤j}(𝐪⁴𝐭²)^i（係数はすべて +1）。m = 2, 3。"""
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "scripts"))
    import rk_k3_torus_family as Tf
    for m in (2, 3):
        Pi, (Q, Tb, A) = Tf.refined_only_part(m)
        guess = sum(A**(2 * j) * Q**(2 * (m - 1 - j)) * Tb**(2 * j) * sum((Q**4 * Tb**2)**i for i in range(j + 1)) for j in range(m))
        assert sp.expand(Pi - guess) == 0


def test_refined_crosscap_c1_blocks():
    """B：T(2,5) の [1] の c = 1 部分（G = z_q z_t、W = (q/t)^{1/4} の多項式）を ν^β ごとに分けると、
    unrefined の台のブロックは係数が同符号（+, −, +）、台の外のブロック（β = −8）は (𝐭² − 1) × 正の多項式。"""
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "scripts"))
    import rk_k3_torus_family as Tf
    assert Tf.c1_blocks(2) == {-10: "coherent+", -8: "tfactor+", -6: "coherent−", -4: "coherent+"}


def test_refined_crosscap_special_points():
    """B の因子の意味：T(2,5) の g^{ref}_[1] は SO(2)（𝐚 = ±𝐪、ν² = t）で (−𝐭)^m = (q/t)^{δ/2}（U(1) の自明さ）、
    𝐭 → −𝐭 で (−1)^m 倍（だから t = q で消える部分は (𝐭² − 1) を因子にもつ）。"""
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "scripts"))
    import rk_k3_torus_family as Tf
    so2, par, (Q, Tb, A) = Tf.special_points(2)
    assert so2 == {Tb**2} and par == 1


def test_trefoil_fixed_N_two_probe():
    """固定 N の 2 箱（probe 2 個、Y の Schur |R| ≤ 2）。f_(1,1)（Sym²V、D_N の (2)）を Laurent にする framing 単項式は
    N = 2 で m_(2) = s^8 = m_(1,1)²（SO(4) = SU(2)² の積の形）、N = 3 で m_(2) = r^4 s^20（scripts/rk_k5_trefoil_2box.py の走査で一意）。
    近くの単項式・m_(2) = 1 では Laurent でない。f は整数係数で、f_(1), f_(2) は 1 probe の計算（別の展開）の f_1, f_2 と一致。"""
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "scripts"))
    import rk_k5_trefoil_2box as T2
    import rk_k5_framing_scan as Fs
    from rcsknots.kauffman_knot import free_energy, r as r_, s as s_
    N = 2
    assert sp.expand(T2.M2[N] - T2.M11[N]**2) == 0
    fu = T2.free_energy_2(N, {lam: T2.unknot_K(N, lam) for lam in ((1,), (1, 1), (2,))})
    assert fu[(2,)] == 0 and fu[(1, 1)] == 0
    f = T2.free_energy_2(N, T2.trefoil_W(N, T2.M2[N]))
    for R in ((1,), (2,), (1, 1)):
        st = T2.stats(f[R])
        assert st["laurent"] and st["integer"]
    for m in (sp.Integer(1), T2.M2[N] * r_**2, T2.M2[N] / r_**2, T2.M2[N] * s_**2, T2.M2[N] / s_**2, -T2.M2[N]):
        assert not T2.laurent(T2.free_energy_2(N, T2.trefoil_W(N, m))[(1, 1)])
    parts = Fs.setup(N)
    W = [sp.Integer(1), sp.cancel(T2.M1 * parts[1][0]), sp.cancel(T2.M11[N] * parts[2][0])]
    f1p = free_energy(N, W, 2)
    assert sp.cancel(f1p[1] - f[(1,)]) == 0 and sp.cancel(f1p[2] - f[(2,)]) == 0
