"""K4'：被覆（composite）の Hopf 絡み目。Pieri 係数、主特殊化、unknot と [1],[1] の解析式との一致。"""
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rcsknots.kauffman_cover import pieri_vertical, principal, cover_Z, log_coeffs  # noqa: E402
from rcsknots.macdonald import macdonald_P, pairing_weight, macdonald_norm  # noqa: E402
from rcsknots.partitions import partitions, z as zl  # noqa: E402

q, s = F(3, 7), F(5, 11)
t = s**2


def test_pieri_against_power_sums():
    conv = lambda c: F(c)
    def e_in_p(r):
        return {rho: F((-1)**(r - len(rho)), zl(rho)) for rho in partitions(r)}
    def pmul(f, g):
        out = {}
        for a, ca in f.items():
            for b, cb in g.items():
                key = tuple(sorted(a + b, reverse=True)); out[key] = out.get(key, F(0)) + ca * cb
        return out
    for n in (2, 3):
        Pn = macdonald_P(n, q, t, conv)
        for mu in partitions(n):
            for r_ in (1, 2):
                prod_ = pmul(Pn[mu], e_in_p(r_)); Pm = macdonald_P(n + r_, q, t, conv)
                for lam in partitions(n + r_):
                    mu2 = tuple(mu) + (0,) * (len(lam) - len(mu))
                    if len(lam) < len(mu) or not all(0 <= a - b <= 1 for a, b in zip(lam, mu2)):
                        continue
                    ip = sum(c * Pm[lam].get(rho, 0) * pairing_weight(rho, q, t, conv) for rho, c in prod_.items())
                    assert pieri_vertical(lam, mu, q, t) == ip / macdonald_norm(lam, q, t)


def test_principal_specialization_small_N():
    """P_λ(t^ρ) を 3 変数の Macdonald 多項式（冪和）で直接評価したものと比べる。"""
    N = 3
    conv = lambda c: F(c)
    xs = [s**(N + 1 - 2 * i) for i in range(1, N + 1)]
    for n in (1, 2, 3):
        for lam, P in macdonald_P(n, q, t, conv).items():
            if len(lam) > N:
                continue
            val = sum(c * __import__("math").prod(sum(x**k for x in xs) for k in rho) for rho, c in P.items())
            assert principal(lam, N, q, t, s) == val


def test_cover_unknot_and_hopf_11():
    """½f^{cov}_{10} = −(ν−1/ν)/(s−1/s)、½f^{cov}_{11} = (ν−1/ν)(ν~−1/ν~) z_q/z_t、ν = s^N、ν~ = ν q^{1/2}/s（q = r²）。"""
    r = F(3, 7); qq = r**2
    for N in (3, 5, 7):
        L = log_coeffs(cover_Z(N, qq, s, 2, 2), 2, 2)
        nu = s**N; nut = nu * r / s
        zt = s - 1 / s; zq = r - 1 / r
        assert L[(1, 0)] / 2 == -(nu - 1 / nu) / zt
        assert L[(0, 1)] / 2 == -(nu - 1 / nu) / zt
        assert L[(1, 1)] / 2 == (nu - 1 / nu) * (nut - 1 / nut) * zq / zt


def test_symbolic_cover_matches_finite_N():
    """A = ν² を形式変数とした Z^{cov} が、有限 N（N = 3, 5, 7）の直接計算と一致する（A = B = 2）。"""
    import sympy as sp
    from rcsknots.kauffman_cover import cover_Z_sym, qs_, ss_, nus_
    Zs = cover_Z_sym(2, 2)
    for N in (3, 5, 7):
        Zn = cover_Z(N, q, s, 2, 2)
        sub = {qs_: sp.Rational(q.numerator, q.denominator), ss_: sp.Rational(s.numerator, s.denominator), nus_: sp.Rational(s.numerator, s.denominator)**N}
        for ab, v in Zn.items():
            assert Zs[ab].as_expr().subs(sub) == sp.Rational(v.numerator, v.denominator), (N, ab)


def test_hopf_crosscap_fundamental_sector():
    """g_{1b} = z_q (ν~ − 1/ν~) χ_{b−1}(ν~/q^{1/2})（b ≤ 3）：t は ν~ = ν (q/t)^{1/2} を通してしか入らない（A = B = 3 の記号計算）。"""
    import pickle, sympy as sp
    from pathlib import Path
    from rcsknots.kauffman_cover import cover_f_sym, qs_, ss_, nus_
    r_, s_, nu_, nut = sp.symbols("r s nu nut", positive=True)
    fc = cover_f_sym(3, 1)
    q_, t_, U_ = sp.symbols("q t U")
    f = {ab: sp.sympify(v, locals={"q": q_, "t": t_, "U": U_}) for ab, v in pickle.load(open(ROOT / "data/refined_kauffman/hopf_f.pkl", "rb")).items()}
    zq = r_ - 1 / r_
    chi = lambda m, x: sum(x**(m - 2 * i) for i in range(m + 1))
    for b in (1, 2, 3):
        fo = f[(b, 1)].subs({q_: r_**2, t_: s_**2, U_: nu_ * s_})
        g = sp.cancel(fo - fc[(b, 1)].subs({qs_: r_**2, ss_: s_, nus_: nu_}) / 2)
        g = sp.cancel(g.subs(nu_, nut * s_ / r_))
        assert sp.cancel(g - zq * (nut - 1 / nut) * chi(b - 1, nut / r_)) == 0, b


def test_hopf_crosscap_all_d_even_sector():
    """全 d の交差帽部分 g_{ab} = f_{ab} − ½f^{cov}_{ab}（記号計算、a, b ≤ 2）：
    (1,2) は整数 Laurent、(2,2) は分母がちょうど 2(1+t)（× 単項式）で 2(1+t)g_{22} は整数 Laurent、
    E_{22}(t=−1) = [2(1+t)g_{22}]_{t=−1} = z_q (z_q^+)² (ν~−ν~⁻¹)(q^{−1/2}ν~ − q^{1/2}ν~⁻¹)。"""
    import sympy as sp
    from rcsknots.kauffman_cover import cover_f_sym, qs_, ss_, nus_
    from rcsknots.kauffman_hopf import hopf_log_coeffs, mobius_f, qs, ts, Us
    r_, s_, nu_ = sp.symbols("r s nu", positive=True)
    nut = sp.Symbol("nut")
    L = hopf_log_coeffs(2, 2)
    fc = cover_f_sym(2, 2)

    def g_of(a, b):
        fo = mobius_f(L, a, b).subs({qs: r_**2, ts: s_**2, Us: nu_ * s_})
        return sp.cancel(fo - fc[(a, b)].subs({qs_: r_**2, ss_: s_, nus_: nu_}) / 2)

    def laurent_int(e):
        n, d = sp.fraction(sp.cancel(e))
        return sp.Poly(d, r_, s_, nu_).is_monomial and all(c.is_integer for c in sp.Poly(sp.expand(n), r_, s_, nu_).coeffs())

    assert laurent_int(g_of(1, 2))
    g22 = g_of(2, 2)
    den = sp.factor(sp.fraction(g22)[1])
    assert sp.cancel(den / (2 * (1 + s_**2))).is_polynomial(s_) and sp.Poly(sp.cancel(den / (2 * (1 + s_**2))), r_, s_, nu_).is_monomial
    R = sp.cancel(2 * (1 + s_**2) * g22)
    assert laurent_int(R)
    E_m1 = sp.cancel(sp.cancel(R.subs(nu_, nut * s_ / r_)).subs(s_, sp.I))
    zq, zp = r_ - 1 / r_, r_ + 1 / r_
    assert sp.simplify(E_m1 - zq * zp**2 * (nut - 1 / nut) * (nut / r_ - r_ / nut)) == 0


def test_sector_resolved_cover_and_fixed_sector_22():
    """セクター別の被覆（probe を x 側と 1/x 側に分けた 4 変数）：セクター和が 1 変数の f^{cov} に等しく、Ĩ 対称 f(k,l,m,n) = f(l,k,n,m)。
    (2,2) の固定セクター (1,1)×(1,1) の全状態数 P = (1−t)f^{cov}_{(11),(11)} = z_q³(ν~−ν~⁻¹)(t q^{−1/2}ν~ − q^{1/2}ν~⁻¹)、
    候補 B の E_22 = −z_q(z_q^+)²(ν~−ν~⁻¹)(t q^{−1/2}ν~ + q^{1/2}ν~⁻¹) は (P − E)/2 が整数 Laurent（固定セクターの部分和）。"""
    import sympy as sp
    from rcsknots.kauffman_cover import cover_f_sym4, cover_f_sym, qs_, ss_, nus_
    r_, s_, nu_ = sp.symbols("r s nu", positive=True)
    nut = sp.Symbol("nut")
    sub = {qs_: r_**2, ss_: s_, nus_: nu_}
    f4 = {k: v.subs(sub) for k, v in cover_f_sym4(2, 2).items()}
    fc = {k: v.subs(sub) for k, v in cover_f_sym(2, 2).items()}
    for (a, b), v in fc.items():
        assert sp.cancel(sum(x for (k, l, m, n), x in f4.items() if k + l == a and m + n == b) - v) == 0
    for (k, l, m, n), v in f4.items():
        assert sp.cancel(v - f4[(l, k, n, m)]) == 0
    t = s_**2; z, zp = r_ - 1 / r_, r_ + 1 / r_
    P = sp.cancel(((1 - t) * f4[(1, 1, 1, 1)]).subs(nu_, nut * s_ / r_))
    assert sp.simplify(P - z**3 * (nut - 1 / nut) * (t / r_ * nut - r_ / nut)) == 0
    E = -z * zp**2 * (nut - 1 / nut) * (t / r_ * nut + r_ / nut)
    Q = sp.cancel((P - E) / 2)
    assert sp.simplify(Q - z * (nut - 1 / nut) * ((r_ + r_**-3) * t * nut + 2 * r_ / nut)) == 0


def test_E22_unique_from_t_minus1_and_t_q():
    """E_22 は固定セクター (1,1)×(1,1) の状態の和なので t（= s²）の 0, 1 次だけをもつ。
    (i) t = −1 で E = 2(1+t)g_22、(ii) t = q で E = 2(1+q)(g_22 − ĝ_22)（ĝ は Mariño の交差帽 = g^{odd} の t = q）
    の 2 条件から E = E_0 + t E_1 は一意に決まり、−z_q(z_q^+)²(ν~−ν~⁻¹)(t q^{−1/2}ν~ + q^{1/2}ν~⁻¹) に等しい。"""
    import sympy as sp
    from rcsknots.kauffman_cover import cover_f_sym, qs_, ss_, nus_
    from rcsknots.kauffman_hopf import hopf_log_coeffs, mobius_f, qs, ts, Us
    r_, s_, nu_ = sp.symbols("r s nu", positive=True)
    T = sp.Symbol("T")
    L = hopf_log_coeffs(2, 2)
    fc = cover_f_sym(2, 2)

    def g_of(a, b):
        fo = mobius_f(L, a, b).subs({qs: r_**2, ts: s_**2, Us: nu_ * s_})
        return sp.cancel(fo - fc[(a, b)].subs({qs_: r_**2, ss_: s_, nus_: nu_}) / 2)

    g22, g11 = g_of(2, 2), g_of(1, 1)
    g11_x2 = g11.subs({r_: r_**2, s_: s_**2, nu_: nu_**2}, simultaneous=True)
    g22_odd_tq = sp.cancel((g22 + g11_x2 / 2).subs(s_, r_))          # g^{odd}_22 = g_22 + ½ g^{odd}_11(x²)、t = q
    R = sp.cancel(2 * (1 + s_**2) * g22)
    Em1 = sp.cancel(R.subs(s_, sp.I))                                   # E_0 − E_1
    Eq = sp.cancel(R.subs(s_, r_) - 2 * (1 + r_**2) * g22_odd_tq)        # E_0 + q E_1
    E1 = sp.cancel((Eq - Em1) / (1 + r_**2)); E0 = sp.cancel(Em1 + E1)
    E = E0 + s_**2 * E1
    nut = nu_ * r_ / s_
    zq, zp = r_ - 1 / r_, r_ + 1 / r_
    assert sp.cancel(E - (-zq * zp**2 * (nut - 1 / nut) * (s_**2 / r_ * nut + r_ / nut))) == 0


def test_multi_probe_cover_and_adams_b2():
    """成分 1 に probe 無限個の被覆（cover_f_multi(2,2)）：1 変数への特殊化が cover_f_sym4 に一致し、
    Ĩ = ρ⊗ρ の入れ替え（Adams p_2 = s_[2] − s_[1,1]）の予言 E_[2],2 = −E_[1,1],2 = E_22 が成り立つ：
    g_[2],2 + g_[1,1],2 と g_[2],2 − ½E_22/(1+t) は (1+t) の分母を持たない（SO 側は hopf_f_multi.pkl）。"""
    import pickle, sympy as sp
    from rcsknots.kauffman_cover import cover_f_multi, cover_f_sym4, qs_, ss_, nus_
    r_, s_, nu_ = sp.symbols("r s nu", positive=True)
    sub = {qs_: r_**2, ss_: s_, nus_: nu_}
    fm = {k: v.subs(sub) for k, v in cover_f_multi(2, 2).items()}
    for (k, l, m, n), v in cover_f_sym4(2, 2).items():
        assert sp.cancel(fm[((k,) if k else (), (l,) if l else (), m, n)] - v.as_expr().subs(sub)) == 0
    half = {}
    for (r1, r2, k, l), v in fm.items():
        if sum(r1) + sum(r2) != 2 or k + l != 2:
            continue
        if r1 and r2:          # s_1 s_1 = s_2 + s_11
            for R in ((2,), (1, 1)):
                half[R] = half.get(R, 0) + v / 2
        else:
            R = r1 or r2
            half[R] = half.get(R, 0) + v / 2
    q_, t_, U_ = sp.symbols("q t U")
    fSO = {k: sp.sympify(v, locals={"q": q_, "t": t_, "U": U_}).subs({q_: r_**2, t_: s_**2, U_: nu_ * s_})
           for k, v in pickle.load(open(ROOT / "data/refined_kauffman/hopf_f_multi.pkl", "rb"))["schur"].items()}
    g2, g11 = sp.cancel(fSO[((2,), 2)] - half[(2,)]), sp.cancel(fSO[((1, 1), 2)] - half[(1, 1)])
    mono = lambda e: sp.Poly(sp.fraction(sp.cancel(e))[1], r_, s_, nu_).is_monomial
    assert not mono(g2) and not mono(g11)
    assert mono(g2 + g11)
    nut = nu_ * r_ / s_
    zq, zp = r_ - 1 / r_, r_ + 1 / r_
    E22 = -zq * zp**2 * (nut - 1 / nut) * (s_**2 / r_ * nut + r_ / nut)
    assert mono(g2 - E22 / (2 * (1 + s_**2)))
