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
