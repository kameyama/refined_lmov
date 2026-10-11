"""Refined reformulated invariants f^q, f^{\\bar t}, \\hat f and the refined BPS numbers (1703.05408 §3–4, App. B).

Symbolic part (sympy, |μ| ≤ 3):
  f^q_μ : from Σ_λ \\bar{rCS}_λ g_λ P_λ(x;q,t) = Exp( Σ_μ f^q_μ/(q^{1/2}−q^{−1/2}) s_μ(x) )            (refined-GT1)
  f^t_μ : from Σ_λ \\bar{rCS}_λ P_{λ^T}(−x;t,q) = Exp( Σ_μ f^t_μ/(t^{−1/2}−t^{1/2}) s_μ(x) )            (refined-GT2)
  in terms of symbols C[λ] = \\bar{rCS}_λ(a,q,t) and C1[d] = \\bar{rCS}_{[1]}(a^d,q^d,t^d).
  M_{μρ}(t) = Σ_σ C_{μσρ} B_σ(t);  f^q = M(t) \\hat f,  f^t = M(q^{−1}) \\hat f.
Variables: Q = q^{1/2}, T = t^{1/2} (sympy symbols), so that all expressions are rational in Q, T.
"""
from fractions import Fraction
from functools import lru_cache

import sympy as sp

from .characters import chi
from .macdonald import macdonald_P, macdonald_norm
from .partitions import partitions, z, conjugate

Q, T = sp.symbols("Q T", positive=True)
q, t = Q**2, T**2
COLORS3 = [(1,), (2,), (1, 1), (3,), (2, 1), (1, 1, 1)]


def C(lam):
    return sp.Symbol("C_" + "".join(map(str, lam)))


def C1(d):
    return sp.Symbol(f"C1_{d}")


def Cd(lam, d):
    """Adams 像 C_λ(x^d)（q, t も d 乗）の記号。d = 1 は C(λ)、λ = [1] は C1(d)（3 箱までの記号と互換）。"""
    if d == 1:
        return C(lam)
    if lam == (1,):
        return C1(d)
    return sp.Symbol("C_" + "".join(map(str, lam)) + f"_{d}")


def _rat(c):
    return sp.Rational(c.numerator, c.denominator)


@lru_cache(maxsize=None)
def _P_sym(n, swap=False):
    """P_λ(x;q,t) (or P_λ(x;t,q) if swap) in the p-basis, sympy."""
    qq, tt = (t, q) if swap else (q, t)
    return macdonald_P(n, qq, tt, _rat)


def g_lambda(lam, qq=q, tt=t):
    """g_λ = (q/t)^{|λ|} / <P_λ,P_λ>_{q,t}."""
    return (qq / tt) ** sum(lam) / macdonald_norm(lam, qq, tt)


def _mul_trunc(u, v, nmax):
    out = {}
    for r1, x in u.items():
        for r2, y in v.items():
            if sum(r1) + sum(r2) <= nmax:
                key = tuple(sorted(r1 + r2, reverse=True))
                out[key] = out.get(key, 0) + x * y
    return out


def _plethystic_log(adams, nmax):
    """Log(Z) = Σ_d μ(d)/d log(ψ_d Z) for Z = 1 + Σ_ρ Z_ρ p_ρ, truncated to degree ≤ nmax.

    adams(d) returns ψ_d Z − 1 as a dict over partitions (p_k → p_{dk}, parameters → d-th powers).
    """
    total = {}
    for d in range(1, nmax + 1):
        mu = sp.mobius(d)
        if mu == 0:
            continue
        X = adams(d)
        power, k = dict(X), 1
        while power:
            coeff = sp.Rational(mu, d) * sp.Rational((-1) ** (k - 1), k)
            for r, c in power.items():
                total[r] = total.get(r, 0) + coeff * c
            power, k = _mul_trunc(power, X, nmax), k + 1
    return total


def _schur_coeffs(F, nmax):
    """Coefficient of s_μ in F = Σ F_ρ p_ρ:  Σ_ρ F_ρ χ_μ(ρ)."""
    return {mu: sum(F.get(rho, 0) * chi(mu, rho) for rho in partitions(sum(mu)))
            for n in range(1, nmax + 1) for mu in partitions(n)}


def _subs_d(expr, d):
    return expr.subs({Q: Q**d, T: T**d}, simultaneous=True)


@lru_cache(maxsize=None)
def f_formulas(nmax=3):
    """{μ: (f^q_μ, f^t_μ)} as sympy expressions in Q, T, C[λ], C1[d]."""
    def zq(d):
        out = {}
        for n in range(1, nmax // d + 1):
            for lam in partitions(n):
                cs = Cd(lam, d)
                g = _subs_d(g_lambda(lam), d)
                for rho, c in _P_sym(n)[lam].items():
                    key = tuple(d * k for k in rho)
                    out[key] = out.get(key, 0) + cs * g * _subs_d(c, d)
        return out

    def zt(d):
        out = {}
        for n in range(1, nmax // d + 1):
            for lam in partitions(n):
                cs = Cd(lam, d)
                for rho, c in _P_sym(n, swap=True)[conjugate(lam)].items():
                    key = tuple(d * k for k in rho)
                    out[key] = out.get(key, 0) + cs * (-1) ** n * _subs_d(c, d)
        return out

    sq = _schur_coeffs(_plethystic_log(zq, nmax), nmax)
    st = _schur_coeffs(_plethystic_log(zt, nmax), nmax)
    # (T/Q)^{|μ|}: App. B and the notebook rescale x by (t/q)^{1/2} relative to the text definition of g_λ
    return {mu: (sq[mu] * (Q - 1 / Q) * (T / Q) ** sum(mu), st[mu] * (1 / T - T)) for mu in sq}


# ---- M matrix ------------------------------------------------------------------------------------

def clebsch(mu, sigma, rho):
    n = sum(mu)
    return sum(Fraction(chi(mu, nu) * chi(sigma, nu) * chi(rho, nu), z(nu)) for nu in partitions(n))


def B(sigma, tt=T):
    """B_σ(t) = (−t)^d t^{−(|σ|−1)/2} for the hook σ = (|σ|−d, 1^d), else 0.  (tt = t^{1/2})"""
    n = sum(sigma)
    if len(sigma) > 1 and sigma[1] > 1:
        return 0
    d = len(sigma) - 1
    return (-tt**2) ** d * tt ** (-(n - 1))


def M_matrix(n, tt=T):
    parts = partitions(n)
    return sp.Matrix([[sum(_rat(clebsch(mu, sg, rho)) * B(sg, tt) for sg in parts) for rho in parts] for mu in parts])


# ---- g-expansion (notebook `gfactor`) ------------------------------------------------------------

def to_bold(poly, a_shift2=0):
    """Laurent poly in (A=a^{1/2}, Q, T) of KN → bold (𝐚, 𝐪, 𝐭) exponents, after multiplying by A^{a_shift2}.

    a = −𝐚²𝐭 (so A² = −𝐚²𝐭; only even total A-exponents are allowed), √q = Q = −𝐪𝐭, √t = T = 𝐪.
    Returns {(e_𝐚, e_𝐪, e_𝐭): c}.
    """
    out = {}
    for (ea, eq, et), c in poly.items():
        ea += a_shift2
        if ea % 2:
            raise ValueError("odd power of a^{1/2} after the a-grading shift")
        m = ea // 2                                   # a^m = (−1)^m 𝐚^{2m} 𝐭^m
        key = (2 * m, eq + et, m + eq)
        sign = (-1) ** ((m + eq) % 2)
        out[key] = out.get(key, 0) + sign * c
    return {k: v for k, v in out.items() if v}


def g_expand(bold):
    """Same result as `g_expand_direct`, computed with polynomials in 𝐭 (python-flint).

    With w = 𝐪²𝐭 and s = 𝐭 + 𝐭^{−1}, G = w + w^{−1} − s.  The 𝐚^{2β} row closes iff its 𝐪-exponents are even and it is
    symmetric under w ↔ w^{−1}; then it is a polynomial in y = w + w^{−1} (peel off y^k = Σ_j C(k,j) w^{k−2j}),
    and Horner's rule in y = G + s gives the coefficients of G^g.  Raises if the expansion does not close.
    """
    from math import comb
    from flint import fmpz_poly

    rows = {}
    for (ea, eq, et), c in bold.items():
        if eq % 2:
            raise ValueError(f"g-expansion does not close: odd power of q in {(ea, eq, et)}")
        rows.setdefault(ea // 2, {})[(eq // 2, et - eq // 2)] = c      # 𝐪^{2k} 𝐭^{e} = w^k 𝐭^{e−k}
    tpoly = fmpz_poly([0, 1])
    result = {}
    for beta, X in rows.items():
        if any(X.get((-k, f)) != c for (k, f), c in X.items()):
            raise ValueError(f"g-expansion does not close at beta={beta}: not symmetric under w <-> 1/w")
        K = max(k for k, _ in X)
        tmin = min(f for _, f in X)
        d = [[0] * (max(f for _, f in X) - tmin + 1) for _ in range(K + 1)]
        for (k, f), c in X.items():
            if k >= 0:
                d[k][f - tmin] = c
        d = [fmpz_poly(v) for v in d]
        e = [None] * (K + 1)                     # X = Σ_k e_k(𝐭) y^k
        for k in range(K, -1, -1):
            e[k] = d[k]
            for j in range(1, k // 2 + 1):
                d[k - 2 * j] -= comb(k, j) * e[k]
        # Q_k = 𝐭^{K−k} P_k with P_K = e_K, P_k = P_{k+1}(G + s) + e_k; (G + s)𝐭 = G𝐭 + 𝐭² + 1
        Qs = [e[K]]                              # Qs[g]: coefficient of G^g (polynomial in 𝐭, offset tmin)
        for k in range(K - 1, -1, -1):
            new = [x * (tpoly * tpoly + 1) for x in Qs] + [fmpz_poly(0)]
            for g in range(len(Qs)):
                new[g + 1] += Qs[g] * tpoly
            new[0] += e[k] * tpoly ** (K - k)
            Qs = new
        for g, x in enumerate(Qs):                # P_0 = Q_0 / 𝐭^K
            for i, c in enumerate(x.coeffs()):
                if c:
                    result[(g, beta, i + tmin - K)] = int(c)
    return result


def g_expand_direct(bold):
    """Notebook `gfactor`: write the coefficient of 𝐚^{2β} as Σ_g c_g(𝐭) G^g with G = (𝐪−𝐪^{−1})(𝐪𝐭−𝐪^{−1}𝐭^{−1}).

    Returns {(g, β, F): N} with N the coefficient of 𝐭^F G^g 𝐚^{2β}.  Raises if the expansion does not close.
    """
    G = {(2, 1): 1, (0, 1): -1, (0, -1): -1, (-2, -1): 1}          # (e_𝐪, e_𝐭) of G = 𝐪²𝐭 − 𝐭 − 𝐭^{−1} + 𝐪^{−2}𝐭^{−1}
    by_beta = {}
    for (ea, eq, et), c in bold.items():
        by_beta.setdefault(ea // 2, {})[(eq, et)] = c

    def mul(u, v):
        out = {}
        for (a1, b1), x in u.items():
            for (a2, b2), y in v.items():
                k = (a1 + a2, b1 + b2)
                out[k] = out.get(k, 0) + x * y
        return {k: w for k, w in out.items() if w}

    result = {}
    for beta, X in by_beta.items():
        X = dict(X)
        top = max(eq for eq, _ in X)
        Gpow = {0: {(0, 0): 1}}
        for n in range(1, top // 2 + 2):
            Gpow[n] = mul(Gpow[n - 1], G)
        for n in range(top // 2 + 1, 0, -1):
            coeffs = {et: c for (eq, et), c in X.items() if eq == 2 * n}
            if not coeffs:
                continue
            # [𝐪^{2n}] G^n = 𝐭^n, so x(𝐭) = coeffs / 𝐭^n
            for et, c in coeffs.items():
                result[(n, beta, et - n)] = c
                for (gq, gt), gc in Gpow[n].items():
                    k = (gq, gt + et - n)
                    X[k] = X.get(k, 0) - c * gc
            X = {k: w for k, w in X.items() if w}
        if any(eq != 0 for eq, _ in X) or any(abs(eq) % 2 for eq, _ in X):
            raise ValueError(f"g-expansion does not close at beta={beta}: remainder {X}")
        for (eq, et), c in X.items():
            result[(0, beta, et)] = c
    return result
