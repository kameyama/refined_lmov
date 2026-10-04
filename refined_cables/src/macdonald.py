"""Macdonald symmetric functions P_λ(x;q,t) in the power-sum basis.

All functions are generic in the coefficient field: q, t are field elements
(sympy expressions for exact work, flint.nmod for arithmetic mod a prime) and
`conv` maps a fractions.Fraction into that field.

Conventions (Macdonald, Ch. VI):
  <p_λ, p_μ>_{q,t} = δ_{λμ} z_λ Π_i (1 - q^{λ_i}) / (1 - t^{λ_i}),
  P_λ = m_λ + Σ_{μ<λ} u_{λμ} m_μ,  orthogonal for < , >_{q,t}.
P_λ is obtained by Gram–Schmidt along any linear extension of dominance (here: lex order).
"""
from fractions import Fraction
from functools import lru_cache

from partitions import partitions, z, boxes, arm, leg, conjugate, n_stat


@lru_cache(maxsize=None)
def _p_in_m(n):
    """R[ν][λ] = coefficient of m_λ in p_ν (integers)."""
    parts = partitions(n)

    def count(nu, lam):
        # number of maps from the parts of ν to the rows of λ with row sums λ
        @lru_cache(maxsize=None)
        def rec(k, rem):
            if k == len(nu):
                return 1 if all(r == 0 for r in rem) else 0
            total = 0
            for j, r in enumerate(rem):
                if r >= nu[k]:
                    total += rec(k + 1, rem[:j] + (r - nu[k],) + rem[j + 1:])
            return total
        return rec(0, lam)

    return {nu: {lam: count(nu, lam) for lam in parts} for nu in parts}


@lru_cache(maxsize=None)
def m_in_p(n):
    """M[λ][ν] = coefficient of p_ν in m_λ (Fractions), inverse of _p_in_m."""
    from flint import fmpq_mat
    parts = partitions(n)
    R = _p_in_m(n)
    A = fmpq_mat([[R[nu][lam] for lam in parts] for nu in parts])  # rows ν, cols λ
    Ainv = A.inv()  # m_λ = Σ_ν Ainv[λ, ν] p_ν
    out = {}
    for i, lam in enumerate(parts):
        out[lam] = {}
        for j, nu in enumerate(parts):
            c = Ainv[i, j]
            if c != 0:
                out[lam][nu] = Fraction(int(c.p), int(c.q))
    return out


def pairing_weight(nu, q, t, conv):
    """<p_ν, p_ν>_{q,t}."""
    w = conv(Fraction(z(nu)))
    for k in nu:
        w = w * (1 - q ** k) / (1 - t ** k)
    return w


def macdonald_P(n, q, t, conv):
    """{λ: {ν: coefficient of p_ν in P_λ}} for all λ ⊢ n."""
    if n == 0:
        return {(): {(): conv(Fraction(1))}}
    parts = list(reversed(partitions(n)))  # lex ascending: lower partitions first
    M = m_in_p(n)
    w = {nu: pairing_weight(nu, q, t, conv) for nu in parts}
    mvec = {lam: {nu: conv(c) for nu, c in M[lam].items()} for lam in parts}

    def ip(f, g):
        s = conv(Fraction(0))
        for nu, c in f.items():
            d = g.get(nu)
            if d is not None:
                s = s + c * d * w[nu]
        return s

    P, norm = {}, {}
    for lam in parts:
        vec = dict(mvec[lam])
        for mu in P:
            c = ip(mvec[lam], P[mu]) / norm[mu]
            for nu, d in P[mu].items():
                vec[nu] = vec.get(nu, 0) - c * d
        P[lam] = vec
        norm[lam] = ip(vec, vec)
    return P


def macdonald_norm(lam, q, t):
    """<P_λ, P_λ>_{q,t} = Π_s (1 - q^{a+1} t^l) / (1 - q^a t^{l+1})."""
    out = 1
    for i, j in boxes(lam):
        a, l = arm(lam, i, j), leg(lam, i, j)
        out = out * (1 - q ** (a + 1) * t ** l) / (1 - q ** a * t ** (l + 1))
    return out


def principal_specialization(lam, a, q, t):
    """P_λ(1, t, ..., t^{N-1}) with a = t^N  (Macdonald (VI.6.11')).

    = Π_s (t^{l'(s)} - a q^{a'(s)}) / (1 - q^{a(s)} t^{l(s)+1}),  the unknot of 1703.05408.
    """
    out = 1
    for i, j in boxes(lam):
        out = out * (t ** i - a * q ** j) / (1 - q ** arm(lam, i, j) * t ** (leg(lam, i, j) + 1))
    return out


def t_eigenvalue_exponents(lam):
    """Exponents (e_q, e_t) of T_λ = q^{Σ λ_i(λ_i-1)/2} t^{-Σ (i-1) λ_i} = q^{n(λ')} t^{-n(λ)}.

    Sign of the t-exponent is a convention to be fixed in Stage 2 (see LOG.md).
    """
    return n_stat(conjugate(lam)), -n_stat(lam)
