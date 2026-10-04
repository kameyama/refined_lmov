"""Refined invariants of torus knots and their cables, evaluated at a numeric point mod a prime.

A point fixes base values x, y, b with q = x^L, t = y^L, a = b^2, so that the
fractional powers q^{k/L}, t^{k/L}, a^{1/2} that appear in the formulas are well defined.
"""
from fractions import Fraction

from flint import nmod

from .partitions import partitions, scale, n_stat, conjugate
from .macdonald import macdonald_P, pairing_weight, principal_specialization, macdonald_norm

PRIME = 2**61 - 1


class Point:
    def __init__(self, x, y, b, L=2, p=PRIME):
        self.p, self.L = p, L
        self.x, self.y, self.b = nmod(x, p), nmod(y, p), nmod(b, p)
        self.q, self.t, self.a = self.x ** L, self.y ** L, self.b ** 2
        self._P = {}

    def conv(self, c):
        return nmod(c.numerator, self.p) / nmod(c.denominator, self.p)

    def mono(self, ea=0, eq=0, et=0):
        """a^{ea} q^{eq} t^{et} for ea ∈ ½Z, eq, et ∈ (1/L)Z."""
        ea2, eqL, etL = Fraction(ea) * 2, Fraction(eq) * self.L, Fraction(et) * self.L
        assert ea2.denominator == eqL.denominator == etL.denominator == 1, (ea, eq, et)
        return self.b ** int(ea2) * self.x ** int(eqL) * self.y ** int(etL)

    def P(self, n):
        if n not in self._P:
            self._P[n] = macdonald_P(n, self.q, self.t, self.conv)
        return self._P[n]

    def norm(self, lam):
        return macdonald_norm(lam, self.q, self.t)

    def weight(self, nu):
        return pairing_weight(nu, self.q, self.t, self.conv)


def plethysm_coeffs(pt, lam, m):
    """c_{λν}: P_λ[p_k -> p_{mk}] = Σ_ν c_{λν} P_ν, with q, t unchanged (hypothesis, checked in Stage 2)."""
    n = sum(lam)
    Plam = pt.P(n)[lam]
    out = {}
    for nu, Pnu in pt.P(m * n).items():
        s = 0
        for rho, c in Plam.items():
            mr = scale(rho, m)
            d = Pnu.get(mr)
            if d is not None:
                s = s + c * d * pt.weight(mr)
        if s != 0:
            out[nu] = s / pt.norm(nu)
    return out


def t_eig(lam):
    """(e_q, e_t) of the T-eigenvalue T_λ = q^{n(λ')} t^{-n(λ)}."""
    return Fraction(n_stat(conjugate(lam))), Fraction(-n_stat(lam))


def unknot(pt, lam):
    """Symmetric principal specialization P_λ(t^ρ) = (t/a)^{|λ|/2} Π (t^{l'} - a q^{a'})/(1 - q^a t^{l+1})."""
    k = sum(lam)
    return pt.mono(ea=Fraction(-k, 2), et=Fraction(k, 2)) * principal_specialization(lam, pt.a, pt.q, pt.t)


def cable_sum(pt, lam, m, n, companion):
    """Σ_ν c^{(m)}_{λν} T_ν^{n/m} companion(ν)  (no framing correction)."""
    s = 0
    for nu, c in plethysm_coeffs(pt, lam, m).items():
        eq, et = t_eig(nu)
        s = s + c * pt.mono(eq=eq * Fraction(n, m), et=et * Fraction(n, m)) * companion(nu)
    return s
