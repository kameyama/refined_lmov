"""Colored HOMFLY-PT of torus knots and cables from the (unrefined) cabling formula, exact in sympy.

Variables: KN's (a, q) at t = q, with Q = q^{1/2} so that fractional twists are polynomial.
  unknot (unreduced, symmetric normalization):  s_λ(q^ρ) = (q/a)^{|λ|/2} Π_□ (q^{l'} − a q^{a'}) / (1 − q^{h})
  twist eigenvalue: θ_λ = q^{κ_λ/?}… see `theta` below (convention fixed by the trefoil test).
Cabling formula (Morton; Rosso–Jones for K = unknot):
  W_λ(C_{p,s}(K)) = θ_λ^{−ps} Σ_{μ ⊢ p|λ|} c^{(p)}_{λμ} θ_μ^{s/p} W_μ(K),   c^{(p)}_{λμ} = Σ_ρ χ_λ(ρ) χ_μ(pρ) / z_ρ,
with K in framing 0 and s measured in the Seifert framing of K.
"""
from fractions import Fraction
from functools import lru_cache

import sympy as sp

from .characters import chi
from .partitions import partitions, boxes, arm, leg, z, scale

a, Q = sp.symbols("a Q", positive=True)
q = Q**2


def content(lam):
    return sum(j - i for i, j in boxes(lam))


@lru_cache(maxsize=None)
def adams_coeffs(lam, p):
    """c^{(p)}_{λμ}: s_λ[p_k ↦ p_{pk}] = Σ_μ c_{λμ} s_μ."""
    n = sum(lam)
    out = {}
    for mu in partitions(p * n):
        c = sum(Fraction(chi(lam, rho) * chi(mu, scale(rho, p)), z(rho)) for rho in partitions(n))
        if c:
            out[mu] = int(c)
    return out


def unknot(lam):
    """Unreduced unknot W_λ(○) = (q/a)^{|λ|/2} Π (q^{l'} − a q^{a'}) / (1 − q^h)."""
    k = sum(lam)
    w = (Q**2 / a) ** sp.Rational(k, 2)
    for i, j in boxes(lam):
        w *= (q**i - a * q**j) / (1 - q ** (arm(lam, i, j) + leg(lam, i, j) + 1))
    return w


def theta(lam, power):
    """θ_λ^{power} with θ_λ = q^{κ_λ} a^{|λ|/2} (= DGR 𝐪^{2κ} 𝐚^{|λ|}, κ_λ = Σ contents); power may be a Fraction."""
    power = sp.Rational(power)
    return Q ** (2 * content(lam) * power) * a ** (sum(lam) * power / 2)


def cable_unreduced(lam, p, s, companion, sign=1):
    """W_λ of the (p, s)-cable of a framing-0 knot whose unreduced invariants are companion(μ)."""
    total = 0
    for mu, c in adams_coeffs(lam, p).items():
        total += c * theta(mu, sp.Rational(sign * s, p)) * companion(mu)
    return sp.expand(theta(lam, -sign * p * s)) * total


def torus_reduced(lam, m, n, sign=1):
    """Reduced colored HOMFLY of T(m, n) = (m, n)-cable of the unknot."""
    val = cable_unreduced(lam, m, n, unknot, sign)
    return sp.factor(sp.simplify(val / unknot(lam)))
