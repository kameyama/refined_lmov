"""Fix the monomial normalization of a refined invariant (docs/conventions.md §5).

Hypothesis (checked on all torus-knot data of the 1703.05408 notebook in Stage 2):
  (N1) at t = q the reduced invariant equals the framing-0 reduced colored HOMFLY-PT (DGR conventions);
  (N2) in the bold (DGR) variables a = −𝐚²𝐭, √q = −𝐪𝐭, √t = 𝐪 the minimal 𝐭-degree is 0.
These two conditions determine the factor ±a^α q^β t^γ uniquely.

Input `raw`: {(e_a, e_q, e_t): c} with integer exponents of a, q, t.
Output: {(2e_a, 2e_q, 2e_t): c}, i.e. exponents of A = a^{1/2}, Q = q^{1/2}, T = t^{1/2} (as in laurent.py).
"""
from fractions import Fraction

import sympy as sp

from .homfly import a as ha, Q as hQ


def t_limit(raw):
    """raw at t = q as a sympy expression in (a, Q)."""
    return sum(c * ha**ea * hQ ** (2 * (eq + et)) for (ea, eq, et), c in raw.items())


def normalize(raw, homfly_reduced):
    ratio = sp.factor(sp.simplify(homfly_reduced / t_limit(raw)))
    coeff, factors = ratio.as_coeff_mul()
    if coeff not in (1, -1):
        raise ValueError(f"t = q limit is not a monomial multiple of the HOMFLY-PT polynomial: {ratio}")
    powers = sp.Mul(*factors).as_powers_dict()
    if set(powers) - {ha, hQ}:
        raise ValueError(f"t = q limit is not a monomial multiple of the HOMFLY-PT polynomial: {ratio}")
    alpha = Fraction(str(powers.get(ha, 0)))
    beta_plus_gamma = Fraction(str(powers.get(hQ, 0))) / 2
    tmin = min(Fraction(ea) + alpha + 2 * eq for (ea, eq, et) in raw)
    beta = -tmin / 2
    gamma = beta_plus_gamma - beta
    out = {}
    for (ea, eq, et), c in raw.items():
        key = (2 * (ea + alpha), 2 * (eq + beta), 2 * (et + gamma))
        assert all(k.denominator == 1 for k in key)
        out[tuple(int(k) for k in key)] = int(coeff) * c
    return out, (int(coeff), alpha, beta, gamma)
