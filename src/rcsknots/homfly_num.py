"""Colored HOMFLY-PT of torus knots and cables from the cabling formula, evaluated at points (A, Q) mod p.

Same formula and conventions as homfly.py (A = a^{1/2}, Q = q^{1/2}, θ_λ = q^{κ_λ} a^{|λ|/2}), but numerical,
so that large colors are fast; exact polynomials are recovered with exact.reconstruct_points.
"""
from fractions import Fraction
from functools import lru_cache

import numpy as np

from .field import V
from .homfly import adams_coeffs, content
from .partitions import boxes, arm, leg


class HomflyNum:
    def __init__(self, Av, Qv):
        self.A, self.Q = Av, Qv
        self.a, self.q = Av * Av, Qv * Qv
        self._cache = {}

    def unknot(self, lam):
        k = sum(lam)
        w = (self.Q / self.A) ** k
        for i, j in boxes(lam):
            w = w * (self.q ** i - self.a * self.q ** j) / (1 - self.q ** (arm(lam, i, j) + leg(lam, i, j) + 1))
        return w

    def theta(self, lam, power):
        power = Fraction(power)
        eq, ea = 2 * content(lam) * power, sum(lam) * power
        assert eq.denominator == 1 and ea.denominator == 1, (lam, power)
        return self.Q ** int(eq) * self.A ** int(ea)

    def cable_unreduced(self, lam, p, s, companion, sign=-1):
        total = 0
        for mu, c in adams_coeffs(lam, p).items():
            total = total + c * self.theta(mu, Fraction(sign * s, p)) * companion(mu)
        return self.theta(lam, -sign * p * s) * total

    def iterated(self, lam, steps, sign=-1):
        """Unreduced invariant of the iterated cable given by steps [(p1, q1), (p2, q2), ...] (innermost = companion first)."""
        key = (lam, tuple(steps))
        if key not in self._cache:
            if not steps:
                val = self.unknot(lam)
            else:
                *inner, (p, s) = steps
                val = self.cable_unreduced(lam, p, s, lambda mu: self.iterated(mu, inner, sign), sign)
            self._cache[key] = val
        return self._cache[key]

    def reduced(self, lam, steps, sign=-1):
        return self.iterated(lam, steps, sign) / self.unknot(lam)


def homfly_evaluator(lam, steps, sign=-1):
    """evaluator(xs=A, ys=Q, p) -> (n, 1) for exact.reconstruct_points."""
    def ev(xs, ys, p):
        H = HomflyNum(V(xs, p), V(ys, p))
        return H.reduced(lam, steps, sign).a.reshape(-1, 1)
    return ev
