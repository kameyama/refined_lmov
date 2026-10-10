"""Refined invariants of torus knots and their cables from the elliptic Hall algebra (Gorsky–Negut 1304.3328).

Physical parameters (q, t) are those of 1703.05408.  GN's "tilde" operators are evaluated in the world
(q, t_w = 1/t):  P̃_{kn,km} = φ ∘ P_{kn,km}(q, t) ∘ φ^{-1},  φ: p_k ↦ p_k / (1 − t^k).
Vectors are kept in the modified Macdonald basis H̃_μ(q, t_w).

Normalization found in Stage 2 (docs/conventions.md §6):
  the SYT formula (fix3) of GN, read in the fixed-point basis I_μ = H̃_μ / ∇_μ, equals c_k · P̃_{kn,km} with
  c_k = (−1)^{k−1} (1 − t_w^{−k}) / [k]_{q,t_w};  then P̃_{k,0} is exactly multiplication by φ(p_k).
Evaluation: ε(f) = f[p_k ↦ (1 − a^k)/(1 − t^k)], i.e. ε(φ^{-1} f̃) = f̃[p_k ↦ 1 − a^k].
"""
from fractions import Fraction

from .ehall import EHA
from .macdonald import macdonald_P
from .modified import modified_H, solve_in_basis, multiply_p, h
from .ehall import boxes_of
from .partitions import n_stat, conjugate, partitions


class GN:
    def __init__(self, q, t, a, conv):
        self.q, self.t, self.a, self.conv = q, t, a, conv
        self.one, self.zero = conv(Fraction(1)), conv(Fraction(0))
        self.tw = 1 / t
        self.eha = EHA(q, self.tw, self.one)
        self._H, self._P, self._nabla, self._nabla_inv, self._cinv = {}, {}, {}, {}, {}

    # -- bases ----------------------------------------------------------------
    def H(self, n):
        if n not in self._H:
            self._H[n] = modified_H(n, self.q, self.tw, self.conv)
        return self._H[n]

    def P(self, n):
        """Physical Macdonald P_λ(q, t) in the p-basis."""
        if n not in self._P:
            self._P[n] = macdonald_P(n, self.q, self.t, self.conv)
        return self._P[n]

    def nabla(self, mu):
        if mu not in self._nabla:
            self._nabla[mu] = self.q ** n_stat(conjugate(mu)) * self.tw ** n_stat(mu) * self.one
        return self._nabla[mu]

    def c(self, k):
        q, tw = self.q, self.tw
        bracket = (q ** k - 1) * (tw ** k - 1) / ((q - 1) * (tw - 1))
        return (-1) ** (k - 1) * (1 - tw ** (-k)) / bracket

    # -- operators --------------------------------------------------------------
    def op(self, k, n, m, vec):
        """P̃_{kn,km} acting on vec (H̃-basis dict)."""
        vI = {mu: c * self.nabla(mu) for mu, c in vec.items()}
        wI = self.eha.apply(k, n, m, vI)
        ck_inv = self._cinv.get(k)
        if ck_inv is None:
            ck_inv = self._cinv[k] = 1 / self.c(k)
        return {mu: c * self.nabla_inv(mu) * ck_inv for mu, c in wI.items()}

    def nabla_inv(self, mu):
        if mu not in self._nabla_inv:
            self._nabla_inv[mu] = 1 / self.nabla(mu)
        return self._nabla_inv[mu]

    def apply_symmetric(self, f, n, m, vec):
        """f[p_k ↦ P̃_{kn,km}] acting on vec, for f = {ρ: coefficient of p_ρ} (physical symmetric function)."""
        out = {}
        for rho, coeff in f.items():
            v = vec
            for k in rho:
                v = self.op(k, n, m, v)
            for mu, c in v.items():
                out[mu] = out[mu] + coeff * c if mu in out else coeff * c
        return out

    # -- conversions ---------------------------------------------------------------
    def to_p(self, vec):
        """H̃-basis vector -> p-basis dict (tilde world)."""
        out = {}
        for mu, c in vec.items():
            for rho, d in self.H(sum(mu))[mu].items():
                out[rho] = out[rho] + c * d if rho in out else c * d
        return out

    def untilde(self, f_tilde):
        """φ^{-1}: p-basis tilde function -> physical symmetric function."""
        out = {}
        for rho, c in f_tilde.items():
            w = c
            for k in rho:
                w = w * (1 - self.t ** k)
            out[rho] = w
        return out

    def evaluate(self, vec):
        """ε of the physical function whose tilde image is vec:  f̃[p_k ↦ 1 − a^k]."""
        s = self.zero
        for rho, c in self.to_p(vec).items():
            w = c
            for k in rho:
                w = w * (1 - self.a ** k)
            s = s + w
        return s

    # -- knots -----------------------------------------------------------------------
    def torus(self, lam, n, m):
        """Unreduced GN invariant ε(P_λ[p_k ↦ P_{kn,km}] · 1) (no framing correction)."""
        return self.evaluate(self.apply_symmetric(self.P(sum(lam))[lam], n, m, {(): self.one}))

    def pattern(self, lam, n, m):
        """Physical symmetric function F_λ = P_λ[p_k ↦ P_{kn,km}] · 1 in the p-basis."""
        return self.untilde(self.to_p(self.apply_symmetric(self.P(sum(lam))[lam], n, m, {(): self.one})))

    def cable(self, lam, inner, outer):
        """Unreduced invariant of the iterated knot: inner = (n2, m2) pattern, outer = (n1, m1) companion.

        F_λ = P_λ[p_k ↦ P_{k·inner}]·1, then ε(F_λ[p_k ↦ P_{k·outer}]·1)  (DAHA iteration, CD14 (2.12)).
        """
        F = self.pattern(lam, *inner)
        return self.evaluate(self.apply_symmetric(F, *outer, {(): self.one}))

    # -- a kept symbolic: values are lists of coefficients of a^0, a^1, ... -----------------
    def evaluate_apoly(self, vec):
        """ε as a polynomial in a, from the H̃-basis vector directly: H̃_μ[1 − a] = Π_{(r,c)∈μ} (1 − a q^c t_w^r).

        Same result as `evaluate_apoly_via_p`, without building the Macdonald basis in the final degree.
        """
        out = []
        for mu, c in vec.items():
            poly = [self.one]
            for r, col in boxes_of(mu):             # multiply by (1 − a q^col t_w^row)
                x = self.q ** col * self.tw ** r * self.one
                new = poly + [self.zero]
                for i, y in enumerate(poly):
                    new[i + 1] = new[i + 1] - y * x
                poly = new
            while len(out) < len(poly):
                out.append(self.zero)
            for i, y in enumerate(poly):
                out[i] = out[i] + c * y
        return out

    def evaluate_apoly_via_p(self, vec):
        """ε as a polynomial in a: f̃[p_k ↦ 1 − a^k] (list of coefficients), through the p-basis."""
        out = []
        for rho, c in self.to_p(vec).items():
            poly = [1]
            for k in rho:                       # multiply by (1 − a^k)
                new = poly + [0] * k
                for i, x in enumerate(poly):
                    new[i + k] -= x
                poly = new
            while len(out) < len(poly):
                out.append(self.zero)
            for i, x in enumerate(poly):
                if x:
                    out[i] = out[i] + c * x
        return out

    def reduce_apoly(self, lam, apoly):
        """Divide by the unknot P_λ(1,…,t^{N−1}) = Π (t^{l'} − a q^{a'}) / h_λ(q,t)  (exact division in a)."""
        num = [self.one]
        for i, j in boxes_of(lam):              # multiply by (t^i − a q^j)
            new = [self.zero] * (len(num) + 1)
            for d, x in enumerate(num):
                new[d] = new[d] + x * self.t ** i
                new[d + 1] = new[d + 1] - x * self.q ** j
            num = new
        rem = list(apoly)
        quot = [self.zero] * max(len(rem) - len(num) + 1, 1)
        lead_inv = 1 / num[-1]
        for d in range(len(rem) - len(num), -1, -1):
            c = rem[d + len(num) - 1] * lead_inv
            quot[d] = c
            for i, x in enumerate(num):
                rem[d + i] = rem[d + i] - c * x
        if any(r != 0 for r in rem):
            raise ArithmeticError("unreduced invariant is not divisible by the unknot")
        h_lam = h(lam, self.q, self.t)
        return [c * h_lam for c in quot]

    def torus_reduced_apoly(self, lam, n, m):
        vec = self.apply_symmetric(self.P(sum(lam))[lam], n, m, {(): self.one})
        return self.reduce_apoly(lam, self.evaluate_apoly(vec))

    def iterated_reduced_apoly(self, lam, slopes):
        """Reduced invariant of an iterated cable; slopes = [outermost pattern, ..., innermost pattern, companion].

        F = P_λ[p_k ↦ P_{k·slopes[0]}]·1, then F ↦ F[p_k ↦ P_{k·slope}]·1 for the middle slopes, then ε at the companion
        (DAHA iteration, CD14 (2.12)). With two slopes this is cable_reduced_apoly.
        """
        F = self.pattern(lam, *slopes[0])
        for sl in slopes[1:-1]:
            F = self.untilde(self.to_p(self.apply_symmetric(F, *sl, {(): self.one})))
        vec = self.apply_symmetric(F, *slopes[-1], {(): self.one})
        return self.reduce_apoly(lam, self.evaluate_apoly(vec))

    def cable_reduced_apoly(self, lam, inner, outer):
        F = self.pattern(lam, *inner)
        vec = self.apply_symmetric(F, *outer, {(): self.one})
        return self.reduce_apoly(lam, self.evaluate_apoly(vec))
