"""Numerical evaluation of the refined reformulated invariants at points (Q, T), with A = a^{1/2} kept symbolic.

A knot is given by its reduced invariants {color: {(e_A, e_Q, e_T): c}} (normalized as in 1703.05408).
Values are Laurent polynomials in A with coefficients `field.V` (class APoly).
"""
import numpy as np
import sympy as sp

from .bps import f_formulas, M_matrix, C, C1, Q as SQ, T as ST, COLORS3
from .field import V
from .partitions import boxes, arm, leg


class APoly:
    """Laurent polynomial in A with field coefficients: {exponent: V}."""

    def __init__(self, terms):
        self.terms = {e: c for e, c in terms.items()}

    def __add__(self, o):
        if not isinstance(o, APoly):
            o = APoly({0: o})
        out = dict(self.terms)
        for e, c in o.terms.items():
            out[e] = out[e] + c if e in out else c
        return APoly(out)

    __radd__ = __add__

    def __neg__(self):
        return APoly({e: -c for e, c in self.terms.items()})

    def __sub__(self, o):
        return self + (-o if isinstance(o, APoly) else -o)

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, APoly):
            return APoly({e: c * o for e, c in self.terms.items()})
        out = {}
        for e1, c1 in self.terms.items():
            for e2, c2 in o.terms.items():
                e = e1 + e2
                out[e] = out[e] + c1 * c2 if e in out else c1 * c2
        return APoly(out)

    __rmul__ = __mul__

    def __truediv__(self, s):
        return self * (1 / s)

    def __pow__(self, k):
        out = APoly({0: 1})
        for _ in range(int(k)):
            out = out * self
        return out


def eval_expr(expr, env):
    """Evaluate a sympy expression; env maps symbols to V or APoly."""
    if expr.is_Symbol:
        return env[expr]
    if expr.is_Integer:
        return int(expr)
    if expr.is_Rational:
        return env["_rat"](expr)
    if expr.is_Add:
        out = 0
        for x in expr.args:
            y = eval_expr(x, env)
            out = y + out if isinstance(y, APoly) else out + y
        return out
    if expr.is_Mul:
        out = 1
        for x in expr.args:
            y = eval_expr(x, env)
            out = y * out if isinstance(y, APoly) else out * y
        return out
    if expr.is_Pow:
        base, e = expr.args
        b = eval_expr(base, env)
        if e.is_Integer:
            if int(e) < 0:
                if isinstance(b, APoly):
                    raise ValueError("negative power of an A-polynomial")
                return b ** int(e)
            return b ** int(e)
        raise ValueError(f"unsupported exponent {e}")
    raise ValueError(f"unsupported expression {expr}")


def _power_table(x, exps):
    """{e: x^e} for the exponents in exps, by successive multiplication (one inverse for negative e)."""
    lo, hi = min(exps + [0]), max(exps + [0])
    tab = {0: x ** 0}
    cur = tab[0]
    for e in range(1, hi + 1):
        cur = cur * x
        tab[e] = cur
    if lo < 0:
        xi, cur = x.inv(), tab[0]
        for e in range(-1, lo - 1, -1):
            cur = cur * xi
            tab[e] = cur
    return tab


def laurent_at(poly, Qv, Tv, d=1):
    """{(eA,eQ,eT): c} evaluated at Q^d, T^d with A ↦ A^d, as an APoly."""
    qt = _power_table(Qv, [d * eq for _, eq, _ in poly])
    tt = _power_table(Tv, [d * et for _, _, et in poly])
    out = {}
    for (ea, eq, et), c in poly.items():
        val = c * qt[d * eq] * tt[d * et]
        e = d * ea
        out[e] = out[e] + val if e in out else val
    return APoly(out)


def unknot_sym(lam, Qv, Tv, d=1):
    """(t/a)^{|λ|/2} Π (t^{l'} − a q^{a'})/(1 − q^a t^{l+1}) at (a^d, q^d, t^d), as an APoly in A."""
    q, t = Qv ** (2 * d), Tv ** (2 * d)
    k = sum(lam)
    out = APoly({-d * k: Tv ** (d * k)})
    for i, j in boxes(lam):
        factor = APoly({0: t ** i, 2 * d: -(q ** j)})
        out = out * factor / (1 - q ** arm(lam, i, j) * t ** (leg(lam, i, j) + 1))
    return out


class Reformulated:
    """f^q, f^t, \\hat f for one knot at points (Qv, Tv)."""

    def __init__(self, reduced, Qv, Tv, conv):
        self.reduced, self.Qv, self.Tv = reduced, Qv, Tv
        self.env = {SQ: Qv, ST: Tv, "_rat": lambda r: conv(__import__("fractions").Fraction(int(r.p), int(r.q)))}
        for lam in reduced:
            self.env[C(lam)] = unknot_sym(lam, Qv, Tv) * laurent_at(reduced[lam], Qv, Tv)
        for d in (2, 3):
            self.env[C1(d)] = unknot_sym((1,), Qv, Tv, d) * laurent_at(reduced[(1,)], Qv, Tv, d)

    def f(self, mu):
        fq, ft = f_formulas(3)[mu]
        return eval_expr(fq, self.env), eval_expr(ft, self.env)

    def hat_f(self, n):
        """(\\hat f_ρ for ρ ⊢ n) = M(t)^{-1} f^q, and the residual f^t − M(q^{-1}) \\hat f (should vanish)."""
        from .partitions import partitions
        parts = partitions(n)
        Minv = M_matrix(n, ST).inv()
        Mq = M_matrix(n, 1 / SQ)               # M(q^{-1}): t^{1/2} ↦ q^{-1/2}
        fq = [self.f(mu)[0] for mu in parts]
        ft = [self.f(mu)[1] for mu in parts]
        hat = []
        for i in range(len(parts)):
            s = 0
            for j in range(len(parts)):
                s = s + fq[j] * eval_expr(sp.together(Minv[i, j]), self.env)
            hat.append(s)
        resid = []
        for i in range(len(parts)):
            s = ft[i]
            for j in range(len(parts)):
                s = s - hat[j] * eval_expr(Mq[i, j], self.env)
            resid.append(s)
        return dict(zip(parts, hat)), dict(zip(parts, resid))
