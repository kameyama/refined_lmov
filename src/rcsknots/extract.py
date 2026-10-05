"""From reduced refined invariants to \\hat f (exact) and the BPS numbers \\hat N_{ρ,g,β,2J_r}."""
from fractions import Fraction

import numpy as np

from .bps import to_bold, g_expand
from .bpsnum import Reformulated
from .exact import reconstruct_points
from .field import V, make_conv
from .partitions import partitions

A_OFFSET = 200          # A-exponents are stored shifted by this amount during interpolation


def _arr(c, k, p):
    return np.broadcast_to(c.a, (k,)) if isinstance(c, V) else np.full(k, int(c) % p, dtype=np.int64)


def hat_f_exact(reduced, n, log=print, check_ft=True):
    """{ρ ⊢ n: {(e_A, e_Q, e_T): c}} exact \\hat f_ρ (no a-grading shift applied)."""
    parts = partitions(n)
    out = {}

    def evaluator_for(rho):
        def ev(xs, ys, p):
            k = len(xs)
            R = Reformulated(reduced, V(xs, p), V(ys, p), make_conv(p, k))
            hat, resid = R.hat_f(n)
            if check_ft:
                for r in resid.values():
                    if any(np.any(_arr(c, k, p) != 0) for c in r.terms.values()):
                        raise ArithmeticError("f^t ≠ M(q^{-1}) \\hat f at some point")
            arr = np.zeros((k, 2 * A_OFFSET + 1), dtype=np.int64)
            for e, c in hat[rho].terms.items():
                arr[:, e + A_OFFSET] = (arr[:, e + A_OFFSET] + _arr(c, k, p)) % p
            return arr
        return ev

    for rho in parts:
        log(f" hat f_{rho}:")
        try:
            raw = reconstruct_points(evaluator_for(rho), seed=hash(rho) % 1000, log=log)
        except RuntimeError:
            # no Laurent polynomial of the probed size fits: \hat f_ρ has denominators (refined LMOV form fails)
            out[rho] = None
            continue
        out[rho] = {(d - A_OFFSET, eq, et): c for (d, eq, et), c in raw.items()}
    return out


# a-grading shift of \\hat f_ρ before the g-expansion, in units of A = a^{1/2} ("up to a-grading shift by ±1/2"):
#   torus notebook (T(2,3) tables etc.):  |ρ|=1: a^{+1/2}, |ρ|=2: 1, |ρ|=3: a^{−1/2}
#   non-torus notebook (4_1 tables etc.): |ρ|=1: a^{−1/2}, |ρ|=2: 1, |ρ|=3: a^{+1/2}
A_SHIFT = {"torus": {1: 1, 2: 0, 3: -1}, "nontorus": {1: -1, 2: 0, 3: 1}}


def bps_numbers(hat_f_rho, rho, convention="torus"):
    """{(g, β, F=2J_r): c} from the g-expansion of the shifted \\hat f_ρ in bold variables."""
    return g_expand(to_bold(hat_f_rho, A_SHIFT[convention][sum(rho)]))


def reformulated_exact(reduced, n, log=print):
    """All reformulated invariants of size n, exactly, in KN's lowercase half-variables (A, Q, T).

    Returns {("fq", μ): poly, ("ft", μ): poly, ("hat", ρ): poly} (no a-shift, no change to bold variables).
    The check f^t = M(q^{-1}) \\hat f is done at every evaluation point.
    """
    parts = partitions(n)
    keys = [("fq", m) for m in parts] + [("ft", m) for m in parts] + [("hat", r) for r in parts]
    span = {}

    def raw_values(xs, ys, p):
        k = len(xs)
        R = Reformulated(reduced, V(xs, p), V(ys, p), make_conv(p, k))
        hat, resid = R.hat_f(n)
        for r in resid.values():
            if any(np.any(_arr(c, k, p) != 0) for c in r.terms.values()):
                raise ArithmeticError("f^t ≠ M(q^{-1}) \\hat f at some point")
        vals = {}
        for m in parts:
            fq, ft = R.f(m)
            vals[("fq", m)], vals[("ft", m)] = fq, ft
        for r in parts:
            vals[("hat", r)] = hat[r]
        return vals, k

    # A-exponent window from one small evaluation
    rng = np.random.default_rng(7)
    xs, ys = rng.integers(2, 10**6, 4), rng.integers(2, 10**6, 4)
    vals, _ = raw_values(xs, ys, 2**31 - 1)
    lo = min(min(v.terms) for v in vals.values()) - 2
    hi = max(max(v.terms) for v in vals.values()) + 2
    width = hi - lo + 1

    def ev(xs, ys, p):
        vals, k = raw_values(xs, ys, p)
        arr = np.zeros((k, width * len(keys)), dtype=np.int64)
        for i, key in enumerate(keys):
            for e, c in vals[key].terms.items():
                arr[:, i * width + e - lo] = (arr[:, i * width + e - lo] + _arr(c, k, p)) % p
        return arr

    raw = reconstruct_points(ev, seed=11 + n, log=log)
    out = {key: {} for key in keys}
    for (d, eq, et), c in raw.items():
        i, e = divmod(d, width)
        out[keys[i]][(e + lo, eq, et)] = c
    return out
