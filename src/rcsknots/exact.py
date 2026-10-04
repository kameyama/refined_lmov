"""Reconstruct exact Laurent polynomials f(a, q, t) from evaluations mod p.

`func(G)` receives a `GN` object whose q, t are vectors of evaluation points and returns the
value as a list of coefficients of a^0, a^1, ... (each a `field.V`).  The result is
{(e_a, e_q, e_t): integer}.  Interpolation uses prime P31; the lifted integer polynomial is then
checked at random points mod a second prime that were not used for interpolation.
"""
from fractions import Fraction
import time

import numpy as np

from .field import V, P31, P31B, make_conv
from .gn import GN
from .interp import newton_to_monomial, symmetric_lift

BATCH = 4096


def gn_evaluator(func):
    """Wrap func(GN) -> list of a-coefficients into an evaluator(xs, ys, p) -> array (n, D)."""
    def ev(qv, tv, p):
        n = len(qv)
        G = GN(V(qv, p), V(tv, p), V(np.full(n, 3), p), make_conv(p, n))
        return np.stack([c.a for c in func(G)], axis=1)
    return ev


def _evaluate(evaluator, qs, ts, p):
    """Values (len(qs), D) of the evaluator at the points (qs[i], ts[i]), in batches."""
    out = [evaluator(qs[s:s + BATCH], ts[s:s + BATCH], p) for s in range(0, len(qs), BATCH)]
    D = max(x.shape[1] for x in out)
    return np.concatenate([np.pad(x, ((0, 0), (0, D - x.shape[1]))) for x in out])


def _distinct(rng, n, p):
    vals = set()
    while len(vals) < n:
        vals.update(int(v) for v in rng.integers(2, p - 1, n))
    return np.array(sorted(vals)[:n], dtype=np.int64)


def _laurent_1d(xs, vals, shift, p):
    """Interpolate x^{shift}·f(x); return lifted coefficients (N, D) indexed by exponent + shift."""
    w = np.array([pow(int(x), shift, p) for x in xs], dtype=np.int64)
    return newton_to_monomial(xs, vals * w[:, None] % p, p)


def _probe_range(func, axis, rng, p, log):
    """Exponent range of f in q (axis=0) or t (axis=1) at a random value of the other variable."""
    other = int(rng.integers(2, p - 1))
    N = 64
    while True:
        xs = _distinct(rng, N + 3, p)
        others = np.full(N + 3, other, dtype=np.int64)
        qs, ts = (xs, others) if axis == 0 else (others, xs)
        vals = _evaluate(func, qs, ts, p)
        shift = N // 2
        c = _laurent_1d(xs[:N], vals[:N], shift, p)
        # check the three extra points
        ok = True
        for i in range(N, N + 3):
            x = int(xs[i])
            pred = np.zeros(vals.shape[1], dtype=np.int64)
            for k in range(N - 1, -1, -1):
                pred = (pred * x + c[k]) % p
            pred = pred * pow(x, p - 1 - shift % (p - 1), p) % p
            ok &= bool(np.all(pred == vals[i]))
        nz = np.nonzero(np.any(c != 0, axis=1))[0]
        if ok and len(nz) and nz[0] > 0 and nz[-1] < N - 1:
            lo, hi = int(nz[0]) - shift, int(nz[-1]) - shift
            log(f"  axis {'qt'[axis]}: exponents {lo}..{hi} (probe N={N})")
            return lo, hi
        N *= 2
        if N > 4096:
            raise RuntimeError("exponent range too large")


def reconstruct(func, seed=0, log=print):
    """Exact polynomial of a GN function func(G) -> a-coefficients; keys (e_a, e_q, e_t)."""
    return reconstruct_points(gn_evaluator(func), seed, log)


def reconstruct_points(func, seed=0, log=print):
    """Exact polynomial from evaluator func(xs, ys, p) -> (n, D); keys (index d, e_x, e_y)."""
    rng = np.random.default_rng(seed)
    t0 = time.time()
    qlo, qhi = _probe_range(func, 0, rng, P31, log)
    tlo, thi = _probe_range(func, 1, rng, P31, log)
    Nq, Nt = qhi - qlo + 1, thi - tlo + 1
    qs, ts = _distinct(rng, Nq, P31), _distinct(rng, Nt, P31)
    Qg, Tg = np.meshgrid(qs, ts, indexing="ij")
    log(f"  grid {Nq} x {Nt} = {Nq * Nt} points")
    vals = _evaluate(func, Qg.ravel(), Tg.ravel(), P31).reshape(Nq, Nt, -1)
    # remove the Laurent shifts, then interpolate in q and in t
    wq = np.array([pow(int(x), (-qlo) % (P31 - 1), P31) for x in qs], dtype=np.int64)
    wt = np.array([pow(int(x), (-tlo) % (P31 - 1), P31) for x in ts], dtype=np.int64)
    vals = vals * wq[:, None, None] % P31 * wt[None, :, None] % P31
    cq = newton_to_monomial(qs, vals, P31)                              # (Nq, Nt, D)
    ct = newton_to_monomial(ts, np.transpose(cq, (1, 0, 2)), P31)       # (Nt, Nq, D)
    ct = symmetric_lift(ct, P31)
    poly = {}
    for it, iq, d in zip(*np.nonzero(ct)):
        poly[(int(d), int(iq) + qlo, int(it) + tlo)] = int(ct[it, iq, d])
    log(f"  {len(poly)} terms, {time.time() - t0:.0f}s")
    verify(func, poly, seed + 1000)
    return poly


def evaluate_poly(poly, qv, tv, p):
    """poly at points (qv, tv): array (n, D+1) of a-coefficients mod p."""
    D = max(k[0] for k in poly) + 1
    out = np.zeros((len(qv), D), dtype=np.int64)
    for (d, eq, et), c in poly.items():
        term = np.array([c % p * pow(int(x), eq % (p - 1), p) % p * pow(int(y), et % (p - 1), p) % p
                         for x, y in zip(qv, tv)], dtype=np.int64)
        out[:, d] = (out[:, d] + term) % p
    return out


def verify(func, poly, seed, npts=8):
    """Check the integer polynomial at random points mod the second prime (not used for interpolation)."""
    rng = np.random.default_rng(seed)
    qv, tv = _distinct(rng, npts, P31B), _distinct(rng, npts, P31B)
    got = _evaluate(func, qv, tv, P31B)
    pred = evaluate_poly(poly, qv, tv, P31B)
    D = max(got.shape[1], pred.shape[1])
    got = np.pad(got, ((0, 0), (0, D - got.shape[1])))
    pred = np.pad(pred, ((0, 0), (0, D - pred.shape[1])))
    if not np.array_equal(got, pred):
        raise ArithmeticError("verification at independent points failed")
