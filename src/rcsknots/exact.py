"""Reconstruct exact Laurent polynomials f(a, q, t) from evaluations mod p.

`func(G)` receives a `GN` object whose q, t are vectors of evaluation points and returns the
value as a list of coefficients of a^0, a^1, ... (each a `field.V`).  The result is
{(e_a, e_q, e_t): integer}.  Interpolation uses prime P31; the lifted integer polynomial is then
checked at random points mod a second prime that were not used for interpolation.
"""
from fractions import Fraction
import time

import numpy as np

from .field import V, P31, P31B, P31C, make_conv
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


_FORK_STATE = {}       # evaluator and points shared with forked workers


def _batch(i):
    st = _FORK_STATE
    bs = st["bs"]
    s = i * bs
    out = st["evaluator"](st["qs"][s:s + bs], st["ts"][s:s + bs], st["p"])
    if st["cache"] is not None:
        np.save(st["cache"] / f"{st['tag']}_{i}.npy", out)
    return out


def _evaluate(evaluator, qs, ts, p, workers=1, cache=None, tag=""):
    """Values (len(qs), D) of the evaluator at the points (qs[i], ts[i]), in batches.

    cache：まとまりごとの結果を保存するディレクトリ（Path）。保存済みのまとまりは読み込む（再開）。
    workers > 1：未計算のまとまりを fork のプロセスプールで並列に計算する。
    """
    # 並列のときは、まとまりを細かくして（ワーカーあたり 4 つ程度、最小 256 点）負荷をならす
    bs = BATCH if workers <= 1 else min(BATCH, max(256, -(-len(qs) // (4 * workers))))
    nb = (len(qs) + bs - 1) // bs
    tag = f"{tag}_b{bs}"
    _FORK_STATE.update(evaluator=evaluator, qs=qs, ts=ts, p=p, cache=cache, tag=tag, bs=bs)
    out = [None] * nb
    if cache is not None:
        cache.mkdir(parents=True, exist_ok=True)
        for i in range(nb):
            f = cache / f"{tag}_{i}.npy"
            if f.exists():
                out[i] = np.load(f)
    todo = [i for i in range(nb) if out[i] is None]
    if workers > 1 and len(todo) > 1:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(min(workers, len(todo))) as pool:
            for i, x in zip(todo, pool.map(_batch, todo)):
                out[i] = x
    else:
        for i in todo:
            out[i] = _batch(i)
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


def _probe_range(func, axis, rng, p, log, workers=1):
    """Exponent range of f in q (axis=0) or t (axis=1) at a random value of the other variable."""
    other = int(rng.integers(2, p - 1))
    N = 64
    while True:
        xs = _distinct(rng, N + 3, p)
        others = np.full(N + 3, other, dtype=np.int64)
        qs, ts = (xs, others) if axis == 0 else (others, xs)
        vals = _evaluate(func, qs, ts, p, workers)
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


def reconstruct(func, seed=0, log=print, workers=1, cache=None):
    """Exact polynomial of a GN function func(G) -> a-coefficients; keys (e_a, e_q, e_t)."""
    return reconstruct_points(gn_evaluator(func), seed, log, workers, cache)


def reconstruct_points(func, seed=0, log=print, workers=1, cache=None):
    """Exact polynomial from evaluator func(xs, ys, p) -> (n, D); keys (index d, e_x, e_y).

    workers, cache：格子点での評価の並列数と、まとまりごとの保存先（再開用、_evaluate を参照）。
    格子点は seed で決まるので、同じ seed・同じ入力なら保存した結果をそのまま再利用できる。
    """
    rng = np.random.default_rng(seed)
    t0 = time.time()
    qlo, qhi = _probe_range(func, 0, rng, P31, log, workers)
    tlo, thi = _probe_range(func, 1, rng, P31, log, workers)
    Nq, Nt = qhi - qlo + 1, thi - tlo + 1
    log(f"  grid {Nq} x {Nt} = {Nq * Nt} points, 2 primes")
    coeffs = []
    for p in (P31, P31B):
        qs, ts = _distinct(rng, Nq, p), _distinct(rng, Nt, p)
        Qg, Tg = np.meshgrid(qs, ts, indexing="ij")
        vals = _evaluate(func, Qg.ravel(), Tg.ravel(), p, workers, cache,
                         f"grid_{p}_{qlo}_{qhi}_{tlo}_{thi}_s{seed}").reshape(Nq, Nt, -1)
        # remove the Laurent shifts, then interpolate in q and in t
        wq = np.array([pow(int(x), (-qlo) % (p - 1), p) for x in qs], dtype=np.int64)
        wt = np.array([pow(int(x), (-tlo) % (p - 1), p) for x in ts], dtype=np.int64)
        vals = vals * wq[:, None, None] % p * wt[None, :, None] % p
        cq = newton_to_monomial(qs, vals, p)                              # (Nq, Nt, D)
        coeffs.append(newton_to_monomial(ts, np.transpose(cq, (1, 0, 2)), p))   # (Nt, Nq, D)
    D = max(c.shape[2] for c in coeffs)
    c1, c2 = (np.pad(c, ((0, 0), (0, 0), (0, D - c.shape[2]))) for c in coeffs)
    # Chinese remainder theorem: symmetric representative modulo P31·P31B (≈ 4.6e18)
    inv = pow(P31, P31B - 2, P31B)
    poly = {}
    for it, iq, d in zip(*np.nonzero((c1 != 0) | (c2 != 0))):
        a1, a2 = int(c1[it, iq, d]), int(c2[it, iq, d])
        x = a1 + P31 * (((a2 - a1) * inv) % P31B)
        M = P31 * P31B
        if x > M // 2:
            x -= M
        if x:
            poly[(int(d), int(iq) + qlo, int(it) + tlo)] = x
    log(f"  {len(poly)} terms, {time.time() - t0:.0f}s")
    verify(func, poly, seed + 1000)
    return poly


def evaluate_poly(poly, qv, tv, p):
    """poly at points (qv, tv): array (n, D+1) of a-coefficients mod p（項についてはベクトル化）。"""
    D = max(k[0] for k in poly) + 1
    keys = np.array(list(poly), dtype=np.int64)
    coef = np.array([c % p for c in poly.values()], dtype=np.int64)
    out = np.zeros((len(qv), D), dtype=np.int64)
    for n, (x, y) in enumerate(zip(qv, tv)):
        xp = np.array([pow(int(x), int(e) % (p - 1), p) for e in range(keys[:, 1].min(), keys[:, 1].max() + 1)], dtype=np.int64)
        yp = np.array([pow(int(y), int(e) % (p - 1), p) for e in range(keys[:, 2].min(), keys[:, 2].max() + 1)], dtype=np.int64)
        term = coef * xp[keys[:, 1] - keys[:, 1].min()] % p * yp[keys[:, 2] - keys[:, 2].min()] % p
        out[n] = _sum_mod(keys[:, 0], term, D, p)
    return out


def _sum_mod(idx, vals, D, p):
    """Σ vals grouped by idx, mod p（int64 のあふれを避けて 2^20 項ずつ足す）。"""
    out = np.zeros(D, dtype=np.int64)
    for s in range(0, len(vals), 1 << 20):
        part = np.zeros(D, dtype=np.int64)
        np.add.at(part, idx[s:s + (1 << 20)], vals[s:s + (1 << 20)])
        out = (out + part % p) % p
    return out


def verify(func, poly, seed, npts=8):
    """Check the integer polynomial at random points mod a third prime (not used for interpolation)."""
    rng = np.random.default_rng(seed)
    qv, tv = _distinct(rng, npts, P31C), _distinct(rng, npts, P31C)
    got = _evaluate(func, qv, tv, P31C)
    pred = evaluate_poly(poly, qv, tv, P31C)
    D = max(got.shape[1], pred.shape[1])
    got = np.pad(got, ((0, 0), (0, D - got.shape[1])))
    pred = np.pad(pred, ((0, 0), (0, D - pred.shape[1])))
    if not np.array_equal(got, pred):
        raise ArithmeticError("verification at independent points failed")
