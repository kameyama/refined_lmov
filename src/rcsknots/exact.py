"""Reconstruct exact Laurent polynomials f(a, q, t) from evaluations mod p.

`func(G)` receives a `GN` object whose q, t are vectors of evaluation points and returns the
value as a list of coefficients of a^0, a^1, ... (each a `field.V`).  The result is
{(e_a, e_q, e_t): integer}.  Interpolation uses prime P31 (and P31B with the Chinese remainder theorem only
when the single-prime lift fails); the lifted integer polynomial is checked at random points mod a third
prime P31C that were not used for interpolation.
"""
from fractions import Fraction
import time

import numpy as np

from .field import V, P31, P31B, P31C, P31D, make_conv
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
                out[i] = np.load(f, mmap_mode="r")      # 読むのは下の詰め直しのときだけ（全体を 2 重に持たない）
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
    res = np.zeros((len(qs), D), dtype=np.int64)
    for i, x in enumerate(out):
        res[i * bs:i * bs + len(x), :x.shape[1]] = x
    return res


INTERP_BLOCK = 16      # 補間を a の係数の列 16 本ずつに分ける（並列の単位。1 まとまりの作業メモリは格子 × 16 × 8 バイトの数倍）
_GRID_STATE = {}       # 補間する値と点（fork したワーカーと共有、読むだけ）


def _interp_block(d0):
    """列 d0..d0+INTERP_BLOCK の 2 変数補間。0 でない係数を (it, iq, d, 値 mod p) の配列で返す。"""
    st = _GRID_STATE
    p = st["p"]
    v = st["vals"][:, :, d0:d0 + INTERP_BLOCK] * st["wq"][:, None, None] % p * st["wt"][None, :, None] % p
    cq = newton_to_monomial(st["qs"], v, p)                                 # (Nq, Nt, b)
    c = newton_to_monomial(st["ts"], np.transpose(cq, (1, 0, 2)), p)       # (Nt, Nq, b)
    it, iq, d = np.nonzero(c)
    return it, iq, d + d0, c[it, iq, d]


def _interp_grid(vals, qs, ts, qlo, tlo, p, workers=1):
    """格子 (Nq, Nt, D) の値から、x^{qlo}·y^{tlo} のずれを戻した係数 {(d, e_q, e_t): c mod p}（0 でないものだけ）。

    列（a の係数）ごとに独立なので、INTERP_BLOCK 本ずつに分けて fork のプロセスプールで並列に計算する。
    """
    wq = np.array([pow(int(x), (-qlo) % (p - 1), p) for x in qs], dtype=np.int64)
    wt = np.array([pow(int(x), (-tlo) % (p - 1), p) for x in ts], dtype=np.int64)
    for xs in (qs, ts):                    # 逆 Vandermonde 行列を先に作り、ワーカーに引き継ぐ
        newton_to_monomial(xs, np.zeros((len(xs), 1), dtype=np.int64), p)
    _GRID_STATE.update(vals=vals, qs=qs, ts=ts, wq=wq, wt=wt, p=p)
    starts = list(range(0, vals.shape[2], INTERP_BLOCK))
    if workers > 1 and len(starts) > 1:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(min(workers, len(starts))) as pool:
            parts = pool.map(_interp_block, starts)
    else:
        parts = [_interp_block(d0) for d0 in starts]
    _GRID_STATE.clear()
    return {(int(d), int(iq) + qlo, int(it) + tlo): int(c)
            for it_, iq_, d_, c_ in parts for it, iq, d, c in zip(it_, iq_, d_, c_)}


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
    log(f"  grid {Nq} x {Nt} = {Nq * Nt} points")
    primes, coeffs = (P31, P31B, P31D), []
    for i, p in enumerate(primes):
        qs, ts = _distinct(rng, Nq, p), _distinct(rng, Nt, p)
        tag = f"grid_{p}_{qlo}_{qhi}_{tlo}_{thi}_s{seed}"
        coeffs.append(_interp_cached(func, qs, ts, qlo, tlo, p, workers, cache, tag))
        # CRT（Garner）：対称な代表を法 Π p で取る。係数が小さければ素数 1 つで足り、独立な点で検証できれば終わり
        M, lifted = primes[0], dict(coeffs[0])
        for pj, cj in zip(primes[1:i + 1], coeffs[1:]):
            inv = pow(M, pj - 2, pj)
            lifted = {k: x + M * (((cj.get(k, 0) - x) * inv) % pj) for k in lifted.keys() | cj.keys()
                      for x in (lifted.get(k, 0),)}
            M *= pj
        poly = {k: (x - M if x > M // 2 else x) for k, x in lifted.items() if x}
        try:
            verify(func, poly, seed + 1000)
            log(f"  {len(poly)} terms, {i + 1} prime{'s' if i else ''}, {time.time() - t0:.0f}s")
            return poly
        except ArithmeticError:
            if i + 1 == len(primes):
                raise
            log(f"  {i + 1} prime{'s' if i else ''} not enough; using prime {i + 2}")


def _interp_cached(func, qs, ts, qlo, tlo, p, workers, cache, tag):
    """格子点で評価して補間した係数 {(d, e_q, e_t): c mod p}。cache があれば補間の結果も保存し、再開のときは読み込む。"""
    f = cache / f"coef_{tag}.npz" if cache is not None else None
    if f is not None and f.exists():
        z = np.load(f)
        return {(int(d), int(eq), int(et)): int(c) for d, eq, et, c in zip(z["d"], z["eq"], z["et"], z["c"])}
    Qg, Tg = np.meshgrid(qs, ts, indexing="ij")
    vals = _evaluate(func, Qg.ravel(), Tg.ravel(), p, workers, cache, tag).reshape(len(qs), len(ts), -1)
    out = _interp_grid(vals, qs, ts, qlo, tlo, p, workers)
    if f is not None:
        k = np.array(list(out), dtype=np.int64).reshape(-1, 3)
        np.savez(f, d=k[:, 0], eq=k[:, 1], et=k[:, 2], c=np.array(list(out.values()), dtype=np.int64))
    return out


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
