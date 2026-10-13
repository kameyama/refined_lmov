"""トーラス結び目の 4 箱以上の refined LMOV 不変量：\\hat f_ρ（|ρ| = N の全 ρ）を厳密に求め、G 展開の前後の正整数性を判定する。

手順（rcsknots.bpsgen）：GN の値 × conventions.md §13 の正規化式で rCS_λ（|λ| ≤ N）を評価点ごとに作り、plethystic log と
M(t)^{-1}（類関数の基底で対角）で \\hat f_ρ を得る。f^{\\bar t} = M(q^{-1}) \\hat f は全ての評価点で検査する。
補間：A は a = A² の多項式（窓 [N·lo₁, N·hi₁]、lo₁, hi₁ は [1] の A 次数）、Q・T は指数の偶奇が一様なので q = Q², t = T² の Laurent 多項式。
法 p の値を素数ごとに補間し CRT で持ち上げ、補間に使っていない素数 P31C の乱数の点で検証する。
使い方：python scripts/framed_knot/highbox_torus.py m n N [workers]   結果は data/framed_knot/highbox/T{m}_{n}_N{N}.json.xz、ログは logs/framed_knot/highbox_torus.log に追記。
"""
import json
import lzma
import resource
import sys
import time
from pathlib import Path

import numpy as np

from rcsknots import exact
from rcsknots.bps import g_expand, to_bold
from rcsknots.bpsgen import Points, torus_rcs_factory
from rcsknots.field import V, P31, P31B, P31C, P31D
from rcsknots.interp import newton_to_monomial
from rcsknots.partitions import partitions

OUT = Path("data/framed_knot/highbox")
LOG = Path("logs/framed_knot/highbox_torus.log")


def a_shift(N):
    """a 次数のずらし（A の単位）。N ≤ 3 は nontorus 規約（extract.A_SHIFT）、奇数の N はその周期的な延長。正値性には影響しない。"""
    return 0 if N % 2 == 0 else (-1 if N % 4 == 1 else 1)


class Problem:
    """|ρ| = N の全ての \\hat f_ρ を補間する問題。rcs_factory() は bpsgen.Points に渡す rCS の関数を返す。

    A の窓 [alo, ahi]（偶奇は alo と同じ）。post(Qs, Ts, p) を与えると、各点の値にその数（長さ len(Qs) の配列）を掛ける。
    """

    def __init__(self, N, rcs_factory, alo, ahi, post=None):
        self.N, self.rcs_factory, self.post = N, rcs_factory, post
        self.rhos = partitions(N)
        self.alo, self.ahi = alo, ahi
        self.Na = (self.ahi - self.alo) // 2 + 1
        self.e0 = (0, 0)                                                     # Q, T の指数の偶奇（probe で決める）
        self.step = (2, 2)                                                   # 偶奇が一様なら Q², T² で補間する（混ざれば 1）

    def raw(self, Qs, Ts, p):
        """(len, len(rhos)·Na)：各点の A^{alo + 2k} の係数。"""
        k = len(Qs)
        rng = np.random.default_rng(int(Qs[0]) % 1000)
        avals = []
        while len(avals) < self.Na + 1:                       # 最後の 1 点は窓の検査用
            x = int(rng.integers(2, p - 1))
            if all((x * x - y * y) % p for y in avals):
                avals.append(x)
        pts = Points(self.rcs_factory(), V(Qs, p), V(Ts, p), avals)
        hat, ok = pts.hat_f(self.N)
        if not ok:
            raise ArithmeticError("f^t ≠ M(q^{-1}) \\hat f at some point")
        A = np.array(avals, dtype=np.int64)
        shift = np.array([pow(int(x), (-self.alo) % (p - 1), p) for x in A], dtype=np.int64)
        cols = []
        for rho in self.rhos:
            vals = hat[self.N][rho].a.reshape(k, self.Na + 1) * shift[None, :] % p     # A^{−alo} \\hat f
            c = newton_to_monomial(A[:-1] * A[:-1] % p, vals[:, :-1].T, p)              # 変数 a = A²、(Na, k)
            x, pred = int(A[-1]) ** 2 % p, np.zeros(k, dtype=np.int64)
            for j in range(self.Na - 1, -1, -1):
                pred = (pred * x + c[j]) % p
            if not np.array_equal(pred, vals[:, -1]):
                raise ArithmeticError("A の窓が足りない")
            cols.append(c.T)
        out = np.concatenate(cols, axis=1)
        if self.post is not None:
            out = out * (self.post(Qs, Ts, p) % p)[:, None] % p
        return out

    def ev_parity(self, Qs, Ts, p):
        """偶奇をそろえた値：Q^{−e0} T^{−f0} × raw。"""
        out = self.raw(Qs, Ts, p)
        w = np.array([pow(int(q), (-self.e0[0]) % (p - 1), p) * pow(int(t), (-self.e0[1]) % (p - 1), p) % p
                      for q, t in zip(Qs, Ts)], dtype=np.int64)
        return out * w[:, None] % p


_STATE = {}


def _chunk(i):
    st = _STATE
    return st["func"](st["Qs"][st["cuts"][i]:st["cuts"][i + 1]], st["Ts"][st["cuts"][i]:st["cuts"][i + 1]], st["p"])


def evaluate(func, Qs, Ts, p, workers):
    """func を点の列で評価する。点を workers 個程度の大きなまとまりに分ける（numpy のベクトル化の効率のため、1 まとまり最大 4096 点）。"""
    n = len(Qs)
    nb = max(workers, -(-n // 4096)) if workers > 1 else -(-n // 4096)
    nb = max(1, min(nb, n))
    cuts = [n * i // nb for i in range(nb + 1)]
    _STATE.update(func=func, Qs=Qs, Ts=Ts, p=p, cuts=cuts)
    if workers > 1 and nb > 1:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(min(workers, nb)) as pool:
            parts = pool.map(_chunk, range(nb))
    else:
        parts = [_chunk(i) for i in range(nb)]
    return np.concatenate(parts, axis=0)


def probe(prob, axis, rng, p, log, workers):
    """Q（axis 0）または T（axis 1）の指数の集合（他方は乱数の値に固定）。"""
    other = int(rng.integers(2, p - 1))
    Np = max(64, 1 << (6 * prob.N * prob.N - 1).bit_length())      # T(2,3) の幅 ≈ 3N² から見積もった初期値（足りなければ倍にする）
    while True:
        xs = exact._distinct(rng, Np + 3, p)
        others = np.full(Np + 3, other, dtype=np.int64)
        Qs, Ts = (xs, others) if axis == 0 else (others, xs)
        try:
            vals = evaluate(prob.raw, Qs, Ts, p, workers)
        except ZeroDivisionError:                    # 分母の零点：他方の値を取り直す
            other = int(rng.integers(2, p - 1))
            continue
        shift = Np // 2
        c = exact._laurent_1d(xs[:Np], vals[:Np], shift, p)
        ok = True
        for i in range(Np, Np + 3):
            x = int(xs[i])
            pred = np.zeros(vals.shape[1], dtype=np.int64)
            for kk in range(Np - 1, -1, -1):
                pred = (pred * x + c[kk]) % p
            pred = pred * pow(x, (-shift) % (p - 1), p) % p
            ok &= bool(np.all(pred == vals[i]))
        nz = np.nonzero(np.any(c != 0, axis=1))[0]
        if ok and len(nz) and nz[0] > 0 and nz[-1] < Np - 1:
            exps = [int(e) - shift for e in nz]
            log(f"  axis {'QT'[axis]}: exponents {exps[0]}..{exps[-1]}, parities {sorted({e % 2 for e in exps})} (probe {Np})")
            return exps, c
        Np *= 2
        if Np > 8192:
            raise RuntimeError("exponent range too large")


def torus_problem(m, n, N):
    d = (m - 1) * (n - 1) // 2
    lo1, hi1 = 2 * d - 1, 2 * d + 2 * (min(m, n) - 1) + 1     # \\hat f_[1] = (A − 1/A) rCS_[1] の A 次数（a の幅は min(m,n) − 1）
    return Problem(N, lambda: torus_rcs_factory(m, n), N * lo1, N * hi1)


def run(m, n, N, workers=8, log=print):
    return run_problem(torus_problem(m, n, N), 1000 * m + 100 * n + N, workers, log)


def run_problem(prob, seed, workers=8, log=print):
    t0 = time.time()
    rng = np.random.default_rng(seed)
    qexps, cq = probe(prob, 0, rng, P31, log, workers)
    texps, ct = probe(prob, 1, rng, P31, log, workers)
    prob.step = tuple(2 if len({e % 2 for e in ex}) == 1 else 1 for ex in (qexps, texps))
    prob.e0 = (qexps[0] % prob.step[0], texps[0] % prob.step[1])
    sq, st = prob.step
    xlo, xhi = (qexps[0] - prob.e0[0]) // sq, (qexps[-1] - prob.e0[0]) // sq
    ylo, yhi = (texps[0] - prob.e0[1]) // st, (texps[-1] - prob.e0[1]) // st
    Nx, Ny = xhi - xlo + 1, yhi - ylo + 1
    log(f"  grid {Nx} x {Ny} = {Nx * Ny} points, A window {prob.alo}..{prob.ahi} ({prob.Na}), {len(prob.rhos)} rho")
    primes, coeffs, M, lifted = (P31, P31B, P31D), [], 1, {}
    for i, p in enumerate(primes):
        tp = time.time()
        for attempt in range(5):                     # 乱数の点が分母の零点に当たったら点を取り直す
            Qs, Ts = exact._distinct(rng, Nx, p), exact._distinct(rng, Ny, p)
            Qg, Tg = np.meshgrid(Qs, Ts, indexing="ij")
            try:
                vals = evaluate(prob.ev_parity, Qg.ravel(), Tg.ravel(), p, workers).reshape(Nx, Ny, -1)
                break
            except ZeroDivisionError:
                log(f"  prime {i + 1}: 評価点が分母の零点に当たったので点を取り直す")
        log(f"  prime {i + 1}: evaluation {time.time() - tp:.0f}s")
        cj = exact._interp_grid(vals, Qs ** sq % p, Ts ** st % p, xlo, ylo, p, workers)
        del vals
        if i == 0:
            lifted, M = dict(cj), p
        else:
            inv = pow(M, p - 2, p)
            lifted = {k: x + M * (((cj.get(k, 0) - x) * inv) % p) for k in lifted.keys() | cj.keys() for x in (lifted.get(k, 0),)}
            M *= p
        poly = {k: (x - M if x > M // 2 else x) for k, x in lifted.items() if x}
        if verify(prob, poly, xlo, ylo):
            log(f"  {len(poly)} terms, {i + 1} prime(s), verified at independent points")
            break
        if i + 1 == len(primes):
            raise ArithmeticError("verification failed with all primes")
        log(f"  {i + 1} prime(s) not enough")
    # 列 d = (ρ の番号)·Na + k → A^{alo + 2k}、x^{ex} y^{ey} → Q^{sq·ex + e0} T^{st·ey + f0}
    hats = {rho: {} for rho in prob.rhos}
    for (d, ex, ey), c in poly.items():
        rho = prob.rhos[d // prob.Na]
        hats[rho][(prob.alo + 2 * (d % prob.Na), sq * ex + prob.e0[0], st * ey + prob.e0[1])] = c
    return hats, time.time() - t0


def verify(prob, poly, xlo, ylo, npts=6):
    p = P31C
    rng = np.random.default_rng(777)
    Qs, Ts = exact._distinct(rng, npts, p), exact._distinct(rng, npts, p)
    got = prob.ev_parity(Qs, Ts, p)
    pred = exact.evaluate_poly(poly, Qs ** prob.step[0] % p, Ts ** prob.step[1] % p, p)
    D = max(got.shape[1], pred.shape[1])
    got = np.pad(got, ((0, 0), (0, D - got.shape[1])))
    pred = np.pad(pred, ((0, 0), (0, D - pred.shape[1])))
    return bool(np.array_equal(got, pred))


def judge(hat, N):
    """(前の項数, 前の負, 後の項数, 後の負, 全体の符号, \\hat N)。符号は前の負が少ない方。G 展開が閉じなければ後は None。"""
    bold = to_bold(hat, a_shift(N))
    pos, neg = sum(v > 0 for v in bold.values()), sum(v < 0 for v in bold.values())
    sign = 1 if neg <= pos else -1
    try:
        Nh = g_expand({k: sign * v for k, v in bold.items()})
        aneg = sum(v < 0 for v in Nh.values())
    except ValueError:
        Nh, aneg = None, None
    return {"before_terms": len(bold), "before_neg": min(pos, neg), "sign": sign,
            "after_terms": None if Nh is None else len(Nh), "after_neg": aneg,
            "max_abs_N": None if Nh is None else max(abs(v) for v in Nh.values())}, Nh


def main():
    m, n, N = map(int, sys.argv[1:4])
    compute(m, n, N, int(sys.argv[4]) if len(sys.argv) > 4 else 8)


def compute(m, n, N, workers=8):
    """1 件（結び目と箱数）を計算して保存する。保存済みなら読まない。戻り値は計算時間（秒、保存済みなら記録の値）。"""
    name = f"T{m}_{n}"
    LOG.parent.mkdir(exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / f"{name}_N{N}.json.xz"
    with LOG.open("a") as lg:
        def log(s):
            print(s, flush=True)
            lg.write(s + "\n")
            lg.flush()
        if out.exists():
            log(f"{name} N={N}: 保存済み（{out}）")
            with lzma.open(out, "rt") as f:
                return json.load(f)["seconds"]
        log(f"== {name} N={N} start {time.strftime('%Y-%m-%d %H:%M:%S')} workers={workers}")
        hats, sec = run(m, n, N, workers, log)
        res = {}
        for rho in partitions(N):
            r, Nh = judge(hats[rho], N)
            res[rho] = (r, Nh)
            log(f"{name} ρ={list(rho)}：前 {r['before_neg']}/{r['before_terms']} 負、"
                f"後 {r['after_neg']}/{r['after_terms']} 負、符号 {r['sign']:+d}、max|N̂| {r['max_abs_N']}")
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30          # macOS はバイト
        rssc = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 2**30
        log(f"{name} N={N}: {sec:.0f}s、最大 RSS 親 {rss:.2f} GB・子 {rssc:.2f} GB")
        data = {"knot": name, "torus": [m, n], "N": N, "a_shift_A_units": a_shift(N),
                "variables_hat_f": "exponents of (A, Q, T) of 1703.05408, before the a-shift",
                "N_keys": "(g, beta, 2J_r) in bold variables, overall sign applied", "seconds": round(sec, 1),
                "rho": {"".join(map(str, rho)): {"judge": r,
                                                  "hat_f": {",".join(map(str, k)): v for k, v in sorted(hats[rho].items())},
                                                  "N": None if Nh is None else {",".join(map(str, k)): v for k, v in sorted(Nh.items())}}
                        for rho, (r, Nh) in res.items()}}
        with lzma.open(out, "wt") as f:
            json.dump(data, f)
        log(f"{name} N={N} done")
        return sec


if __name__ == "__main__":
    main()
