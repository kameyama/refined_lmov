"""[3,1] の t = q の Z 因子展開を「シフトは箱の組ごとの寄与の和」という仮定で総当たりする（方針 2 の段階 A の別法）。

各部分集合 S の各箱 b のシフト k_b(S) = Σ_{b'∈S∖b} e(b, b')、e ∈ {−1, 0, 1}（12 個）。連結な S の係数は 1、
連結でない S の係数は ε² の 2 次までの未知数。乱数の実数点 (A, q) で一次方程式を最小二乗で解き、残差 0 のものを解とする。
[2,1] では e(C,R) = −1, e(C,K) = +1, 他 0 で AMMM の展開になる（同じスクリプトで確かめる）。
"""
import itertools
import sys

import numpy as np

sys.path.insert(0, "scripts")
import f41_zexp as Z  # noqa: E402


def run(name, vals=(-1, 0, 1), npts=60, seed=0, all_coef=False):
    base = Z.BASES[name]
    boxes = list(base)
    H = Z.target(name)
    rng = np.random.default_rng(seed)
    A = rng.uniform(0.6, 1.6, npts)
    q = rng.uniform(0.6, 1.6, npts)
    Hv = sum(c * A ** a * q ** e for (a, e), c in H.items())

    def Zv(b, k):
        I, J = base[b][0] + k, base[b][1] - k
        return (A * q ** I - 1 / (A * q ** I)) * (A * q ** -J - 1 / (A * q ** -J))
    zt = {(b, k): Zv(b, k) for b in boxes for k in range(-4, 5)}
    e2 = (q - 1 / q) ** 2
    subsets = [S for r in range(1, len(boxes) + 1) for S in itertools.combinations(boxes, r)]
    disc = [S for S in subsets if len(S) > 1 and (all_coef or not Z.connected(name, S))]
    pairs = [(a, b) for a in boxes for b in boxes if a != b]
    sols = []
    for es in itertools.product(vals, repeat=len(pairs)):
        e = dict(zip(pairs, es))
        rhs = Hv - 1
        cols = []
        for S in subsets:
            pv = np.ones(npts)
            for b in S:
                k = sum(e[(b, b2)] for b2 in S if b2 != b)
                pv = pv * zt[(b, k)]
            if S in disc:
                cols += [pv, pv * e2, pv * e2 * e2]
            else:
                rhs = rhs - pv
        M = np.stack(cols, axis=1)
        c, *_ = np.linalg.lstsq(M, rhs, rcond=None)
        if np.max(np.abs(M @ c - rhs)) < 1e-7 * (1 + np.max(np.abs(rhs))):
            sols.append((e, np.round(c, 6)))
    return sols, disc


if __name__ == "__main__":
    name = sys.argv[1]
    allc = len(sys.argv) > 2 and sys.argv[2] == "all"
    sols, disc = run(name, npts=120 if allc else 60, all_coef=allc)
    print(f"{name}：解 {len(sols)} 個、連結でない部分集合 {disc}")
    for e, c in sols[:20]:
        print({k: v for k, v in e.items() if v}, c)
