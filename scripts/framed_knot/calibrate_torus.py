"""Stage 2：GN の構成によるトーラス結び目の reduced 値と、1703.05408 ノートブックの値の比（単項式）を求める。"""
import sys
from fractions import Fraction

import numpy as np

from rcsknots.field import V, P31, make_conv
from rcsknots.gn import GN
from rcsknots.laurent import load_table, evaluate
from rcsknots.macdonald import principal_specialization

COLORS = {"1": (1,), "2": (2,), "11": (1, 1), "3": (3,), "111": (1, 1, 1), "21": (2, 1)}


def setup(seed, npts=3):
    rng = np.random.default_rng(seed)
    x, y, b = (V(rng.integers(2, P31, npts), P31) for _ in range(3))
    conv = make_conv(P31, npts)
    return x, y, b, GN(x * x, y * y, b * b, conv)


def monomial_table(x, y, b, R=130):
    tab = {}
    for i in range(-30, 31):
        for j in range(-R, R + 1):
            for k in range(-R, R + 1):
                tab[int((b ** i * x ** j * y ** k).a[0])] = (i, j, k)
    return tab


def main(knots):
    ref = load_table("data/framed_knot/reference/notebook_torus.txt")
    x, y, b, G = setup(11)
    tab = monomial_table(x, y, b)
    for knot, (n, m) in knots:
        for cname, lam in COLORS.items():
            if (knot, cname) not in ref:
                continue
            red = G.torus(lam, n, m) / principal_specialization(lam, G.a, G.q, G.t)
            nb = evaluate(ref[(knot, cname)], b, x, y)
            found = None
            for sign in (1, -1):
                r = sign * nb / red
                e = tab.get(int(r.a[0]))
                if e and r == b ** e[0] * x ** e[1] * y ** e[2]:
                    found = (sign, e)
            print(knot, (n, m), cname, "ratio = sign, (2·exp_a, 2·exp_q, 2·exp_t):", found, flush=True)


if __name__ == "__main__":
    which = sys.argv[1:] or ["T2_3"]
    spec = {"T2_3": (2, -3), "T2_5": (2, -5), "T2_7": (2, -7), "T2_9": (2, -9), "T3_4": (3, -4)}
    main([(k, spec[k]) for k in which])
