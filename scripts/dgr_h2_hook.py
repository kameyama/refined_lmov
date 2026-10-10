"""H2（docs/dgr_lmov/open_problems.md2）：hook 性質を、正規化によらない形で GN の値から直接確かめる。

hook λ = (a+1, 1^l)、n = |λ| について、予想は
  raw_λ(q, t)|_{a=t/q} / raw_{[1]}(q', t')|_{a'=t'/q'} = ± q^α t^β（単項式）、　q' = q^{a+1}t^l、t' = q^a t^{l+1}
（raw は GN の reduced の値 = torus(λ)/torus(λ; 自明な結び目)、framing などの正規化の前）。
単項式であることは、4 組の (q, t)（有理数の平方）で、2 組から α, β を決め、残り 2 組で厳密に確かめる。
使い方：PYTHONPATH=src python scripts/dgr_h2_hook.py
"""
import math
import sys
from fractions import Fraction as Fr

from rcsknots.gn import GN

POINTS = [(Fr(3, 7) ** 2, Fr(5, 11) ** 2), (Fr(2, 9) ** 2, Fr(7, 13) ** 2), (Fr(5, 3) ** 2, Fr(4, 7) ** 2), (Fr(6, 5) ** 2, Fr(3, 2) ** 2)]


def raw(q, t, lam, slope):
    G = GN(q, t, t / q, lambda x: x)
    return G.torus(lam, *slope) / G.torus(lam, 1, 0)


def ratio(q, t, lam, slope):
    a, l = lam[0] - 1, len(lam) - 1
    q2, t2 = q ** (a + 1) * t ** l, q ** a * t ** (l + 1)
    return raw(q, t, lam, slope) / raw(q2, t2, (1,), slope)


def exact_power(x, al):
    """x^al（al は半整数、x は有理数の平方）。"""
    num, den = al.numerator, al.denominator
    if den == 1:
        return x ** num
    r = Fr(math.isqrt(x.numerator), math.isqrt(x.denominator))
    assert r * r == x
    return r ** num


def check(lam, slope):
    vals = [ratio(q, t, lam, slope) for q, t in POINTS]
    if any(v == 0 for v in vals):
        return "0 を含む"
    # α, β を最初の 2 点の log から半整数に丸めて決める
    (q1, t1), (q2, t2) = POINTS[:2]
    import numpy as np
    A = np.array([[math.log(q1), math.log(t1)], [math.log(q2), math.log(t2)]], dtype=float)
    b = np.array([math.log(abs(vals[0])), math.log(abs(vals[1]))])
    al, be = np.linalg.solve(A, b)
    al, be = Fr(round(2 * al), 2), Fr(round(2 * be), 2)
    sign = 1 if vals[0] > 0 else -1
    ok = all(v == sign * exact_power(q, al) * exact_power(t, be) for v, (q, t) in zip(vals, POINTS))
    return f"{'単項式' if ok else '単項式でない'}（± = {sign}, α = {al}, β = {be}）"


def main():
    hooks = [(2,), (1, 1), (3,), (2, 1), (1, 1, 1), (4,), (3, 1), (2, 1, 1), (1, 1, 1, 1), (5,), (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)]
    nmax = {(2, 3): 5, (2, 5): 4, (3, 2): 4, (3, 4): 3, (2, -3): 4}
    for slope, nm in nmax.items():
        for lam in hooks:
            if sum(lam) <= nm:
                print(f"T{slope} {lam}: {check(lam, slope)}", flush=True)


if __name__ == "__main__":
    main()
