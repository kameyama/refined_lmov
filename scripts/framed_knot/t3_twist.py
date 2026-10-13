"""ステージ T3：ツイスト結び目 K_p の 2 箱の refined LMOV 不変量。

入力：Nawata–Ramadevi–Zodinmawia–Sun（arXiv:1209.1409）の対称彩色 superpolynomial の二重和の式 (superpoly3), (superpoly4)
（TeX references/1209.1409/src/superpolynomials-final.tex の 921-946 行）。予想に基づく式である。
P_n は S^{n−1} 彩色、変数は GS の old 規約（a, q は DGR の半分、t は t_r）。K_1 = 3_1, K_{−1} = 4_1, K_2 = 5_2, K_{−2} = 6_1, …
t_c 次数の [2], [1,1] は scripts/framed_knot/t4_colored.py の変換（thin なので δ' は一定、EGP を満たす値に決める）。
使い方：python scripts/framed_knot/t3_twist.py P_MIN P_MAX
"""
import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from t4_colored import judge2, tc_colored  # noqa: E402
from t4_egp import egp_solve  # noqa: E402

a, q, t = sp.symbols("a q t")


def qp(x, k):
    out = 1
    for i in range(k):
        out *= 1 - x * q ** i
    return out


def qbinom(n, k):
    return sp.cancel(qp(q, n) / (qp(q, k) * qp(q, n - k)))


def superpoly(p, n):
    """NRZS の P_n(K_p; a, q, t)（p ≠ 0）。"""
    total = 0
    for k in range(n):                       # (q^{1−n}; q)_k = 0 for k ≥ n
        head = q ** k * qp(-a * t / q, k) / qp(q, k) * qp(q ** (1 - n), k) * qp(-a * t ** 3 * q ** (n - 1), k)
        for l in range(k + 1):
            total += head * (-1) ** l * (a * t ** 2) ** (p * l) * q ** sp.Rational((2 * p + 1) * l * (l - 1), 2) \
                * (1 - a * t ** 2 * q ** (2 * l - 1)) / qp(a * t ** 2 * q ** (l - 1), k + 1) * qbinom(k, l)
    if p > 0:
        total *= (-t) ** (1 - n)
    return sp.factor(sp.cancel(sp.together(total)))


def terms(expr):
    e = sp.expand(expr * (a * q * t) ** 200)
    return {(i - 200, j - 200, k - 200): int(c) for (i, j, k), c in sp.Poly(e, a, q, t).terms()}


def inputs(p):
    """[1]（太字）と GS 規約の S² の生成元。"""
    P1 = {(2 * i, 2 * j, k): c for (i, j, k), c in terms(superpoly(p, 2)).items()}
    return P1, terms(superpoly(p, 3))


def run(p):
    P1, gens = inputs(p)
    if min(P1.values()) < 0 or min(gens.values()) < 0:
        return f"K_{p}: 入力に負の係数"
    sol = None
    for d in [x / 2 for x in range(-24, 25)]:
        if egp_solve(gens, P1, [d]) is not None:
            sol = d
            break
    if sol is None:
        return f"K_{p}: 一定の δ' で EGP を満たすものが無い"
    P2, P11 = tc_colored(gens, lambda x: sol)
    return f"K_{p}: rank {sum(P1.values())}, {sum(gens.values())}, δ' = {sol}: {judge2(P1, P2, P11)}"


if __name__ == "__main__":
    lo, hi = int(sys.argv[1]), int(sys.argv[2])
    for p in range(lo, hi + 1):
        if p:
            print(run(p), flush=True)
