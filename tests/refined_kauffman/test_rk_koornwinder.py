"""refined Kauffman K2'：Koornwinder 多項式（q-差分作用素の固有関数）の検査。

- n = 2 の主特殊化が van Diejen–Sahi の評価式と一致（パラメータ (1, −1, q^{1/2}, −q^{1/2}) = AS の D_n）。
- Mimachi の双対 Cauchy 恒等式（n = 2, m = 1 と n = 1, m = 2）。
"""
import pathlib
import sys

import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from rcsknots.refined_kauffman.koornwinder import koornwinder  # noqa: E402
from rcsknots.refined_kauffman.kauffman_refined import koornwinder_eval, q as qs, t as ts, U as Us  # noqa: E402
from rcsknots.partitions import partitions, conjugate  # noqa: E402

q = sp.Rational(9, 49); sq = sp.Rational(3, 7); t = sp.Rational(5, 11)
tbar = (1, -1, sq, -sq)


def test_evaluation_n2():
    xs = sp.symbols("x1:3")
    for lam in ((1,), (2,), (1, 1), (2, 1)):
        K = koornwinder(lam, 2, q, t, tbar, xs)
        val = K.subs({xs[0]: t, xs[1]: 1})
        pred = koornwinder_eval(lam, tbar).subs({qs: q, ts: t, Us: t**2})
        assert sp.nsimplify(val - pred) == 0


def _mimachi(n, m):
    xs = sp.symbols(f"x1:{n+1}"); ys = sp.symbols(f"y1:{m+1}")
    lhs = 0
    for k in range(0, n*m + 1):
        for lam in partitions(k):
            if len(lam) > n or (lam and lam[0] > m):
                continue
            full = tuple(lam) + (0,)*(n - len(lam))
            comp = tuple(m - full[n-1-i] for i in range(n))
            lamT = conjugate(lam) if lam else ()
            lhs += (-1)**k * koornwinder(comp, n, q, t, tbar, xs) * koornwinder(lamT, m, t, q, tbar, ys)
    rhs = sp.prod([(1 - x*yy)*(1 - x/yy)/x for x in xs for yy in ys])
    return sp.simplify(sp.expand(lhs) - sp.expand(rhs))


def test_mimachi():
    assert _mimachi(2, 1) == 0
    assert _mimachi(1, 2) == 0
