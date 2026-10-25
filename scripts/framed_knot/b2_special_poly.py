"""open_problems B2：GNSSS（1512.07883）の四重次数の式の S² が、彩色 HOMFLY-PT の定理（special polynomial の性質 P_[2](q=1) = P_[1](q=1)²）を満たすかを、
pretzel K_{s,t}（s, t ∈ [−5, 5]）の全てで調べ、2 箱の判定（次元が二乗か、logs/framed_knot/t3_pretzel.log）と比べる。結果は logs/framed_knot/b2_special_poly.log。
"""
import re

import sympy as sp

from t3_pretzel import Q, a, q, quadruple, t, tc, tr


def euler_q1(e):
    return sp.cancel(e.subs({Q: q, tr: t / q, tc: q}, simultaneous=True)).subs({t: -1}).subs({q: 1})


squares = {}
for line in open("logs/framed_knot/t3_pretzel.log"):
    m = re.match(r"K_\{(-?\d+),(-?\d+)\}: rank (\d+), (\d+)", line)
    if m:
        squares[(int(m[1]), int(m[2]))] = int(m[4]) == int(m[3]) ** 2
table = {}
for st, sq in sorted(squares.items()):
    ok = sp.simplify(euler_q1(quadruple(*st, 2)) - euler_q1(quadruple(*st, 1)) ** 2) == 0
    table[(sq, ok)] = table.get((sq, ok), []) + [st]
    print(f"K_{{{st[0]},{st[1]}}}: 次元が二乗 {sq}、special polynomial の性質 {ok}", flush=True)
for (sq, ok), v in sorted(table.items()):
    print(f"次元が二乗 {sq}・special polynomial {ok}：{len(v)} 個", flush=True)
