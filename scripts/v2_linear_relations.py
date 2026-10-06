"""第 4 期 V2-2：特殊化 A = s·Q^i·T^j で、rCS_[2,1] = α·X + β·Y（α, β は結び目によらない (Q, T) の有理関数）となる関係を探す。

α, β を T(2,3), T(2,5) の 2 結び目で解き、T(2,7), T(3,4), 4_1 で確かめる。X, Y は低い彩色の積から選ぶ。結果は logs/v2_linear_relations.log。
"""
import sys
from itertools import combinations, product
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from v2_specializations import A, PRODUCTS, Q, T, ex, load  # noqa: E402

SOLVE, CHECK = ["T2_3", "T2_5"], ["T2_7", "T3_4", "4_1"]


def main():
    raw = load()
    data = {k: {c: ex(p) for c, p in raw[k].items()} for k in SOLVE + CHECK}
    names = [n for n in PRODUCTS if n != "1"] + ["1"]
    a_, b_ = sp.symbols("alpha beta")
    for s, i, j in product((1, -1), range(-3, 4), range(-3, 4)):
        if (i, j) == (0, 0):
            continue
        sub = {A: s * Q ** i * T ** j}
        vals = {k: {n: sp.expand(f(data[k]).subs(sub)) for n, f in PRODUCTS.items()} | {"21": sp.expand(data[k]["21"].subs(sub))}
                for k in data}
        for X, Y in combinations(names, 2):
            eqs = [sp.Eq(vals[k]["21"], a_ * vals[k][X] + b_ * vals[k][Y]) for k in SOLVE]
            try:
                sol = sp.solve(eqs, [a_, b_], dict=True)
            except Exception:
                continue
            if len(sol) != 1:
                continue
            al, be = sp.cancel(sol[0][a_]), sp.cancel(sol[0][b_])
            if all(sp.cancel(vals[k]["21"] - al * vals[k][X] - be * vals[k][Y]) == 0 for k in CHECK):
                print(f"A = {'+' if s > 0 else '-'}Q^{i} T^{j}: r21 = ({sp.factor(al)})·{X} + ({sp.factor(be)})·{Y}", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
