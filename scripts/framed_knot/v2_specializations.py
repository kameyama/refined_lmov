"""第 4 期 V2-2：rCS_[2,1] の特殊化 A = s·Q^i·T^j（s = ±1）で、低い彩色の積の単項式倍になるものを、正解のデータで探す。

正解：KN のノートブック（T(2,3), T(2,5), T(2,7), T(2,9), 4_1）と第 3 期の楕円 Hall 代数（T(3,4), T(3,5), T(2,11), T(3,7), T(4,5)）。
ある特殊化と積の組が、全ての結び目で「比が単項式」になれば恒等式の候補として出力する。結果は logs/framed_knot/v2_specializations.log。
"""
import json
import sys
from itertools import product
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rcsknots.laurent import load_table  # noqa: E402

A, Q, T = sp.symbols("A Q T")
NB = ["T2_3", "T2_5", "T2_7", "T2_9", "4_1"]
EHA = ["T3_4", "T3_5", "T2_11", "T3_7", "T4_5"]
COLORS = ["1", "2", "11", "3", "111", "21"]


def load():
    R = load_table("data/framed_knot/reference/notebook_torus.txt")
    out = {k: {c: R[(k, c)] for c in COLORS} for k in NB}
    for k in EHA:
        out[k] = {c: {tuple(map(int, e.split(","))): v for e, v in json.loads(Path(f"data/framed_knot/torus/{k}_{c}.json").read_text())["poly"].items()}
                  for c in COLORS}
    return out


def ex(P):
    return sum(c * A ** i * Q ** j * T ** k for (i, j, k), c in P.items())


PRODUCTS = {"r1·r11": lambda r: r["1"] * r["11"], "r1·r2": lambda r: r["1"] * r["2"], "r3": lambda r: r["3"],
            "r111": lambda r: r["111"], "r1^3": lambda r: r["1"] ** 3, "r1": lambda r: r["1"], "r2": lambda r: r["2"],
            "r11": lambda r: r["11"], "1": lambda r: sp.Integer(1)}


def is_monomial(e):
    e = sp.factor(sp.cancel(e))
    if e == 0:
        return False
    num, den = sp.fraction(e)
    return len(sp.Add.make_args(sp.expand(num))) == 1 and len(sp.Add.make_args(sp.expand(den))) == 1


def main():
    data = {k: {c: ex(p) for c, p in v.items()} for k, v in load().items()}
    knots = list(data)
    found = []
    for s, i, j in product((1, -1), range(-3, 4), range(-3, 4)):
        if (i, j) == (0, 0):
            continue
        sub = {A: s * Q ** i * T ** j}
        for name, f in PRODUCTS.items():
            ok = True
            for k in knots:
                r21 = sp.expand(data[k]["21"].subs(sub))
                rhs = sp.expand(f(data[k]).subs(sub))
                if r21 == 0 and rhs == 0:
                    continue
                if r21 == 0 or rhs == 0 or not is_monomial(r21 / rhs):
                    ok = False
                    break
            if ok:
                line = f"A = {'+' if s > 0 else '-'}Q^{i} T^{j}: rCS_[2,1] / {name} は全 {len(knots)} 結び目で単項式"
                print(line, flush=True)
                found.append(line)
    print(f"候補 {len(found)} 個", flush=True)


if __name__ == "__main__":
    main()
