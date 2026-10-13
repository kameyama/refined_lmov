"""4_1 の hook 型の彩色の Z 因子展開（AMMM 1211.6375 の形）を t = q で求める（方針 2 の段階 A）。

AMMM の [2,1]（:396-403）：H/S = 1 + (Z_{3|3} + Z_{2|0} + Z_{0|2}) + (Z_{4|2}Z_{2|0} + Z_{2|4}Z_{0|2} + (1−ε²)Z_{2|0}Z_{0|2}) + Z_{3|3}Z_{2|0}Z_{0|2}、
Z_{I|J} = {Aq^I}{Aq^{−J}} = D_I D_{−J}、ε = q − 1/q。箱の部分集合ごとに 1 項で、各箱の Z 因子は Z_{I+k|J−k}（k は部分集合ごとのシフト）、
係数 c_S は ε² の多項式。ここでは部分集合ごとに (シフト, 係数) を 1 つ選ぶ 0-1 整数計画で H を再現する分解を探す。
使い方：python scripts/framed_knot/f41_zexp.py 21|31
"""
import itertools
import sys
from collections import defaultdict

import numpy as np
import scipy.sparse as sps
from scipy.optimize import Bounds, LinearConstraint, milp


def mul(*ps):
    out = {(0, 0): 1}
    for p in ps:
        nx = defaultdict(int)
        for k1, c1 in out.items():
            for k2, c2 in p.items():
                nx[(k1[0] + k2[0], k1[1] + k2[1])] += c1 * c2
        out = {k: v for k, v in nx.items() if v}
    return out


def D(k):
    """{A q^k}（A, q の指数の辞書）"""
    return {(1, k): 1, (-1, -k): -1}


def Z(I, J):
    return mul(D(I), D(-J))


EPS2 = {(0, 2): 1, (0, 0): -2, (0, -2): 1}             # ε² = q² − 2 + q^{−2}
COEFS = {"1": {(0, 0): 1}, "1-e2": None, "1-2e2": None, "1+e2": None, "1-e2+e4": None, "1-3e2": None}


def coef(name):
    e4 = mul(EPS2, EPS2)
    table = {"1": {(0, 0): 1},
             "1-e2": {k: v for k, v in {**{(0, 0): 1}}.items()},
             }
    out = defaultdict(int)
    out[(0, 0)] += 1
    n2 = {"1": 0, "1-e2": -1, "1-2e2": -2, "1-3e2": -3, "1+e2": 1, "1-e2+e4": -1}[name]
    n4 = 1 if name == "1-e2+e4" else 0
    for k, v in EPS2.items():
        out[k] += n2 * v
    for k, v in e4.items():
        out[k] += n4 * v
    del table
    return {k: v for k, v in out.items() if v}


BASES = {"21": {"C": (3, 3), "R": (0, 2), "K": (2, 0)},
         "31": {"C": (5, 3), "R1": (2, 2), "R2": (0, 2), "K": (3, -1)}}


def target(name):
    import sympy as sp
    from rcsknots.double_twist import H21_1612, A as LA, q
    if name == "21":
        H = H21_1612(-1, 1)
    else:
        sys.path.insert(0, "scripts/framed_knot")
        from f41_box4 import homfly31
        h = homfly31()                                    # 既に (e_A, e_q) の辞書
        return {k: v for k, v in h.items()}
    e = sp.expand(sp.cancel(H) * LA ** 40 * q ** 80)
    return {(i - 40, j - 80): int(c) for (i, j), c in sp.Poly(e, LA, q).terms()}


ADJ = {"21": {("C", "R"), ("C", "K")}, "31": {("C", "R1"), ("R1", "R2"), ("C", "K")}}


def connected(name, S):
    """箱の部分集合が辺でつながっているか（AMMM：係数の補正は連結でない部分集合だけ）。"""
    S = set(S)
    seen, todo = set(), [next(iter(S))]
    while todo:
        b = todo.pop()
        seen.add(b)
        todo += [y for x, y in ADJ[name] if x == b and y in S and y not in seen]
        todo += [x for x, y in ADJ[name] if y == b and x in S and x not in seen]
    return seen == S


def build(name, kmax=3, coefs=("1", "1-e2", "1-2e2", "1+e2", "1-e2+e4"), conn_rule=False):
    base = BASES[name]
    boxes = list(base)
    H = target(name)
    rest = dict(H)
    rest[(0, 0)] = rest.get((0, 0), 0) - 1
    cols, polys = [], []
    groups = []
    for r in range(1, len(boxes) + 1):
        for S in itertools.combinations(boxes, r):
            g = []
            shifts = [(0,) * r] if r == 1 else itertools.product(range(-kmax, kmax + 1), repeat=r)
            cs = ("1",) if r == 1 or (conn_rule and connected(name, S)) else coefs
            for ks in shifts:
                prod = mul(*[Z(base[b][0] + k, base[b][1] - k) for b, k in zip(S, ks)])
                for c in cs:
                    g.append(len(cols))
                    cols.append((S, ks, c))
                    polys.append(mul(prod, coef(c)))
            groups.append(g)
    keys = sorted(set(rest) | {k for p in polys for k in p})
    kidx = {k: i for i, k in enumerate(keys)}
    ri, ci, data = [], [], []
    for j, p in enumerate(polys):
        for k, v in p.items():
            ri.append(kidx[k])
            ci.append(j)
            data.append(v)
    Aeq = sps.csr_matrix((data, (ri, ci)), shape=(len(keys), len(cols)), dtype=float)
    beq = np.array([rest.get(k, 0) for k in keys], dtype=float)
    G = sps.csr_matrix(([1.0] * sum(len(g) for g in groups),
                        ([i for i, g in enumerate(groups) for _ in g], [j for g in groups for j in g])),
                       shape=(len(groups), len(cols)))
    return cols, Aeq, beq, G


def solve(name, exclude=(), log=print, kw_time=600, **kw):
    cols, Aeq, beq, G = build(name, **kw)
    cons = [LinearConstraint(Aeq, beq, beq), LinearConstraint(G, 1, 1)]
    for sol in exclude:                                   # 既知の解を除く（no-good 制約）
        row = np.zeros(len(cols))
        row[[cols.index(c) for c in sol]] = 1
        cons.append(LinearConstraint(row, -np.inf, len(sol) - 1))
    # 目的：シフトの絶対値と係数の補正を小さく
    cost = np.array([sum(abs(k) for k in ks) + (0 if c == "1" else 1) for _, ks, c in cols], dtype=float)
    res = milp(cost, constraints=cons, integrality=np.ones(len(cols)), bounds=Bounds(0, 1),
               options={"time_limit": kw_time})
    log(f"  {name}：変数 {len(cols)}、{res.message}")
    if res.x is None:
        return None
    return [cols[j] for j in np.nonzero(res.x > 0.5)[0]]


if __name__ == "__main__":
    name = sys.argv[1]
    nsol = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    kmax = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    tl = float(sys.argv[4]) if len(sys.argv) > 4 else 600
    found = []
    for _ in range(nsol):
        sol = solve(name, exclude=found, kmax=kmax, conn_rule=True, kw_time=tl,
                    coefs=("1", "1-e2", "1-2e2", "1+e2", "1-e2+e4", "1-3e2"))
        if sol is None:
            break
        found.append(sol)
        print("解", len(found))
        for S, ks, c in sol:
            print("  ", S, ks, c)
