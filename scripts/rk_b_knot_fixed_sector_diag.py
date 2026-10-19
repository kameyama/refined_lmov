r"""結び目の固定セクターの整数計画の切り分け：条件 (i)（t = −1）と (ii)（t = q）を片方ずつ外し、0 ≤ Q/P ≤ 1 の範囲の有無で解を探す（progress.md 続き 34）。使い方：python scripts/rk_b_knot_fixed_sector_diag.py m"""
import sys, io, contextlib; sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import numpy as np, sympy as sp
from collections import defaultdict
from scipy.optimize import milp, LinearConstraint, Bounds
import rk_k4_fixed_sector as H
import rk_b_knot_fixed_sector as K
import rk_b_cover_split as C
r, s, nu = H.r, H.s, H.nu
u = sp.Symbol("u", positive=True); t = s**2
m = int(sys.argv[1])
g = K.load_g(m)
with contextlib.redirect_stdout(io.StringIO()):
    sec = C.cover_sectors(m)
x = sp.cancel(sec["x"].to_sympy().subs(u, nu * s))
P = H.lp(sp.cancel((1 - t) * x)); R = H.lp(sp.cancel(2 * (1 + t) * g[(2,)]))
g1x2 = g[(1,)].subs({r: r**2, s: s**2, nu: nu**2}, simultaneous=True)
Eq = H.lp(sp.cancel((-(1 + r**2) * (g1x2 + 1)).subs(s, r)), (r, nu))
mons = sorted(P); idx = {mm: i for i, mm in enumerate(mons)}
print("P の係数の符号：正", sum(1 for v in P.values() if v > 0), "負", sum(1 for v in P.values() if v < 0))
# R の単項式が P の台に入っているか（t=−1 での比較の前に）
Ri = defaultdict(lambda: 0)
for (a, b, c), v in R.items(): Ri[(a, c)] += v * sp.I**b
Pi = defaultdict(lambda: 0)
for (a, b, c), v in P.items(): Pi[(a, c)] += v * sp.I**b
print("t=−1（s=i）での R と P の (r,ν) の台：R だけにある", len(set(k for k, v in Ri.items() if sp.expand(v) != 0) - set(Pi)), "個")
def build(use_i, use_ii):
    S = defaultdict(lambda: sp.Integer(0))
    for mm, c in R.items(): S[mm] += c / 2
    for mm, c in P.items(): S[mm] -= c / 2
    rows, rhs = [], []
    if use_i:
        Si = defaultdict(lambda: 0)
        for (a, b, c), v in S.items(): Si[(a, c)] += v * sp.I**b
        for key in set(Si) | {(a, c) for (a, b, c) in mons}:
            target = -sp.expand(Si.get(key, 0))
            for f in (sp.re, sp.im):
                row = np.zeros(len(mons))
                for mm in mons:
                    if (mm[0], mm[2]) == key: row[idx[mm]] = float(f(sp.I**mm[1]))
                tv = float(f(target))
                if np.any(row) or tv != 0: rows.append(row); rhs.append(tv)
    if use_ii:
        Pq = defaultdict(lambda: sp.Integer(0))
        for (a, b, c), v in P.items(): Pq[(a + b, c)] += v
        Qq = {k: (Pq.get(k, 0) - Eq.get(k, 0)) / 2 for k in set(Pq) | set(Eq)}
        print("   (ii) の Q|_{t=q} が整数:", all(v.is_integer for v in Qq.values()))
        for key, v in Qq.items():
            row = np.zeros(len(mons))
            for mm in mons:
                if (mm[0] + mm[1], mm[2]) == key: row[idx[mm]] = 1
            if np.any(row) or v != 0: rows.append(row); rhs.append(float(v))
    A = np.array(rows); b = np.array(rhs)
    lo = np.array([min(0, float(P[mm])) for mm in mons]); hi = np.array([max(0, float(P[mm])) for mm in mons])
    res = milp(np.zeros(len(mons)), constraints=LinearConstraint(A, b, b), bounds=Bounds(lo, hi), integrality=np.ones(len(mons)))
    # 範囲の制約なしでの実行可能性
    res2 = milp(np.zeros(len(mons)), constraints=LinearConstraint(A, b, b), bounds=Bounds(-50 * np.ones(len(mons)), 50 * np.ones(len(mons))), integrality=np.ones(len(mons)))
    return res.status == 0, res2.status == 0
for ui, uii in ((True, False), (False, True), (True, True)):
    print(f"(i)={ui}、(ii)={uii}：0≤Q/P≤1 で解あり {build(ui, uii)[0]}、範囲なし（整数）で解あり {build(ui, uii)[1]}", flush=True)
