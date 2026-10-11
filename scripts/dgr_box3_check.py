"""3 箱：答えが既知の結び目で、rCS_{[2,1]} の構成（scripts/dgr_box3.py の build21）が正解を一意に再現するかを確かめる。

入力は load_reduced の rCS（[2,1] 以外）と、正解の [2,1] の t = q（非 refined のデータの代わり）。全体の符号は引数で与える（これまで (1,1,1) だけが実行可能）。
使い方：PYTHONPATH=src python scripts/dgr_box3_check.py 結び目 1,1,1 [条件（既定 U,T,H,E,U2,L。B, R は仮説・予想の条件）] [方向の数（既定 1）]
"""
import sys, time; sys.path.insert(0, "scripts")
import numpy as np
from dgr_box3 import build21
import dgr_d7_construct as DC
from framing_check import load_reduced
from rcsknots.bps import to_bold
K = sys.argv[1]; sig = tuple(int(x) for x in sys.argv[2].split(","))
conds = tuple(sys.argv[3].split(",")) if len(sys.argv) > 3 else ("U", "T", "H", "E", "U2", "L")
R = load_reduced(K, 3)
true21 = {k: v for k, v in to_bold(R[(2, 1)], 1 if any(k[0] % 2 for k in R[(2, 1)]) else 0).items() if v}
H21u = {}
for (a, q, t), c in true21.items(): H21u[(a, q)] = H21u.get((a, q), 0) + c * (-1) ** (t % 2)
H21u = {k: v for k, v in H21u.items() if v}
small = {l: P for l, P in R.items() if l != (2, 1)}
ts = [k[2] for k in true21]; w = (min(ts) - 3, max(ts) + 4)
t0 = time.time(); lin, mons = build21(small, H21u, w, sig, conds)
print(K, sig, "".join(conds), "変数", len(lin.col), "式", len(lin.rows), f"{time.time()-t0:.0f}s", flush=True)
fix = DC.Lin(); fix.col, fix.rows, fix.rhs = dict(lin.col), list(lin.rows), list(lin.rhs)
for m in mons:
    fix.eq({("x", m): 1}, true21.get(m, 0))
missing = [m for m in true21 if m not in set(mons)]
rt = DC.solve(fix, xbound=400, integral=False)
print("  正解が全ての条件を満たす", rt.status == 0 and not missing, "（窓の外の項", len(missing), "）", flush=True)
r = DC.solve(lin, xbound=400); print("  実行可能", r.x is not None, f"{time.time()-t0:.0f}s", flush=True)
if r.x is None:
    print(K, sig, "実行不能", f"{time.time()-t0:.0f}s", flush=True); sys.exit()
sol = DC.solution(lin, r, mons)
print("  最初の解が正解", sol == true21, flush=True)


def other(cvec, bound):
    """x ≠ sol の整数解を探す：c·x ≥ c·sol + 1 または c·x ≤ c·sol − 1（c は整数ベクトル）を足した実行可能性問題。"""
    from scipy.optimize import Bounds, LinearConstraint, milp
    import scipy.sparse as sps
    n = len(lin.col)
    Aeq = sps.lil_matrix((len(lin.rows), n))
    for i, row in enumerate(lin.rows):
        for j, v in row.items():
            Aeq[i, j] = v
    lb = np.array([-400 if nm[0] == "x" else 0 for nm in lin.col], float)
    ub = np.array([400 if nm[0] == "x" else np.inf for nm in lin.col], float)
    cv = sum(cvec[lin.col[("x", m)]] * v for m, v in sol.items())
    rhs = np.array(lin.rhs, float)
    cons = [LinearConstraint(Aeq.tocsr(), rhs, rhs), LinearConstraint(cvec.reshape(1, -1), *bound(cv))]
    rr = milp(np.zeros(n), constraints=cons, bounds=Bounds(lb, ub), integrality=np.ones(n), options={"time_limit": 3000})
    return rr


# 一意性：整数の方向 c をとり、c·x > c·sol と c·x < c·sol のそれぞれで別の整数解があるかを調べる。
# どちらも実行不能なら、c·(x − sol) = 0 となる解を除いて一意（c は乱数なので、残る可能性は小さい）。
rng = np.random.default_rng(0)
sols = [sol]
for _ in range(int(sys.argv[4]) if len(sys.argv) > 4 and sys.argv[4].isdigit() else 1):
    c = np.zeros(len(lin.col))
    for m in mons: c[lin.col[("x", m)]] = rng.integers(1, 10 ** 6)
    for name, bound in (("上", lambda v: (v + 1, np.inf)), ("下", lambda v: (-np.inf, v - 1))):
        rr = other(c, bound)
        st = "実行不能" if rr.status == 2 else ("解あり" if rr.x is not None else f"打ち切り（status {rr.status}）")
        if rr.x is not None:
            x = DC.solution(lin, rr, mons)
            if x not in sols: sols.append(x)
        print("  別の解（" + name + "）", st, f"{time.time()-t0:.0f}s", flush=True)


def l1_min():
    """【仮説】rCS_{[2,1]} の係数の絶対値の和が最小の解（打ち消しの最小）。x = p − n、p, n ≥ 0 として Σ(p + n) を最小化する。"""
    from scipy.optimize import Bounds, LinearConstraint, milp
    import scipy.sparse as sps
    n = len(lin.col); xs = [lin.col[("x", m)] for m in mons]; k = len(xs)
    A = sps.lil_matrix((len(lin.rows) + k, n + 2 * k))
    for i, row in enumerate(lin.rows):
        for j, v in row.items():
            A[i, j] = v
    for i, j in enumerate(xs):
        A[len(lin.rows) + i, j] = 1; A[len(lin.rows) + i, n + i] = -1; A[len(lin.rows) + i, n + k + i] = 1
    rhs = np.array(list(lin.rhs) + [0] * k, float)
    lb = np.array([-400 if nm[0] == "x" else 0 for nm in lin.col] + [0] * (2 * k), float)
    ub = np.array([400 if nm[0] == "x" else np.inf for nm in lin.col] + [400] * (2 * k), float)
    obj = np.array([0.0] * n + [1.0] * (2 * k))
    rr = milp(obj, constraints=[LinearConstraint(A.tocsr(), rhs, rhs)], bounds=Bounds(lb, ub), integrality=np.ones(n + 2 * k),
              options={"time_limit": 3000})
    return {m: int(round(rr.x[lin.col[("x", m)]])) for m in mons if round(rr.x[lin.col[("x", m)]]) != 0} if rr.x is not None else None


if "--l1" in sys.argv:
    xl = l1_min()
    print("  L1 最小の解が正解", xl == true21, "Σ|係数|", sum(abs(v) for v in xl.values()) if xl else None,
          "（正解", sum(abs(v) for v in true21.values()), "）", f"{time.time()-t0:.0f}s", flush=True)
import json
from pathlib import Path
out = Path(f"data/dgr_lmov/box3/check_{K}_{''.join(conds)}.json")
out.write_text(json.dumps({"sigma": sig, "true": [[list(m), v] for m, v in true21.items()],
                           "solutions": [[[list(m), v] for m, v in x.items()] for x in sols]}))
print(K, sig, "".join(conds), "異なる解（下界）", len(sols), "正解を含む", true21 in sols, f"{time.time()-t0:.0f}s", flush=True)
