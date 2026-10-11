"""3 箱：答えが既知の結び目で、rCS_{[2,1]} の構成（scripts/dgr_box3.py の build21）が正解を一意に再現するかを確かめる。

入力は load_reduced の rCS（[2,1] 以外）と、正解の [2,1] の t = q（非 refined のデータの代わり）。全体の符号は引数で与える（これまで (1,1,1) だけが実行可能）。
使い方：PYTHONPATH=src python scripts/dgr_box3_check.py 結び目 1,1,1 [乱数の試行回数]
"""
import sys, time; sys.path.insert(0, "scripts")
import numpy as np
from dgr_box3 import build21
import dgr_d7_construct as DC
from framing_check import load_reduced
from rcsknots.bps import to_bold
K = sys.argv[1]; sig = tuple(int(x) for x in sys.argv[2].split(","))
R = load_reduced(K, 3)
true21 = {k: v for k, v in to_bold(R[(2, 1)], 1 if any(k[0] % 2 for k in R[(2, 1)]) else 0).items() if v}
H21u = {}
for (a, q, t), c in true21.items(): H21u[(a, q)] = H21u.get((a, q), 0) + c * (-1) ** (t % 2)
H21u = {k: v for k, v in H21u.items() if v}
small = {l: P for l, P in R.items() if l != (2, 1)}
ts = [k[2] for k in true21]; w = (min(ts) - 3, max(ts) + 4)
t0 = time.time(); lin, mons = build21(small, H21u, w, sig)
print(K, sig, "変数", len(lin.col), "式", len(lin.rows), f"{time.time()-t0:.0f}s", flush=True)
r = DC.solve(lin, xbound=400); print("  実行可能", r.x is not None, f"{time.time()-t0:.0f}s", flush=True)
sols = []
if r.x is not None:
    sols.append(DC.solution(lin, r, mons))
    rng = np.random.default_rng(0)
    for _ in range(int(sys.argv[3]) if len(sys.argv) > 3 else 5):
        obj = np.zeros(len(lin.col))
        for m in mons: obj[lin.col[("x", m)]] = rng.normal()
        rr = DC.solve(lin, objective=obj, xbound=400)
        if rr.x is not None:
            s = DC.solution(lin, rr, mons)
            if s not in sols: sols.append(s)
        print("  試行", len(sols), f"{time.time()-t0:.0f}s", flush=True)
print(K, sig, "異なる解", len(sols), "正解を含む", true21 in sols, f"{time.time()-t0:.0f}s", flush=True)
