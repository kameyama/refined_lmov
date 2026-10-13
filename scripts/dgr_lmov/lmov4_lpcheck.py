"""保存した正値性の不等式（lmov4_pos_K_<p>.pkl、S なし）の多面体の有界性を、HiGHS の複数の解法で判定し直す（Solve error の切り分け）。
使い方：PYTHONPATH=src:scripts/framed_knot:scripts/dgr_lmov python scripts/dgr_lmov/lmov4_lpcheck.py p"""
import pickle
import sys
import time

import numpy as np
import scipy.sparse as sps
from scipy.optimize import linprog

p = int(sys.argv[1])
S = pickle.loads(open(f"data/dgr_lmov/lmovstruct/lmov4_pos_K_{p}.pkl", "rb").read())
ineqs = S["ineqs"][:S.get("nhat", 5461)]  # K_{−1} の古い保存は nhat なし（5461）
d = len(S["V"])
A = sps.lil_matrix((len(ineqs), d))
for i, (r, _) in enumerate(ineqs):
    for j, v in r.items():
        A[i, j] = v
A = A.tocsr(); b = np.array([float(rr) for _, rr in ineqs])
scale = np.array([max(1.0, max(abs(float(v)) for v in r.values())) for r, _ in ineqs])
An = sps.diags(1.0 / scale) @ A; bn = b / scale
print(f"K_{p}: 不等式 {len(ineqs)}、自由度 {d}、係数の最大 {abs(An).max():.3g}、最小の非零 {min(abs(v) for v in An.data):.3g}", flush=True)
rng = np.random.default_rng(0)
for k in range(3):
    c = rng.standard_normal(d)
    for method in ("highs-ds", "highs-ipm"):
        for sg in (1, -1):
            t0 = time.time()
            r = linprog(sg * c, A_ub=-An, b_ub=-bn, bounds=[(None, None)] * d, method=method,
                        options={"presolve": True, "time_limit": 600})
            print(f"  方向 {k} 符号 {sg:+d} {method}: status {r.status} ({r.message[:50]}) fun {None if r.fun is None else round(sg * r.fun, 2)} ({time.time() - t0:.0f}s)", flush=True)
# 有界な箱の中で解けるか（|w| ≤ 10⁶）：非有界なら箱の端に当たる
c = rng.standard_normal(d)
r = linprog(c, A_ub=-An, b_ub=-bn, bounds=[(-1e6, 1e6)] * d, method="highs")
print(f"  箱 ±1e6 の中：status {r.status} fun {r.fun} 端に当たる変数 {int(np.sum(np.abs(r.x) > 1e6 - 1)) if r.x is not None else None}", flush=True)
