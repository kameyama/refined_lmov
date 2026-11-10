"""第 3 期 U3：正解の分かる結び目で、条件から rCS_[2,1] が一意に決まるかを調べる（scripts/framed_knot/u3_21.py）。

有界の範囲（|x| ≤ BOUND）で乱数の目的関数の整数計画を何度か解き、得られた解と正解の差の階数を見る。
使い方：python scripts/framed_knot/u3_calibrate.py 名前 KN名 [amph] [u2only]   （例：3_1 T2_3、4_1 4_1 amph）
"""
import pickle
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

import u3_21 as U
from colored_from_cables import colored_p
from rcsknots.laurent import load_table

BOUND = 20
name, kn = sys.argv[1], sys.argv[2]
amph = "amph" in sys.argv[3:]
u2_cases = (True,) if "u2only" in sys.argv[3:] else (False, True)
R = load_table("data/framed_knot/reference/notebook_torus.txt")
C = {"1": (1,), "2": (2,), "11": (1, 1), "3": (3,), "111": (1, 1, 1)}
red = {lam: R[(kn, c)] for c, lam in C.items()}
truth = R.get((kn, "21"))
cache = Path(f"data/cache/u3_{kn}.pkl")
if cache.exists():
    S = pickle.load(open(cache, "rb"))
else:
    S = (U.alpha_unk_cached(), U.dK(red))
    pickle.dump(S, open(cache, "wb"))
H21 = U.homfly_t_eq_q(colored_p(name, 3)[(2, 1)])
for u2 in u2_cases:
    for positive in (False, True):
        xs, A, b, lb, ub = U.build(red, H21, structure=S, log=lambda s: None, positive=positive, amphichiral=amph, u2=u2)
        nx = len(xs)
        lb, ub = lb.copy(), ub.copy()
        lb[:nx], ub[:nx] = -BOUND, BOUND
        tv = np.array([truth.get(e, 0) for e in xs]) if truth else None
        sols = []
        for seed in range(6):
            c = np.zeros(A.shape[1])
            c[:nx] = np.random.default_rng(seed).uniform(-1, 1, nx)
            res = milp(c, constraints=LinearConstraint(A, b, b), integrality=np.ones(A.shape[1]),
                       bounds=Bounds(lb, ub), options={"time_limit": 300})
            if res.x is not None:
                sols.append(np.round(res.x[:nx]))
        if not sols:
            print(f"{name} U(2)={u2} 正値性={positive}: 解なし", flush=True)
            continue
        ref = tv if tv is not None else sols[0]
        D = np.array([s - ref for s in sols])
        same = sum(int(np.all(s == ref)) for s in sols)
        print(f"{name} U(2)={u2} 正値性={positive}: 解 {len(sols)} 個、正解（または最初の解）と一致 {same}、差の階数 {np.linalg.matrix_rank(D)}", flush=True)
