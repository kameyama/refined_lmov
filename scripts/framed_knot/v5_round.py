"""第 4 期 V5：5_2（σ = +1）の未決を詰める試み。LP 緩和の解を丸めた点 r の近くを、目的関数 Σ|x − r| の整数計画で探す（制限時間 30 分）。

結果は logs/framed_knot/v5_weak.log に追記。解が見つかれば data/cache/v5_5_2_sol_sigma1.pkl に保存する。
"""
import pickle

import numpy as np
import scipy.sparse as sps
from scipy.optimize import Bounds, LinearConstraint, linprog, milp

import u3_21 as U
from v5_weak import AMPHICHIRAL, dK_cached, inputs
from colored_from_cables import colored_p

NAME, SIGMA, SECONDS = "5_2", 1, 1800
red, _ = inputs(NAME)
H21 = U.homfly_t_eq_q(colored_p(NAME, 3)[(2, 1)])        # 5_2 は鏡像なし（v5_weak.log）
au = U.alpha_unk_cached()
dk = dK_cached(NAME, red)
xs, A, b, lb, ub = U.build(red, H21, structure=(au, dk), positive=True, sigma=SIGMA, amphichiral=NAME in AMPHICHIRAL, log=lambda s: None)
nx, nc = len(xs), A.shape[1]
lp = linprog(np.zeros(nc), A_eq=A, b_eq=b, bounds=list(zip(lb, ub)), method="highs")
r = np.round(lp.x[:nx])
# 変数 (x, h, N, d)：d_i ≥ x_i − r_i、d_i ≥ r_i − x_i、目的は Σ d_i
I = sps.identity(nx, format="csr")
Z = sps.csr_matrix((nx, nc - nx))
D1 = sps.hstack([I, Z, -I])           # x − d ≤ r
D2 = sps.hstack([-I, Z, -I])          # −x − d ≤ −r
Aeq = sps.hstack([A, sps.csr_matrix((A.shape[0], nx))])
cons = [LinearConstraint(Aeq, b, b), LinearConstraint(sps.vstack([D1, D2]), -np.inf, np.concatenate([r, -r]))]
c = np.concatenate([np.zeros(nc), np.ones(nx)])
res = milp(c, constraints=cons, integrality=np.concatenate([np.ones(nc), np.zeros(nx)]),
           bounds=Bounds(np.concatenate([lb, np.zeros(nx)]), np.concatenate([ub, np.full(nx, np.inf)])),
           options={"time_limit": SECONDS})
with open("logs/framed_knot/v5_weak.log", "a") as log:
    if res.x is not None:
        sol = {e: int(round(v)) for e, v in zip(xs, res.x[:nx]) if round(v)}
        pickle.dump(sol, open(f"data/cache/v5_{NAME}_sol_sigma{SIGMA}.pkl", "wb"))
        msg = f"{NAME} σ=+1（LP 解の丸めの近くを探索）：正値性を満たす整数解あり（[2,1] の項数 {len(sol)}、Σ|x − r| = {res.fun:.0f}）"
    else:
        msg = f"{NAME} σ=+1（LP 解の丸めの近くを探索）：{SECONDS} 秒で整数解は見つからず（未決：{res.message[:50]}）"
    print(msg, flush=True)
    log.write(msg + "\n")
