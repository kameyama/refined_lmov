"""ステージ T4：refined 指数増大性 rCS_[2](q=1) = rCS_[1](q=1)² を満たす δ' の割り当てを整数計画で探す。

太字で q = 1 は 𝐭 = −𝐪⁻¹。GS の生成元 x = (i, j, k) と δ' から [2] の項は (−1)^k 𝐚^{2i} 𝐪^{k − 4i − 4δ'}（scripts/framed_knot/t4_colored.py の変換）。
使い方：from t4_egp import egp_solutions
"""
from collections import Counter

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp


def square_at_q1(P1):
    """rCS_[1](𝐭 = −𝐪⁻¹)² を {(e_𝐚, e_𝐪): c} で。"""
    s = Counter()
    for (a, q, t), c in P1.items():
        s[(a, q - t)] += c * (-1) ** (t % 2)
    out = Counter()
    for k1, c1 in s.items():
        for k2, c2 in s.items():
            out[(k1[0] + k2[0], k1[1] + k2[1])] += c1 * c2
    return {k: v for k, v in out.items() if v}


def egp_problem(gens, P1, choices):
    """変数 z_{x,d}（x：生成元、d ∈ choices、4d は整数）。各 x で Σ_d z = 重複度、EGP の等式。"""
    xs = sorted(gens)
    var = [(x, d) for x in xs for d in choices]
    rows, rhs = {}, square_at_q1(P1)
    for n, ((i, j, k), d) in enumerate(var):
        key = (2 * i, k - 4 * i - round(4 * d))
        rows.setdefault(key, {})[n] = (-1) ** (k % 2)
    keys = sorted(set(rows) | set(rhs))
    A = np.zeros((len(keys) + len(xs), len(var)))
    b = np.zeros(len(keys) + len(xs))
    for r, key in enumerate(keys):
        for n, v in rows.get(key, {}).items():
            A[r, n] = v
        b[r] = rhs.get(key, 0)
    for r, x in enumerate(xs):
        for n, (y, d) in enumerate(var):
            if y == x:
                A[len(keys) + r, n] = 1
        b[len(keys) + r] = gens[x]
    return var, A, b


def egp_solve(gens, P1, choices, cost=None):
    var, A, b = egp_problem(gens, P1, choices)
    c = np.zeros(len(var)) if cost is None else np.array([cost(x, d) for x, d in var])
    res = milp(c, constraints=LinearConstraint(A, b, b), integrality=np.ones(len(var)),
               bounds=Bounds(0, np.array([gens[x] for x, _ in var], dtype=float)))
    if res.x is None:
        return None
    return {(x, d): int(round(v)) for (x, d), v in zip(var, res.x) if round(v)}
