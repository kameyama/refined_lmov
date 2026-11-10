"""第 4 期 V5：thin な双曲結び目の 3 箱の弱い判定（正値性を満たす rCS_[2,1] が存在するか）。

入力：[1], [2], [3]（GNSSS 1512.07883 の四重次数の予想式、t_c 次数 \\tilde P(Q=q, t_r=1/q, t_c=tq)。第 2・3 期で検証済み）、
[1,1], [1,1,1]（ミラー (e_a, e_q, e_t) ↦ (e_a, −e_q, e_t − e_q)）、[2,1] は未知数（scripts/framed_knot/u3_21.py の整数計画）。
[2,1] の t = q の値は 3 本ケーブルから（colored_from_cables.colored_p）。カイラリティは [3] の t = q で決める。
使い方：python scripts/framed_knot/v5_weak.py 名前 ...（4_1, 5_2, 6_1, 7_2, 8_1）。結果は logs/framed_knot/v5_weak.log に追記。
判定：4_1 は KN の候補を代入した LP。他は LP 緩和（不能なら解なし確定）→ 制限時間付きの整数計画（見つからなければ未決）。
"""
import pickle
import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import u3_21 as U  # noqa: E402
from colored_from_cables import colored_p  # noqa: E402
from t3_pretzel import tc_graded  # noqa: E402
from t4_colored import bold_to_small  # noqa: E402

PRETZEL = {"4_1": (1, -1), "5_2": (1, 2), "6_1": (1, -2), "7_2": (1, 3), "8_1": (1, -3)}
AMPHICHIRAL = {"4_1"}


def mir(P):
    return {(a, -q, t - q): c for (a, q, t), c in P.items()}


def t_eq_q(small):
    """小文字 (A, Q, T) の多項式の T = Q での値 {(e_A, e_Q): c}"""
    out = {}
    for (a, q, t), c in small.items():
        out[(a, q + t)] = out.get((a, q + t), 0) + c
    return {k: v for k, v in out.items() if v}


def inputs(name):
    s, t = PRETZEL[name]
    P1, P2, P3 = (tc_graded(s, t, r) for r in (1, 2, 3))
    red = {(1,): P1, (2,): P2, (1, 1): mir(P2), (3,): P3, (1, 1, 1): mir(P3)}
    return {lam: bold_to_small(P) for lam, P in red.items()}, (sum(P1.values()), sum(P3.values()))


def dK_cached(name, red):
    """U.dK(red)。保存値（data/cache/v5_{name}.pkl）があれば読み、無ければ計算して保存する。"""
    cache = Path(f"data/cache/v5_{name}.pkl")
    if cache.exists():
        return pickle.load(open(cache, "rb"))
    dk = U.dK(red)
    cache.parent.mkdir(parents=True, exist_ok=True)
    pickle.dump(dk, open(cache, "wb"))
    return dk


def main(names):
    au = U.alpha_unk_cached()
    log = open("logs/framed_knot/v5_weak.log", "a")

    def out(s):
        print(s, flush=True)
        log.write(s + "\n")
        log.flush()

    for name in names:
        red, (d1, d3) = inputs(name)
        out(f"{name}: dim H^[1] = {d1}, dim H^[3] = {d3}, 三乗？ {d3 == d1 ** 3}")
        C = colored_p(name, 3)
        Hm = {lam: U.homfly_t_eq_q(C[lam]) for lam in C}
        h3 = t_eq_q(red[(3,)])
        flip = {(-a, -q): c for (a, q), c in Hm[(3,)].items()}
        if h3 == Hm[(3,)]:
            mirror_cable = False
        elif h3 == flip:
            mirror_cable = True
        else:
            out(f"{name}: [3] の t = q が 3 本ケーブルの値と一致しない（鏡像でも）。中止")
            continue
        out(f"{name}: [3] の t = q が 3 本ケーブルの値と一致（鏡像 {mirror_cable}）")
        H21 = Hm[(2, 1)] if not mirror_cable else {(-a, -q): c for (a, q), c in Hm[(2, 1)].items()}
        dk = dK_cached(name, red)
        amph = name in AMPHICHIRAL
        for sigma in (1, -1):
            out(f"{name} σ={sigma:+d}: {judge(red, H21, au, dk, amph, sigma, name)}")


def judge(red, H21, au, dk, amph, sigma, name):
    """正値性を課した系の判定。KN の候補があれば代入して LP、無ければ LP 緩和 → 制限時間付きの整数計画。"""
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, linprog, milp
    xs, A, b, lb, ub = U.build(red, H21, structure=(au, dk), positive=True, sigma=sigma, amphichiral=amph, log=lambda s: None)
    nx = len(xs)
    if name in KN_CANDIDATE:
        from rcsknots.laurent import load_table
        cand = load_table("data/framed_knot/reference/notebook_torus.txt")[KN_CANDIDATE[name]]
        if any(e not in set(xs) for e in cand):
            return "KN の候補が x の箱の外"
        lb, ub = lb.copy(), ub.copy()
        for n, e in enumerate(xs):
            lb[n] = ub[n] = cand.get(e, 0)
        res = linprog(np.zeros(A.shape[1]), A_eq=A, b_eq=b, bounds=list(zip(lb, ub)), method="highs")
        return "KN の候補で全条件（正値性を含む）を満たす" if res.status == 0 else f"KN の候補は条件を満たさない（{res.message[:60]}）"
    res = linprog(np.zeros(A.shape[1]), A_eq=A, b_eq=b, bounds=list(zip(lb, ub)), method="highs")
    if res.status == 2:
        return "LP 緩和が実行不能 → 正値性を満たす [2,1] は（x の箱の中に）存在しない"
    if res.status != 0:
        return f"LP 緩和が終わらない（{res.message[:60]}）"
    r = milp(np.zeros(A.shape[1]), constraints=LinearConstraint(A, b, b), integrality=np.ones(A.shape[1]),
             bounds=Bounds(lb, ub), options={"time_limit": MILP_SECONDS})
    if r.x is not None:
        sol = {e: int(round(v)) for e, v in zip(xs, r.x[:nx]) if round(v)}
        pickle.dump(sol, open(f"data/cache/v5_{name}_sol_sigma{sigma}.pkl", "wb"))
        return f"正値性を満たす整数解あり（[2,1] の項数 {len(sol)}）"
    return f"LP 緩和は実行可能、整数解は {MILP_SECONDS} 秒で見つからず（未決：{r.message[:50]}）"


KN_CANDIDATE = {"4_1": ("4_1", "21")}      # KN が構成した [2,1]（ノートブック）
MILP_SECONDS = 1800


if __name__ == "__main__":
    main(sys.argv[1:])
