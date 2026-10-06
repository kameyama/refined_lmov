"""第 4 期 V5：thin な双曲結び目の 3 箱の弱い判定（正値性を満たす rCS_[2,1] が存在するか）。

入力：[1], [2], [3]（GNSSS 1512.07883 の四重次数の予想式、t_c 次数 \\tilde P(Q=q, t_r=1/q, t_c=tq)。第 2・3 期で検証済み）、
[1,1], [1,1,1]（ミラー (e_a, e_q, e_t) ↦ (e_a, −e_q, e_t − e_q)）、[2,1] は未知数（scripts/u3_21.py の整数計画）。
[2,1] の t = q の値は 3 本ケーブルから（colored_from_cables.colored_p）。カイラリティは [3] の t = q で決める。
使い方：python scripts/v5_weak.py 名前 ...（4_1, 5_2, 6_1, 7_2, 8_1）。結果は logs/v5_weak.log に追記。
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


def main(names):
    au = pickle.load(open("data/cache/u3_alpha_unk.pkl", "rb"))
    log = open("logs/v5_weak.log", "a")

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
        cache = Path(f"data/cache/v5_{name}.pkl")
        if cache.exists():
            dk = pickle.load(open(cache, "rb"))
        else:
            dk = U.dK(red)
            pickle.dump(dk, open(cache, "wb"))
        amph = name in AMPHICHIRAL
        for positive in (False, True):
            res = []
            for sigma in (1, -1):
                sol = U.solve(red, H21, structure=(au, dk), positive=positive, sigma=sigma, amphichiral=amph, log=lambda s: None)
                res.append(sol is not None)
                if sol is not None and positive:
                    pickle.dump(sol, open(f"data/cache/v5_{name}_sol_sigma{sigma}.pkl", "wb"))
            out(f"{name}: 正値性{'あり' if positive else 'なし'}（両手型性 {amph}）：σ=+1 {'解あり' if res[0] else '解なし'}、σ=−1 {'解あり' if res[1] else '解なし'}")


if __name__ == "__main__":
    main(sys.argv[1:])
