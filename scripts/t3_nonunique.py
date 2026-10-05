"""ステージ T3：2 箱の正値性は自動的には従わないことの確認（5_2）。

正値性を課さず、Laurent 性と G 展開が閉じることだけを満たす δ' の割り当てを、乱数の目的関数で 12 個探し、それぞれの正値性を判定する。
結果は logs/t3_nonunique.log。
"""
import random
from collections import Counter

from t3_twist import inputs
from t4_colored import judge2, tc_colored
from t4_ilp import solve, structure

S = structure()
P1, G = inputs(2)
for seed in range(12):
    rnd, w = random.Random(seed), {}
    sol = solve(P1, G, [-2, -1.5, -1, -0.5, 0], -1, S, log=lambda s: None, positive=False,
                cost=lambda x, d: w.setdefault((x, d), rnd.random()))
    two, one_one = {}, {}
    for (x, d), m in sol.items():
        a2, b2 = tc_colored({x: m}, lambda y, d=d: d)
        for D, E in ((two, a2), (one_one, b2)):
            for k, v in E.items():
                D[k] = D.get(k, 0) + v
    counts = dict(Counter(d for (x, d), m in sol.items() for _ in range(m)))
    print(f"seed {seed}: δ' の個数 {counts}: {judge2(P1, two, one_one)}", flush=True)
