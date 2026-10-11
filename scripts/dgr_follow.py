"""DGR の superpolynomial の構成を例で追い、条件ごとに一意性を調べる（ユーザーの指示「DGR や quadruply の計算をいくつかの例でフォローし…」）。

正解：yui-kr の reduced HOMFLY-PT ホモロジー（Nakagane–Sano 2111.00388、data/reference/yui_kr）の Poincaré 多項式 P（太字）。
入力：P の 𝐭 = −1（HOMFLY-PT）、同じブレイド（KnotInfo）から計算した reduced Khovanov（khoca）と knot Floer（snappy）。
未知数：P の係数 c ≥ 0 と、各微分の対の始点 Q_d ≥ 0・残り r_d ≥ 0（scripts/construct_superpolynomial.py と同じ定式化）。
条件の段階：
  C0 HOMFLY-PT ＋ 非負性
  C1 ＋ DGR の対称性 (e_a, e_q, e_t) ↔ (e_a, −e_q, e_t − e_q)
  C2 ＋ d_{±1}（次数 (−2, ±2, −1/−3)、残りはちょうど 1 項）
  C3 ＋ 残りの位置 (Σ, −Σ, 0), (Σ, Σ, Σ)
  C4 ＋ d_2 の残り = Khovanov（𝐚 = 𝐪² での (2e_a + e_q, e_t)）
  C5 ＋ d_0 の残り = knot Floer（(A, M) = (e_q/2, e_t − e_a)）
各段階で：最小の総階数、正解の総階数、最小解が正解か、最小解の数（乱数の費用で 8 回解いて異なるものを数える）。
使い方：PYTHONPATH=src python scripts/dgr_follow.py 3_1 4_1 ...
"""
import csv
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import lil_matrix

from rcsknots.knot_table import homfly, load_homology, superpoly

ROOT = Path(__file__).resolve().parents[1]
DIFFS = {"d1": (-2, 2, -1), "d-1": (-2, -2, -3), "d2": (-2, 4, -1), "d0": (-2, 0, -3)}
STAGES = {"C0": dict(sym=False, diffs=(), sigma=False), "C1": dict(sym=True, diffs=(), sigma=False),
          "C2": dict(sym=True, diffs=("d1", "d-1"), sigma=False), "C3": dict(sym=True, diffs=("d1", "d-1"), sigma=True),
          "C4": dict(sym=True, diffs=("d1", "d-1", "d2"), sigma=True), "C5": dict(sym=True, diffs=("d1", "d-1", "d2", "d0"), sigma=True)}


def braid(name):
    with open(ROOT / "data/reference/yui_kr/braid.csv") as f:
        for row in csv.DictReader(f):
            if row["name"] == name:
                return json.loads(row["braid"])
    raise KeyError(name)


def khovanov(word):
    from khoca import InteractiveCalculator
    KH = InteractiveCalculator(101, (0, 0), 0)
    _, mess = KH(tuple(word), print_messages=True)
    line = mess[mess.index("Reduced Homology:") + 3]
    return {(int(q), -int(t)): (int(c) if c else 1) for c, t, q in re.findall(r"(\d*)t\^(-?\d+)q\^(-?\d+)", line)}


def knot_floer(word):
    import snappy
    return snappy.Link(braid_closure=word).knot_floer_homology()["ranks"]


def model(H, kh, hfk, grid_a, grid_q, grid_t, sym, diffs, sigma):
    grid = [(a, q, t) for a in grid_a for q in grid_q for t in grid_t]
    gset, idx, rows, lo, hi = set(grid), {}, [], [], []

    def var(k):
        if k not in idx:
            idx[k] = len(idx)
        return idx[k]

    def add(co, l, h):
        rows.append(co); lo.append(l); hi.append(h)

    for e in grid:
        var(("c", e))
    for a in grid_a:
        for q in grid_q:
            add({var(("c", (a, q, t))): (-1) ** (t % 2) for t in grid_t}, H.get((a, q), 0), H.get((a, q), 0))
    for d in diffs:
        da, dq, dt = DIFFS[d]
        for e in grid:
            co = {var(("c", e)): 1, var((d + "_r", e)): -1, var((d, e)): -1}
            src = (e[0] - da, e[1] - dq, e[2] - dt)
            if src in gset:
                co[var((d, src))] = co.get(var((d, src)), 0) - 1
            add(co, 0, 0)
            if (e[0] + da, e[1] + dq, e[2] + dt) not in gset:
                add({var((d, e)): 1}, 0, 0)
    if sym:
        for e in grid:
            m = (e[0], -e[1], e[2] - e[1])
            if m not in gset:
                add({var(("c", e)): 1}, 0, 0)
            elif m > e:
                add({var(("c", e)): 1, var(("c", m)): -1}, 0, 0)
    for d in ("d1", "d-1"):
        if d in diffs:
            add({var((d + "_r", e)): 1 for e in grid}, 1, 1)
            if sigma:
                for e in grid:
                    ok = (e[1] == -e[0] and e[2] == 0) if d == "d1" else (e[1] == e[0] and e[2] == e[0])
                    if not ok:
                        add({var((d + "_r", e)): 1}, 0, 0)
    for d, target, key in (("d2", kh, lambda e: (2 * e[0] + e[1], e[2])), ("d0", hfk, lambda e: (e[1] // 2, e[2] - e[0]) if e[1] % 2 == 0 else None)):
        if d not in diffs:
            continue
        spec = {}
        for e in grid:
            k = key(e)
            if k is None:
                add({var((d + "_r", e)): 1}, 0, 0)
            else:
                spec.setdefault(k, []).append(e)
        for k in set(spec) | set(target):
            if k not in spec:
                return None                       # 目標の生成元が格子の外
            add({var((d + "_r", e)): 1 for e in spec[k]}, target.get(k, 0), target.get(k, 0))
    M = lil_matrix((len(rows), len(idx)))
    for i, co in enumerate(rows):
        for j, v in co.items():
            M[i, j] = v
    return M.tocsr(), np.array(lo, float), np.array(hi, float), idx, grid


def solve(mod, cost, extra=None, time_limit=300):
    M, lo, hi, idx, grid = mod
    cons = [LinearConstraint(M, lo, hi)]
    if extra is not None:
        cons.append(extra)
    res = milp(cost, constraints=cons, integrality=np.ones(len(idx)), bounds=Bounds(0, np.inf), options={"time_limit": time_limit})
    if res.x is None:
        return None
    return {e: int(round(res.x[idx[("c", e)]])) for e in grid if round(res.x[idx[("c", e)]])}


def follow(name, margin_t=4, samples=8, log=print):
    P = superpoly(load_homology()[name])
    H = homfly(P)
    word = braid(name)
    kh = khovanov(word)
    hfk = knot_floer(word)
    A = sorted({a for a, _ in H}); Q = sorted({q for _, q in H})
    ga = list(range(A[0] - 2, A[-1] + 3, 2)); gq = list(range(Q[0] - 2, Q[-1] + 3, 2))
    ts = [t for (_, _, t) in P] + [t for _, t in kh]
    gt = list(range(min(ts) - margin_t, max(ts) + margin_t + 1))
    rank = sum(P.values())
    out = {"knot": name, "rank": rank, "thin": len({e[0] // 2 + e[1] // 2 - e[2] for e in P}) == 1, "stages": {}}
    for st, kw in STAGES.items():
        t0 = time.time()
        mod = model(H, kh, hfk, ga, gq, gt, **kw)
        if mod is None:
            out["stages"][st] = "目標が格子の外"
            continue
        n = len(mod[3])
        cost = np.zeros(n)
        for e in mod[4]:
            cost[mod[3][("c", e)]] = 1
        # 正解が実行可能か（規約の確認）
        fix = np.zeros((1, n))
        best = solve(mod, cost)
        if best is None:
            out["stages"][st] = "解なし"
            log(f"{name} {st}: 解なし")
            continue
        mn = sum(best.values())
        row = np.zeros(n)
        for e in mod[4]:
            row[mod[3][("c", e)]] = 1
        same_rank = LinearConstraint(row.reshape(1, -1), mn, mn)
        rng = np.random.default_rng(0)
        sols = [best]
        for _ in range(samples):
            c2 = np.zeros(n)
            for e in mod[4]:
                c2[mod[3][("c", e)]] = rng.integers(1, 100)
            s = solve(mod, c2, same_rank)
            if s is not None and s not in sols:
                sols.append(s)
        out["stages"][st] = {"vars": n, "rows": mod[0].shape[0], "min_rank": mn, "min_is_true": best == P,
                             "true_among_min": P in sols, "n_min_solutions_found": len(sols), "seconds": round(time.time() - t0, 1)}
        log(f"{name}（総階数 {rank}、thin {out['thin']}） {st}: 最小 {mn}、正解を含む {P in sols}、最小解 {len(sols)} 個以上"
            f"（変数 {n}、式 {mod[0].shape[0]}、{time.time() - t0:.0f} 秒）")
    return out


def main(names):
    outdir = ROOT / "data/dgr_lmov/follow"
    outdir.mkdir(parents=True, exist_ok=True)
    logf = open(ROOT / "logs/dgr_follow.log", "a")

    def log(s):
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()
    for n in names:
        res = follow(n, log=log)
        (outdir / f"{n}.json").write_text(json.dumps(res, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main(sys.argv[1:])
