"""四重次数（GGS）の S² 側の構成を例で追い、条件ごとに一意性を調べる（scripts/dgr_lmov/dgr_follow.py の S² 版）。

正解：S² の t_c 次数の Poincaré 多項式 P₂（太字、≥ 0）。ツイスト結び目は GNSSS の予想式、10_124/10_139/10_152 は quiver と KRSS の予想式。
未知数：P₂ の係数 c ≥ 0、打ち消す微分の対の始点 Q ≥ 0 と残り r ≥ 0。
条件の段階（scripts/dgr_lmov/dgr_s2_discover.py で正解から見つけた構造。全ての例で成立）：
  S0  𝐭 = −1 で色付き HOMFLY-PT [2] ＋ 非負性
  S1  ＋ 打ち消す微分 d_{1|0}, d_{0|2}（次数 −u_{n|m}、残りはちょうど 1 項で、その特殊化での値は正解と同じ単項式。位置は直線上で自由）
  S2  ＋ 色付き微分の特殊化：P₂|_{d_{1|1}} = c·P₁|_{d_{1|1}}(𝐪²𝐭, 𝐭)、P₂|_{d_{0|3}} = c'·P₁|_{d_{0|3}}（係数 c, c' は正解から）
  S3  ＋ refined の指数増大性 P₂(𝐪 = 𝐭⁻¹) = P₁(𝐪 = 𝐭⁻¹)²
  S4  ＋ sl(2) の像（Khovanov の類似）：P₂|_{d_{2|0}} を正解の値に固定（特殊化のみ。ホモロジーとしての残りの情報は使わない）
各段階で：最小の総階数、正解を含むか、最小解の数（乱数の費用 8 回）。
使い方：PYTHONPATH=src python scripts/dgr_lmov/dgr_s2_follow.py [--spec]（--spec は打ち消す微分も特殊化だけで課す）
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import lil_matrix

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "framed_knot"))
from v4_closed import inputs  # noqa: E402
from rcsknots.dgr_lmov.dgr import specialize, u  # noqa: E402
from rcsknots.knot_table import load_homology, superpoly  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
STAGES = ["S0", "S1", "S2", "S3", "S4"]
HOMOLOGY = True                     # False：打ち消す微分もホモロジーではなく特殊化（単項式）としてだけ課す（離散トモグラフィーの比較）


def ratio(A, B):
    ka, kb = min(A), min(B)
    return (ka[0] - kb[0], ka[1] - kb[1]), A[ka] // B[kb]


def build(P1, P2, stage, margin_t=3):
    P2 = {k: int(v) for k, v in P2.items() if v}
    P1 = {k: int(v) for k, v in P1.items() if v}
    H = {}
    for (a, q, t), c in P2.items():
        H[(a, q)] = H.get((a, q), 0) + c * (-1) ** (t % 2)
    H = {k: v for k, v in H.items() if v}
    As = sorted({e[0] for e in P2}); Qs = sorted({e[1] for e in P2}); Ts = sorted({e[2] for e in P2})
    ga = list(range(As[0] - 2, As[-1] + 3, 2)); gq = list(range(Qs[0] - 2, Qs[-1] + 3, 2))
    gt = list(range(Ts[0] - margin_t, Ts[-1] + margin_t + 1))
    grid = [(a, q, t) for a in ga for q in gq for t in gt]
    gset, idx, rows, lo, hi = set(grid), {}, [], [], []

    def var(k):
        if k not in idx:
            idx[k] = len(idx)
        return idx[k]

    def add(co, l, h=None):
        rows.append(co); lo.append(l); hi.append(l if h is None else h)
    for e in grid:
        var(("c", e))
    for a in ga:
        for q in gq:
            add({var(("c", (a, q, t))): (-1) ** (t % 2) for t in gt}, H.get((a, q), 0))
    k = STAGES.index(stage)

    def special(n, m, target):
        img = {}
        for e in grid:
            for key, v in specialize({e: 1}, n, m).items():
                img.setdefault(key, {})[var(("c", e))] = v
        for key in set(img) | set(target):
            if key not in img:
                return False
            add(img[key], target.get(key, 0))
        return True
    if k >= 1 and not HOMOLOGY:
        for n, m in ((1, 0), (0, 2)):
            if not special(n, m, specialize(P2, n, m)):
                return None
    if k >= 1 and HOMOLOGY:
        for n, m in ((1, 0), (0, 2)):
            ea, eq, et = u(n, m)
            mono = specialize(P2, n, m)
            assert len(mono) == 1
            d = f"d{n}{m}"
            for e in grid:
                co = {var(("c", e)): 1, var((d + "_r", e)): -1, var((d, e)): -1}
                src = (e[0] - ea, e[1] - eq, e[2] - et)
                if src in gset:
                    co[var((d, src))] = co.get(var((d, src)), 0) - 1
                add(co, 0)
                if (e[0] + ea, e[1] + eq, e[2] + et) not in gset:
                    add({var((d, e)): 1}, 0)
                if specialize({e: 1}, n, m) != {next(iter(mono)): 1} and specialize({e: 1}, n, m) != {next(iter(mono)): -1}:
                    add({var((d + "_r", e)): 1}, 0)          # 残りの位置は特殊化で正解の単項式に写る直線上
            add({var((d + "_r", e)): 1 for e in grid}, 1)
    if k >= 2:
        for (n, m), (al, be) in (((1, 1), (2, 1)), ((0, 3), (1, 0))):
            S2 = specialize(P2, n, m)
            S1 = {(al * q, t + be * q): c for (q, t), c in specialize(P1, n, m).items()}
            (s, c) = ratio(S2, S1)
            target = {(q + s[0], t + s[1]): c * v for (q, t), v in S1.items()}
            assert target == S2
            if not special(n, m, target):
                return None
    if k >= 3:
        sq = {}
        for (a1, q1, t1), v1 in P1.items():
            for (a2, q2, t2), v2 in P1.items():
                kk = (a1 + a2, t1 + t2 - q1 - q2)
                sq[kk] = sq.get(kk, 0) + v1 * v2
        img = {}
        for e in grid:
            img.setdefault((e[0], e[2] - e[1]), {})[var(("c", e))] = 1
        for kk in set(img) | set(sq):
            if kk not in img:
                return None
            add(img[kk], sq.get(kk, 0))
    if k >= 4:
        if not special(2, 0, specialize(P2, 2, 0)):
            return None
    M = lil_matrix((len(rows), len(idx)))
    for i, co in enumerate(rows):
        for j, v in co.items():
            M[i, j] = v
    return M.tocsr(), np.array(lo, float), np.array(hi, float), idx, grid, P2


def solve(mod, cost, extra=None):
    M, lo, hi, idx, grid, _ = mod
    cons = [LinearConstraint(M, lo, hi)] + ([extra] if extra is not None else [])
    res = milp(cost, constraints=cons, integrality=np.ones(len(idx)), bounds=Bounds(0, np.inf), options={"time_limit": 300})
    if res.x is None:
        return None
    return {e: int(round(res.x[idx[("c", e)]])) for e in grid if round(res.x[idx[("c", e)]])}


def follow(name, P1, P2, log, samples=8):
    out = {}
    for st in STAGES:
        mod = build(P1, P2, st)
        if mod is None:
            log(f"{name} {st}: 目標が格子の外"); continue
        n = len(mod[3]); cost = np.zeros(n); row = np.zeros(n)
        for e in mod[4]:
            cost[mod[3][("c", e)]] = 1; row[mod[3][("c", e)]] = 1
        best = solve(mod, cost)
        if best is None:
            log(f"{name} {st}: 解なし"); continue
        mn = sum(best.values()); truth = mod[5]
        rng = np.random.default_rng(0); sols = [best]
        for _ in range(samples):
            c2 = np.zeros(n)
            for e in mod[4]:
                c2[mod[3][("c", e)]] = rng.integers(1, 100)
            s = solve(mod, c2, LinearConstraint(row.reshape(1, -1), mn, mn))
            if s is not None and s not in sols:
                sols.append(s)
        out[st] = {"min_rank": mn, "true_rank": sum(truth.values()), "true_among_min": truth in sols, "n_min_found": len(sols)}
        log(f"{name}（S² の総階数 {sum(truth.values())}） {st}: 最小 {mn}、正解を含む {truth in sols}、最小解 {len(sols)} 個以上（変数 {n}）")
    return out


def main():
    global HOMOLOGY
    if "--spec" in sys.argv:
        HOMOLOGY = False
    logf = open(ROOT / ("logs/dgr_lmov/dgr_s2_follow.log" if HOMOLOGY else "logs/dgr_lmov/dgr_s2_follow_spec.log"), "w")

    def log(s):
        print(s, flush=True); logf.write(s + "\n"); logf.flush()
    res = {}
    for p in (1, -1, 2, -2, 3, -3):
        P1, P2 = inputs(p)
        res[f"K_{p}"] = follow(f"K_{p}", P1, P2, log)
    H = load_homology()
    for K in ("10_124", "10_139", "10_152"):
        d = json.loads((ROOT / f"data/framed_knot/b4/{K}_s2_tc.json").read_text())
        res[K] = follow(K, superpoly(H[K]), {tuple(map(int, k.split(","))): v for k, v in d["[2]"].items()}, log)
    (ROOT / ("data/dgr_lmov/s2_follow.json" if HOMOLOGY else "data/dgr_lmov/s2_follow_spec.json")).write_text(json.dumps(res, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
