"""非彩色 [1] の superpolynomial を、Gukov らの原典（DGR, math/0505662）の手法で組み合わせ的に構成する（整数計画）。

未知数：太字の変数での P(𝐚,𝐪,𝐭) = Σ c[e_a,e_q,e_t] 𝐚^{e_a} 𝐪^{e_q} 𝐭^{e_t}、c ≥ 0 の整数。条件：
  (H)  P(𝐭 = −1) = reduced HOMFLY-PT（DGR の規約；本研究の rCS_□ の 𝐭 = −1 の値。ブレイドの HOMFLY-PT と照合済み）
  (d1) 次数 (−2, 2, −1) の微分：c = s1 + Q1 + Q1[ずらし]、s1 はちょうど 1 項
  (d−1) 次数 (−2, −2, −3)：同様にちょうど 1 項が残る
  (d2) 次数 (−2, 4, −1)：残り K ≥ 0 を 𝐚 = 𝐪^2 に特殊化すると reduced Khovanov ホモロジー（data/homology、khoca）
  (d0) 次数 (−2, 0, −3)：残り F ≥ 0 を (A, M) = (e_q/2, e_t − e_a) に移すと knot Floer ホモロジー（snappy）
  (Σ)  d_1 の残りは (Σ, −Σ, 0)、d_{−1} の残りは (Σ, Σ, Σ)（S 不変量 Σ。DGR、Gorsky–Gukov–Stošić 1304.3481 §3.6）
  (S)  DGR の対称性 (e_a, e_q, e_t) ↔ (e_a, −e_q, e_t − e_q)（三葉結び目で成立。KN の q ↔ t^{-1} の対称性の [1] の場合）
目的：総階数 Σc の最小化。三葉結び目では DGR の P = 𝐚²𝐪⁻² + 𝐚²𝐪²𝐭² + 𝐚⁴𝐭³ を再現する。
これは DGR の構造的な仮定（各微分が単独でホモロジーを与える）に基づく候補であり、HOMFLY-PT ホモロジーの直接計算ではない。
出力：data/superpoly/{K}.json
使い方：python scripts/construct_superpolynomial.py [--margin-t=10 --margin-a=4 --time=600] T2_3 K1 K2 ...
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import lil_matrix

from rcsknots.bps import to_bold
from rcsknots.laurent import load_table

DIFFS = {"d1": (-2, 2, -1), "d-1": (-2, -2, -3), "d2": (-2, 4, -1), "d0": (-2, 0, -3)}


def rcs_box(name):
    if name.startswith("K"):
        d = json.loads(Path(f"data/cables/{name}_1.json").read_text())["poly"]
        return {tuple(map(int, k.split(","))): v for k, v in d.items()}
    return load_table("data/reference/notebook_torus.txt")[(name, "1")]


def daha_bold(name):
    """rCS_□ in bold variables with the absolute sign: its 𝐭 = −1, 𝐚 = 𝐪^2 value equals the Euler characteristic
    of the reduced Khovanov homology (the Jones polynomial), as for the true superpolynomial."""
    P = rcs_box(name)
    b = to_bold(P, next(iter(P))[0] % 2)
    kh = json.loads(Path(f"data/homology/{name}.json").read_text())["kh"]
    jones, euler = {}, {}
    for k, r in kh.items():
        q, t = map(int, k.split(","))
        jones[q] = jones.get(q, 0) + (-1) ** t * r
    for (ea, eq, et), c in b.items():
        euler[2 * ea + eq] = euler.get(2 * ea + eq, 0) + (-1) ** et * c
    jones = {k: v for k, v in jones.items() if v}
    euler = {k: v for k, v in euler.items() if v}
    if euler == jones:
        return b
    assert {k: -v for k, v in euler.items()} == jones, "Euler characteristic does not match the Jones polynomial"
    return {k: -v for k, v in b.items()}


def homfly(bold):
    H = {}
    for (ea, eq, et), c in bold.items():
        H[(ea, eq)] = H.get((ea, eq), 0) + c * (-1) ** et
    return {k: v for k, v in H.items() if v}


def construct(name, margin_t=10, margin_a=4, diffs=("d1", "d-1", "d2", "d0"), symmetric=True, time_limit=600, fix_rank=None, seed=None, sigma=True):
    hom = json.loads(Path(f"data/homology/{name}.json").read_text())
    kh = {tuple(map(int, k.split(","))): v for k, v in hom["kh"].items()}
    hfk = {tuple(map(int, k.split(","))): v for k, v in hom["hfk"].items()}
    daha = daha_bold(name)
    H = homfly(daha)
    A = sorted({ea for ea, _ in H})
    A = list(range(A[0] - margin_a, A[-1] + margin_a + 1, 2))
    qs = sorted({eq for _, eq in H})
    Qs = list(range(qs[0] - 2, qs[-1] + 3, 2))
    ts = [t for _, t in kh] + [et for (_, _, et) in daha]
    Ts = list(range(min(ts) - margin_t, max(ts) + margin_t + 1))
    grid = [(a, q, t) for a in A for q in Qs for t in Ts]
    idx = {}

    def var(key):
        if key not in idx:
            idx[key] = len(idx)
        return idx[key]

    for e in grid:
        var(("c", e))
        for d in DIFFS:
            var((d, e))                     # pair sources Q_d
            var((d + "_rest", e))           # survivors
    rows, rhs_lo, rhs_hi = [], [], []

    def add(coeffs, lo, hi):
        rows.append(coeffs)
        rhs_lo.append(lo)
        rhs_hi.append(hi)

    gset = set(grid)
    # (H)
    for a in A:
        for q in Qs:
            add({var(("c", (a, q, t))): (-1) ** t for t in Ts}, H.get((a, q), 0), H.get((a, q), 0))
    # differentials: c[e] − rest[e] − Q[e] − Q[e − deg] = 0 ; Q[e] = 0 if e + deg outside the grid
    for d, (da, dq, dt) in DIFFS.items():
        if d not in diffs:
            continue
        for e in grid:
            co = {var(("c", e)): 1, var((d + "_rest", e)): -1, var((d, e)): -1}
            src = (e[0] - da, e[1] - dq, e[2] - dt)
            if src in gset:
                co[var((d, src))] = co.get(var((d, src)), 0) - 1
            add(co, 0, 0)
            tgt = (e[0] + da, e[1] + dq, e[2] + dt)
            if tgt not in gset:
                add({var((d, e)): 1}, 0, 0)
    # (S) DGR の対称性（KN の rCS_{λ^T}(a,q,t) = rCS_λ(a,t^{-1},q^{-1}) の [1] の場合）：(e_a, e_q, e_t) ↔ (e_a, −e_q, e_t − e_q)
    if symmetric:
        for e in grid:
            m = (e[0], -e[1], e[2] - e[1])
            if m not in gset:
                add({var(("c", e)): 1}, 0, 0)
            elif m > e:
                add({var(("c", e)): 1, var(("c", m)): -1}, 0, 0)
    for d in ("d1", "d-1"):
        if d not in diffs:
            continue
        add({var((d + "_rest", e)): 1 for e in grid}, 1, 1)
        # (Σ) GGS（1304.3481 §3.6）・DGR：d_1 の残りは (Σ, −Σ, 0)、対称性により d_{−1} の残りは (Σ, Σ, Σ)
        if sigma:
            for e in grid:
                ok = (e[1] == -e[0] and e[2] == 0) if d == "d1" else (e[1] == e[0] and e[2] == e[0])
                if not ok:
                    add({var((d + "_rest", e)): 1}, 0, 0)
    # (d2) → Khovanov
    spec = {} if "d2" in diffs else None
    for e in grid if spec is not None else []:
        spec.setdefault((2 * e[0] + e[1], e[2]), []).append(e)
    for k in (set(spec) | set(kh)) if spec is not None else []:
        es = spec.get(k, [])
        if not es and kh.get(k, 0):
            raise ValueError(f"Kh generator {k} outside the grid")
        if es:
            add({var(("d2_rest", e)): 1 for e in es}, kh.get(k, 0), kh.get(k, 0))
    # (d0) → knot Floer
    spec = {} if "d0" in diffs else None
    for e in grid if spec is not None else []:
        if e[1] % 2 == 0:
            spec.setdefault((e[1] // 2, e[2] - e[0]), []).append(e)
    for k in (set(spec) | set(hfk)) if spec is not None else []:
        es = spec.get(k, [])
        if not es and hfk.get(k, 0):
            raise ValueError(f"HFK generator {k} outside the grid")
        if es:
            add({var(("d0_rest", e)): 1 for e in es}, hfk.get(k, 0), hfk.get(k, 0))

    n = len(idx)
    M = lil_matrix((len(rows), n))
    for i, co in enumerate(rows):
        for j, v in co.items():
            M[i, j] = v
    cost = np.zeros(n)
    rng = np.random.default_rng(seed)
    for e in grid:
        cost[idx[("c", e)]] = 1 if seed is None else rng.integers(1, 100)
    if fix_rank is not None:                       # explore other solutions of the same total rank
        M.resize((len(rows) + 1, n))
        for e in grid:
            M[len(rows), idx[("c", e)]] = 1
        rhs_lo, rhs_hi = rhs_lo + [fix_rank], rhs_hi + [fix_rank]
    res = milp(cost, constraints=LinearConstraint(M.tocsr(), rhs_lo, rhs_hi), integrality=np.ones(n),
               bounds=Bounds(0, np.inf), options={"time_limit": time_limit})
    if res.x is None:
        return None, daha, res.message
    P = {e: int(round(res.x[idx[("c", e)]])) for e in grid if round(res.x[idx[("c", e)]])}
    return P, daha, res.message


def main(args):
    """args: knot names, optionally preceded by --margin-t=, --margin-a=, --time= (seconds)."""
    opts = dict(a.lstrip("-").split("=") for a in args if a.startswith("--"))
    kw = {"margin_t": int(opts.get("margin-t", 10)), "margin_a": int(opts.get("margin-a", 4)), "time_limit": int(opts.get("time", 600))}
    out = Path("data/superpoly")
    out.mkdir(parents=True, exist_ok=True)
    for name in [a for a in args if not a.startswith("--")]:
        P, daha, msg = construct(name, **kw)
        if P is None:
            print(name, "解なし：", msg, flush=True)
            (out / f"{name}.json").write_text(json.dumps({"knot": name, "solution": None, "message": msg}))
            continue
        same = P == {k: v for k, v in daha.items()}
        print(name, f"総階数 {sum(P.values())}（DAHA の |係数| の総和 {sum(map(abs, daha.values()))}）",
              "DAHA と一致" if same else "DAHA と不一致", flush=True)
        (out / f"{name}.json").write_text(json.dumps(
            {"knot": name, "solution": {",".join(map(str, k)): v for k, v in sorted(P.items())},
             "rank": sum(P.values()), "equals_daha": same, "message": msg, "options": kw}))


if __name__ == "__main__":
    main(sys.argv[1:])
