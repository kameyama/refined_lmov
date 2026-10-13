"""D12 段階 2：仮説 H_u″（P 正規化の hook の Schur 係数の冪級数の非負性）を rCS_{[2,1]} の構成の条件に課す。

未知数：x = x_p + Φ·g（x_p は正解または条件の等式を満たす一つの解、Φ は dgr_box3.phi21、g は窓に収まる多項式）。
条件：箱（(𝐚, 𝐪) は非 refined の台、𝐭 は r₂r₁ の範囲）、転置、(R)、
  S：N(x) = x·(𝐪² + 𝐚²𝐭)(1 + 𝐚²𝐪²𝐭³) を Π_□(1 − q^{a(□)} t^{l(□)+1}) = (1 − 𝐪²)²(1 − 𝐪⁶𝐭²) で割った 𝐪 → 0 の冪級数の係数 ≥ 0
     （𝐪 の次数 ≤ N の 𝐪 の最大次数 + ORDER まで。N は g の一次式なので係数も一次式）。全体の符号は x_p の N の多数派。
  P（任意）：\\hat f_ρ の spin 表示の非負性（framing 0）。
判定：多面体の有界性（乱数方向の LP）、整数解の列挙（上限 100）、正解を含むか。
使い方：PYTHONPATH=src:scripts/framed_knot:scripts/dgr_lmov python scripts/dgr_lmov/lmov21_hu.py T2_5 [--R] [--P] [--order=40]
記録：logs/dgr_lmov/lmov21_hu.log、data/dgr_lmov/lmovstruct/lmov21_hu.json。
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sps
from scipy.optimize import Bounds, LinearConstraint, linprog, milp

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "framed_knot"))
from dgr_box3 import _mul, formula3, hat3, phi21, rcs_window21  # noqa: E402
from dgr_hu import macdonald_denominators_bold, unknot_numerator_factors  # noqa: E402
from framing_check import load_reduced  # noqa: E402
from t4_colored import bold_to_small  # noqa: E402
from v4_closed import L  # noqa: E402
from rcsknots.bps import to_bold  # noqa: E402

LOG = Path("logs/dgr_lmov/lmov21_hu.log")
OUT = Path("data/dgr_lmov/lmovstruct/lmov21_hu.json")


def log(s):
    print(s, flush=True)
    with LOG.open("a") as f:
        f.write(s + "\n")


def bold(P):
    return {k: v for k, v in to_bold(P, 1 if any(k[0] % 2 for k in P) else 0).items() if v}


def series_rows(poly_cols, lam, order):
    """poly_cols：{列: {(a, q, t): 係数}}（列 None は定数項）。各列の N·Π(Σ m^k) を 𝐪 の次数 ≤ top まで展開し、
    {(a, q, t): {列: 係数}} を返す。top = 定数項の N の 𝐪 の最大次数 + order。"""
    dens = macdonald_denominators_bold(lam)
    const = poly_cols[None]
    top = max(k[1] for k in const) + order
    out = {}
    for col, P in poly_cols.items():
        cur = dict(P)
        for kq, kt in dens:
            acc = {}
            for (a, q, t), c in cur.items():
                k = 0
                while q + k * kq <= top:
                    key = (a, q + k * kq, t + k * kt)
                    acc[key] = acc.get(key, 0) + c
                    k += 1
            cur = acc
        for key, c in cur.items():
            if c and key[1] <= top:
                out.setdefault(key, {})[col] = c
    return out


def run(K, use_R=False, use_P=False, order=40, log=log):
    t0 = time.time()
    R = load_reduced(K, 3)
    small = {l: P for l, P in R.items() if l != (2, 1)}
    xp = bold(R[(2, 1)])
    H21u = {}
    for (a, q, t), c in xp.items():
        H21u[(a, q)] = H21u.get((a, q), 0) + c * (-1) ** (t % 2)
    H21u = {k: v for k, v in H21u.items() if v}
    b1, b2 = bold(small[(1,)]), bold(small[(2,)])
    Phi = phi21()
    win = rcs_window21(b1, b2)
    prod = _mul(b2, b1)
    ts = [m[2] for m in prod]
    box = [(min(k[0] for k in H21u), max(k[0] for k in H21u)), (min(k[1] for k in H21u), max(k[1] for k in H21u)), (min(ts), max(ts))]
    ext = [(min(k[i] for k in Phi), max(k[i] for k in Phi)) for i in range(3)]
    G = [(a, q, t) for a in range(box[0][0] - ext[0][0], box[0][1] - ext[0][1] + 1, 2)
         for q in range(box[1][0] - ext[1][0], box[1][1] - ext[1][1] + 1, 2)
         for t in range(box[2][0] - ext[2][0], box[2][1] - ext[2][1] + 1)]
    gi = {e: i for i, e in enumerate(G)}
    n = len(G)
    log(f"{K}: g の変数 {n}（箱 {box}、R {use_R}、P {use_P}、order {order}）")
    if n == 0:
        return {"n_g": 0, "unique": True}
    xrow = {}
    for e in G:
        for f, c in Phi.items():
            m = (e[0] + f[0], e[1] + f[1], e[2] + f[2])
            xrow.setdefault(m, {})[gi[e]] = c
    eqs, ineqs = [], []
    for m, row in xrow.items():
        if not box[0][0] <= m[0] <= box[0][1] or not box[1][0] <= m[1] <= box[1][1] or not box[2][0] <= m[2] <= box[2][1]:
            eqs.append((row, -xp.get(m, 0)))
        elif use_R and (m[0] not in win or not win[m[0]][0] <= m[2] <= win[m[0]][1]):
            eqs.append((row, -xp.get(m, 0)))
    for m, row in xrow.items():
        s = (m[0], -m[1], m[2] - m[1])
        if s > m:
            r2 = dict(row)
            for j, v in xrow.get(s, {}).items():
                r2[j] = r2.get(j, 0) - v
            r2 = {j: v for j, v in r2.items() if v}
            if r2:
                eqs.append((r2, xp.get(s, 0) - xp.get(m, 0)))
    # S：N(x) の冪級数。列ごとの多項式 N(Φ·e)、定数 N(x_p)
    F = {(0, 0, 0): 1}
    for f in unknot_numerator_factors((2, 1)):
        F = _mul(F, f)
    cols = {None: _mul(xp, F)}
    for e in G:
        cols[gi[e]] = _mul(_mul({e: 1}, Phi), F)
    Np = cols[None]
    pos, neg = sum(1 for v in Np.values() if v > 0), sum(1 for v in Np.values() if v < 0)
    sgn = 1 if pos >= neg else -1
    srows = series_rows(cols, (2, 1), order)
    nS = 0
    for key, r in srows.items():
        c0 = r.pop(None, 0)
        if r:
            ineqs.append(({j: sgn * v for j, v in r.items()}, -sgn * c0))
            nS += 1
        elif sgn * c0 < 0:
            log(f"  警告：x_p の冪級数に負の係数 {key}（g で直せない）")
    # P：\hat f の spin 表示の非負性（framing 0）
    nP = 0
    if use_P:
        phis = L.from_dict(bold_to_small(Phi))
        hp = hat3({**small, (2, 1): bold_to_small(xp)})
        for rho, (ns, den) in formula3().items():
            Qr = (ns[1] * phis / den).to_dict()
            shift = 0 if all(k[0] % 2 == 0 for k in hp[rho]) else 1
            hb = to_bold(hp[rho], shift)
            p_, n_ = sum(v > 0 for v in hb.values()), sum(v < 0 for v in hb.values())
            s = 1 if n_ <= p_ else -1
            rows = {}
            for e in G:
                (se, ce), = bold_to_small({e: 1}).items()
                for k, v in Qr.items():
                    (kb, cb), = to_bold({(k[0] + se[0], k[1] + se[1], k[2] + se[2]): v * ce}, shift).items()
                    rows.setdefault(kb, {})[gi[e]] = rows.get(kb, {}).get(gi[e], 0) + cb
            for kb, r in rows.items():
                r = {j: v for j, v in r.items() if v}
                if r:
                    ineqs.append(({j: s * v for j, v in r.items()}, -s * hb.get(kb, 0)))
                    nP += 1
    log(f"  等式 {len(eqs)}、不等式 S {nS}、P {nP} ({time.time() - t0:.0f}s)")

    def mat(rows):
        A = sps.lil_matrix((len(rows), n))
        for i, (row, _) in enumerate(rows):
            for j, v in row.items():
                A[i, j] = v
        return A.tocsr()
    Ae, Ai = mat(eqs), mat(ineqs)
    be = np.array([r for _, r in eqs], float); bi = np.array([r for _, r in ineqs], float)
    out = {"n_g": n, "n_eq": len(eqs), "n_S": nS, "n_P": nP, "order": order}
    # 正解（g = 0）は実行可能か
    ok = (np.all(np.abs(Ae @ np.zeros(n) - be) < 1e-9) if len(eqs) else True) and (np.all(Ai @ np.zeros(n) - bi >= -1e-9) if len(ineqs) else True)
    out["true_feasible"] = bool(ok)
    # 有界性：乱数方向の LP
    rng = np.random.default_rng(0)
    widths = []
    for _ in range(6):
        c = rng.standard_normal(n)
        w = []
        for sg in (1, -1):
            r = linprog(sg * c, A_ub=-Ai if len(ineqs) else None, b_ub=-bi if len(ineqs) else None, A_eq=Ae if len(eqs) else None,
                        b_eq=be if len(eqs) else None, bounds=[(None, None)] * n, method="highs")
            w.append(None if r.status == 3 else float(sg * r.fun))
        widths.append(w)
    out["bounded"] = all(a is not None and b is not None for a, b in widths)
    out["widths"] = [(None if a is None else round(a, 3), None if b is None else round(b, 3)) for a, b in widths]
    log(f"  正解は実行可能 {ok}、多面体は有界 {out['bounded']}（幅の例 {out['widths'][:3]}）({time.time() - t0:.0f}s)")
    # 整数解の列挙（有界なら）
    cons = ([LinearConstraint(Ae, be, be)] if len(eqs) else []) + ([LinearConstraint(Ai, bi, np.inf)] if len(ineqs) else [])
    bnd = Bounds(np.full(n, -np.inf), np.full(n, np.inf))
    if out["bounded"]:
        c = rng.integers(-30, 31, size=n).astype(float)
        sols, lo = [], -np.inf
        while len(sols) < 100:
            rr = milp(c, constraints=cons + [LinearConstraint(c.reshape(1, -1), lo, np.inf)], bounds=bnd, integrality=np.ones(n),
                      options={"time_limit": 1800})
            if rr.x is None:
                break
            g = tuple(int(round(v)) for v in rr.x)
            if g in sols:
                break
            sols.append(g)
            lo = float(c @ np.array(g)) + 1
        out["count"] = len(sols) if len(sols) < 100 else ">=100"
        out["contains_true"] = tuple([0] * n) in sols
        log(f"  整数解 {out['count']} 個（正解を含む {out['contains_true']}）({time.time() - t0:.0f}s)")
    return out


if __name__ == "__main__":
    K = sys.argv[1]
    use_R = "--R" in sys.argv
    use_P = "--P" in sys.argv
    order = next((int(a.split("=")[1]) for a in sys.argv if a.startswith("--order=")), 40)
    log(f"=== lmov21_hu {time.strftime('%Y-%m-%d %H:%M:%S')} {' '.join(sys.argv[1:])} ===")
    res = run(K, use_R, use_P, order)
    prev = json.loads(OUT.read_text()) if OUT.exists() else {}
    prev[f"{K} R={use_R} P={use_P} order={order}"] = json.loads(json.dumps(res, default=str))
    OUT.write_text(json.dumps(prev, ensure_ascii=False, indent=1))
