"""H11：unrefined から refined の [2,1] を一意（高々 O(1) 個）に決める追加の仮定の探索。

サブコマンド
  support K ...   既知の結び目で、rCS_{[2,1]}（太字）の台が入力（[1], [2], [1,1]）の積の台から決まる窓に収まるかを調べ、
                  収まる規則の窓で一次の条件（U, T, E, U2, H, L）の零空間の次元を数える。
                    col   ：(𝐚, 𝐪) の列ごとの 𝐭 の範囲 ⊂ r₂r₁ の同じ列の範囲（col2：r₂r₁ ∪ r₁₁r₁）
                    hull  ：凸包 conv(supp r₂r₁)（hull2：conv(supp r₂r₁ ∪ supp r₁₁r₁)）
                    delta ：δ = e_𝐭 − e_𝐪/2 − e_𝐚 の範囲
  framing K       既知の結び目で、正解と整数計画の別解（data/dgr_lmov/box3/check_<K>_UTHEU2L.json）の framed \\hat f^{(x)}_ρ の
                  符号（全体の符号と負の項の数、G 展開の後も）を x の範囲で比べる。別解の差 Δ の像は 𝐭^{5x}·Δ\\hat f^{(0)} なので、
                  正解の全体の符号が x で変わるなら framing の正値性は上界になる。
使い方：PYTHONPATH=src:scripts python scripts/lmov21_h11.py support T2_3 4_1 T2_5 T2_7 T3_4 K1
        PYTHONPATH=src:scripts python scripts/lmov21_h11.py framing T2_5 [XMIN XMAX]
出力：logs/lmov21_h11.log、data/lmovstruct/lmov21_h11_*.json。
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lmov21_linear as LL  # noqa: E402
from dgr_box3 import _mul, build21, hat3  # noqa: E402
from framing_check import framed, framing_monomial, load_reduced  # noqa: E402
from t4_colored import bold_to_small  # noqa: E402
from rcsknots.bps import g_expand, to_bold  # noqa: E402

OUT = Path("data/lmovstruct")
NMAX_ENUM = 3000
LOG = Path("logs/lmov21_h11.log")


def log(s):
    print(s, flush=True)
    with LOG.open("a") as f:
        f.write(s + "\n")


def bold(P):
    return {k: v for k, v in to_bold(P, 1 if any(k[0] % 2 for k in P) else 0).items() if v}


def col_ranges(S):
    out = {}
    for a, q, t in S:
        lo, hi = out.get((a, q), (10 ** 9, -10 ** 9))
        out[(a, q)] = (min(lo, t), max(hi, t))
    return out


def in_hull(points, m):
    """m が points の凸包に入るか（LP の実行可能性）。"""
    P = np.array(points, dtype=float).T          # 3 × n
    n = P.shape[1]
    A_eq = np.vstack([P, np.ones((1, n))])
    b_eq = np.array(list(m) + [1.0])
    r = linprog(np.zeros(n), A_eq=A_eq, b_eq=b_eq, bounds=[(0, None)] * n, method="highs")
    return r.status == 0


def delta(m):
    a, q, t = m
    return t - q / 2 - a


def support(K):
    R = load_reduced(K, 3)
    b = {l: bold(P) for l, P in R.items()}
    t21 = b[(2, 1)]
    prodA, prodB, prodC = _mul(b[(2,)], b[(1,)]), _mul(b[(1, 1)], b[(1,)]), _mul(_mul(b[(1,)], b[(1,)]), b[(1,)])
    res = {"knot": K, "terms21": len(t21)}
    # 列ごとの 𝐭 の範囲
    for name, S in (("col", set(prodA)), ("col2", set(prodA) | set(prodB)), ("col3", set(prodA) | set(prodB) | set(prodC))):
        cr = col_ranges(S)
        bad = [m for m in t21 if (m[0], m[1]) not in cr or not cr[(m[0], m[1])][0] <= m[2] <= cr[(m[0], m[1])][1]]
        res[name] = {"violations": len(bad), "window_size": sum(hi - lo + 1 for lo, hi in cr.values()), "examples": bad[:5]}
    # 凸包
    for name, S in (("hull", list(prodA)), ("hull2", list(set(prodA) | set(prodB)))):
        bad = [m for m in t21 if not in_hull(S, m)]
        res[name] = {"violations": len(bad), "examples": bad[:5]}
    # δ の範囲
    d21 = sorted({delta(m) for m in t21})
    res["delta"] = {"rCS21": (d21[0], d21[-1]),
                    **{nm: (min(delta(m) for m in S), max(delta(m) for m in S)) for nm, S in
                       (("r1", b[(1,)]), ("r2", b[(2,)]), ("r11", b[(1, 1)]), ("r3", b[(3,)]), ("r111", b[(1, 1, 1)]),
                        ("r2r1", prodA), ("r11r1", prodB), ("r1^3", prodC))}}
    log(f"{K}: [2,1] {len(t21)} 項。列の範囲の違反 col {res['col']['violations']} / col2 {res['col2']['violations']} / col3 {res['col3']['violations']}"
        f"（窓の大きさ {res['col']['window_size']} / {res['col2']['window_size']} / {res['col3']['window_size']}）、"
        f"凸包の違反 hull {res['hull']['violations']} / hull2 {res['hull2']['violations']}")
    log(f"  δ の範囲：{res['delta']}")
    return res, b, t21


def nullspace_with_mask(K, b, t21, mask_name, log=log):
    """窓を mask（許す単項式の集合）に制限して、U+T+E+U2(+H)+L の零空間の次元。"""
    R = load_reduced(K, 3)
    small = {l: P for l, P in R.items() if l != (2, 1)}
    H21u = {}
    for (a, q, t), c in t21.items():
        H21u[(a, q)] = H21u.get((a, q), 0) + c * (-1) ** (t % 2)
    H21u = {k: v for k, v in H21u.items() if v}
    prodA, prodB = _mul(b[(2,)], b[(1,)]), _mul(b[(1, 1)], b[(1,)])
    S = set(prodA) if mask_name == "col" else set(prodA) | set(prodB)
    cr = col_ranges(S)
    ts = [m[2] for m in S]
    window = (min(ts), max(ts) + 1)
    out = {}
    for conds in (("U", "T", "E", "U2", "L"), ("U", "T", "E", "U2", "H", "L")):
        t0 = time.time()
        lin, mons = build21(small, H21u, window, (1, 1, 1), conds=tuple(c for c in conds if c != "L"))
        rows = [(dict(r), bb) for r, bb in zip(lin.rows, lin.rhs)] + LL.laurent_rows(small, mons, lin.col)
        nmask = 0
        for m in mons:
            c = cr.get((m[0], m[1]))
            if c is None or not c[0] <= m[2] <= c[1]:
                rows.append(({lin.col[("x", m)]: 1}, 0))
                nmask += 1
        xv = {lin.col[("x", m)]: t21.get(m, 0) for m in mons}
        ok = all(sum(v * xv.get(c, 0) for c, v in row.items()) == bb for row, bb in rows) and all(m in set(mons) for m in t21)
        rank, feas = LL.rank_nullity(rows, len(mons))
        key = "+".join(conds)
        out[key] = {"unknowns": len(mons), "masked": nmask, "nullity": len(mons) - rank, "feasible": feas, "true_satisfies": ok}
        log(f"  {K:6s} 窓 {mask_name:5s} {key:16s} 未知数 {len(mons)}（窓の外 {nmask}）零空間 {len(mons) - rank} 実行可能 {feas} 正解 {ok} ({time.time() - t0:.0f}s)")
        if len(mons) - rank == 0 and feas:
            x = LL.solve_unique(rows, len(mons))
            inv = {v: k for k, v in lin.col.items()}
            match = x is not None and {inv[c][1]: v for c, v in x.items() if v} == t21
            out[key]["unique_matches_true"] = match
            log(f"  {K:6s} 一意解、既知と一致 {match}")
    return out


def framing_compare(K, xs):
    R = load_reduced(K, 3)
    small = {l: P for l, P in R.items() if l != (2, 1)}
    true21 = bold(R[(2, 1)])
    path = Path(f"data/dgr_lmov/box3/check_{K}_UTHEU2L.json")
    alts = []
    if path.exists():
        d = json.loads(path.read_text())
        alts = [{tuple(m): v for m, v in sol} for sol in d["solutions"]]
    cands = [("true", true21)] + [(f"alt{i}", a) for i, a in enumerate(alts)]
    res = {}
    for x in xs:
        fr = framed(small, x)
        s21, (da, dq, dt) = framing_monomial((2, 1), x)
        for name, c21 in cands:
            sm = bold_to_small(c21)
            fr21 = {(ea + da, eq + dq, et + dt): s21 * v for (ea, eq, et), v in sm.items()}
            hats = hat3({**fr, (2, 1): fr21})
            line = []
            for rho, h in hats.items():
                hb = bold(h)
                pos, neg = sum(v > 0 for v in hb.values()), sum(v < 0 for v in hb.values())
                sign = 1 if neg <= pos else -1
                try:
                    N = g_expand({k: sign * v for k, v in hb.items()})
                    gneg = sum(v < 0 for v in N.values())
                except ValueError:
                    gneg = None
                res[(x, name, rho)] = (sign, min(pos, neg), len(hb), gneg)
                line.append(f"{rho}:{'+' if sign > 0 else '-'} 前{min(pos, neg)}/{len(hb)} 後{gneg}")
            log(f"{K} x={x:+d} {name:5s} " + "  ".join(line))
    return res


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    log(f"=== lmov21_h11 {time.strftime('%Y-%m-%d %H:%M:%S')} {' '.join(sys.argv[1:])} ===")
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "support":
        out = {}
        for K in args:
            res, b, t21 = support(K)
            for mask in ("col", "col2"):
                if res[mask]["violations"] == 0:
                    res[f"nullspace_{mask}"] = nullspace_with_mask(K, b, t21, mask)
            out[K] = res
            path = OUT / "lmov21_h11_support.json"
            prev = json.loads(path.read_text()) if path.exists() else {}
            prev.update({K: json.loads(json.dumps(res, default=str))})
            path.write_text(json.dumps(prev, ensure_ascii=False, indent=1))
    elif cmd == "framing":
        K = args[0]
        x0, x1 = (int(args[1]), int(args[2])) if len(args) > 2 else (-6, 6)
        res = framing_compare(K, range(x0, x1 + 1))
        path = OUT / f"lmov21_h11_framing_{K}.json"
        path.write_text(json.dumps({f"{x},{n},{r}": v for (x, n, r), v in res.items()}, ensure_ascii=False, indent=1))


if __name__ == "__main__" and not (len(sys.argv) > 1 and sys.argv[1] in ("face", "predict")):
    main()


# ---- face：全条件 ＋ framing の正値性 ＋ BPS 状態の総数の最小化で、残る整数解を数える ----------------

def face(K, xs, use_R=False, genus=False, irr=False, seeds=6, log=log):
    """既知の結び目：x = x_true + Φ·g。face_core を参照。"""
    R = load_reduced(K, 3)
    small = {l: P for l, P in R.items() if l != (2, 1)}
    xp = bold(R[(2, 1)])
    H21u = {}
    for (a, q, t), c in xp.items():
        H21u[(a, q)] = H21u.get((a, q), 0) + c * (-1) ** (t % 2)
    H21u = {k: v for k, v in H21u.items() if v}
    return face_core(K, small, xp, H21u, xs, use_R, genus, irr, seeds, known=True, log=log)


def particular_solution(small, H21u, window, conds=("U", "T", "H", "E", "U2", "L", "R"), sigmas=(1, 1, 1), log=log):
    """条件を全て満たす rCS_{[2,1]}（太字）を一つ、整数計画の実行可能性で求める（未知の結び目の出発点）。"""
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp
    from dgr_box3 import _matrix
    t0 = time.time()
    lin, mons = build21(small, H21u, window, sigmas, conds=conds)
    n = len(lin.col)
    lb = np.array([-400 if nm[0] == "x" else 0 for nm in lin.col], float)
    ub = np.array([400 if nm[0] == "x" else np.inf for nm in lin.col], float)
    rhs = np.array(lin.rhs, float)
    r = milp(np.zeros(n), constraints=[LinearConstraint(_matrix(lin).tocsr(), rhs, rhs)], bounds=Bounds(lb, ub),
             integrality=np.ones(n), options={"time_limit": 3000})
    log(f"  実行可能解の探索：{r.message}（変数 {n}、{time.time() - t0:.0f}s）")
    if r.x is None:
        return None
    return {m: int(round(r.x[lin.col[("x", m)]])) for m in mons if round(r.x[lin.col[("x", m)]]) != 0}


def face_core(K, small, xp, H21u, xs, use_R=False, genus=False, irr=False, seeds=6, known=False, log=log):
    """x = xp + Φ·g（xp は条件の等式を満たす一つの解。既知の結び目では正解）。
    条件：転置、箱（(𝐚, 𝐪) は非 refined の台、𝐭 は r₂r₁ の範囲）、(R)、各 framing x ∈ xs での \\hat f_ρ^{(x)} の
    spin 表示（genus, irr も選べる）の非負性。全体の符号は既知なら xp（正解）から、未知なら ＋ を仮定する【仮定】。
    目的：framing 0 の BPS 状態の総数 Σ_ρ σ_ρ Σ \\hat𝐍 の最小化。最適面の上の整数解を列挙する。"""
    import numpy as np
    import scipy.sparse as sps
    from scipy.optimize import Bounds, LinearConstraint, milp
    from dgr_box3 import formula3, phi21, rcs_window21
    from v4_closed import L
    b1, b2 = bold(small[(1,)]), bold(small[(2,)])
    Phi = phi21()
    win = rcs_window21(b1, b2)
    prod = _mul(b2, b1)
    ts = [m[2] for m in prod]
    box = [(min(k[0] for k in H21u), max(k[0] for k in H21u)), (min(k[1] for k in H21u), max(k[1] for k in H21u)), (min(ts), max(ts))]
    ext = [(min(k[i] for k in Phi), max(k[i] for k in Phi)) for i in range(3)]
    assert all(k[1] % 2 == 0 for k in Phi), "Φ の 𝐪 の指数は偶数のはず"
    G = [(a, q, t) for a in range(box[0][0] - ext[0][0], box[0][1] - ext[0][1] + 1, 2)
         for q in range(box[1][0] - ext[1][0], box[1][1] - ext[1][1] + 1, 2)     # 𝐪 の偶奇は非 refined の台に合わせる
         for t in range(box[2][0] - ext[2][0], box[2][1] - ext[2][1] + 1)]
    gi = {e: i for i, e in enumerate(G)}
    n = len(G)
    log(f"{K}: g の変数 {n}（箱 {box}、R {use_R}、genus {genus}、irr {irr}、framing {list(xs)}）")
    if n == 0:
        return {"n_g": 0, "unique": True, "face_count": 1,
                "face_solutions_rcs21": [{",".join(map(str, m)): v for m, v in sorted(xp.items()) if v}]}
    xrow = {}
    for e in G:
        for f, c in Phi.items():
            m = (e[0] + f[0], e[1] + f[1], e[2] + f[2])
            xrow.setdefault(m, {})[gi[e]] = c
    eqs, ineqs = [], []
    for m, row in xrow.items():
        if not box[0][0] <= m[0] <= box[0][1] or not box[1][0] <= m[1] <= box[1][1] or not box[2][0] <= m[2] <= box[2][1]:
            eqs.append((row, -xp.get(m, 0)))                   # 箱の外は 0
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
    # framing 0 の像の行（g_e → \hat f_ρ の太字の係数）
    phis = L.from_dict(bold_to_small(Phi))
    rows0, shifts, sig = {}, {}, {}
    hp0 = hat3({**small, (2, 1): bold_to_small(xp)})
    for rho, (ns, den) in formula3().items():
        Qr = (ns[1] * phis / den).to_dict()
        shift = 0 if all(k[0] % 2 == 0 for k in hp0[rho]) else 1
        shifts[rho] = shift
        rows = {}
        for e in G:
            (se, ce), = bold_to_small({e: 1}).items()
            for k, v in Qr.items():
                (kb, cb), = to_bold({(k[0] + se[0], k[1] + se[1], k[2] + se[2]): v * ce}, shift).items()
                rows.setdefault(kb, {})[gi[e]] = rows.get(kb, {}).get(gi[e], 0) + cb
        rows0[rho] = {kb: {j: v for j, v in r.items() if v} for kb, r in rows.items()}
    obj = np.zeros(n)
    for x in xs:
        fr = framed(small, x)
        s21, (da, dq, dt) = framing_monomial((2, 1), x)
        sm = bold_to_small(xp)
        fr21 = {(ea + da, eq + dq, et + dt): s21 * v for (ea, eq, et), v in sm.items()}
        hx = hat3({**fr, (2, 1): fr21})
        for rho in rows0:
            hb = to_bold(hx[rho], shifts[rho])
            pos, neg = sum(v > 0 for v in hb.values()), sum(v < 0 for v in hb.values())
            s = (1 if neg <= pos else -1) if known else 1
            sig[(x, rho)] = s
            if known:
                assert min(pos, neg) == 0, f"正解が framing {x} で正でない ρ={rho}"
            for kb, r in rows0[rho].items():
                kx = (kb[0], kb[1], kb[2] + 5 * x)
                ineqs.append(({j: s * v for j, v in r.items()}, -s * hb.get(kx, 0)))
                if x == 0:
                    for j, v in r.items():
                        obj[j] += s * v
            if genus or irr:
                # 列ごとの像を σ 対称化（2 倍して整数に）してから G 展開・既約表示にする。合計が σ 不変なら線形性で正しい
                hbx = {k: 2 * v for k, v in hb.items()}
                cols = {}
                for kb, r in rows0[rho].items():
                    kx = (kb[0], kb[1], kb[2] + 5 * x)
                    ks = (kx[0], -kx[1], kx[2] - kx[1])
                    for j, v in r.items():
                        d = cols.setdefault(j, {})
                        d[kx] = d.get(kx, 0) + v
                        d[ks] = d.get(ks, 0) + v
                if genus:
                    Nx = g_expand(hbx)
                    grows = {}
                    for j, d in cols.items():
                        for k, v in g_expand({kk: vv for kk, vv in d.items() if vv}).items():
                            grows.setdefault(k, {})[j] = v
                    for k in set(grows) | set(Nx):
                        r = grows.get(k, {})
                        if r:
                            ineqs.append(({j: s * v for j, v in r.items()}, -s * Nx.get(k, 0)))
                if irr:
                    from rcsknots.lmovstruct import rows as w_rows
                    def irrc(bd):
                        out = {}
                        for (bb, ff), rr in w_rows(bd).items():
                            Kmax = max(abs(k) for k in rr)
                            for k in range(0, Kmax + 1):
                                v = rr.get(k, 0) - rr.get(k + 2, 0)
                                if v:
                                    out[(bb, k, ff)] = v
                        return out
                    Ix = irrc(hbx)
                    irows = {}
                    for j, d in cols.items():
                        for k, v in irrc({kk: vv for kk, vv in d.items() if vv}).items():
                            irows.setdefault(k, {})[j] = v
                    for k in set(irows) | set(Ix):
                        r = irows.get(k, {})
                        if r:
                            ineqs.append(({j: s * v for j, v in r.items()}, -s * Ix.get(k, 0)))
    log(f"  等式 {len(eqs)}、不等式 {len(ineqs)}、符号 {sorted(set(sig.values()))}")

    def mat(rows):
        A = sps.lil_matrix((len(rows), n))
        for i, (row, _) in enumerate(rows):
            for j, v in row.items():
                A[i, j] = v
        return A.tocsr()
    rhs_e = np.array([r for _, r in eqs], float); rhs_i = np.array([r for _, r in ineqs], float)
    cons = [LinearConstraint(mat(eqs), rhs_e, rhs_e), LinearConstraint(mat(ineqs), rhs_i, np.inf)]
    bnd = Bounds(np.full(n, -np.inf), np.full(n, np.inf))

    def solve(c, extra=()):
        return milp(c, constraints=cons + list(extra), bounds=bnd, integrality=np.ones(n), options={"time_limit": 1800})
    t0 = time.time()
    r0 = solve(obj)
    out = {"n_g": n, "n_eq": len(eqs), "n_ineq": len(ineqs)}
    if r0.x is None:
        out["status"] = r0.message
        log(f"  総数の最小化：{r0.message}")
        return out
    opt = round(float(obj @ r0.x))
    out["min_total"] = opt
    if known:
        log(f"  総数の最小値 {opt}（正解は 0 = 最小か {opt == 0}）({time.time() - t0:.0f}s)")
        if opt < 0:
            log("  正解より総数の小さい解がある → 総数の最小化は正解を選ばない")
    else:
        log(f"  総数の最小値（出発点との差）{opt} ({time.time() - t0:.0f}s)")
    face_cons = [LinearConstraint(obj.reshape(1, -1), -np.inf, opt + 1e-6)]
    rng = np.random.default_rng(0)
    sols = set()
    for s in range(seeds):
        c = rng.integers(-10 ** 4, 10 ** 4, size=n).astype(float)
        rr = solve(c, face_cons)
        if rr.x is not None:
            sols.add(tuple(int(round(v)) for v in rr.x))
    out["face_solutions_found"] = len(sols)
    out["face_contains_true"] = tuple([0] * n) in sols
    # 正確な列挙：整数の乱数 c で c·g の小さい順に、c·g ≥ 前の値 + 1 を足して解き直す（c·g の同値は乱数で起きにくい）
    c = rng.integers(-10 ** 6, 10 ** 6, size=n).astype(float)
    allsols, lo, nmax = [], -np.inf, NMAX_ENUM
    while len(allsols) < nmax:
        rr = solve(c, face_cons + [LinearConstraint(c.reshape(1, -1), lo, np.inf)])
        if rr.x is None:
            break
        g = tuple(int(round(v)) for v in rr.x)
        allsols.append(g)
        lo = float(c @ np.array(g)) + 1
    out["face_count"] = len(allsols) if len(allsols) < nmax else f">={nmax}"
    out["face_count_contains_true"] = tuple([0] * n) in allsols
    if len(allsols) <= 5:
        sols_x = []
        for g in allsols:
            x = dict(xp)
            for e, ge in zip(G, g):
                if ge:
                    for f, c in Phi.items():
                        m = (e[0] + f[0], e[1] + f[1], e[2] + f[2])
                        x[m] = x.get(m, 0) + ge * c
            sols_x.append({",".join(map(str, m)): v for m, v in sorted(x.items()) if v})
        out["face_solutions_rcs21"] = sols_x
    sols = set(allsols) | sols
    log(f"  面上の整数解の列挙：{out['face_count']} 個（正解を含む {out['face_count_contains_true']}）({time.time() - t0:.0f}s)")
    # 列挙した解の特徴：正解が各特徴で何位か（選択の規則の候補の検査）
    def feats(g):
        dx = {}
        for e, ge in zip(G, g):
            if ge:
                for f, c in Phi.items():
                    m = (e[0] + f[0], e[1] + f[1], e[2] + f[2])
                    dx[m] = dx.get(m, 0) + ge * c
        x = dict(xp)
        for m, v in dx.items():
            x[m] = x.get(m, 0) + v
        x = {m: v for m, v in x.items() if v}
        hf = {}
        for rho, rows in rows0.items():
            hb = to_bold(hp0[rho], shifts[rho])
            d = dict(hb)
            for kb, r in rows.items():
                v = sum(r.get(j, 0) * gj for j, gj in enumerate(g) if gj)
                if v:
                    d[kb] = d.get(kb, 0) + v
            hf[rho] = {k: v for k, v in d.items() if v}
        allh = [k for d in hf.values() for k in d]
        gmax = 0
        for rho, d in hf.items():
            try:
                N = g_expand(d)
                gmax = max(gmax, max(k[0] for k in N))
            except ValueError:
                gmax = None
        return {"hat_terms": sum(len(d) for d in hf.values()), "rcs_terms": len(x), "rcs_l1": sum(abs(v) for v in x.values()),
                "rcs_neg": sum(1 for v in x.values() if v < 0), "hat_t_width": max(k[2] for k in allh) - min(k[2] for k in allh),
                "hat_q_width": max(k[1] for k in allh) - min(k[1] for k in allh), "hat_gmax": gmax,
                "hat_max_coeff": max(abs(v) for d in hf.values() for v in d.values())}
    if isinstance(out["face_count"], int) and allsols:
        F = [feats(g) for g in allsols]
        ft = feats(tuple([0] * n))
        rank = {}
        for key in ft:
            vals = [f[key] for f in F if f[key] is not None]
            rank[key] = {"true": ft[key], "n_le": sum(1 for v in vals if v <= ft[key]), "n_lt": sum(1 for v in vals if v < ft[key]),
                         "min": min(vals) if vals else None, "max": max(vals) if vals else None}
        out["feature_rank"] = rank
        for key, r in rank.items():
            log(f"    特徴 {key:14s} 正解 {r['true']}、より小さい解 {r['n_lt']}、以下の解 {r['n_le']} / {len(allsols)}（範囲 {r['min']}〜{r['max']}）")
    # 面の有界性（LP 緩和で乱数方向の最大・最小が有限か）
    from scipy.optimize import linprog
    Ae, Ai = mat(eqs), mat(ineqs)
    Ai2 = sps.vstack([-Ai, obj.reshape(1, -1)]).tocsr(); bi2 = np.concatenate([-rhs_i, [opt + 1e-6]])
    widths = []
    for _ in range(4):
        c = rng.standard_normal(n)
        w = []
        for sgn in (1, -1):
            r = linprog(sgn * c, A_ub=Ai2, b_ub=bi2, A_eq=Ae, b_eq=rhs_e, bounds=[(None, None)] * n, method="highs")
            w.append(None if r.status == 3 else float(sgn * r.fun))
        widths.append(w)
    out["face_bounded"] = all(a is not None and b is not None for a, b in widths)
    # 別解の要約：Δx（rCS）と Δ\hat f の項数
    summ = []
    for g in list(sols)[:3]:
        if not any(g):
            continue
        dx = {}
        for e, ge in zip(G, g):
            if ge:
                for f, c in Phi.items():
                    m = (e[0] + f[0], e[1] + f[1], e[2] + f[2])
                    dx[m] = dx.get(m, 0) + ge * c
        dx = {m: v for m, v in dx.items() if v}
        dh = {}
        for rho, rows in rows0.items():
            d = {}
            for kb, r in rows.items():
                v = sum(r.get(j, 0) * gj for j, gj in enumerate(g) if gj)
                if v:
                    d[kb] = v
            dh[str(rho)] = (len(d), sum(1 for v in d.values() if v > 0), sum(1 for v in d.values() if v < 0),
                            sum(dx_ for dx_ in d.values()))
        summ.append({"d_rcs_terms": len(dx), "d_hat": dh})
    out["alt_summary"] = summ
    log(f"  面は有界 {out['face_bounded']}。別解の例：{summ[:2]}")
    # 一意性：c·g ≥ 1 または ≤ −1（正解 g = 0 からの距離）
    c = rng.integers(1, 10 ** 6, size=n).astype(float)
    st = [solve(np.zeros(n), face_cons + [LinearConstraint(c.reshape(1, -1), lo, hi)]).status for lo, hi in ((1, np.inf), (-np.inf, -1))]
    out["face_unique"] = all(x == 2 for x in st)
    log(f"  最適面の上の整数解：乱数 {seeds} 回で {len(sols)} 個（正解を含む {out['face_contains_true']}）、一意 {out['face_unique']} ({time.time() - t0:.0f}s)")
    return out


def predict_twist(p, xs, use_R=True, genus=True, irr=True, log=log):
    """ツイスト結び目 K_p の [2,1] の予言：GNSSS の長方形と X0 の非 refined の [2,1] を入力に、実行可能解 → face_core。"""
    small, H21u, window, known21 = LL.twist_inputs(p)
    LL.MARGIN = 0
    small, H21u, window, known21 = LL.twist_inputs(p)         # 窓の規則（margin 0）
    xp = particular_solution(small, H21u, window, log=log)
    if xp is None:
        return {"status": "実行可能解なし"}
    res = face_core(f"K_{p}", small, xp, H21u, xs, use_R, genus, irr, known=False, log=log)
    if known21 is not None and res.get("face_solutions_rcs21"):
        kb = {",".join(map(str, m)): v for m, v in sorted(known21.items())}
        res["matches_known"] = [sol == kb for sol in res["face_solutions_rcs21"]]
        log(f"  既知の値との一致：{res['matches_known']}")
    return res


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "predict":
    p_ = int(sys.argv[2])
    xs = [int(v) for v in sys.argv[3].split(",")] if len(sys.argv) > 3 else [0]
    OUT.mkdir(parents=True, exist_ok=True)
    log(f"=== lmov21_h11 predict {time.strftime('%Y-%m-%d %H:%M:%S')} {' '.join(sys.argv[2:])} ===")
    res = predict_twist(p_, xs, "--noR" not in sys.argv, "--nogenus" not in sys.argv, "--noirr" not in sys.argv)
    path = OUT / "lmov21_h11_predict.json"
    prev = json.loads(path.read_text()) if path.exists() else {}
    prev[f"K_{p_} x={sys.argv[3] if len(sys.argv) > 3 else '0'}"] = json.loads(json.dumps(res, default=str))
    path.write_text(json.dumps(prev, ensure_ascii=False, indent=1))

if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "face":
    K = sys.argv[2]
    xs = [int(v) for v in sys.argv[3].split(",")] if len(sys.argv) > 3 else [0]
    use_R = "--R" in sys.argv
    genus = "--genus" in sys.argv
    irr = "--irr" in sys.argv
    OUT.mkdir(parents=True, exist_ok=True)
    log(f"=== lmov21_h11 face {time.strftime('%Y-%m-%d %H:%M:%S')} {' '.join(sys.argv[2:])} ===")
    res = face(K, xs, use_R, genus, irr)
    path = OUT / "lmov21_h11_face.json"
    prev = json.loads(path.read_text()) if path.exists() else {}
    prev[f"{K} x={sys.argv[3] if len(sys.argv) > 3 else '0'} R={use_R} genus={genus} irr={irr}"] = json.loads(json.dumps(res, default=str))
    path.write_text(json.dumps(prev, ensure_ascii=False, indent=1))
