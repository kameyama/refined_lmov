"""[2,1] の構成の線形解析：各仮定（一次の条件）が rCS_{[2,1]} の自由度をどれだけ削るかを、正値性なしで数える。

未知数：rCS_{[2,1]} の太字の係数 x_m（窓は build21 と同じ）。条件は scripts/dgr_box3.py の build21 の
  U：𝐭 = −1 で非 refined の [2,1]（t = q と同じ特殊化）、T：転置で不変（W）、E：EGP（q = 1 で r₁r₁₁）、U2：a = t² で r₁r₁₁、H：hook（a = t/q で c²A(𝐪³𝐭, 𝐭)）
と、本スクリプトの
  L：\\hat f_{[3]}, \\hat f_{[2,1]}, \\hat f_{[1,1,1]} の Laurent 性。分子 N_ρ(x) が分母 D_ρ の各既約因子で割り切れることを、因子ごとの代入（単項式の径数づけ、
     円分多項式は剰余）で一次式にする。重複度 2 の (Q ∓ 1)², (T ∓ 1)² は微分の条件も課す。
零空間の次元は法 p = 2³¹ − 1 の疎な消去で数える（有理数の次元と一致する確率が高い）。既知の結び目では正解が全ての式を満たすことを整数で確かめる。
使い方：PYTHONPATH=src python scripts/lmov21_linear.py T2_3 4_1 T2_5 T3_4 K1 tw:2 tw:-2 ...
出力：data/lmovstruct/lmov21_linear.json、logs/lmov21_linear.log。
"""
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dgr_box3 import KEYS, build21, formula3, inputs_R  # noqa: E402
from t4_colored import bold_to_small  # noqa: E402
from v4_closed import L  # noqa: E402
from rcsknots.bps import to_bold  # noqa: E402

P = 2 ** 31 - 1
OUT = Path("data/lmovstruct")


# ---- Laurent 性の一次式 ---------------------------------------------------------------------------

def _images(e):
    """小文字の単項式 (ea, eq, et) の、各因子での像 {(因子名, key): 係数}。"""
    ea, eq, et = e
    out = {}
    # (Q − 1)², (Q + 1)²：値と微分
    out[("Q=1", (ea, et))] = 1
    out[("dQ=1", (ea, et))] = eq
    out[("Q=-1", (ea, et))] = (-1) ** (eq % 2)
    out[("dQ=-1", (ea, et))] = eq * (-1) ** ((eq - 1) % 2)
    out[("T=1", (ea, eq))] = 1
    out[("dT=1", (ea, eq))] = et
    out[("T=-1", (ea, eq))] = (-1) ** (et % 2)
    out[("dT=-1", (ea, eq))] = et * (-1) ** ((et - 1) % 2)
    # Q² + 1：Q^e = i^e → (成分 e mod 2, 符号)
    s = (1, 1, -1, -1)[eq % 4]
    out[("Q^2+1", (ea, et, eq % 2))] = s
    s = (1, 1, -1, -1)[et % 4]
    out[("T^2+1", (ea, eq, et % 2))] = s
    # QT ∓ 1：T = ±1/Q
    out[("QT-1", (ea, eq - et))] = 1
    out[("QT+1", (ea, eq - et))] = (-1) ** (et % 2)
    # QT² ∓ 1：Q = ±T^{−2}
    out[("QT2-1", (ea, et - 2 * eq))] = 1
    out[("QT2+1", (ea, et - 2 * eq))] = (-1) ** (eq % 2)
    # Q²T ∓ 1：T = ±Q^{−2}
    out[("Q2T-1", (ea, eq - 2 * et))] = 1
    out[("Q2T+1", (ea, eq - 2 * et))] = (-1) ** (et % 2)
    # Q² + Q + 1：Q³ = 1、Q² = −1 − Q。成分 (1, Q)
    c3 = {0: ((1, 0),), 1: ((0, 1),), 2: ((-1, -1),)}
    for comp, v in enumerate(c3[eq % 3][0]):
        if v:
            out[("Q^2+Q+1", (ea, et, comp))] = v
    for comp, v in enumerate(c3[et % 3][0]):
        if v:
            out[("T^2+T+1", (ea, eq, comp))] = v
    # Q² − Q + 1：Q³ = −1、Q² = Q − 1。Q^0..Q^5 = 1, Q, Q−1, −1, −Q, 1−Q
    c6 = {0: (1, 0), 1: (0, 1), 2: (-1, 1), 3: (-1, 0), 4: (0, -1), 5: (1, -1)}
    for comp, v in enumerate(c6[eq % 6]):
        if v:
            out[("Q^2-Q+1", (ea, et, comp))] = v
    for comp, v in enumerate(c6[et % 6]):
        if v:
            out[("T^2-T+1", (ea, eq, comp))] = v
    return out


FACTORS = {"(3,)": None, "(2, 1)": {"Q^2+1", "T^2+1"}, "(1, 1, 1)": None}   # ρ = [2,1] の分母には Q²+1, T²+1 が無い


def laurent_rows(small, mons, lin_col):
    """{ρ: [(row {列: 係数}, rhs)]}。"""
    R = inputs_R({**small, (2, 1): {(0, 0, 0): 0}})
    rows_all = []
    for rho, (ns, den) in formula3().items():
        skip = FACTORS[str(rho)] or set()
        const = None
        for k, n in ns.items():
            if KEYS[k] == "r21":
                continue
            term = n * R[KEYS[k]]
            const = term if const is None else const + term
        const = const.to_dict()
        eqs = defaultdict(dict)          # (因子, key) → {列: 係数}
        rhs = defaultdict(int)
        for m in mons:
            poly = (ns[1] * L.from_dict(bold_to_small({m: 1}))).to_dict()
            col = lin_col[("x", m)]
            for e, c in poly.items():
                for key, v in _images(e).items():
                    if key[0] in skip:
                        continue
                    eqs[key][col] = eqs[key].get(col, 0) + c * v
        for e, c in const.items():
            for key, v in _images(e).items():
                if key[0] in skip:
                    continue
                rhs[key] -= c * v
                eqs.setdefault(key, {})
        for key, row in eqs.items():
            row = {c: v for c, v in row.items() if v}
            if row or rhs[key]:
                rows_all.append((row, rhs[key]))
    return rows_all


# ---- 法 p の疎な消去 -------------------------------------------------------------------------------

def rank_nullity(rows, ncols, p=P):
    """rows：[(row {列: 係数}, rhs)]。(rank, feasible)。"""
    pivots = []                                   # (列, 行) 挿入順
    infeasible = False
    for row, b in rows:
        r = {c: v % p for c, v in row.items() if v % p}
        rb = b % p
        for pc, prow, pb in pivots:
            f = r.get(pc)
            if f:
                for c, v in prow.items():
                    nv = (r.get(c, 0) - f * v) % p
                    if nv:
                        r[c] = nv
                    else:
                        r.pop(c, None)
                rb = (rb - f * pb) % p
        if not r:
            if rb:
                infeasible = True
            continue
        pc = min(r)
        inv = pow(r[pc], p - 2, p)
        r = {c: v * inv % p for c, v in r.items()}
        rb = rb * inv % p
        pivots.append((pc, r, rb))
    return len(pivots), not infeasible


def nullspace_basis(rows, ncols, p=P):
    """零空間（斉次）の基底（法 p）。rank_nullity と同じ消去で、自由な列ごとに後退代入。"""
    pivots = []
    for row, _ in rows:
        r = {c: v % p for c, v in row.items() if v % p}
        for pc, prow in pivots:
            f = r.get(pc)
            if f:
                for c, v in prow.items():
                    nv = (r.get(c, 0) - f * v) % p
                    if nv:
                        r[c] = nv
                    else:
                        r.pop(c, None)
        if not r:
            continue
        pc = min(r)
        inv = pow(r[pc], p - 2, p)
        pivots.append((pc, {c: v * inv % p for c, v in r.items()}))
    # 完全簡約：後ろの pivot で前の行を掃除
    for i in range(len(pivots) - 1, -1, -1):
        pc, prow = pivots[i]
        for j in range(i):
            qc, qrow = pivots[j]
            f = qrow.get(pc)
            if f:
                for c, v in prow.items():
                    nv = (qrow.get(c, 0) - f * v) % p
                    if nv:
                        qrow[c] = nv
                    else:
                        qrow.pop(c, None)
    pivcols = {pc for pc, _ in pivots}
    free = [c for c in range(ncols) if c not in pivcols]
    basis = []
    for fc in free:
        vec = {fc: 1}
        for pc, prow in pivots:
            v = prow.get(fc)
            if v:
                vec[pc] = (-v) % p
        basis.append(vec)
    return basis


def solve_unique(rows, ncols, p=P):
    """零空間が 0 のとき、法 p の完全簡約で解を求め、対称な剰余で整数に戻し、全ての式を整数で検算する。失敗なら None。"""
    pivots = []
    for row, b in rows:
        r = {c: v % p for c, v in row.items() if v % p}
        rb = b % p
        for pc, prow, pb in pivots:
            f = r.get(pc)
            if f:
                for c, v in prow.items():
                    nv = (r.get(c, 0) - f * v) % p
                    if nv:
                        r[c] = nv
                    else:
                        r.pop(c, None)
                rb = (rb - f * pb) % p
        if not r:
            continue
        pc = min(r)
        inv = pow(r[pc], p - 2, p)
        pivots.append((pc, {c: v * inv % p for c, v in r.items()}, rb * inv % p))
    if len(pivots) < ncols:
        return None
    for i in range(len(pivots) - 1, -1, -1):
        pc, prow, pb = pivots[i]
        for j in range(i):
            qc, qrow, qb = pivots[j]
            f = qrow.get(pc)
            if f:
                for c, v in prow.items():
                    nv = (qrow.get(c, 0) - f * v) % p
                    if nv:
                        qrow[c] = nv
                    else:
                        qrow.pop(c, None)
                pivots[j] = (qc, qrow, (qb - f * pb) % p)
    x = {}
    for pc, prow, pb in pivots:
        assert len(prow) == 1
        x[pc] = pb if pb <= p // 2 else pb - p
    for row, b in rows:
        if sum(v * x.get(c, 0) for c, v in row.items()) != b:
            return None
    return x


# ---- 入力 ------------------------------------------------------------------------------------------

def t_window_rule(small, margin):
    """𝐭 の窓の規則：rCS_{[2,1]} の 𝐭 の範囲は rCS_{[2]}·rCS_{[1]} の範囲に収まる（既知の 5 結び目で等しい）。±margin。"""
    from dgr_box3 import _mul
    b = {l: {k: v for k, v in to_bold(P, 1 if any(k[0] % 2 for k in P) else 0).items() if v} for l, P in small.items()}
    prod = _mul(b[(2,)], b[(1,)])
    ts = [m[2] for m in prod]
    return (min(ts) - margin, max(ts) + 1 + margin)


MARGIN = None      # None：従来の窓。整数：𝐭 の窓の規則 ±MARGIN


def known_inputs(K):
    from framing_check import load_reduced
    R = load_reduced(K, 3)
    r21 = R[(2, 1)]
    sh = 1 if any(k[0] % 2 for k in r21) else 0
    true21 = {k: v for k, v in to_bold(r21, sh).items() if v}
    H21u = defaultdict(int)
    for (a, q, t), c in true21.items():
        H21u[(a, q)] += c * (-1) ** (t % 2)
    H21u = {k: v for k, v in H21u.items() if v}
    small = {l: Pp for l, Pp in R.items() if l != (2, 1)}
    ts = [k[2] for k in true21]
    window = (min(ts) - 3, max(ts) + 4) if MARGIN is None else t_window_rule(small, MARGIN)
    return small, H21u, window, true21


def twist_inputs(p):
    from dgr_box3_predict import inputs3, x0, known21
    bold, small = inputs3(p)
    H21u = x0(p, (2, 1))
    ts = [e[2] for e in bold[(1,)]]
    window = (min(6 * min(ts), 0) - 4, max(6 * max(ts), 0) + 6) if MARGIN is None else t_window_rule(small, MARGIN)
    return small, H21u, window, known21(p)


def analyze(name, small, H21u, window, true21, log, scale=1, full=True):
    t0 = time.time()
    if scale != 1:
        c = (window[0] + window[1]) // 2
        half = (window[1] - window[0]) * scale // 2
        window = (c - half, c + half + 1)
    res = {"knot": name, "window_t": window, "scale": scale}
    steps = [("U",), ("U", "T"), ("U", "T", "E"), ("U", "T", "E", "U2"), ("U", "T", "E", "U2", "H"),
             ("U", "T", "E", "U2", "H", "L"), ("U", "T", "L"), ("U", "T", "H"), ("U", "T", "H", "L"), ("U", "T", "E", "U2", "L")]
    if not full:
        steps = [("U", "T", "E", "U2", "L"), ("U", "T", "E", "U2", "H", "L")]
    lin0, mons = build21(small, H21u, window, (1, 1, 1), conds=("U",))
    ncols = len(mons)
    res["n_unknowns"] = ncols
    lrows = None
    for conds in steps:
        lin, _ = build21(small, H21u, window, (1, 1, 1), conds=tuple(c for c in conds if c != "L"))
        rows = [(dict(r), b) for r, b in zip(lin.rows, lin.rhs)]
        if "L" in conds:
            if lrows is None:
                lrows = laurent_rows(small, mons, lin.col)
            rows += lrows
        # 既知の結び目：正解が式を満たすか（整数）
        ok = None
        if true21 is not None:
            xv = {lin.col[("x", m)]: true21.get(m, 0) for m in mons}
            ok = all((sum(v * xv.get(c, 0) for c, v in row.items()) - b) == 0 for row, b in rows)
            missing = [m for m in true21 if m not in set(mons)]
            if missing:
                ok = f"窓の外に正解の項 {len(missing)} 個"
        rank, feas = rank_nullity(rows, ncols)
        null = ncols - rank
        key = "+".join(conds)
        res[key] = {"rows": len(rows), "rank": rank, "nullity": null, "feasible": feas, "true_satisfies": ok}
        log(f"  {name:6s} {key:16s} 式 {len(rows):6d} 階数 {rank:5d} 零空間 {null:4d} 実行可能 {feas} 正解 {ok} ({time.time() - t0:.0f}s)")
    # 全条件での零空間の基底の台（太字の単項式）
    lin, _ = build21(small, H21u, window, (1, 1, 1), conds=("U", "T", "E", "U2", "H"))
    rows = [(dict(r), b) for r, b in zip(lin.rows, lin.rhs)] + lrows
    basis = nullspace_basis(rows, ncols)
    inv = {v: k for k, v in lin.col.items()}
    res["kernel_dim"] = len(basis)
    res["kernel_support"] = [sorted(inv[c][1] for c in vec) for vec in basis[:5]]
    if not basis:
        x = solve_unique(rows, ncols)
        if x is not None:
            sol = {",".join(map(str, inv[c][1])): v for c, v in x.items() if v}
            res["solution_bold"] = sol
            res["solution_terms"] = len(sol)
            res["solution_matches_known"] = (true21 is not None) and ({inv[c][1]: v for c, v in x.items() if v} == true21)
            log(f"  {name:6s} 一意解 {len(sol)} 項、既知と一致 {res['solution_matches_known']}")
        else:
            res["solution_bold"] = None
            log(f"  {name:6s} 一意解の整数での検算に失敗")
    return res


def main():
    global MARGIN
    OUT.mkdir(parents=True, exist_ok=True)
    logf = open("logs/lmov21_linear.log", "a")

    def log(s):
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    log(f"=== lmov21_linear {time.strftime('%Y-%m-%d %H:%M:%S')} ===")
    out = {}
    scale = 1
    full = True
    args = []
    for arg in sys.argv[1:]:
        if arg.startswith("--scale="):
            scale = int(arg.split("=")[1])
        elif arg.startswith("--margin="):
            MARGIN = int(arg.split("=")[1])
        elif arg == "--brief":
            full = False
        else:
            args.append(arg)
    tag = (f"_s{scale}" if scale != 1 else "") + (f"_m{MARGIN}" if MARGIN is not None else "")
    for arg in args:
        if arg.startswith("tw:"):
            p = int(arg[3:])
            small, H21u, window, true21 = twist_inputs(p)
            name = f"K_{p}"
        else:
            small, H21u, window, true21 = known_inputs(arg)
            name = arg
        out[name] = analyze(name, small, H21u, window, true21, log, scale, full)
        path = OUT / f"lmov21_linear{tag}.json"
        prev = json.loads(path.read_text()) if path.exists() else {}
        prev.update(out)
        path.write_text(json.dumps(prev, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
