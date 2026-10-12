"""ツイスト結び目 K_p の refined [2,1]・[3,1]（[2,1,1]）を、一次の条件だけで求める（D12、4 箱の Laurent 性による構成）。

入力：refined の長方形（[1], [2], [1,1], [3], [1,1,1] は GNSSS、[4], [2,2], [1⁴] は KNTZ の式 x1_check.triply）、
      非 refined の [2,1], [3,1]（3 本・4 本ケーブルの逆解き data/x0）。
未知数：rCS_{[2,1]}, rCS_{[3,1]} の太字の係数。窓：(𝐚, 𝐪) は非 refined の台の箱、𝐭 は行ごとの両端の線（予想 R：[λ_i] の行の min-plus・max-plus）。
条件（全て一次）：
  [2,1]：t = q、転置、EGP、U(2)、hook、3 箱の \\hat f の Laurent 性、R。
  [3,1]：t = q、EGP（q = 1：r₁₁r₁r₁、t = 1：r₃r₁）、U(2)（r₂r₁₁）、hook（c³·A(𝐪⁴𝐭², 𝐭)）、R。[2,1,1] は転置。
  4 箱：\\hat f_ρ（ρ ⊢ 4）の Laurent 性（分子の D_ρ による剰余 = 0）。
解き方：[2,1], [3,1] それぞれの一次の条件の解空間 x = x₀ + B y（有理数の特解、整数の基底）を求め、4 箱の Laurent 性を y の一次方程式にして
      法 p の密な消去で解く。零空間が 0 なら一意。x を 2 つの素数の CRT と対称な剰余で整数に戻し、全ての条件を整数で検算する。
使い方：PYTHONPATH=src:scripts python scripts/lmov4_predict.py p   （p = ±1 は既知との照合、±2 は予言）
記録：logs/lmov4_predict.log、data/lmovstruct/lmov4_predict_K_<p>.json。
"""
import json
import sys
import time
from fractions import Fraction
from math import gcd
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lmov21_linear as LL  # noqa: E402
from dgr_box3 import _mul, build21, hat3, rcs_window21  # noqa: E402
from dgr_box3_predict import inputs3, x0 as unref  # noqa: E402
from dgr_d7_construct import Lin  # noqa: E402
from lmov21_polytope import exact_affine  # noqa: E402
from lmov4_formula import formula_n, hat_n, inputs_n  # noqa: E402
from lmov4_joint import remainder_rows, transpose_bold  # noqa: E402
from t4_colored import bold_to_small  # noqa: E402
from v4_closed import L  # noqa: E402
from rcsknots.bps import g_expand, to_bold  # noqa: E402
from rcsknots.dgr import specialize  # noqa: E402

LOG = Path("logs/lmov4_predict.log")
PRIMES = (1048573, 1048571)


def log(s):
    print(s, flush=True)
    with LOG.open("a") as f:
        f.write(s + "\n")


def bold(P):
    return {k: v for k, v in to_bold(P, 1 if any(k[0] % 2 for k in P) else 0).items() if v}


def ends(b):
    o = {}
    for a, _, t in b:
        lo, hi = o.get(a, (10 ** 9, -10 ** 9))
        o[a] = (min(lo, t), max(hi, t))
    return o


def minplus(e1, e2):
    out = {}
    for a1, (l1, h1) in e1.items():
        for a2, (l2, h2) in e2.items():
            lo, hi = out.get(a1 + a2, (10 ** 9, -10 ** 9))
            out[a1 + a2] = (min(lo, l1 + l2), max(hi, h1 + h2))
    return out


def twist_inputs(p):
    """小文字の reduced（|λ| ≤ 4 の長方形）と、非 refined の [2,1], [3,1], [2,1,1]（太字 {(e_𝐚, e_𝐪): c}）。"""
    from x1_check import bold_dict, triply
    boldd, small = inputs3(p)
    for lam in ((4,), (2, 2), (1, 1, 1, 1)):
        B = {k: v for k, v in bold_dict(triply(lam[0], len(lam), p, 1)).items() if v}
        boldd[lam] = B
        small[lam] = bold_to_small(B)
    H = {lam: unref(p, lam) for lam in ((2, 1), (3, 1), (2, 1, 1))}
    return small, boldd, H


def affine21(small, H21u, bold_in):
    win = rcs_window21(bold_in[(1,)], bold_in[(2,)])
    ts = [v for lo, hi in win.values() for v in (lo, hi)]
    window = (min(ts), max(ts) + 1)
    lin, mons = build21(small, H21u, window, (1, 1, 1), conds=("U", "T", "E", "U2", "H", "R"))
    rows = [(dict(r), b) for r, b in zip(lin.rows, lin.rhs)] + LL.laurent_rows(small, mons, lin.col)
    x0, basis = exact_affine(rows, len(mons))
    inv = {v: k for k, v in lin.col.items()}
    return ({inv[c][1]: f for c, f in x0.items()}, [{inv[c][1]: z for c, z in vec.items()} for vec in basis], mons)


def affine31(small, H31u, bold_in):
    b = bold_in
    win = minplus(ends(b[(3,)]), ends(b[(1,)]))
    As = sorted({a for a, _ in H31u}); Qs = sorted({q for _, q in H31u})
    tlo = min(lo for lo, _ in win.values()); thi = max(hi for _, hi in win.values())
    mons = [(a, q, t) for a in range(As[0], As[-1] + 1, 2) for q in range(Qs[0], Qs[-1] + 1, 2) for t in range(tlo, thi + 1)]
    lin = Lin()
    for m in mons:
        lin.var(("x", m))
    for m in mons:                                            # R：行ごとの両端の線の外は 0
        if m[0] not in win or not win[m[0]][0] <= m[2] <= win[m[0]][1]:
            lin.eq({("x", m): 1}, 0)
    for a in range(As[0], As[-1] + 1, 2):                     # U
        for q in range(Qs[0], Qs[-1] + 1, 2):
            lin.eq({("x", (a, q, t)): (-1) ** (t % 2) for t in range(tlo, thi + 1)}, H31u.get((a, q), 0))

    def special(n, m, target):
        img = {}
        for mm in mons:
            for k, v in specialize({mm: 1}, n, m).items():
                img.setdefault(k, {})[("x", mm)] = v
        for k in set(img) | set(target):
            lin.eq(img.get(k, {}), target.get(k, 0))
    prod = _mul(_mul(b[(1, 1)], b[(1,)]), b[(1,)])           # E（q = 1）
    target = {}
    for (a, q, t), v in prod.items():
        target[(a, t - q)] = target.get((a, t - q), 0) + v
    img = {}
    for m in mons:
        img.setdefault((m[0], m[2] - m[1]), {})[("x", m)] = 1
    for k in set(img) | set(target):
        lin.eq(img.get(k, {}), target.get(k, 0))
    prod = _mul(b[(3,)], b[(1,)])                             # E'（t = 1 ⇔ 𝐪 = ±1）
    for sgn in (1, -1):
        target = {}
        for (a, q, t), v in prod.items():
            target[(a, t)] = target.get((a, t), 0) + v * sgn ** (q % 2)
        img = {}
        for m in mons:
            img.setdefault((m[0], m[2]), {})[("x", m)] = sgn ** (m[1] % 2)
        for k in set(img) | set(target):
            lin.eq(img.get(k, {}), target.get(k, 0))
    special(2, 0, _mul(specialize(b[(2,)], 2, 0), specialize(b[(1, 1)], 2, 0)))      # U2
    (cq, ct), cs = next(iter(specialize(b[(1,)], 1, 0).items()))                     # hook
    A = specialize(b[(1,)], 1, 1)
    special(1, 1, {(4 * eq, et + 2 * eq + 3 * ct): cs ** 3 * v for (eq, et), v in A.items()})
    rows = [(dict(r), bb) for r, bb in zip(lin.rows, lin.rhs)]
    x0, basis = exact_affine(rows, len(mons))
    inv = {v: k for k, v in lin.col.items()}
    return ({inv[c][1]: f for c, f in x0.items()}, [{inv[c][1]: z for c, z in vec.items()} for vec in basis], mons)


def scale_int(x0):
    d = 1
    for f in x0.values():
        d = d * f.denominator // gcd(d, f.denominator)
    return d, {m: int(f * d) for m, f in x0.items()}


def new_state(ncols):
    """法 p の RREF の状態：ピボット列（挿入順）、RREF の行（最後の列は右辺）、矛盾の有無。"""
    return {"piv": [], "B": np.zeros((0, ncols + 1)), "bad": False}


def eliminate(st, rows, rhs, p, batch=1500):
    """A y = rhs（法 p）の行を状態 st の RREF に取り込む（密な消去、まとまりごと）。rows：[{列: 係数}]。"""
    ncols = st["B"].shape[1] - 1
    piv_cols, B = st["piv"], st["B"]
    for start in range(0, len(rows), batch):
        chunk = rows[start:start + batch]
        M = np.zeros((len(chunk), ncols + 1))
        for i, row in enumerate(chunk):
            for c, v in row.items():
                M[i, c] = v % p
            M[i, ncols] = rhs[start + i] % p
        if len(piv_cols):
            M = (M - (M[:, piv_cols] @ B) % p) % p
        for i in range(len(chunk)):
            nz = np.nonzero(M[i, :ncols])[0]
            if len(nz) == 0:
                if M[i, ncols] % p:
                    st["bad"] = True
                continue
            c = int(nz[0])
            inv = pow(int(M[i, c]), p - 2, p)
            M[i] = (M[i] * inv) % p
            col = M[:, c].copy(); col[i] = 0
            M = (M - np.outer(col, M[i])) % p
            if len(piv_cols):
                B = (B - np.outer(B[:, c], M[i])) % p
            B = np.vstack([B, M[i:i + 1]])
            piv_cols.append(c)
    st["B"] = B
    return st


def solution_of(st):
    """RREF の状態から (特解 y（自由な変数は 0）, 零空間の次元, 矛盾か)。"""
    ncols = st["B"].shape[1] - 1
    if st["bad"]:
        return None, None, True
    y = np.zeros(ncols, dtype=object)
    for r, pc in enumerate(st["piv"]):
        y[pc] = int(st["B"][r, ncols])
    return y, ncols - len(st["piv"]), False


def nullspace_of(st, p):
    """RREF の状態から零空間の基底（{列: 係数} 法 p、自由な列ごとに 1 本）。"""
    ncols = st["B"].shape[1] - 1
    pivset = set(st["piv"])
    basis = []
    for fc in range(ncols):
        if fc in pivset:
            continue
        vec = {fc: 1}
        for r, pc in enumerate(st["piv"]):
            v = int(st["B"][r, fc])
            if v:
                vec[pc] = (-v) % p
        basis.append(vec)
    return basis


def solve_modp(rows, ncols, rhs, p, batch=1500):
    """A y = rhs（法 p）を密な消去で解く。返り値 (解 y または None, 零空間の次元, 矛盾か)。"""
    return solution_of(eliminate(new_state(ncols), rows, rhs, p, batch))


def rref_path(p):
    return Path(f"data/lmovstruct/lmov4_rref_K_{p}.pkl")


def rref_rho(p, rho):
    """ρ の剰余の行（predict_rows_K_<p>_<ρ>.pkl）を 2 つの素数の RREF に取り込み、状態を保存する。
    行の辞書は巨大（1 ρ で数 GB）なので、ρ ごとに別プロセスで呼び、終わればメモリを OS に返す。"""
    import pickle
    S = pickle.loads(Path(f"data/lmovstruct/lmov4_sys_K_{p}.pkl").read_bytes())
    ncols = len(S["B21"]) + len(S["B31"])
    del S
    st = pickle.loads(rref_path(p).read_bytes()) if rref_path(p).exists() else {"done": [], "states": {pp: new_state(ncols) for pp in PRIMES}}
    if rho in st["done"]:
        return st
    t0 = time.time()
    rows, rhs = pickle.loads((Path("data/lmovstruct") / f"predict_rows_K_{p}_{''.join(map(str, rho))}.pkl").read_bytes())
    for pp in PRIMES:
        eliminate(st["states"][pp], rows, rhs, pp)
    st["done"].append(rho)
    rref_path(p).write_bytes(pickle.dumps(st))
    log(f"  RREF に取り込み ρ={rho}：式 {len(rows)}、ピボット " + ", ".join(f"{len(st['states'][pp]['piv'])}" for pp in PRIMES) + f" ({time.time() - t0:.0f}s)")
    return st


def rref_stage(p):
    """全ての ρ を（別プロセスで 1 つずつ）RREF に取り込み、状態を返す。保存済みの ρ は飛ばす。"""
    import pickle
    import subprocess
    from rcsknots.partitions import partitions
    for rho in partitions(4):
        if rref_path(p).exists() and rho in pickle.loads(rref_path(p).read_bytes())["done"]:
            continue
        r = subprocess.run([sys.executable, __file__, str(p), "rref", "".join(map(str, rho))], cwd=Path.cwd())
        if r.returncode:
            raise RuntimeError(f"rref ρ={rho} が失敗（{r.returncode}）")
    st = pickle.loads(rref_path(p).read_bytes())
    for pp in PRIMES:
        _, nullity, bad = solution_of(st["states"][pp])
        log(f"  法 {pp}：矛盾 {bad}、零空間 {nullity}")
    return st


def crt_sym(r1, r2):
    p1, p2 = PRIMES
    inv = pow(p1, -1, p2)
    x = (r1 + p1 * (((r2 - r1) * inv) % p2)) % (p1 * p2)
    return x - p1 * p2 if x > p1 * p2 // 2 else x


def run(p):
    t0 = time.time()
    small, boldd, H = twist_inputs(p)
    log(f"K_{p}: 入力 長方形 {sorted(small)}、非 refined [2,1] {len(H[(2, 1)])} 項、[3,1] {len(H[(3, 1)])} 項")
    # 非 refined の [2,1,1] は [3,1] の転置と整合するか
    tr = {(a, -q): c for (a, q), c in H[(3, 1)].items()}
    log(f"  非 refined の [2,1,1] = [3,1] の転置：{tr == H[(2, 1, 1)]}")
    x021, B21, mons21 = affine21(small, H[(2, 1)], boldd)
    log(f"  [2,1]：未知数 {len(mons21)}、零空間 {len(B21)} ({time.time() - t0:.0f}s)")
    x031, B31, mons31 = affine31(small, H[(3, 1)], boldd)
    log(f"  [3,1]：未知数 {len(mons31)}、零空間 {len(B31)} ({time.time() - t0:.0f}s)")
    d21, x021i = scale_int(x021); d31, x031i = scale_int(x031)
    dd = d21 * d31 // gcd(d21, d31)
    x021i = {m: v * (dd // d21) for m, v in x021i.items()}; x031i = {m: v * (dd // d31) for m, v in x031i.items()}
    fo, keys = formula_n(4)
    R = inputs_n({**small, (2, 1): {(0, 0, 0): 0}, (3, 1): {(0, 0, 0): 0}, (2, 1, 1): {(0, 0, 0): 0}}, 4)
    idx = {k: i for i, k in enumerate(keys)}
    e21 = tuple(1 if i in (idx[((2, 1), 1)], idx[((1,), 1)]) else 0 for i in range(len(keys)))
    e31 = tuple(1 if i == idx[((3, 1), 1)] else 0 for i in range(len(keys)))
    e211 = tuple(1 if i == idx[((2, 1, 1), 1)] else 0 for i in range(len(keys)))
    ncols = len(B21) + len(B31)
    import pickle
    cache_dir = Path("data/lmovstruct"); cache_dir.mkdir(parents=True, exist_ok=True)
    Path(f"data/lmovstruct/lmov4_sys_K_{p}.pkl").write_bytes(pickle.dumps(
        {"mons21": mons21, "mons31": mons31, "B21": B21, "B31": B31, "x021i": x021i, "x031i": x031i, "dd": dd, "rows": None, "rhs": None,
         "small": small, "H": H, "boldd": boldd}))
    for rho, (ns, D) in fo.items():
        cache = cache_dir / f"predict_rows_K_{p}_{''.join(map(str, rho))}.pkl"      # 再開用：ρ ごとの剰余の行
        if cache.exists():
            log(f"  ρ={rho}: 剰余の式（保存済み）")
            continue
        # 定数項：既知の単項式 ＋ x₀ の寄与（dd 倍して整数）
        const = None
        for expo, n in ns.items():
            if expo in (e21, e31, e211):
                continue
            term = n
            for (lam, dgr), k in zip(keys, expo):
                for _ in range(k):
                    term = term * R[(lam, dgr)]
            term = term * L.from_dict({(0, 0, 0): dd})
            const = term if const is None else const + term
        const = const + ns[e21] * R[((1,), 1)] * L.from_dict(bold_to_small(x021i)) \
            + ns[e31] * L.from_dict(bold_to_small(x031i)) + ns[e211] * L.from_dict(bold_to_small(transpose_bold(x031i)))
        polys = [const.to_dict()]
        base21 = ns[e21] * R[((1,), 1)]
        for v in B21:
            polys.append((base21 * L.from_dict(bold_to_small(v)) * L.from_dict({(0, 0, 0): dd})).to_dict())
        for v in B31:
            Pv = ns[e31] * L.from_dict(bold_to_small(v)) + ns[e211] * L.from_dict(bold_to_small(transpose_bold(v)))
            polys.append((Pv * L.from_dict({(0, 0, 0): dd})).to_dict())
        rr = remainder_rows(polys, D)
        del polys
        rr_rows, rr_rhs = [], []
        for key, r in rr.items():
            c0 = r.pop(0, 0)
            r2 = {j - 1: v for j, v in r.items()}
            if r2 or c0:
                rr_rows.append(r2)
                rr_rhs.append(-c0)
        del rr
        cache.write_bytes(pickle.dumps((rr_rows, rr_rhs)))
        del rr_rows, rr_rhs
        log(f"  ρ={rho}: 剰余の式 ({time.time() - t0:.0f}s)")
    log(f"  消去（ρ ごとに別プロセスで RREF へ）：列 {ncols}")
    st = rref_stage(p)
    sols = []
    for pp in PRIMES:
        y, nullity, bad = solution_of(st["states"][pp])
        if bad or nullity:
            return {"p": p, "inconsistent": bad, "nullity": nullity}
        sols.append(y)
    # x = x₀ + B y を 2 つの素数で復元（dd で割る）
    def recon(x0i, Bv, mons, off):
        out = {}
        for m in mons:
            vals = []
            for yv, pp in zip(sols, PRIMES):
                v = x0i.get(m, 0)
                for j, vec in enumerate(Bv):
                    z = vec.get(m, 0)
                    if z:
                        v += z * int(yv[off + j]) * dd
                vals.append(v % pp)
            vals = [v % pp for v, pp in zip(vals, PRIMES)]
            xi = crt_sym(vals[0], vals[1])
            if xi % dd:
                raise ValueError(f"整数に戻らない {m}")
            if xi:
                out[m] = xi // dd
        return out
    # 注意：x₀ の寄与は dd 倍済み、基底の係数 y は dd 倍の系の解なので x = (x₀·dd + dd·B y)/dd
    x21 = recon(x021i, B21, mons21, 0)
    x31 = recon(x031i, B31, mons31, len(B21))
    log(f"  解：[2,1] {len(x21)} 項（負 {sum(1 for v in x21.values() if v < 0)}）、[3,1] {len(x31)} 項（負 {sum(1 for v in x31.values() if v < 0)}）({time.time() - t0:.0f}s)")
    # 検算：3 箱・4 箱の \hat f が Laurent（公式の割り算が割り切れる）か、正値性、t = q
    full = {**small, (2, 1): bold_to_small(x21), (3, 1): bold_to_small(x31), (2, 1, 1): bold_to_small(transpose_bold(x31))}
    checks = {}
    try:
        h3 = hat3(full); h4 = hat_n(full, 4)
        for name, hs in (("3box", h3), ("4box", h4)):
            for rho, h in hs.items():
                hb = bold(h)
                pos, neg = sum(1 for v in hb.values() if v > 0), sum(1 for v in hb.values() if v < 0)
                try:
                    N = g_expand({k: (1 if pos >= neg else -1) * v for k, v in hb.items()})
                    gneg = sum(1 for v in N.values() if v < 0)
                except ValueError:
                    gneg = None
                checks[f"{name} {rho}"] = {"terms": len(hb), "spin_neg": min(pos, neg), "genus_neg": gneg}
        laurent = True
    except Exception as e:  # noqa: BLE001
        laurent = f"割り算に失敗: {e}"
    u21 = {}
    for (a, q, t), c in x21.items():
        u21[(a, q)] = u21.get((a, q), 0) + c * (-1) ** (t % 2)
    checks["t=q [2,1]"] = {k: v for k, v in u21.items() if v} == H[(2, 1)]
    log(f"  検算：Laurent {laurent}；t = q の [2,1] {checks['t=q [2,1]']}；\\hat f の負（spin, genus）: " +
        ", ".join(f"{k}:({v['spin_neg']},{v['genus_neg']})" for k, v in checks.items() if isinstance(v, dict)))
    return {"p": p, "nullity": 0, "x21": {",".join(map(str, m)): v for m, v in sorted(x21.items())},
            "x31": {",".join(map(str, m)): v for m, v in sorted(x31.items())}, "checks": checks, "seconds": time.time() - t0}


def nullspace_two_primes(st):
    """法 p₁, p₂ の RREF（状態 st）の零空間の基底を、CRT と有理数の復元で整数ベクトルにする（次元が一致することを確かめる）。"""
    from lmov21_polytope import ratrec
    bases = [nullspace_of(st["states"][pp], pp) for pp in PRIMES]
    assert len(bases[0]) == len(bases[1]), "零空間の次元が素数で違う"
    out = []
    m = PRIMES[0] * PRIMES[1]
    for v1, v2 in zip(bases[0], bases[1]):
        assert set(v1) == set(v2) or True
        vec = {}
        for c in set(v1) | set(v2):
            x = crt_sym(v1.get(c, 0) % PRIMES[0], v2.get(c, 0) % PRIMES[1])
            f = ratrec(x % m, m)
            if f is None:
                raise ValueError("有理数の復元に失敗")
            if f:
                vec[c] = f
        den = 1
        for f in vec.values():
            den = den * f.denominator // gcd(den, f.denominator)
        out.append({c: int(f * den) for c, f in vec.items()})
    return out


def positivity_stage(p, series=True, hat4=True, log=log):
    """4 箱の Laurent 性までの一次の系（lmov4_sys_K_<p>.pkl）の解空間 x = x_p + V w の上で、
    3 箱・4 箱の \\hat f の spin 表示の非負性と、S（[2,1] と [3,1] の P 正規化の冪級数の非負性）を課し、多面体の有界性と整数解を調べる。"""
    import pickle
    import scipy.sparse as sps
    from scipy.optimize import Bounds, LinearConstraint, linprog, milp
    from dgr_box3 import formula3
    from dgr_hu import unknot_numerator_factors
    from lmov21_hu import series_rows
    t0 = time.time()
    S = pickle.loads(Path(f"data/lmovstruct/lmov4_sys_K_{p}.pkl").read_bytes())
    mons21, mons31, B21, B31, x021i, x031i, dd, rows, rhs, small, boldd = (S[k] for k in ("mons21", "mons31", "B21", "B31", "x021i", "x031i", "dd", "rows", "rhs", "small", "boldd"))
    del rows, rhs
    ncols = len(B21) + len(B31)
    # 特解 y_p（RREF から。法 p で解き、2 素数で復元）
    st = rref_stage(p)
    ys = []
    for pp in PRIMES:
        y, nullity, bad = solution_of(st["states"][pp])
        assert not bad
        ys.append(y)
    from lmov21_polytope import ratrec
    m = PRIMES[0] * PRIMES[1]
    yp = {}
    for j in range(ncols):
        f = ratrec(crt_sym(int(ys[0][j]) % PRIMES[0], int(ys[1][j]) % PRIMES[1]) % m, m)
        if f is None:
            raise ValueError("特解の復元に失敗")
        yp[j] = f
    N = nullspace_two_primes(st)
    del st
    M21 = np.array([[float(v.get(j, 0)) for j in range(len(B21))] for v in N]) if N else np.zeros((0, len(B21)))
    rank21 = int(np.linalg.matrix_rank(M21)) if N and len(B21) else 0
    log(f"K_{p}: 連立の零空間 {len(N)}、[2,1] への射影の次元 {rank21}（[2,1] 単独 {len(B21)}）、特解の分母の最大 {max(f.denominator for f in yp.values())} ({time.time() - t0:.0f}s)")
    # x 空間へ：x = (x0i + dd·Σ y B)/dd
    def to_x(yvec, scale):
        x21, x31 = {}, {}
        for j, vec in enumerate(B21):
            c = yvec.get(j, 0) if isinstance(yvec, dict) else yvec[j]
            if c:
                for mm, z in vec.items():
                    x21[mm] = x21.get(mm, 0) + c * z * scale
        for j, vec in enumerate(B31):
            c = yvec.get(len(B21) + j, 0) if isinstance(yvec, dict) else yvec[len(B21) + j]
            if c:
                for mm, z in vec.items():
                    x31[mm] = x31.get(mm, 0) + c * z * scale
        return {k: v for k, v in x21.items() if v}, {k: v for k, v in x31.items() if v}
    xp21, xp31 = to_x(yp, Fraction(1))
    xp21 = {mm: Fraction(x021i.get(mm, 0), dd) + v for mm, v in {**{mm: Fraction(0) for mm in x021i}, **xp21}.items()}
    xp31 = {mm: Fraction(x031i.get(mm, 0), dd) + v for mm, v in {**{mm: Fraction(0) for mm in x031i}, **xp31}.items()}
    xp21 = {k: v for k, v in xp21.items() if v}; xp31 = {k: v for k, v in xp31.items() if v}
    V = [to_x(vec, 1) for vec in N]                            # 整数の基底（[2,1] 側, [3,1] 側）
    d = len(N)
    if d == 0:
        log("  自由度なし")
        return {"nullity": 0}
    # 線形な像：\hat f（3 箱・4 箱）と S を、特解と各基底ベクトルについて計算する
    def small_of(bd):
        den = 1
        for v in bd.values():
            den = den * Fraction(v).denominator // gcd(den, Fraction(v).denominator)
        return bold_to_small({k: int(Fraction(v) * den) for k, v in bd.items()}), den
    fo3 = formula3(); fo4, keys = formula_n(4)
    idx = {k: i for i, k in enumerate(keys)}
    e21 = tuple(1 if i in (idx[((2, 1), 1)], idx[((1,), 1)]) else 0 for i in range(len(keys)))
    e31 = tuple(1 if i == idx[((3, 1), 1)] else 0 for i in range(len(keys)))
    e211 = tuple(1 if i == idx[((2, 1, 1), 1)] else 0 for i in range(len(keys)))
    Rk = inputs_n({**small, (2, 1): {(0, 0, 0): 0}, (3, 1): {(0, 0, 0): 0}, (2, 1, 1): {(0, 0, 0): 0}}, 4)
    R3 = {}
    from dgr_box3 import KEYS, inputs_R
    R3 = inputs_R({**small, (2, 1): {(0, 0, 0): 0}})

    def lin_nums(x21, x31):
        """x に対する 3 箱・4 箱の \\hat f の分子の一次の部分（L）。{(箱, ρ): (分子, 分母)}。"""
        out = {}
        s21 = L.from_dict(bold_to_small(x21)); s31 = L.from_dict(bold_to_small(x31)); s211 = L.from_dict(bold_to_small(transpose_bold(x31)))
        for rho, (ns, den) in fo3.items():
            out[(3, rho)] = (ns[1] * s21, den)
        for rho, (ns, D) in fo4.items():
            out[(4, rho)] = (ns[e21] * Rk[((1,), 1)] * s21 + ns[e31] * s31 + ns[e211] * s211, D)
        return out

    def const_nums():
        out = {}
        for rho, (ns, den) in fo3.items():
            num = None
            for k, n in ns.items():
                if KEYS[k] == "r21":
                    continue
                term = n * R3[KEYS[k]]
                num = term if num is None else num + term
            out[(3, rho)] = num
        for rho, (ns, D) in fo4.items():
            num = None
            for expo, n in ns.items():
                if expo in (e21, e31, e211):
                    continue
                term = n
                for (lam, dg), k in zip(keys, expo):
                    for _ in range(k):
                        term = term * Rk[(lam, dg)]
                num = term if num is None else num + term
            out[(4, rho)] = num
        return out

    def hat_images(x21, x31):
        """基底ベクトル（解の差）の像：一次の部分だけで割り切れる。"""
        return {key: (num / den).to_dict() for key, (num, den) in lin_nums(x21, x31).items()}

    def hat_full(x21i, x31i, scale):
        """特解（scale 倍して整数）の像：scale·定数項 ＋ 一次の部分 を割る。"""
        cn = const_nums()
        return {key: ((cn[key] * L.from_dict({(0, 0, 0): scale}) + num) / den).to_dict() for key, (num, den) in lin_nums(x21i, x31i).items()}
    # 特解の像（有理数）：整数化して den_p 倍で扱う
    sp21, den21 = small_of(xp21); sp31, den31 = small_of(xp31)
    denp = den21 * den31 // gcd(den21, den31)
    xp21i = {k: int(v * denp) for k, v in xp21.items()}; xp31i = {k: int(v * denp) for k, v in xp31.items()}
    img_p = hat_full(xp21i, xp31i, denp)                     # denp·\\hat f(x_p)
    img_V = [hat_images(v21, v31) for v21, v31 in V]
    log(f"  像の計算 ({time.time() - t0:.0f}s)")
    # 太字の係数ごとの不等式：denp·const + img_p + Σ w_j denp·img_Vj ≥ 0（符号は定数側の多数派）
    ineqs = []
    shift_of = {}
    for key in img_p:
        allk = list(img_p[key]) + [k for im in img_V for k in im[key]]
        shift = 1 if any(k[0] % 2 for k in allk) else 0
        shift_of[key] = shift
        base = {}
        for k, v in to_bold(img_p[key], shift).items():
            base[k] = base.get(k, 0) + v
        sgn = 1 if sum(1 for v in base.values() if v > 0) >= sum(1 for v in base.values() if v < 0) else -1
        cols = {}
        for j, im in enumerate(img_V):
            for k, v in to_bold(im[key], shift).items():
                cols.setdefault(k, {})[j] = cols.get(k, {}).get(j, 0) + v * denp
        for k in set(cols) | set(base):
            r = {j: sgn * v for j, v in cols.get(k, {}).items() if v}
            if r:
                ineqs.append((r, -sgn * base.get(k, 0)))
    nhat = len(ineqs)
    # S：[2,1] と [3,1] の P 正規化の冪級数
    nS = 0
    if series:
        for lam, xpi, Vside in (((2, 1), xp21i, [v21 for v21, _ in V]), ((3, 1), xp31i, [v31 for _, v31 in V])):
            F = {(0, 0, 0): 1}
            for f in unknot_numerator_factors(lam):
                F = _mul(F, f)
            cols = {None: _mul(xpi, F)}
            for j, v in enumerate(Vside):
                cols[j] = _mul({k: val * denp for k, val in v.items()}, F)
            Np = cols[None]
            sgn = 1 if sum(1 for v in Np.values() if v > 0) >= sum(1 for v in Np.values() if v < 0) else -1
            for key, r in series_rows(cols, lam, 40).items():
                c0 = r.pop(None, 0)
                if r:
                    ineqs.append(({j: sgn * v for j, v in r.items()}, -sgn * c0))
                    nS += 1
    log(f"  不等式：\\hat f {nhat}、S {nS} ({time.time() - t0:.0f}s)")
    A = sps.lil_matrix((len(ineqs), d))
    for i, (r, _) in enumerate(ineqs):
        for j, v in r.items():
            A[i, j] = v
    A = A.tocsr(); b = np.array([float(rr) for _, rr in ineqs])
    Path(f"data/lmovstruct/lmov4_pos_K_{p}.pkl").write_bytes(pickle.dumps({"ineqs": ineqs, "V": V, "xp21i": xp21i, "xp31i": xp31i, "denp": denp, "nhat": nhat, "nS": nS}))
    # 実行可能性と有界性（w は実数）。不等式は整数係数で大きいので、行を最大係数で正規化する
    scale = np.array([max(1.0, max(abs(float(v)) for v in r.values())) for r, _ in ineqs])
    An = sps.diags(1.0 / scale) @ A; bn = b / scale
    r0 = linprog(np.zeros(d), A_ub=-An, b_ub=-bn, bounds=[(None, None)] * d, method="highs")
    feas = r0.status == 0
    log(f"  実行可能性の LP：status {r0.status}（{r0.message[:60]}）")
    rng = np.random.default_rng(0)
    widths, bounded = [], feas
    for _ in range(6 if feas else 0):
        c = rng.standard_normal(d)
        w = []
        for sg in (1, -1):
            r = linprog(sg * c, A_ub=-An, b_ub=-bn, bounds=[(None, None)] * d, method="highs")
            w.append(float(sg * r.fun) if r.status == 0 else None)
            if r.status not in (0, 3):
                log(f"  LP status {r.status}: {r.message[:60]}")
        widths.append(w)
        bounded = bounded and all(v is not None for v in w)
    log(f"  正値性の多面体：実行可能 {feas}、有界 {bounded}、幅の例 {[(None if a is None else round(a, 2), None if b_ is None else round(b_, 2)) for a, b_ in widths[:3]]} ({time.time() - t0:.0f}s)")
    out = {"nullity": d, "feasible": feas, "bounded": bounded, "n_hat_rows": nhat, "n_S_rows": nS}
    if feas and bounded:
        # 整数解の列挙：denp·x = denp·x_p + Σ w_j V_j が整数 ⇔ w_j ∈ Z（V は整数で原始的とは限らないが、w 整数なら十分）。w を整数に限って列挙
        c = rng.integers(-30, 31, size=d).astype(float)
        sols, lo = [], -np.inf
        while len(sols) < 50:
            rr = milp(c, constraints=[LinearConstraint(An, bn, np.inf), LinearConstraint(c.reshape(1, -1), lo, np.inf)],
                      bounds=Bounds(np.full(d, -np.inf), np.full(d, np.inf)), integrality=np.ones(d), options={"time_limit": 1800})
            if rr.x is None:
                break
            g = tuple(int(round(v)) for v in rr.x)
            if g in sols:
                break
            sols.append(g)
            lo = float(c @ np.array(g)) + 1
        out["count_w_integer"] = len(sols) if len(sols) < 50 else ">=50"
        log(f"  w が整数の解 {out['count_w_integer']} 個 ({time.time() - t0:.0f}s)")
        xs = []
        for g in sols[:5]:
            x21 = dict(xp21i); x31 = dict(xp31i)
            for j, (v21, v31) in enumerate(V):
                if g[j]:
                    for k, v in v21.items():
                        x21[k] = x21.get(k, 0) + g[j] * v * denp
                    for k, v in v31.items():
                        x31[k] = x31.get(k, 0) + g[j] * v * denp
            if all(v % denp == 0 for v in list(x21.values()) + list(x31.values())):
                xs.append({"x21": {",".join(map(str, k)): v // denp for k, v in sorted(x21.items()) if v}, "x31": {",".join(map(str, k)): v // denp for k, v in sorted(x31.items()) if v}})
        out["solutions"] = xs
    return out


def positivity_cached(p, use_S, log=log):
    """保存した不等式（lmov4_pos_K_<p>.pkl）で、S を外した／入れた多面体の実行可能性・有界性・整数解を調べる。"""
    import pickle
    import scipy.sparse as sps
    from scipy.optimize import Bounds, LinearConstraint, linprog, milp
    t0 = time.time()
    S = pickle.loads(Path(f"data/lmovstruct/lmov4_pos_K_{p}.pkl").read_bytes())
    ineqs = S["ineqs"]
    nhat = S.get("nhat", 5461)
    if not use_S:
        ineqs = ineqs[:nhat]
    V, xp21i, xp31i, denp = S["V"], S["xp21i"], S["xp31i"], S["denp"]
    d = len(V)
    A = sps.lil_matrix((len(ineqs), d))
    for i, (r, _) in enumerate(ineqs):
        for j, v in r.items():
            A[i, j] = v
    A = A.tocsr(); b = np.array([float(rr) for _, rr in ineqs])
    scale = np.array([max(1.0, max(abs(float(v)) for v in r.values())) for r, _ in ineqs])
    An = sps.diags(1.0 / scale) @ A; bn = b / scale
    r0 = linprog(np.zeros(d), A_ub=-An, b_ub=-bn, bounds=[(None, None)] * d, method="highs")
    feas = r0.status == 0
    log(f"K_{p} pos（S {use_S}、不等式 {len(ineqs)}、自由度 {d}）：実行可能 {feas}（status {r0.status}）")
    out = {"nullity": d, "use_S": use_S, "feasible": feas}
    if not feas:
        return out
    rng = np.random.default_rng(0)
    bounded, widths = True, []
    for _ in range(6):
        c = rng.standard_normal(d)
        w = []
        for sg in (1, -1):
            r = linprog(sg * c, A_ub=-An, b_ub=-bn, bounds=[(None, None)] * d, method="highs")
            w.append(float(sg * r.fun) if r.status == 0 else None)
        widths.append(w); bounded = bounded and all(v is not None for v in w)
    out["bounded"] = bounded; out["widths"] = widths[:3]
    log(f"  有界 {bounded}、幅の例 {[(None if a is None else round(a, 2), None if b_ is None else round(b_, 2)) for a, b_ in widths[:3]]} ({time.time() - t0:.0f}s)")
    if bounded:
        c = rng.integers(-30, 31, size=d).astype(float)
        sols, lo = [], -np.inf
        while len(sols) < 50:
            rr = milp(c, constraints=[LinearConstraint(An, bn, np.inf), LinearConstraint(c.reshape(1, -1), lo, np.inf)],
                      bounds=Bounds(np.full(d, -np.inf), np.full(d, np.inf)), integrality=np.ones(d), options={"time_limit": 1800})
            if rr.x is None:
                break
            g = tuple(int(round(v)) for v in rr.x)
            if g in sols:
                break
            sols.append(g); lo = float(c @ np.array(g)) + 1
        out["count_w_integer"] = len(sols) if len(sols) < 50 else ">=50"
        log(f"  w が整数の解 {out['count_w_integer']} 個 ({time.time() - t0:.0f}s)")
        xs = []
        for g in sols[:5]:
            x21 = dict(xp21i); x31 = dict(xp31i)
            for j, (v21, v31) in enumerate(V):
                if g[j]:
                    for k, v in v21.items():
                        x21[k] = x21.get(k, 0) + g[j] * v * denp
                    for k, v in v31.items():
                        x31[k] = x31.get(k, 0) + g[j] * v * denp
            ok = all(v % denp == 0 for v in list(x21.values()) + list(x31.values()))
            xs.append({"integral": ok, "x21": {",".join(map(str, k)): v // denp for k, v in sorted(x21.items()) if v} if ok else None,
                       "x31": {",".join(map(str, k)): v // denp for k, v in sorted(x31.items()) if v} if ok else None,
                       "x31_terms": sum(1 for v in x31.values() if v), "x31_neg": sum(1 for v in x31.values() if v < 0)})
        out["solutions"] = xs
    return out


if __name__ == "__main__" and len(sys.argv) > 2 and sys.argv[2] == "pos-cached":
    p = int(sys.argv[1])
    use_S = "--S" in sys.argv
    log(f"=== lmov4_predict pos-cached {time.strftime('%Y-%m-%d %H:%M:%S')} K_{p} S={use_S} ===")
    res = positivity_cached(p, use_S)
    Path(f"data/lmovstruct/lmov4_predict_pos_K_{p}_S{int(use_S)}.json").write_text(json.dumps(res, ensure_ascii=False, indent=1, default=str))

if __name__ == "__main__" and len(sys.argv) > 2 and sys.argv[2] == "pos":
    p = int(sys.argv[1])
    log(f"=== lmov4_predict pos {time.strftime('%Y-%m-%d %H:%M:%S')} K_{p} ===")
    res = positivity_stage(p)
    Path(f"data/lmovstruct/lmov4_predict_pos_K_{p}.json").write_text(json.dumps(res, ensure_ascii=False, indent=1, default=str))

if __name__ == "__main__" and len(sys.argv) > 3 and sys.argv[2] == "rref":
    rref_rho(int(sys.argv[1]), tuple(int(ch) for ch in sys.argv[3]))

if __name__ == "__main__" and not (len(sys.argv) > 2 and sys.argv[2] in ("pos", "pos-cached", "fix21", "rref")):
    p = int(sys.argv[1])
    log(f"=== lmov4_predict {time.strftime('%Y-%m-%d %H:%M:%S')} K_{p} ===")
    res = run(p)
    Path(f"data/lmovstruct/lmov4_predict_K_{p}.json").write_text(json.dumps(res, ensure_ascii=False, indent=1, default=str))
    if p in (1, -1):
        from dgr_box3_predict import known21
        k21 = known21(p)
        got = {tuple(map(int, m.split(","))): v for m, v in res.get("x21", {}).items()}
        log(f"  既知の [2,1] と一致：{got == k21}")


# ---- fix21：[2,1] を H11 の方法の予言に固定し、[3,1] の系が無矛盾か（4 箱の整数性の独立な検査） --------------------

def fix21_check(p, x21_pred, log=log):
    """保存した一次の系（lmov4_sys_K_<p>.pkl）で [2,1] を x21_pred に固定し、[3,1] の未知数だけの系の無矛盾性と零空間の次元を調べる。
    無矛盾なら、x21_pred は（ある [3,1] と組んで）4 箱の \\hat f の Laurent 性を満たす。"""
    import pickle
    from rcsknots.partitions import partitions
    t0 = time.time()
    S = pickle.loads(Path(f"data/lmovstruct/lmov4_sys_K_{p}.pkl").read_bytes())
    mons21, B21, x021i, dd, B31 = (S[k] for k in ("mons21", "B21", "x021i", "dd", "B31"))
    del S
    # y21：x21_pred = x021 + Σ y_j B21_j を有理数で解く（過剰決定、整合性を確かめる）
    n21 = len(B21)
    eqs = []
    for m in mons21:
        eqs.append(([Fraction(B21[j].get(m, 0)) for j in range(n21)], Fraction(x21_pred.get(m, 0)) - Fraction(x021i.get(m, 0), dd)))
    # ガウス消去
    piv = []
    M = [row[:] + [b] for row, b in eqs]
    r = 0
    for c in range(n21):
        pr = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        M[r] = [v / M[r][c] for v in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
    inconsistent = any(all(v == 0 for v in row[:n21]) and row[n21] != 0 for row in M)
    if inconsistent or len(piv) < n21:
        log(f"K_{p} fix21：予言が [2,1] の一次の解空間に無い（矛盾 {inconsistent}、階数 {len(piv)}/{n21}）")
        return {"in_affine_space": False}
    y21 = [M[i][n21] for i in range(n21)]
    log(f"K_{p} fix21：y21 = {[str(v) for v in y21]}（dd = {dd}）")
    # 4 箱の行：列 0..n21−1 を固定して右辺へ（行は dd 倍の系なので y21 もそのまま）。y21 の分母 den を掛けて整数にし、ρ ごとに流して消去する
    den = 1
    for v in y21:
        den = den * v.denominator // gcd(den, v.denominator)
    y21i = [int(v * den) for v in y21]
    ncols = len(B31)
    sts = {pp: new_state(ncols) for pp in PRIMES}
    nrows = 0
    for rho in partitions(4):
        rows, rhs = pickle.loads((Path("data/lmovstruct") / f"predict_rows_K_{p}_{''.join(map(str, rho))}.pkl").read_bytes())
        rows31, rhs31 = [], []
        for row, b in zip(rows, rhs):
            r31 = {c - n21: v * den for c, v in row.items() if c >= n21}
            b31 = b * den - sum(v * y21i[c] for c, v in row.items() if c < n21)
            if r31 or b31:
                rows31.append(r31); rhs31.append(b31)
        del rows, rhs
        nrows += len(rows31)
        for pp in PRIMES:
            eliminate(sts[pp], rows31, rhs31, pp)
        del rows31, rhs31
        log(f"  ρ={rho} を取り込み ({time.time() - t0:.0f}s)")
    out = {"in_affine_space": True}
    for pp in PRIMES:
        y, nullity, bad = solution_of(sts[pp])
        log(f"  [3,1] の系（式 {nrows}、列 {ncols}）法 {pp}：矛盾 {bad}、零空間 {nullity} ({time.time() - t0:.0f}s)")
        out[f"p{pp}"] = {"inconsistent": bad, "nullity": nullity}
    return out


if __name__ == "__main__" and len(sys.argv) > 2 and sys.argv[2] == "fix21":
    p = int(sys.argv[1])
    pred = json.load(open("data/lmovstruct/lmov21_h11_predict.json"))
    key = [k for k in pred if k.startswith(f"K_{p} x=-10")][0]
    x21 = {tuple(map(int, m.split(","))): v for m, v in pred[key]["face_solutions_rcs21"][0].items()}
    log(f"=== lmov4_predict fix21 {time.strftime('%Y-%m-%d %H:%M:%S')} K_{p}（H11 の予言 {len(x21)} 項）===")
    res = fix21_check(p, x21)
    Path(f"data/lmovstruct/lmov4_predict_fix21_K_{p}.json").write_text(json.dumps(res, ensure_ascii=False, indent=1, default=str))
