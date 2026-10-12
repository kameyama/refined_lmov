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


def solve_modp(rows, ncols, rhs, p, batch=1500):
    """A y = rhs（法 p）を密な消去で解く。返り値 (解 y または None, 零空間の次元, 矛盾か)。rows：[{列: 係数}]。"""
    piv_cols, B = [], np.zeros((0, ncols + 1))
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
                    return None, None, True
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
    nullity = ncols - len(piv_cols)
    y = np.zeros(ncols, dtype=object)
    for r, pc in enumerate(piv_cols):
        y[pc] = int(B[r, ncols])                               # 自由な変数は 0
    return y, nullity, False


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
    rows, rhs = [], []
    for rho, (ns, D) in fo.items():
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
        n_before = len(rows)
        for key, r in rr.items():
            c0 = r.pop(0, 0)
            r2 = {j - 1: v for j, v in r.items()}
            if r2 or c0:
                rows.append(r2)
                rhs.append(-c0)
        log(f"  ρ={rho}: 剰余の式 {len(rows) - n_before} ({time.time() - t0:.0f}s)")
    log(f"  消去：式 {len(rows)}、列 {ncols}")
    sols = []
    for pp in PRIMES:
        y, nullity, bad = solve_modp(rows, ncols, rhs, pp)
        log(f"  法 {pp}：矛盾 {bad}、零空間 {nullity}")
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


if __name__ == "__main__":
    p = int(sys.argv[1])
    log(f"=== lmov4_predict {time.strftime('%Y-%m-%d %H:%M:%S')} K_{p} ===")
    res = run(p)
    Path(f"data/lmovstruct/lmov4_predict_K_{p}.json").write_text(json.dumps(res, ensure_ascii=False, indent=1, default=str))
    if p in (1, -1):
        from dgr_box3_predict import known21
        k21 = known21(p)
        got = {tuple(map(int, m.split(","))): v for m, v in res.get("x21", {}).items()}
        log(f"  既知の [2,1] と一致：{got == k21}")
