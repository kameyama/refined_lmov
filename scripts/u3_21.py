"""第 3 期 U3：rCS_[2,1] を構造的な条件から決める整数計画。

未知数 x：rCS_[2,1] の係数（KN の小文字 (A, Q, T) = (a, q, t)^{1/2} の指数、整数、符号は自由）。条件（全て一次）：
  (1) t = q で [2,1] 彩色 HOMFLY-PT（scripts/colored_from_cables.py の colored_p）
  (2) ミラー対称性 x(e_A, e_Q, e_T) = x(e_A, −e_T, −e_Q)
  (3) refined 指数増大性：Q = ±1 で r21 = r1·r11、T = ±1 で r21 = r1·r2
  (4) 両手型性（指定時）：x(e_A, e_Q, e_T) = x(−e_A, −e_Q, −e_T)
  (8) U(2) への簡約【仮説】：A = ±T²（a = t²）で r21 = r1·r11（トーラス結び目と 4_1 の正解で厳密に成立）
  (5) \\hat f_ρ（ρ ⊢ 3）が Laurent：\\hat f_ρ = α_ρ·unk_21·r21 + K_ρ。D·h_ρ = (D α_ρ unk_21)·x + D K_ρ（h は \\hat f の係数）
  (6) G 展開：太字（a 次数シフトは nontorus の +1）の各 𝐚 行 = Σ N 𝐭^F G^g
  (7) 正値性（指定時）：σ·h（太字）≥ 0、σ·N ≥ 0
K_ρ は r21 = 0 とおいた \\hat f_ρ。D·K_ρ を mod p の評価と補間で厳密に求める。
"""
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import scipy.sparse as sps
import sympy as sp
from scipy.optimize import Bounds, LinearConstraint, milp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rcsknots.bps import C, M_matrix, Q, T, f_formulas  # noqa: E402
from rcsknots.bpsnum import Reformulated  # noqa: E402
from rcsknots.exact import reconstruct_points  # noqa: E402
from rcsknots.extract import A_OFFSET, _arr  # noqa: E402
from rcsknots.field import V, make_conv  # noqa: E402
from rcsknots.partitions import arm, boxes, leg  # noqa: E402

A = sp.symbols("A")
PARTS = [(3,), (2, 1), (1, 1, 1)]
# D：α_ρ の分母 (Q² − 1)(Q⁴T² − 1) と unk_21 の分母 (1 − Q²T⁴)(1 − T²)²。D·K_ρ が多項式になることは dK の補間で確かめる
FACTORS = [(1, 0), (2, 1), (1, 2), (0, 1), (0, 1)]                            # 1 − q^i t^j（(i, j) は q, t の指数）


def den_poly():
    D = A ** 4
    for i, j in FACTORS:
        D *= (1 - Q ** (2 * i) * T ** (2 * j))
    return sp.expand(D)


def _poly(e):
    e = sp.expand(e * (A * Q * T) ** 120)
    return {(i - 120, j - 120, k - 120): int(c) for (i, j, k), c in sp.Poly(e, A, Q, T).terms()}


def alpha_unk():
    """{ρ: D·α_ρ·unk_21}（多項式）。α_ρ = ∂\\hat f_ρ/∂C_21。"""
    F = f_formulas(3)
    d = [sp.diff(F[m][0], C((2, 1))) for m in PARTS]
    Minv = M_matrix(3).inv()
    unk = (T / A) ** 3
    for i, j in boxes((2, 1)):
        unk *= (T ** (2 * i) - A ** 2 * Q ** (2 * j)) / (1 - Q ** (2 * arm((2, 1), i, j)) * T ** (2 * (leg((2, 1), i, j) + 1)))
    D = den_poly()
    out = {}
    for r, rho in enumerate(PARTS):
        al = sum(sp.together(Minv[r, j]) * d[j] for j in range(3))
        out[rho] = _poly(sp.cancel(al * unk * D))
    return out


def dK(reduced):
    """{ρ: D·K_ρ}（r21 = 0 の \\hat f_ρ に D を掛けた多項式）。"""
    red = dict(reduced)
    red[(2, 1)] = {}
    Dp = _poly(den_poly())
    out = {}
    for rho in PARTS:
        def ev(xs, ys, p, rho=rho):
            k = len(xs)
            Qv, Tv = V(xs, p), V(ys, p)
            R = Reformulated(red, Qv, Tv, make_conv(p, k))
            hat, _ = R.hat_f(3)
            Dv = np.zeros(k, dtype=np.int64)
            for (ea, eq, et), c in Dp.items():
                Dv = (Dv + c * (Qv ** eq * Tv ** et).a * 1) % p if ea == 0 else Dv
            arr = np.zeros((k, 2 * A_OFFSET + 1), dtype=np.int64)
            for (ea, eq, et), c in Dp.items():
                mono = (Qv ** eq * Tv ** et).a * (c % p) % p
                for e, v in hat[rho].terms.items():
                    arr[:, e + ea + A_OFFSET] = (arr[:, e + ea + A_OFFSET] + mono * _arr(v, k, p)) % p
            return arr
        raw = reconstruct_points(ev, seed=7, log=lambda s: None)
        out[rho] = {(d - A_OFFSET, eq, et): c for (d, eq, et), c in raw.items()}
    return out


def homfly_t_eq_q(H21):
    """sympy の (a, Q) の式（rcsknots.homfly の変数）→ {(e_A, e_Q): c}"""
    from rcsknots.homfly import Q as HQ, a as Ha
    e = sp.expand(H21.subs(Ha, A ** 2).subs(HQ, Q) * (A * Q) ** 120)
    return {(i - 120, j - 120): int(c) for (i, j), c in sp.Poly(e, A, Q).terms()}


def _mul(p1, p2):
    out = defaultdict(int)
    for k1, c1 in p1.items():
        for k2, c2 in p2.items():
            out[tuple(x + y for x, y in zip(k1, k2))] += c1 * c2
    return {k: v for k, v in out.items() if v}


def build(reduced, H21, amphichiral=False, positive=True, sigma=1, margin=2, structure=None, log=print, u2=True):
    """整数計画を作って解く。戻り値：(解の x の辞書 or None, 等式系の情報)。"""
    au, dk = structure or (alpha_unk(), dK(reduced))
    r1, r2, r11 = reduced[(1,)], reduced[(2,)], reduced[(1, 1)]
    egp = list(_mul(r1, r11)) + list(_mul(r1, r2))
    lo = [min(e[i] for e in egp) - margin for i in range(3)]
    hi = [max(e[i] for e in egp) + margin for i in range(3)]
    # x の箱：ミラー対称性で閉じるように Q と T の範囲をそろえる
    qt_lo, qt_hi = min(lo[1], -hi[2]), max(hi[1], -lo[2])
    tt_lo, tt_hi = min(lo[2], -hi[1]), max(hi[2], -lo[1])
    a_par = next(iter(r1))[0] % 2 + 1        # [2,1] の A の偶奇は 3 箱ぶん（r1 の偶奇 × 3）
    xs = [(ea, eq, et) for ea in range(lo[0], hi[0] + 1) if (ea - 3 * next(iter(r1))[0]) % 2 == 0
          for eq in range(qt_lo, qt_hi + 1) for et in range(tt_lo, tt_hi + 1)
          if (eq + et - sum(next(iter(r1))[1:]) * 3) % 2 == 0]
    del a_par
    xcol = {e: n for n, e in enumerate(xs)}
    ncols = len(xs)
    lb, ub = [-np.inf] * ncols, [np.inf] * ncols
    rows, rhs = [], []

    def eq_row(row, b):
        rows.append(row)
        rhs.append(b)

    # (1) t = q
    tq = defaultdict(dict)
    for e, n in xcol.items():
        tq[(e[0], e[1] + e[2])][n] = 1
    for key in set(tq) | set(H21):
        eq_row(tq.get(key, {}), H21.get(key, 0))
    # (2) ミラー、(4) 両手型
    for e, n in xcol.items():
        m = (e[0], -e[2], -e[1])
        if m not in xcol:
            eq_row({n: 1}, 0)
        elif m > e:
            eq_row({n: 1, xcol[m]: -1}, 0)
        if amphichiral:
            m2 = (-e[0], -e[1], -e[2])
            if m2 not in xcol:
                eq_row({n: 1}, 0)
            elif m2 > e:
                eq_row({n: 1, xcol[m2]: -1}, 0)
    # (3) EGP：Q = s で r21 = r1 r11、T = s で r21 = r1 r2（s = ±1）
    for which, prod in (("Q", _mul(r1, r11)), ("T", _mul(r1, r2))):
        for s in (1, -1):
            ev = defaultdict(dict)
            for e, n in xcol.items():
                key = (e[0], e[2]) if which == "Q" else (e[0], e[1])
                sg = s ** (e[1] if which == "Q" else e[2])
                ev[key][n] = ev[key].get(n, 0) + sg
            tgt = defaultdict(int)
            for e, c in prod.items():
                key = (e[0], e[2]) if which == "Q" else (e[0], e[1])
                tgt[key] += c * s ** (e[1] if which == "Q" else e[2])
            for key in set(ev) | set(tgt):
                eq_row(ev.get(key, {}), tgt.get(key, 0))
    # (8) A = s T² で r21 = r1·r11（s = ±1）
    if u2:
        prod = _mul(r1, r11)
        for s_ in (1, -1):
            ev = defaultdict(dict)
            for e, n in xcol.items():
                key = (e[1], e[2] + 2 * e[0])
                ev[key][n] = ev[key].get(n, 0) + s_ ** e[0]
            tgt = defaultdict(int)
            for e, c in prod.items():
                tgt[(e[1], e[2] + 2 * e[0])] += c * s_ ** e[0]
            for key in set(ev) | set(tgt):
                eq_row(ev.get(key, {}), tgt.get(key, 0))
    # (5)(6)(7)
    Dp = _poly(den_poly())
    for rho in PARTS:
        num = defaultdict(dict)
        for e, n in xcol.items():
            for k, c in au[rho].items():
                kk = tuple(x + y for x, y in zip(e, k))
                num[kk][n] = num[kk].get(n, 0) + c
        keys = set(num) | set(dk[rho])
        dlo = [min(e[i] for e in Dp) for i in range(3)]
        dhi = [max(e[i] for e in Dp) for i in range(3)]
        lo_ = [min(e[i] for e in keys) - dlo[i] for i in range(3)]
        hi_ = [max(e[i] for e in keys) - dhi[i] for i in range(3)]
        hbox = [(a_, q_, t_) for a_ in range(lo_[0], hi_[0] + 1) for q_ in range(lo_[1], hi_[1] + 1)
                for t_ in range(lo_[2], hi_[2] + 1) if (a_ + 1) % 2 == 0]
        hcol = {e: ncols + n for n, e in enumerate(hbox)}
        ncols += len(hbox)
        lb += [0.0 if positive else -np.inf] * len(hbox)
        ub += [np.inf] * len(hbox)
        # 太字の係数 H = sgn·h（a 次数シフト +1、to_bold と同じ符号）
        sgn = {e: (-1) ** ((((e[0] + 1) // 2) + e[1]) % 2) for e in hbox}
        conv = defaultdict(dict)
        for e in hbox:
            for de, dc in Dp.items():
                conv[tuple(x + y for x, y in zip(e, de))][hcol[e]] = dc * sgn[e]
        for key in set(conv) | keys:
            row = dict(conv.get(key, {}))
            for n, c in num.get(key, {}).items():
                row[n] = row.get(n, 0) - sigma * c
            eq_row(row, sigma * dk[rho].get(key, 0))
        bold = {}
        for e in hbox:
            m = (e[0] + 1) // 2
            bold[(2 * m, e[1] + e[2], m + e[1])] = hcol[e]
        betas = sorted({b[0] // 2 for b in bold})
        gmax = max(abs(b[1]) for b in bold) // 2 + 1
        tlo, thi = min(b[2] for b in bold), max(b[2] for b in bold)
        G = {(2, 1): 1, (0, 1): -1, (0, -1): -1, (-2, -1): 1}
        Gp = [{(0, 0): 1}]
        for _ in range(gmax):
            nxt = defaultdict(int)
            for (a1, b1), x in Gp[-1].items():
                for (a2, b2), y in G.items():
                    nxt[(a1 + a2, b1 + b2)] += x * y
            Gp.append({k: v for k, v in nxt.items() if v})
        gcons = defaultdict(dict)
        for beta in betas:
            for g in range(gmax + 1):
                for F in range(tlo - g, thi + g + 1):
                    col = ncols
                    ncols += 1
                    lb.append(0.0 if positive else -np.inf)
                    ub.append(np.inf)
                    for (gq, gt), gc in Gp[g].items():
                        gcons[(2 * beta, gq, gt + F)][col] = gc
        for key in set(gcons) | set(bold):
            row = {c: -v for c, v in gcons.get(key, {}).items()}
            if key in bold:
                row[bold[key]] = row.get(bold[key], 0) + 1
            eq_row(row, 0)
    data, ri, ci = [], [], []
    for r, row in enumerate(rows):
        for c, v in row.items():
            if v:
                ri.append(r)
                ci.append(c)
                data.append(v)
    Amat = sps.csr_matrix((data, (ri, ci)), shape=(len(rows), ncols))
    b = np.array(rhs, dtype=float)
    log(f"  x {len(xs)}、行 {len(rows)}、列 {ncols}")
    return xs, Amat, b, np.array(lb), np.array(ub)


def solve(reduced, H21, cost=None, structure=None, log=print, **kw):
    xs, Amat, b, lb, ub = build(reduced, H21, structure=structure, log=log, **kw)
    c = np.zeros(Amat.shape[1])
    if cost is not None:
        c[:len(xs)] = [cost(e) for e in xs]
    res = milp(c, constraints=LinearConstraint(Amat, b, b), integrality=np.ones(Amat.shape[1]),
               bounds=Bounds(lb, ub), options={"time_limit": 3600})
    log(f"  {res.message}")
    if res.x is None:
        return None
    return {e: int(round(v)) for e, v in zip(xs, res.x[:len(xs)]) if round(v)}
