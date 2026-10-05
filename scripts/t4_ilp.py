"""ステージ T4：GS の S² 彩色ホモロジー（t_r）に δ' を割り当てて t_c 次数の [2], [1,1] を作るとき、
\\hat f_{[2]}, \\hat f_{[1,1]} が Laurent 多項式で、G 展開の前後とも非負になる割り当てが存在するかを整数計画で判定する。

\\hat f_ρ = num_ρ / den（den = 2A²(Q⁴−1)(T⁴−1)(Q²T²−1)、(A, Q, T) = (a, q, t)^{1/2}）。num_ρ は reduced の r₂, r₁₁ について一次
（係数は f_formulas と unknot の因子から sympy で求める）。変数：
  z_{x,d} ∈ ℤ_{≥0}：生成元 x に δ' = d を割り当てる個数（Σ_d z_{x,d} = 重複度）
  H_e ∈ ℤ_{≥0}   ：σ·\\hat f_ρ の太字の係数（σ = ±1 は全体の符号）
  N_{g,β,F} ∈ ℤ_{≥0}：σ·\\hat N
条件：den·\\hat f = num（係数ごと）、太字の 𝐚^{2β} 行 = Σ N 𝐭^F G^g。
"""
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import scipy.sparse as sps
import sympy as sp
from scipy.optimize import Bounds, LinearConstraint, milp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from t4_colored import bold_to_small  # noqa: E402
from rcsknots.bps import C, C1, M_matrix, Q, T, f_formulas  # noqa: E402
from rcsknots.partitions import arm, boxes, leg  # noqa: E402

A = sp.symbols("A")
r1, r2, r11, r1d = sp.symbols("r1 r2 r11 r1d")


def _unk(lam, d=1):
    q, t, a = Q ** (2 * d), T ** (2 * d), A ** (2 * d)
    out = (T ** d / A ** d) ** sum(lam)
    for i, j in boxes(lam):
        out *= (t ** i - a * q ** j) / (1 - q ** arm(lam, i, j) * t ** (leg(lam, i, j) + 1))
    return out


def _poly(e):
    """sympy の Laurent 多項式 → {(e_A, e_Q, e_T): 整数}"""
    e = sp.expand(e * (A * Q * T) ** 80)
    return {(i - 80, j - 80, k - 80): int(c) for (i, j, k), c in sp.Poly(e, A, Q, T).terms()}


def structure():
    """{ρ: (c2, c11, rest_expr)} と den（num_ρ = c2·r2 + c11·r11 + rest）。"""
    F, M = f_formulas(2), M_matrix(2)
    parts = [(2,), (1, 1)]
    hat = M.inv() * sp.Matrix([F[m][0] for m in parts])
    sub = {C((1,)): _unk((1,)) * r1, C((2,)): _unk((2,)) * r2, C((1, 1)): _unk((1, 1)) * r11, C1(2): _unk((1,), 2) * r1d}
    den = 2 * A ** 2 * (Q ** 4 - 1) * (T ** 4 - 1) * (Q ** 2 * T ** 2 - 1)
    out = {}
    for i, rho in enumerate(parts):
        num = sp.expand(sp.cancel(sp.together(hat[i].subs(sub)) * den))
        out[rho] = (_poly(num.coeff(r2)), _poly(num.coeff(r11)), num.subs({r2: 0, r11: 0}))
    return out, _poly(den)


def mul(p, m):
    return {tuple(x + y for x, y in zip(k, m)): c for k, c in p.items()}


def solve(P1, gens, choices, sigma=1, struct=None, log=print, positive=True, cost=None):
    """P1：[1] の太字、gens：GS の {(i,j,k): 重複度}。実行可能なら {(x, d): 個数}、不能なら None。

    positive=False なら正値性を課さない（Laurent 性と G 展開が閉じることだけ）。cost は目的関数（最小化）。
    """
    struct, den = struct or structure()
    r1p = bold_to_small(P1)
    r1e = sum(c * A ** i * Q ** j * T ** k for (i, j, k), c in r1p.items())
    r1de = r1e.subs({A: A ** 2, Q: Q ** 2, T: T ** 2}, simultaneous=True)
    xs = sorted(gens)
    var = [(x, d) for x in xs for d in choices]
    nz = len(var)
    cand = []
    for (i, j, k), d in var:
        d4 = round(4 * d)
        m2 = next(iter(bold_to_small({(2 * i, 2 * j, 4 * i - k + 2 * j + d4): 1}).items()))
        m11 = next(iter(bold_to_small({(2 * i, -2 * j, 4 * i - k + d4): 1}).items()))
        cand.append((m2, m11))
    rows_eq, rhs_eq = [], []          # list of {col: coef}
    ncols = nz
    bounds_lo, bounds_hi = [0.0] * nz, [float(gens[x]) for x, _ in var]
    for x in xs:                       # Σ_d z = multiplicity
        rows_eq.append({n: 1 for n, (y, _) in enumerate(var) if y == x})
        rhs_eq.append(gens[x])
    for rho, (c2, c11, rest) in struct.items():
        # num = Σ_n z_n (c2·m2_n + c11·m11_n) + rest
        num = defaultdict(dict)
        for n, ((e2, s2), (e11, s11)) in enumerate(cand):
            for e, c in mul(c2, e2).items():
                num[e][n] = num[e].get(n, 0) + s2 * c
            for e, c in mul(c11, e11).items():
                num[e][n] = num[e].get(n, 0) + s11 * c
        restp = _poly(rest.subs({r1: r1e, r1d: r1de}))
        keys = set(num) | set(restp)
        # \hat f の箱（太字の a 次数シフトは |ρ|=2 で 0）
        lo = [min(e[i] for e in keys) for i in range(3)]
        hi = [max(e[i] for e in keys) for i in range(3)]
        dlo = [min(e[i] for e in den) for i in range(3)]
        dhi = [max(e[i] for e in den) for i in range(3)]
        hbox = [(a_, q_, t_) for a_ in range(lo[0] - dlo[0], hi[0] - dhi[0] + 1, 2)
                for q_ in range(lo[1] - dlo[1], hi[1] - dhi[1] + 1)
                for t_ in range(lo[2] - dlo[2], hi[2] - dhi[2] + 1)]
        hcol = {e: ncols + n for n, e in enumerate(hbox)}
        ncols += len(hbox)
        # small h_e = sgn_e · H_bold,  sgn from to_bold: (−1)^{m+eq}, m = e_A/2
        sgn = {e: (-1) ** ((e[0] // 2 + e[1]) % 2) for e in hbox}
        bounds_lo += [0.0 if positive else -np.inf] * len(hbox)
        bounds_hi += [np.inf] * len(hbox)
        # den · h = σ num
        conv = defaultdict(dict)
        for e in hbox:
            for de, dc in den.items():
                conv[tuple(x + y for x, y in zip(e, de))][hcol[e]] = dc * sgn[e]
        for key in set(conv) | keys:
            row = dict(conv.get(key, {}))
            for n, c in num.get(key, {}).items():
                row[n] = row.get(n, 0) - sigma * c
            rows_eq.append(row)
            rhs_eq.append(sigma * restp.get(key, 0))
        # G 展開：太字の (𝐚, 𝐪, 𝐭) = (2m, e_Q + e_T, m + e_Q)
        bold = {}
        for e in hbox:
            m = e[0] // 2
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
                    bounds_lo.append(0.0 if positive else -np.inf)
                    bounds_hi.append(np.inf)
                    for (gq, gt), gc in Gp[g].items():
                        gcons[(2 * beta, gq, gt + F)][col] = gc
        for key in set(gcons) | set(bold):
            row = {c: -v for c, v in gcons.get(key, {}).items()}
            if key in bold:
                row[bold[key]] = row.get(bold[key], 0) + 1
            rows_eq.append(row)
            rhs_eq.append(0)
        log(f"  ρ={rho}: num 項 {len(keys)}, h 箱 {len(hbox)}, 列 {ncols}")
    data, ri, ci = [], [], []
    for r, row in enumerate(rows_eq):
        for c, v in row.items():
            if v:
                ri.append(r)
                ci.append(c)
                data.append(v)
    Amat = sps.csr_matrix((data, (ri, ci)), shape=(len(rows_eq), ncols))
    b = np.array(rhs_eq, dtype=float)
    log(f"  行 {len(rows_eq)}, 列 {ncols}")
    c = np.zeros(ncols)
    if cost is not None:                # cost(x, d)：δ' の選択に付ける重み
        c[:nz] = [cost(x, d) for x, d in var]
    res = milp(c, constraints=LinearConstraint(Amat, b, b), integrality=np.ones(ncols),
               bounds=Bounds(np.array(bounds_lo), np.array(bounds_hi)), options={"time_limit": 3600})
    log(f"  {res.message}")
    if res.x is None:
        return None
    return {var[n]: int(round(res.x[n])) for n in range(nz) if round(res.x[n])}
