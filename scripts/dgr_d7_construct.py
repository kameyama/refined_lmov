"""段階 D7（docs/research_plan_dgr_lmov.md）：非 refined のデータと見つけた構造から、2 箱の refined LMOV（と rCS_{[2]}）を構成する。

未知数（全て太字 (𝐚, 𝐪, 𝐭) の係数、整数）：
  x_m：rCS_{[2]} の係数（窓の中の単項式 m）。rCS に正値性は課さない（open_problems §F）。rCS_{[1,1]} は転置対称性（太字の σ）で rCS_{[2]} から決める。
  h_{ρ,e} ≥ 0：σ_ρ·\\hat f_ρ の係数（ρ = [2], [1,1]、σ_ρ = ±1 は全体の符号）。この計画では LMOV の正値性を条件に使ってよい（ユーザーの指示）。
条件：
  (U) 𝐭 = −1 で非 refined の [2]（t = q の色付き HOMFLY-PT、太字）に戻る。
  (S) (S_{1|0})：rCS_{[2]}|_{d_{1|0}} = c²（c = rCS_{[1]}|_{d_{1|0}}、H1 で導出）。
  (H) 変種 (b) のみ：hook 性質 rCS_{[2]}|_{d_{1|1}} = c·A(𝐪²𝐭, 𝐭)、A = rCS_{[1]}|_{d_{1|1}}。
  (L) Laurent 性と \\hat f の定義：den·\\hat f_ρ = num_ρ（2 箱の公式、scripts/v4_closed.py の formula。rCS について一次）。
  (E) 任意：refined の指数増大性 rCS_{[2]}(q=1) = rCS_{[1]}(q=1)²（KN (EGP)）。
  (P) h ≥ 0（spin 表示の正値性）。任意で genus 表示：σ 不変性と、genus 展開の係数 N ≥ 0 を等式で課す。
  (Π) 任意：parity knot の 2J_r の偶奇（F − β − g ≡ c mod 2、c は [1] の genus 表示から読む）。
一意性：LP 緩和で各 x_m の最小と最大を求め、全て一致すれば（実数の範囲でも）一意。
入力 [1] は変種 (a) でも使う（2 箱の公式に [1] が入るため）。(a) と (b) の違いは条件 (H) の有無。
"""
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sps
from scipy.optimize import Bounds, LinearConstraint, linprog, milp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from t4_colored import bold_to_small  # noqa: E402
from v4_closed import RHOS, L, formula, mir  # noqa: E402
from rcsknots.bps import g_expand, to_bold  # noqa: E402
from rcsknots.dgr import specialize  # noqa: E402


def unrefined(P):
    out = {}
    for (ea, eq, et), c in P.items():
        out[(ea, eq)] = out.get((ea, eq), 0) + c * (-1) ** (et % 2)
    return {k: v for k, v in out.items() if v}


class Lin:
    """変数の名前 → 列番号と、疎な等式・不等式の組み立て。"""

    def __init__(self):
        self.col, self.rows, self.rhs = {}, [], []

    def var(self, name):
        if name not in self.col:
            self.col[name] = len(self.col)
        return self.col[name]

    def eq(self, coeffs, rhs):
        self.rows.append({self.var(n): c for n, c in coeffs.items() if c})
        self.rhs.append(rhs)


def build(P1, H2u, window_t, variant="b", sigmas=(1, 1), genus=False, egp=False, parity=None):
    """P1：rCS_{[1]}（太字）、H2u：非 refined の [2]（{(e_𝐚, e_𝐪): c}）、window_t：rCS_{[2]} の 𝐭 の窓。"""
    lin = Lin()
    mons = [(a, q, t) for (a, q) in H2u for t in range(*window_t)]
    # 𝐭 = −1 で 0 になる (a, q) の組も、(U) の右辺 0 として未知数を置く必要がある：H2u の (a, q) 範囲の箱をとる
    As = sorted({a for a, _ in H2u}); Qs = sorted({q for _, q in H2u})
    mons = [(a, q, t) for a in range(As[0], As[-1] + 1, 2) for q in range(Qs[0], Qs[-1] + 1, 2) for t in range(*window_t)]
    for m in mons:
        lin.var(("x", m))
    # (U)
    for a in range(As[0], As[-1] + 1, 2):
        for q in range(Qs[0], Qs[-1] + 1, 2):
            lin.eq({("x", (a, q, t)): (-1) ** (t % 2) for t in range(*window_t)}, H2u.get((a, q), 0))
    # (S)：specialize は線形。単項式ごとの像を集める
    c = specialize(P1, 1, 0)
    (cq, ct), cs = next(iter(c.items()))
    target = {(0, 2 * ct): cs * cs}
    img = {}
    for m in mons:
        for k, v in specialize({m: 1}, 1, 0).items():
            img.setdefault(k, {})[("x", m)] = v
    for k in set(img) | set(target):
        lin.eq(img.get(k, {}), target.get(k, 0))
    # (H)
    if variant == "b":
        A = specialize(P1, 1, 1)
        tgt = {(2 * eq, et + eq + ct): cs * v for (eq, et), v in A.items()}       # c·A(𝐪²𝐭, 𝐭)（c = cs·𝐭^{ct}）
        img = {}
        for m in mons:
            for k, v in specialize({m: 1}, 1, 1).items():
                img.setdefault(k, {})[("x", m)] = v
        for k in set(img) | set(tgt):
            lin.eq(img.get(k, {}), tgt.get(k, 0))
    # (E) refined の指数増大性（KN (EGP)）：q = 1（太字で 𝐪 = 𝐭^{−1}、𝐪 の指数は偶数）で rCS_{[2]} = rCS_{[1]}²
    if egp:
        sq = {}
        for (a1, q1, t1), v1 in P1.items():
            for (a2, q2, t2), v2 in P1.items():
                k = (a1 + a2, t1 + t2 - q1 - q2)
                sq[k] = sq.get(k, 0) + v1 * v2
        img = {}
        for m in mons:
            img.setdefault((m[0], m[2] - m[1]), {})[("x", m)] = 1
        for k in set(img) | set(sq):
            lin.eq(img.get(k, {}), sq.get(k, 0))
    # (L)：num_ρ = ns0·r2 + ns1·r11 + ns2·r1² + ns3·r1(x²)。r2, r11 の部分は x について一次
    r1 = L.from_dict(bold_to_small(P1))
    r1d = L.from_dict({(2 * i, 2 * j, 2 * k): v for (i, j, k), v in bold_to_small(P1).items()})
    for idx, (rho, (ns, den)) in enumerate(formula().items()):
        const = (ns[2] * (r1 * r1) + ns[3] * r1d).to_dict()
        lin_part = {}
        for m in mons:
            v = (ns[0] * L.from_dict(bold_to_small({m: 1})) + ns[1] * L.from_dict(bold_to_small(mir({m: 1})))).to_dict()
            for k, w in v.items():
                lin_part.setdefault(k, {})[("x", m)] = w
        dd = den.to_dict()
        # \\hat f の窓（小文字）：num の箱 − den の箱
        keys = set(lin_part) | set(const)
        lo = [min(k[i] for k in keys) - min(k[i] for k in dd) for i in range(3)]
        hi = [max(k[i] for k in keys) - max(k[i] for k in dd) for i in range(3)]
        hmons = [(a, q, t) for a in range(lo[0], hi[0] + 1) for q in range(lo[1], hi[1] + 1) for t in range(lo[2], hi[2] + 1)
                 if a % 2 == 0]
        # 小文字の単項式 e の太字での符号（to_bold の符号）
        sgn = {e: next(iter(to_bold({e: 1}, 0).values())) for e in hmons}
        conv = {}
        for e in hmons:
            for k, w in dd.items():
                key = (e[0] + k[0], e[1] + k[1], e[2] + k[2])
                conv.setdefault(key, {})[("h", rho, e)] = w * sgn[e] * sigmas[idx]       # \\hat f[e] = σ·sgn·h
        for key in set(conv) | keys:
            row = dict(conv.get(key, {}))
            for n_, w in lin_part.get(key, {}).items():
                row[n_] = row.get(n_, 0) - w
            lin.eq(row, const.get(key, 0))
        if genus:
            add_genus(lin, rho, hmons, sgn, parity=parity)
    return lin, mons


def add_genus(lin, rho, hmons, sgn, gmax=None, parity=None):
    """σ 不変性と genus 表示の正値性：太字の h（σ_ρ \\hat f の係数）= Σ_{g,F} N_{g,β,F} 𝐭^F G^g（N ≥ 0）を等式で課す。"""
    bold = {}
    for e in hmons:
        b = next(iter(to_bold({e: 1}, 0)))
        bold[b] = ("h", rho, e)
    rows = {}
    for (ea, eq, et), name in bold.items():
        rows.setdefault(ea // 2, {})[(eq, et)] = name
    Gp = {(2, 1): 1, (0, 1): -1, (0, -1): -1, (-2, -1): 1}

    def mul(u, v):
        out = {}
        for (a1, b1), x in u.items():
            for (a2, b2), y in v.items():
                k = (a1 + a2, b1 + b2)
                out[k] = out.get(k, 0) + x * y
        return {k: w for k, w in out.items() if w}
    for beta, X in rows.items():
        qmax = max(eq for eq, _ in X)
        g_top = qmax // 2 if gmax is None else gmax
        ts = [et for _, et in X]
        powers = [{(0, 0): 1}]
        for _ in range(g_top):
            powers.append(mul(powers[-1], Gp))
        expr = {}
        for g in range(g_top + 1):
            for F in range(min(ts) - g, max(ts) + 1):
                if parity is not None and (F - beta - g - parity) % 2:
                    continue                      # parity knot：F − β − g の偶奇は [1] と同じ
                name = ("N", rho, g, beta, F)
                for (dq, dt), c in powers[g].items():
                    expr.setdefault((dq, F + dt), {})[name] = c
        for key in set(expr) | set(X):
            row = {n_: -c for n_, c in expr.get(key, {}).items()}
            if key in X:
                row[X[key]] = row.get(X[key], 0) + 1
            lin.eq(row, 0)


def solve(lin, xbound=60, objective=None, integral=True):
    n = len(lin.col)
    Aeq = sps.lil_matrix((len(lin.rows), n))
    for i, row in enumerate(lin.rows):
        for j, v in row.items():
            Aeq[i, j] = v
    lb = np.array([(-xbound if name[0] == "x" else 0) for name in lin.col], dtype=float)
    # h と N は 0 以上（正値性）。h は太字の係数で、σ_ρ を掛けた後の値
    ub = np.array([(xbound if name[0] == "x" else np.inf) for name in lin.col], dtype=float)
    c = np.zeros(n) if objective is None else objective
    cons = LinearConstraint(Aeq.tocsr(), np.array(lin.rhs, float), np.array(lin.rhs, float))
    if integral:
        integ = np.ones(n)
        return milp(c, constraints=[cons], bounds=Bounds(lb, ub), integrality=integ, options={"time_limit": 600})
    return linprog(c, A_eq=Aeq.tocsr(), b_eq=np.array(lin.rhs, float), bounds=list(zip(lb, [None if np.isinf(u) else u for u in ub])),
                   method="highs")


def ranges(lin, mons):
    """LP 緩和での各 x_m の最小・最大。"""
    out = {}
    for m in mons:
        j = lin.col[("x", m)]
        obj = np.zeros(len(lin.col)); obj[j] = 1
        lo = solve(lin, objective=obj, integral=False)
        hi = solve(lin, objective=-obj, integral=False)
        out[m] = (lo.fun if lo.status == 0 else None, -hi.fun if hi.status == 0 else None)
    return out


def solution(lin, res, mons):
    x = res.x
    return {m: int(round(x[lin.col[("x", m)]])) for m in mons if round(x[lin.col[("x", m)]]) != 0}


def parity_class(P1):
    """[1] の genus 表示の F − β − g の偶奇（一様でなければ None）。\\hat f_{[1]} = (1+𝐚²𝐭)P1 の a 次数のずらしは nontorus 規約。"""
    from rcsknots.dgr import quotient  # noqa: F401
    N = g_expand(P1)
    cls = {(F - b - g) % 2 for (g, b, F) in N}
    return cls.pop() if len(cls) == 1 else None
