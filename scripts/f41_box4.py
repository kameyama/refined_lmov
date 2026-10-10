"""4_1 の 4 箱：未知の rCS_[3,1]（[2,1,1] は転置で決まる）に一次の条件を課し、4 箱の refined LMOV 不変量が正になる解が存在するかを判定する（弱い判定）。

入力（KN の小文字 (A, Q, T) = (a, q, t)^{1/2} の指数、reduced、framing 0）：
  [1], [2], [1,1], [3], [1,1,1], [4], [2,2], [1⁴]：KNTZ 1902.02275 の refined 版（scripts/x1_check.triply、4_1 = K_{−1,1}）。
    [1]〜[1,1,1] は KN のノートブックと、[r] は GNSSS の式と全係数一致（このスクリプトで確かめる）。
  [2,1]：KN の候補（KN のノートブック、未確立）。4 箱の \\hat f は非線形の項を通じてこれに依存する。
  [3,1] の t = q：Morozov 1612.00422 の [3,1] の微分展開（:1987-2043）。4_1 = K_{−1,1} では F の比が 1 なので Z 因子の和。
未知数 x = rCS_[3,1] の係数。条件（一次）：
  (1) t = q で [3,1] 彩色 HOMFLY-PT、(2) [2,1,1] = 転置 (e_A, e_Q, e_T) ↦ (e_A, −e_T, −e_Q)、(3) 両手型 x(e) = x(−e)、
  (4) 指数増大：Q = ±1 で r31 = r11·r1·r1、T = ±1 で r31 = r3·r1、(5) U(2)：A = ±T² で r31 = r2·r11【仮説、T(2,3), T(2,5), T(3,4) で厳密に成立】、
  (6) \\hat f_ρ（ρ ⊢ 4）が Laurent、(7) G 展開が閉じる、(8) 正値性（G 展開の前後、全体の符号 σ）。
\\hat f_ρ = Σ_{λ∈{[3,1],[2,1,1]}} α_{λρ}·unk_λ·r_λ + K_ρ。D·α·unk は sympy で、D·K_ρ は評価点で求めて補間する。
使い方：python scripts/f41_box4.py [check|lp0|lp|milp] [σ]（lp0 は G 展開の前の正値性だけ）
"""
import json
import pickle
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import scipy.sparse as sps
import sympy as sp
from scipy.optimize import Bounds, LinearConstraint, linprog, milp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import highbox_torus as H  # noqa: E402
from rcsknots.bps import g_lambda, _P_sym, Q, T  # noqa: E402
from rcsknots.bpsgen import B_class  # noqa: E402
from rcsknots.characters import chi  # noqa: E402
from rcsknots.field import V  # noqa: E402
from rcsknots.laurent import load_table  # noqa: E402
from rcsknots.partitions import arm, boxes, leg, partitions  # noqa: E402

A = sp.symbols("A")
CACHE = Path("data/f41_box4")
LOG = Path("logs/f41_box4.log")
UNK = [(3, 1), (2, 1, 1)]
N = 4


def tr(P):
    return {(a, -t, -q): c for (a, q, t), c in P.items()}


def mul(*ps):
    out = {(0, 0, 0): 1}
    for p in ps:
        nx = defaultdict(int)
        for k1, c1 in out.items():
            for k2, c2 in p.items():
                nx[tuple(x + y for x, y in zip(k1, k2))] += c1 * c2
        out = {k: v for k, v in nx.items() if v}
    return out


def inputs(log=print):
    path = CACHE / "inputs.pkl"
    if path.exists():
        return pickle.loads(path.read_bytes())
    from x1_check import triply, bold_dict
    from t4_colored import bold_to_small
    from t3_pretzel import tc_graded
    nb = load_table("data/reference/notebook_torus.txt")
    red = {}
    for lam in [(1,), (2,), (1, 1), (3,), (1, 1, 1), (4,), (2, 2), (1, 1, 1, 1)]:
        red[lam] = bold_to_small(bold_dict(triply(lam[0], len(lam), -1, 1)))
        c = "".join(map(str, lam))
        if ("4_1", c) in nb:
            assert nb[("4_1", c)] == red[lam], lam
            log(f"  {lam}：KNTZ = KN のノートブック")
        if len(lam) == 1:
            assert bold_to_small(tc_graded(1, -1, lam[0])) == red[lam], lam
            log(f"  {lam}：KNTZ = GNSSS")
    red[(2, 1)] = nb[("4_1", "21")]
    h31 = homfly31()
    out = (red, h31)
    CACHE.mkdir(parents=True, exist_ok=True)
    path.write_bytes(pickle.dumps(out))
    return out


def homfly31():
    """1612.00422 の [3,1] の微分展開の 4_1 での値 {(e_A, e_Q): c}（文献の (A, q) → KN の A² = a, q = Q）。"""
    from rcsknots.double_twist import Dn as D, qnum, br, q, A as LA
    Z1 = D(5) * D(-3) + D(3) * D(1) + D(2) * D(-2) + D(0) * D(-2)
    Z2 = D(3) * D(0) * (D(5) * D(-2) + qnum(6) / qnum(3) * D(4) * D(-2) + D(2) * D(-1))
    Z11 = D(3) * D(0) * D(-2) * D(-3)
    Z3 = D(5) * D(4) * D(3) * D(1) * D(0) * D(-2)
    Z21 = D(4) * D(3) * D(1) * D(-1) * D(-2) * D(-3)
    Z31 = D(5) * D(4) * D(3) * D(2) * D(0) * D(-1) * D(-2) * D(-3)
    tZ21 = -qnum(4) ** 2 * qnum(2) * br(q) ** 4 * D(3) * D(-2)
    tZ31 = -qnum(4) ** 2 * qnum(2) * br(q) ** 4 * D(4) * D(3) * D(0) * D(-2)
    Hs = 1 + Z1 + Z2 + qnum(4) / qnum(2) * Z11 + qnum(4) / qnum(3) * Z3 + qnum(4) * qnum(2) / qnum(3) * Z21 + Z31 + tZ21 + tZ31
    e = sp.expand(sp.cancel(sp.together(Hs)) * LA ** 40 * q ** 80)
    return {(i - 40, j - 80): int(c) for (i, j), c in sp.Poly(e, LA, q).terms()}


def t_eq_q(P):
    out = defaultdict(int)
    for (a, q, t), c in P.items():
        out[(a, q + t)] += c
    return {k: v for k, v in out.items() if v}


def spec(P, which, s):
    out = defaultdict(int)
    for (a, q, t), c in P.items():
        if which == "Q":
            out[(a, t)] += c * s ** q
        elif which == "T":
            out[(a, q)] += c * s ** t
        else:
            out[(q, t + 2 * a)] += c * s ** a
    return {k: v for k, v in out.items() if v}


# ---- \hat f の一次の部分 -------------------------------------------------------------------------

def _poly(e):
    e = sp.expand(e * (A * Q * T) ** 200)
    return {(i - 200, j - 200, k - 200): int(c) for (i, j, k), c in sp.Poly(e, A, Q, T).terms()}


def alpha():
    """({(λ, ρ): D·α_{λρ}·unk_λ}, D, D の文字列)。α_{λρ} = ∂\\hat f_ρ/∂\\bar{rCS}_λ = Σ_κ χ_ρ(κ) (Q − 1/Q)(T/Q)^4 g_λ [p_κ]P_λ / B(κ)。

    有理関数体 Q(A, Q, T)（sympy の sparse な field）で計算する。D は分母の最小公倍数から単項式の因子を除いたもの。
    """
    path = CACHE / "alpha.pkl"
    if path.exists():
        return pickle.loads(path.read_bytes())
    from fractions import Fraction
    from rcsknots.macdonald import macdonald_P, macdonald_norm
    K, a_, Q_, T_ = sp.field("A,Q,T", sp.QQ)
    conv = lambda c: K(sp.Rational(c.numerator, c.denominator))  # noqa: E731
    q_, t_ = Q_ ** 2, T_ ** 2
    P = macdonald_P(N, q_, t_, conv)
    raw = {}
    for lam in UNK:
        unk = (T_ / a_) ** N
        for i, j in boxes(lam):
            unk = unk * (T_ ** (2 * i) - a_ ** 2 * Q_ ** (2 * j)) / (1 - Q_ ** (2 * arm(lam, i, j)) * T_ ** (2 * (leg(lam, i, j) + 1)))
        g = (q_ / t_) ** N / macdonald_norm(lam, q_, t_)
        for rho in partitions(N):
            s = K(0)
            for k in partitions(N):
                c = P[lam].get(k)
                if c is not None and chi(rho, k):
                    s = s + c * chi(rho, k) / B_class(k, T_)
            raw[(lam, rho)] = s * (Q_ - 1 / Q_) * (T_ / Q_) ** N * g * unk
    den = None
    for v in raw.values():
        den = v.denom if den is None else den.lcm(v.denom)
    _, fl = den.factor_list()
    den = K.ring(1)
    for f, k in fl:
        if len(f.terms()) > 1:                     # 単項式の因子は Laurent 多項式の中で吸収する
            den = den * f ** k
    out = {}
    for key, v in raw.items():
        w = v * den
        num, dd = w.numer, w.denom
        assert len(dd.terms()) == 1, "D で割り切れない"
        (mono, c), = dd.terms()
        out[key] = {tuple(e - m for e, m in zip(ex, mono)): int(sp.Rational(cf) / c) for ex, cf in num.terms()}
    Dp = {ex: int(c) for ex, c in den.terms()}
    res = (out, Dp, str(den.as_expr().factor()))
    path.write_bytes(pickle.dumps(res))
    return res


def dK(red, Dp, log=print, workers=8, path=CACHE / "dK.pkl"):
    """{ρ: D·K_ρ}（r31 = r211 = 0 とした \\hat f_ρ に D を掛けたもの、厳密）。"""
    if path is not None and path.exists():
        return pickle.loads(path.read_bytes())
    red0 = dict(red)
    for lam in UNK:
        red0[lam] = {}

    def factory():
        def rcs(lam, Qv, Tv):
            P = red0[lam]
            if not P:
                return 0, [Qv * 0]
            lo, hi = min(k[0] for k in P), max(k[0] for k in P)
            cs = [Qv * 0 for _ in range(hi - lo + 1)]
            for (ea, eq, et), c in P.items():
                cs[ea - lo] = cs[ea - lo] + Qv ** eq * Tv ** et * c
            return lo, cs
        return rcs

    def post(Qs, Ts, p):
        Qv, Tv = V(Qs, p), V(Ts, p)
        out = Qv * 0
        for (ea, eq, et), c in Dp.items():
            assert ea == 0
            out = out + Qv ** eq * Tv ** et * c
        return out.a

    prob = H.Problem(N, factory, -24, 24, post=post)
    hats, sec = H.run_problem(prob, 41, workers, log)
    log(f"  D·K：{sum(len(h) for h in hats.values())} 項、{sec:.0f}s")
    if path is not None:
        path.write_bytes(pickle.dumps(hats))
    return hats


def window(prods, h31, margin):
    """x = rCS_[3,1] の窓。各軸の範囲を特殊化の値から決める：e_A は 3 つの積の A の範囲の共通部分、
    e_Q は T = ±1 の積 r3·r1 の Q の範囲、e_T は Q = ±1 の積 r11·r1·r1 の T の範囲、e_Q + e_T は t = q の値の範囲（各 ±margin）。
    （特殊化で打ち消し合う項が窓の外にある解は除かれる。T(2,3), T(2,5), T(3,4) の真の値が窓に入ることはテストで確かめる。）"""
    pQ, pT, pU = prods
    arng = (max(min(e[0] for e in P) for P in prods) - margin, min(max(e[0] for e in P) for P in prods) + margin)
    qrng = (min(e[1] for e in pT) - margin, max(e[1] for e in pT) + margin)
    trng = (min(e[2] for e in pQ) - margin, max(e[2] for e in pQ) + margin)
    srng = (min(k[1] for k in h31) - margin, max(k[1] for k in h31) + margin)
    par = next(iter(pQ))
    return [(ea, eq, et) for ea in range(arng[0], arng[1] + 1) if (ea - par[0]) % 2 == 0
            for eq in range(qrng[0], qrng[1] + 1) for et in range(trng[0], trng[1] + 1)
            if srng[0] <= eq + et <= srng[1] and (eq + et - par[1] - par[2]) % 2 == 0]


# ---- 一次の系 -------------------------------------------------------------------------------------

def build(red, h31, al, Dp, dk, sigma=1, positive=True, margin=2, log=print, amphichiral=True, with_G=True):
    r1, r2, r11, r3 = red[(1,)], red[(2,)], red[(1, 1)], red[(3,)]
    prods = [mul(r11, r1, r1), mul(r3, r1), mul(r2, r11)]
    xs = window(prods, h31, margin)
    xcol = {e: n for n, e in enumerate(xs)}
    ncols = len(xs)
    lb, ub = [-np.inf] * ncols, [np.inf] * ncols
    rows, rhs = [], []

    def eq_row(row, v):
        rows.append(row)
        rhs.append(v)

    # (1) t = q
    tq = defaultdict(dict)
    for e, n in xcol.items():
        tq[(e[0], e[1] + e[2])][n] = 1
    for key in set(tq) | set(h31):
        eq_row(tq.get(key, {}), h31.get(key, 0))
    # (3) 両手型
    for e, n in (xcol.items() if amphichiral else []):
        m = (-e[0], -e[1], -e[2])
        if m > e:
            eq_row({n: 1, xcol[m]: -1}, 0)
    # (4)(5) 特殊化
    for which, prod in (("Q", prods[0]), ("T", prods[1]), ("U", prods[2])):
        for s in (1, -1):
            ev = defaultdict(dict)
            for e, n in xcol.items():
                key = spec({e: 1}, which, s)
                (k, c), = key.items()
                ev[k][n] = ev[k].get(n, 0) + c
            tgt = spec(prod, which, s)
            for key in set(ev) | set(tgt):
                eq_row(ev.get(key, {}), tgt.get(key, 0))
    # (6)(7)(8)：D·h_ρ = Σ_λ (D α unk)_λρ · x_λ + D·K_ρ、h_ρ の太字の係数 ≥ 0、G 展開
    dlo = [min(e[i] for e in Dp) for i in range(3)]
    dhi = [max(e[i] for e in Dp) for i in range(3)]
    for rho in partitions(N):
        num = defaultdict(dict)
        for e, n in xcol.items():
            for lam, ee in ((UNK[0], e), (UNK[1], (e[0], -e[2], -e[1]))):
                for k, c in al[(lam, rho)].items():
                    kk = tuple(x + y for x, y in zip(ee, k))
                    num[kk][n] = num[kk].get(n, 0) + c
        keys = set(num) | set(dk[rho])
        lo_ = [min(e[i] for e in keys) - dlo[i] for i in range(3)]
        hi_ = [max(e[i] for e in keys) - dhi[i] for i in range(3)]
        hbox = [(a_, q_, t_) for a_ in range(lo_[0], hi_[0] + 1) for q_ in range(lo_[1], hi_[1] + 1)
                for t_ in range(lo_[2], hi_[2] + 1) if a_ % 2 == 0 and (q_ + t_) % 2 == 0]
        hcol = {e: ncols + n for n, e in enumerate(hbox)}
        ncols += len(hbox)
        lb += [0.0 if positive else -np.inf] * len(hbox)
        ub += [np.inf] * len(hbox)
        sgn = {e: (-1) ** (((e[0] // 2) + e[1]) % 2) for e in hbox}      # to_bold の符号（a 次数のずらし 0）
        conv = defaultdict(dict)
        for e in hbox:
            for de, dc in Dp.items():
                conv[tuple(x + y for x, y in zip(e, de))][hcol[e]] = dc * sgn[e]
        for key in set(conv) | keys:
            row = dict(conv.get(key, {}))
            for n, c in num.get(key, {}).items():
                row[n] = row.get(n, 0) - sigma * c
            eq_row(row, sigma * dk[rho].get(key, 0))
        if not with_G:
            continue
        bold = {}
        for e in hbox:
            m = e[0] // 2
            bold[(2 * m, e[1] + e[2], m + e[1])] = hcol[e]
        betas = sorted({bb[0] // 2 for bb in bold})
        gmax = max(abs(bb[1]) for bb in bold) // 2 + 1
        tlo, thi = min(bb[2] for bb in bold), max(bb[2] for bb in bold)
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
    Amat = sps.csr_matrix((data, (ri, ci)), shape=(len(rows), ncols), dtype=float)
    log(f"  x {len(xs)}、行 {len(rows)}、列 {ncols}（G 展開 {'あり' if with_G else 'なし'}）")
    return xs, Amat, np.array(rhs, dtype=float), np.array(lb), np.array(ub)


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "lp"
    LOG.parent.mkdir(exist_ok=True)
    lg = LOG.open("a")

    def log(s):
        print(s, flush=True)
        lg.write(s + "\n")
        lg.flush()

    log(f"== f41_box4 {mode} {time.strftime('%Y-%m-%d %H:%M:%S')}")
    t0 = time.time()
    red, h31 = inputs(log)
    # 入力の検査：[3,1] の t = q の値の性質、既知の彩色の指数増大
    r1 = red[(1,)]
    q1 = defaultdict(int)
    for (a, q), c in h31.items():
        q1[a] += c
    p1 = defaultdict(int)
    for (a, q, t), c in mul(r1, r1, r1, r1).items():
        p1[a] += c
    log(f"  [3,1] の t = q：{len(h31)} 項、q = 1 で r1⁴ と一致 {dict(q1) == {k: v for k, v in p1.items() if v}}、"
        f"両手型 {h31 == {(-a, -q): c for (a, q), c in h31.items()}}")
    for lam in [(4,), (2, 2), (1, 1, 1, 1)]:
        log(f"  {lam}：転置で閉じる {tr(red[lam]) == red[{(4,): (1, 1, 1, 1), (2, 2): (2, 2), (1, 1, 1, 1): (4,)}[lam]]}、"
            f"両手型 {red[lam] == {(-a, -q, -t): c for (a, q, t), c in red[lam].items()}}")
    log(f"  [2,1]（KN の候補）：U(2) r21 = r1·r11 {all(spec(red[(2, 1)], 'U', s) == spec(mul(r1, red[(1, 1)]), 'U', s) for s in (1, -1))}")
    al, Dp, dstr = alpha()
    log(f"  D = {dstr}（{time.time() - t0:.0f}s）")
    dk = dK(red, Dp, log)
    check_linear(red, al, Dp, dk, log)
    if mode == "check":
        return
    sigmas = [int(sys.argv[2])] if len(sys.argv) > 2 else [1, -1]
    for sigma in sigmas:
        xs, Amat, b, lb, ub = build(red, h31, al, Dp, dk, sigma=sigma, log=log, margin=0, with_G=(mode not in ("lp0", "eq", "milp0")),
                                    positive=(mode != "eq"))
        t1 = time.time()
        if mode in ("lp", "lp0", "eq"):
            res = linprog(np.zeros(Amat.shape[1]), A_eq=Amat, b_eq=b, bounds=list(zip(lb, ub)), method="highs",
                          options={"time_limit": 3000})
            log(f"  σ = {sigma:+d}：{mode} {res.status} {res.message}（{time.time() - t1:.0f}s）")
            if res.x is not None:
                sol = {e: float(v) for e, v in zip(xs, res.x[:len(xs)]) if abs(v) > 1e-9}
                (CACHE / f"lp_{mode}_sigma{sigma:+d}.json").write_text(json.dumps({",".join(map(str, k)): v for k, v in sol.items()}))
        else:
            res = milp(np.zeros(Amat.shape[1]), constraints=LinearConstraint(Amat, b, b), integrality=np.ones(Amat.shape[1]),
                       bounds=Bounds(lb, ub), options={"time_limit": 3000})
            log(f"  σ = {sigma:+d}：MILP {res.status} {res.message}（{time.time() - t1:.0f}s）")
            if res.x is not None:
                sol = {e: int(round(v)) for e, v in zip(xs, res.x[:len(xs)]) if round(v)}
                (CACHE / f"sol_{mode}_sigma{sigma:+d}.json").write_text(json.dumps({",".join(map(str, k)): v for k, v in sol.items()}))
                log(f"  解：{len(sol)} 項")


def check_linear(red, al, Dp, dk, log=print):
    """分解 D·\\hat f_ρ = Σ_λ (D α unk)_λρ·r_λ + D·K_ρ を、乱数の r31（r211 はその転置）について評価点で確かめる。"""
    from rcsknots.bpsgen import Points
    from rcsknots.field import P31
    rng = np.random.default_rng(5)
    X = {(2 * int(rng.integers(-3, 4)), int(e), int(e) + 2 * int(rng.integers(-2, 3))): int(rng.integers(-5, 6)) for e in rng.integers(-4, 5, 6)}
    X = {k: v for k, v in X.items() if v}
    full = dict(red)
    full[UNK[0]], full[UNK[1]] = X, tr(X)
    p, k = P31, 3
    Qv, Tv = V(rng.integers(2, 10**6, k), p), V(rng.integers(2, 10**6, k), p)
    Av = [int(x) for x in rng.integers(2, 10**6, 2)]

    def rcs(lam, Q_, T_):
        P = full[lam]
        lo, hi = min(e[0] for e in P), max(e[0] for e in P)
        cs = [Q_ * 0 for _ in range(hi - lo + 1)]
        for (ea, eq, et), c in P.items():
            cs[ea - lo] = cs[ea - lo] + Q_ ** eq * T_ ** et * c
        return lo, cs
    hat, ok = Points(rcs, Qv, Tv, Av).hat_f(N)
    Aa = V(np.tile(Av, k), p)
    Qa, Ta = V(np.repeat(Qv.a, 2), p), V(np.repeat(Tv.a, 2), p)

    def ev(P):
        out = Qa * 0
        for (ea, eq, et), c in P.items():
            out = out + Aa ** ea * Qa ** eq * Ta ** et * c
        return out
    D = ev(Dp)
    good = ok
    for rho in partitions(N):
        rhs = ev(dk[rho]) + ev(mul(al[(UNK[0], rho)], X)) + ev(mul(al[(UNK[1], rho)], tr(X)))
        good &= (hat[N][rho] * D == rhs)
    log(f"  一次の分解の検査（乱数の r31、評価点 {k}×2）：{good}")
    return good


if __name__ == "__main__":
    main()
