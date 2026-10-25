r"""K4'(a)：E（Ĩ で固定される oriented 状態の符号付き数）を固定セクターの状態の部分和として決める。
  g_{ab} = ½E_{ab}/(1+t) + C_{ab}（全 d、report.md §7）。被覆をセクター別に分けると（kauffman_cover.cover_f_sym4）、
  Ĩ で固定される状態は (k,l,m,n) = (a/2,a/2,b/2,b/2) のセクターにしかないので、P := (1−t) f^{cov}_{固定セクター} を状態の和とみて
  E = P − 2Q（Q：Ĩ の固有値が −1 の状態の和、0 ≤ Q_m/P_m ≤ 1 の整数）と書ける。条件：
   (i)  (1+t) | (R − E)、R = 2(1+t)g_{ab}（⇔ t = −1 で E = R）、
   (ii) t = q で C = ĝ_{ab}（Mariño の unrefined の交差帽、= Σ_{d 奇数} の定義 g^{odd} の t = q）、⇔ E|_{t=q} = R|_{t=q} − 2(1+q) g^{odd}_{ab}|_{t=q}。
  Q を整数線形計画（scipy.optimize.milp）で求め、各変数の最小・最大で一意性を判定する。
入力：data/refined_kauffman/hopf_fcov_sector_A4.pkl（cover_f_sym4(4,4)、約 6 分）、hopf_g_alld.pkl、hopf_g_odd.pkl。
出力：data/refined_kauffman/hopf_E_fixed.pkl。"""
import sys, pickle
from collections import defaultdict
import numpy as np
import sympy as sp
from scipy.optimize import milp, LinearConstraint, Bounds

r, s, nu = sp.symbols("r s nu", positive=True)
nut = sp.Symbol("nut")
t = s**2
loc = {"r": r, "s": s, "nu": nu}
f4 = {k: sp.sympify(v, locals=loc) for k, v in pickle.load(open("data/refined_kauffman/hopf_fcov_sector_A4.pkl", "rb")).items()}
g = {k: sp.sympify(v, locals=loc) for k, v in pickle.load(open("data/refined_kauffman/hopf_g_alld.pkl", "rb")).items()}
god = {k: sp.sympify(v, locals=loc) for k, v in pickle.load(open("data/refined_kauffman/hopf_g_odd.pkl", "rb")).items()}
SH = r**40 * s**60 * nu**16


def lp(e, gens=(r, s, nu)):
    n, d = sp.fraction(sp.cancel(sp.together(e)))
    sh = SH if len(gens) == 3 else r**40 * nu**16
    Pn = sp.Poly(sp.expand(n * sh), *gens); Pd = sp.Poly(sp.expand(d * sh), *gens)
    assert Pd.is_monomial, d
    (md, cd), = Pd.terms()
    return {tuple(mi - di for mi, di in zip(m, md)): sp.Rational(c, cd) for m, c in Pn.terms()}


def solve(a, b):
    sec = (a // 2, a // 2, b // 2, b // 2)
    Eq_expr = sp.cancel((2 * (1 + t) * g[(a, b)] - 2 * (1 + t) * god[(a, b)]).subs(s, r))
    out = solve_generic((1 - t) * f4[sec], 2 * (1 + t) * g[(a, b)], Eq_expr)
    out["sec"] = sec
    return out


def solve_generic(P_expr, R_expr, Eq_expr):
    """E = P − 2Q（Q は P の部分和）、t = −1 で E = R、t = q で E = Eq を整数線形計画で解き、一意性を調べる。"""
    P = lp(P_expr); R = lp(R_expr)
    mons = sorted(P)
    idx = {m: i for i, m in enumerate(mons)}
    S = defaultdict(lambda: sp.Integer(0))
    for m, c in R.items(): S[m] += c / 2
    for m, c in P.items(): S[m] -= c / 2
    rows, rhs = [], []
    # (i) s = i：(a, c) ごとに Σ_b Q_{a,b,c} i^b = −S(i)_{a,c}（実部・虚部）
    Si = defaultdict(lambda: 0)
    for (x, y, z), v in S.items(): Si[(x, z)] += v * sp.I**y
    keys = set(Si) | {(x, z) for (x, y, z) in mons}
    for key in keys:
        target = -sp.expand(Si.get(key, 0))
        for part, f in (("re", sp.re), ("im", sp.im)):
            row = np.zeros(len(mons))
            for m in mons:
                if (m[0], m[2]) == key:
                    row[idx[m]] = float(f(sp.I**m[1]))
            tv = float(f(target))
            if np.any(row) or tv != 0:
                rows.append(row); rhs.append(tv)
    # (ii) s = r：E|_{t=q} = R|_{t=q} − 2(1+q) g^{odd}|_{t=q} ⇔ Q|_{t=q} = (P − E)|_{t=q}/2
    Eq = lp(Eq_expr, (r, nu))
    Pq = defaultdict(lambda: sp.Integer(0))
    for (x, y, z), v in P.items(): Pq[(x + y, z)] += v
    Qq = {k: (Pq.get(k, 0) - Eq.get(k, 0)) / 2 for k in set(Pq) | set(Eq)}
    assert all(v.is_integer for v in Qq.values()), "t = q の条件が整数にならない"
    for key, v in Qq.items():
        row = np.zeros(len(mons))
        for m in mons:
            if (m[0] + m[1], m[2]) == key:
                row[idx[m]] = 1
        if np.any(row) or v != 0:
            rows.append(row); rhs.append(float(v))
    A = np.array(rows); bvec = np.array(rhs)
    lo = np.array([min(0, float(P[m])) for m in mons]); hi = np.array([max(0, float(P[m])) for m in mons])
    cons = LinearConstraint(A, bvec, bvec)
    integ = np.ones(len(mons))
    res = milp(np.zeros(len(mons)), constraints=cons, bounds=Bounds(lo, hi), integrality=integ)
    if res.status != 0:
        return dict(feasible=False, n=len(mons))
    # 一意性：各変数の最小と最大
    free = []
    for i in range(len(mons)):
        c = np.zeros(len(mons)); c[i] = 1
        mn = milp(c, constraints=cons, bounds=Bounds(lo, hi), integrality=integ).fun
        mx = -milp(-c, constraints=cons, bounds=Bounds(lo, hi), integrality=integ).fun
        if round(mx) != round(mn):
            free.append((mons[i], round(mn), round(mx)))
    Qsol = {m: int(round(res.x[idx[m]])) for m in mons}
    E = sum(P[m] * r**m[0] * s**m[1] * nu**m[2] for m in mons) - 2 * sum(Qsol[m] * r**m[0] * s**m[1] * nu**m[2] for m in mons)
    return dict(feasible=True, n=len(mons), free=free, E=sp.expand(E))


if __name__ == "__main__":
    out = {}
    z, zp = r - 1 / r, r + 1 / r
    for ab in ((2, 2), (2, 4), (4, 4)):
        res = solve(*ab)
        print(f"{ab}: 固定セクター {res['sec']}、状態の単項式 {res['n']} 個、解あり = {res['feasible']}", flush=True)
        if not res["feasible"]:
            continue
        print(f"   一意でない変数: {len(res['free'])} 個 {res['free'][:6]}")
        E = res["E"]; out[ab] = str(E)
        Ew = sp.factor(sp.cancel(E.subs(nu, nut * s / r)))
        C = sp.cancel(g[ab] - E / (2 * (1 + t)))
        Cw = sp.cancel(C.subs(nu, nut * s / r))
        print("   E (ν~) =", Ew if len(str(Ew)) < 600 else str(Ew)[:600] + " …")
        print("   C は ν~ 固定で t 非依存:", not Cw.has(s), " C の項数:", len(sp.Poly(sp.expand(sp.fraction(Cw)[0]), r, s, nut).terms()))
        if ab == (2, 2):
            EB = -z * zp**2 * (nut - 1 / nut) * (t / r * nut + r / nut)
            print("   E_22 == 候補 B −z(z+)²(ν~−ν~⁻¹)(t q^{−1/2}ν~ + q^{1/2}ν~⁻¹):", sp.simplify(Ew - EB) == 0)
            E22 = E
        if ab == (2, 4):
            chi = (nu * r / s)**2 / r**2 + r**2 / (nu * r / s)**2
            print("   E_24 == E_22 · χ_1(ν~²/q):", sp.expand(E - E22 * chi) == 0)
    pickle.dump(out, open("data/refined_kauffman/hopf_E_fixed.pkl", "wb"))
