r"""固定セクターの多重項ごとの Ĩ の捻ったトレースで E を作る（progress.md 続き 39、theory_orientifold_index.md §9）。
  固定セクター（被覆の brane–image の状態）P = (1 − t)x = s z_q Pk、Pk（σ 対称）を G 展開して
  Pk = Σ_m c_m b_m、b_m = ν^β G^g W^j（多重項。G = z_q z_t、W = (q/t)^{1/4}）。
  Ĩ = I_X × (平面の π 回転) の捻ったトレースは、多重項の指標の変数の代入で与えられるとする【仮説】：
    捻り A（q 平面）：r → i r（z_q → i z_q^+）、捻り B（t 平面）：s → i s（z_t → i z_t^+）。
  b_m の捻った指標 b̃_m（位相を除いた実の形）をとり、E = s z_q Σ_m e_m b̃_m とする。
  e_m は不変な多重項の数に符号 ε を付けたもの：|e_m| ≤ |c_m|、e_m ≡ c_m (mod 2)（対は 2 つずつ）。
  条件：Hopf では E が既知（scripts/rk_k4_fixed_sector.py で一意に決まった値）と一致すること。
       結び目では (i) t = −1 で E = R = 2(1 + t)g_(2)、(ii) t = q で E = −(1 + q)(g_(1)(x²) + 1)。
  これを整数計画（scipy.optimize.milp）で解き、各 e_m の範囲で一意性を調べる。"""
import sys, io, contextlib, pickle
from collections import defaultdict
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import numpy as np
import sympy as sp
from scipy.optimize import milp, LinearConstraint, Bounds
from rk_k4_crosscap import V, W, G, in_G

r, s = sp.symbols("r s", positive=True)
nu = sp.Symbol("nu", positive=True)
nuk = sp.Symbol("nu")
I = sp.I
t = s**2


def multiplets(Pk):
    """Pk を G 展開し、{(β, g, j): c}。"""
    ev = sp.expand(sp.cancel(Pk.subs(nu, nuk)).subs({r: V * W, s: V / W}))
    P1 = in_G(ev)
    assert P1 is not None, "Pk は σ 対称でない"
    out = defaultdict(int)
    for tm in sp.Add.make_args(sp.expand(P1)):
        c, pw = tm.as_coeff_mul()
        d = sp.Mul(*pw).as_powers_dict()
        key = (int(d.get(nuk, 0)), int(d.get(G, 0)), int(d.get(W, 0)))
        out[key] += c
    assert all(k[2] % 2 == 0 for k in out), "W の奇数冪"
    return {k: v for k, v in out.items() if v != 0}


def monomials(Pk):
    """Pk を単項式に分け、{("mono", a, b, c): 係数}（r^a s^b ν^c）。"""
    return {("mono",) + k: v for k, v in lp(Pk, (r, s, nu)).items()}


def twisted(key, twist):
    if key[0] == "mono":
        return r**key[1] * s**key[2] * nu**key[3]
    beta, g, j = key
    zq, zt, zqp, ztp = r - 1 / r, s - 1 / s, r + 1 / r, s + 1 / s
    Gt = zqp * zt if twist == "A" else zq * ztp if twist == "B" else zq * zt
    return nu**beta * Gt**g * (r / s)**(j // 2)


def lp(e, gens):
    n, d = sp.fraction(sp.cancel(sp.together(e)))
    Pn, Pd = sp.Poly(sp.expand(n), *gens), sp.Poly(sp.expand(d), *gens)
    assert Pd.is_monomial, d
    (md, cd), = Pd.terms()
    return {tuple(a - b for a, b in zip(m, md)): sp.nsimplify(c / cd) for m, c in Pn.terms()}


def rows_for(basis, target, gens_sub):
    """basis の各関数と target を、代入 gens_sub の後の単項式ごとの線形条件（実部・虚部）にする。"""
    cols = [lp(sp.expand(b.subs(gens_sub)), (r, s, nu)) for b in basis]
    T = lp(sp.expand(target.subs(gens_sub)), (r, s, nu)) if target != 0 else {}
    mons = set(T).union(*cols)
    A, bvec = [], []
    for mo in mons:
        for f in (sp.re, sp.im):
            row = [float(f(c.get(mo, 0))) for c in cols]
            v = float(f(T.get(mo, 0)))
            if any(row) or v:
                A.append(row); bvec.append(v)
    return A, bvec


def solve(mult, basis_fn, conds, uniq=True):
    keys = sorted(mult)
    basis = [basis_fn(k) for k in keys]
    A, bv = [], []
    for target, sub in conds:
        a, b = rows_for(basis, target, sub)
        A += a; bv += b
    n = len(keys)
    # e_m = c_m − 2 y_m、y_m ∈ [min(0, c_m), max(0, c_m)]（|e_m| ≤ |c_m| と偶奇）
    c = np.array([float(mult[k]) for k in keys])
    A = np.array(A); bv = np.array(bv)
    A2, b2 = -2 * A, bv - A @ c
    lo = np.minimum(0, c); hi = np.maximum(0, c)
    cons = LinearConstraint(A2, b2 - 1e-9, b2 + 1e-9)
    res = milp(np.zeros(n), constraints=cons, bounds=Bounds(lo, hi), integrality=np.ones(n))
    if res.status != 0:
        return dict(feasible=False, n=n, rows=len(bv))
    y = np.round(res.x)
    free = []
    if uniq:
        for i in range(n):
            cc = np.zeros(n); cc[i] = 1
            mn = milp(cc, constraints=cons, bounds=Bounds(lo, hi), integrality=np.ones(n)).fun
            mx = -milp(-cc, constraints=cons, bounds=Bounds(lo, hi), integrality=np.ones(n)).fun
            if round(mx) != round(mn):
                free.append(keys[i])
    e = {k: int(c[i] - 2 * y[i]) for i, k in enumerate(keys)}
    return dict(feasible=True, n=n, rows=len(bv), free=free, e=e)


def E_of(e, twist):
    return sp.expand(s * (r - 1 / r) * sum(v * twisted(k, twist) for k, v in e.items()))


def hopf():
    import rk_k4_fixed_sector as H
    Eknown = {eval(k) if isinstance(k, str) else k: sp.sympify(v, locals={"r": r, "s": s, "nu": nu})
              for k, v in pickle.load(open("data/refined_kauffman/hopf_E_fixed.pkl", "rb")).items()}
    for ab in ((2, 2), (2, 4), (4, 4)):
        sec = (ab[0] // 2, ab[0] // 2, ab[1] // 2, ab[1] // 2)
        P = sp.cancel((1 - t) * H.f4[sec])
        Pk = sp.cancel(P / (s * (r - 1 / r)))
        mult = monomials(Pk)
        for tw in ("mono",):
            res = solve(mult, lambda k: s * (r - 1 / r) * twisted(k, tw), [(Eknown[ab], {})])
            msg = f"解あり、一意でない多重項 {len(res['free'])}" if res["feasible"] else "解なし"
            if res["feasible"]:
                allinv = all(abs(v) == abs(mult[k]) for k, v in res["e"].items())
                msg += f"、全多重項が不変（|e| = |c|）{allinv}"
            print(f"Hopf {ab}：多重項 {len(mult)}、捻り {tw}：{msg}", flush=True)
        R = sp.cancel(2 * (1 + t) * H.g[ab])
        Eq = sp.cancel((2 * (1 + t) * H.g[ab] - 2 * (1 + t) * H.god[ab]).subs(s, r))
        res = solve(mult, lambda k: s * (r - 1 / r) * twisted(k, "mono"), [(R, {s: I}), (Eq, {s: r})])
        msg = f"解あり、一意でない多重項 {len(res['free'])}/{res['n']}" if res["feasible"] else "解なし"
        print(f"   条件 (i)+(ii)：{msg}", flush=True)


def knot(m, twists=("A", "B", "mono")):
    import rk_b_knot_fixed_sector as K
    import rk_b_cover_split as C
    u = sp.Symbol("u", positive=True)
    g = K.load_g(m)
    with contextlib.redirect_stdout(io.StringIO()):
        sec = C.cover_sectors(m)
    x = sp.cancel(sec["x"].to_sympy().subs(u, nu * s))
    P = sp.cancel((1 - t) * x)
    Pk = sp.cancel(P / (s * (r - 1 / r)))
    multG = multiplets(Pk)
    multM = monomials(Pk)
    R = sp.cancel(2 * (1 + t) * g[(2,)])
    g1x2 = g[(1,)].subs({r: r**2, s: s**2, nu: nu**2}, simultaneous=True)
    Eq = sp.cancel(-(1 + r**2) * (g1x2 + 1))
    print(f"T(2,{2*m+1})：G の多重項 {len(multG)}（係数が負 {sum(1 for v in multG.values() if v < 0)}）、Pk の単項式 {len(multM)}", flush=True)
    out = {}
    for tw, pre in [(tw, pre) for tw in twists for pre in ("z_q", "z_q^+")]:
        mult = multM if tw == "mono" else multG
        pf = (r - 1 / r) if pre == "z_q" else (r + 1 / r)
        for lab, conds in (("(i)", [(R, {s: I})]), ("(ii)", [(Eq, {s: r})]), ("(i)+(ii)", [(R, {s: I}), (Eq, {s: r})])):
            res = solve(mult, lambda k: s * pf * twisted(k, tw), conds, uniq=(lab == "(i)+(ii)"))
            msg = (f"解あり" + (f"、一意でない多重項 {len(res['free'])}/{res['n']}" if lab == "(i)+(ii)" else "")) if res["feasible"] else "解なし"
            print(f"   捻り {tw}・前の因子 {pre} {lab}：{msg}（条件 {res['rows']} 行）", flush=True)
            if res["feasible"] and lab == "(i)+(ii)":
                out[(tw, pre)] = res
    return out


if __name__ == "__main__":
    if sys.argv[1] == "hopf":
        hopf()
    else:
        res = knot(int(sys.argv[1]))
        pickle.dump({k: {"e": {str(kk): vv for kk, vv in v["e"].items()}, "free": [str(f) for f in v["free"]]} for k, v in res.items()},
                    open(f"data/refined_kauffman/daha/twisted_multiplet_m{sys.argv[1]}.pkl", "wb"))


def hidden_pairs(m, tw="mono", pre="z_q^+"):
    """|e_m| ≤ |c_m| + 2h_m（h_m は指数の打ち消しで見えない Ĩ 不変な対の数）に緩め、(i)+(ii) のもとで Σ h_m の最小値を求める。"""
    import rk_b_knot_fixed_sector as K
    import rk_b_cover_split as C
    u = sp.Symbol("u", positive=True)
    g = K.load_g(m)
    with contextlib.redirect_stdout(io.StringIO()):
        sec = C.cover_sectors(m)
    x = sp.cancel(sec["x"].to_sympy().subs(u, nu * s))
    Pk = sp.cancel((1 - t) * x / (s * (r - 1 / r)))
    mult = monomials(Pk) if tw == "mono" else multiplets(Pk)
    pf = (r - 1 / r) if pre == "z_q" else (r + 1 / r)
    R = sp.cancel(2 * (1 + t) * g[(2,)])
    g1x2 = g[(1,)].subs({r: r**2, s: s**2, nu: nu**2}, simultaneous=True)
    Eq = sp.cancel(-(1 + r**2) * (g1x2 + 1))
    keys = sorted(mult)
    basis = [s * pf * twisted(k, tw) for k in keys]
    A, bv = [], []
    for target, sub in ((R, {s: I}), (Eq, {s: r})):
        a, b = rows_for(basis, target, sub)
        A += a; bv += b
    n = len(keys)
    c = np.array([float(mult[k]) for k in keys]); ac = np.abs(c)
    A = np.array(A); bv = np.array(bv)
    # 変数 [y (n), h (n)]、e = c − 2y、−|c| − 2h ≤ e ≤ |c| + 2h
    Aeq = np.hstack([-2 * A, np.zeros_like(A)])
    beq = bv - A @ c
    I_ = np.eye(n)
    A1 = np.hstack([-2 * I_, -2 * I_]); u1 = ac - c          # e ≤ |c| + 2h
    A2 = np.hstack([2 * I_, -2 * I_]); u2 = ac + c           # −e ≤ |c| + 2h
    cons = [LinearConstraint(Aeq, beq - 1e-9, beq + 1e-9), LinearConstraint(A1, -np.inf, u1), LinearConstraint(A2, -np.inf, u2)]
    cost = np.concatenate([np.zeros(n), np.ones(n)])
    res = milp(cost, constraints=cons, bounds=Bounds(np.concatenate([-np.inf * np.ones(n), np.zeros(n)]), np.inf), integrality=np.ones(2 * n))
    if res.status != 0:
        print(f"T(2,{2*m+1}) 捻り {tw}・{pre}：緩めても解なし", flush=True); return None
    h = np.round(res.x[n:])
    print(f"T(2,{2*m+1}) 捻り {tw}・{pre}：隠れた対の最小数 Σh = {int(h.sum())}（変数 {n}、Σ|c| = {int(ac.sum())}）", flush=True)
    return int(h.sum())
