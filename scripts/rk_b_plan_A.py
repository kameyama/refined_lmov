r"""案 A（docs/refined_kauffman/proposals_SR_action.md）：固定セクター x = R_adj (qn² − κ) + κ − W₁² の Ĩ 捻りを、
unknot で決まった規則で組み立て、未知を R_adj の単項式の符号 ε_m だけにする（progress.md 続き 49）。
  x̃ = R̃_adj (ψ²qn − κ̃) + κ̃ − s_W ψ²(W₁)、E = (1 + t) x̃、
  κ̃ = (1 + A)(1 + q)/((1 + t)(1 + qA/t))（κ の三因子と z₂ の塔を捻る）、ψ² は Adams（brane と image の入れ替えのトレース）、
  R̃_adj = Σ_m ε_m c_m m（ε_m ∈ {±1}）、s_W ∈ {±1}。
  条件：(0) E が ν について Laurent ⟺ R̃_adj(ν² = −t/q) = 1、(0′) t について Laurent（t = 1 の留数が 0）、
        (i) (1 + t) | 2(1 + t) g_(2) − E、(ii) t = q で E = −(1 + q)(g_(1)(q², ν²) + 1)。
  線形条件は無作為な有理点での評価で作り、ε_m を整数計画で解く。"""
import sys, io, contextlib, pickle, random
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import numpy as np
import sympy as sp
from scipy.optimize import milp, LinearConstraint, Bounds

r, s = sp.symbols("r s", positive=True)
nu = sp.Symbol("nu", positive=True)
u = sp.Symbol("u", positive=True)
I = sp.I
t, q = s**2, r**2


def load(m):
    import rk_b_knot_fixed_sector as K
    import rk_b_cover_split as C
    import rk_c_2box_largeN as L
    g = {k: v.subs({K.r: r, K.s: s, K.nu: nu}) for k, v in K.load_g(m).items()}
    with contextlib.redirect_stdout(io.StringIO()):
        sec = C.cover_sectors(m)
    x = sp.cancel(sec["x"].to_sympy().subs(u, nu * s))
    d = pickle.load(open(f"data/refined_kauffman/daha/Radj_A_m{m}.pkl", "rb"))
    Radj = sp.expand(sum(sp.sympify(v, locals={"r": r, "s": s}) * nu**j for j, v in d.items()))
    W1 = sp.cancel(((s / r)**m * L.Rdata_nu("1", f"T2_{2*m+1}") * L.Psym((1,))).to_sympy().subs(u, nu * s))
    A = nu**2
    qn = (nu - 1 / nu) / (s - 1 / s)
    kappa = (1 - A) * (1 - q) / ((1 - t) * (1 - q * A / t))
    chk = sp.cancel(Radj * (qn**2 - kappa) + kappa - W1**2 - x)
    assert chk == 0, "x の分解が合わない"
    return g, x, Radj, W1, qn, kappa


def adams(e):
    return e.subs({r: r**2, s: s**2, nu: nu**2}, simultaneous=True)


def build(m, sW):
    g, x, Radj, W1, qn, kappa = load(m)
    A = nu**2
    kt = (1 + A) * (1 + q) / ((1 + t) * (1 + q * A / t))
    terms = sp.Poly(sp.expand(sp.cancel(Radj * (r * s * nu)**60)), r, s, nu).terms()
    mons = [(c, r**(i - 60) * s**(j - 60) * nu**(l - 60)) for (i, j, l), c in terms]
    # E = (1+t)[ Σ ε_m c_m m (ψ²qn − κ̃) + κ̃ − sW ψ²(W1) ]
    base = (1 + t) * (kt - sW * adams(W1))
    cols = [(1 + t) * c * mo * (adams(qn) - kt) for c, mo in mons]
    return g, base, cols, mons, Radj, kt


def eqs_at(points, fn):
    """fn(sub) → (base_val, [col_vals], target)。実部・虚部の線形方程式を返す。"""
    rows, rhs = [], []
    for pt in points:
        b, cs, tg = fn(pt)
        for f in (sp.re, sp.im):
            row = [float(f(c)) for c in cs]
            v = float(f(tg - b))
            if any(abs(x) > 1e-12 for x in row) or abs(v) > 1e-12:
                rows.append(row); rhs.append(v)
    return rows, rhs


def solve(m, sW=1, npts=14, seed=1):
    g, base, cols, mons, Radj, kt = build(m, sW)
    random.seed(seed)
    pts = [(sp.Rational(random.randint(2, 9), random.randint(2, 9)), sp.Rational(random.randint(2, 9), random.randint(2, 9))) for _ in range(npts)]
    R_ = sp.cancel(2 * (1 + t) * g[(2,)])
    g1x2 = adams(g[(1,)])
    Eq = sp.cancel((-(1 + r**2) * (g1x2 + 1)))
    A_rows, b_rows = [], []
    # (0) R̃_adj(ν² = −t/q) = 1：ν = i s/r
    def f0(pt):
        rv, sv = pt
        sub = {r: rv, s: sv, nu: I * sv / rv}
        cs = [sp.nsimplify(c * mo.subs(sub)) for c, mo in mons]
        return sp.Integer(0), cs, sp.Integer(1)
    a, b = eqs_at(pts, f0); A_rows += a; b_rows += b; n0 = len(a)
    # (0′) t = 1 の留数：(1 − t)E を s = 1 で評価して 0
    def f0p(pt):
        rv, nv = pt
        sub = {r: rv, s: 1, nu: nv}
        bb = sp.nsimplify(sp.cancel((1 - t) * base).subs(sub))
        cs = [sp.nsimplify(sp.cancel((1 - t) * c).subs(sub)) for c in cols]
        return bb, cs, sp.Integer(0)
    a, b = eqs_at(pts, f0p); A_rows += a; b_rows += b; n0p = len(a)
    # (i) t = −1：(R − E)|_{s=i} = 0。E の (1 − t)(1 + qA/t) の極は (0)(0′) で消えるので、評価は (1−t)(1+qA/t) を掛けた形で行う
    def fi(pt):
        rv, nv = pt
        sub = {r: rv, s: I, nu: nv}
        fac = (1 - t) * (1 + q * nu**2 / t)
        bb = sp.nsimplify(sp.cancel(fac * base).subs(sub))
        cs = [sp.nsimplify(sp.cancel(fac * c).subs(sub)) for c in cols]
        tg = sp.nsimplify(sp.cancel(fac * R_).subs(sub))
        return bb, cs, tg
    a, b = eqs_at(pts, fi); A_rows += a; b_rows += b; ni = len(a)
    # (ii) t = q
    def fii(pt):
        rv, nv = pt
        sub = {r: rv, s: rv, nu: nv}
        fac = (1 - t) * (1 + q * nu**2 / t)
        bb = sp.nsimplify(sp.cancel(fac * base).subs(sub))
        cs = [sp.nsimplify(sp.cancel(fac * c).subs(sub)) for c in cols]
        tg = sp.nsimplify(sp.cancel(fac * Eq).subs(sub))
        return bb, cs, tg
    a, b = eqs_at(pts, fii); A_rows += a; b_rows += b; nii = len(a)
    A = np.array(A_rows); bvec = np.array(b_rows)
    n = len(mons)
    # 列は c_m を含む（cols = c_m · m · …）ので、ε_m c_m = e_m と書き直す：e_m = c_m − 2 y_m、y_m ∈ [min(0,c), max(0,c)]（e = 0 を許す）
    cvec = np.array([float(c) for c, _ in mons])
    Acol = A / cvec                       # e_m の係数
    A2, b2 = -2 * Acol, bvec - Acol @ cvec
    lo, hi = np.minimum(0, cvec), np.maximum(0, cvec)
    out = {}
    idx = {"(0)": list(range(0, n0)), "(0′)": list(range(n0, n0 + n0p)), "(i)": list(range(n0 + n0p, n0 + n0p + ni)), "(ii)": list(range(n0 + n0p + ni, len(b_rows)))}
    combos = (("(0)", ["(0)"]), ("(0′)", ["(0′)"]), ("(0)+(0′)", ["(0)", "(0′)"]), ("(0′)+(i)", ["(0′)", "(i)"]), ("(0′)+(ii)", ["(0′)", "(ii)"]),
              ("(0′)+(i)+(ii)", ["(0′)", "(i)", "(ii)"]), ("(0)+(0′)+(i)", ["(0)", "(0′)", "(i)"]), ("all", ["(0)", "(0′)", "(i)", "(ii)"]))
    for lab, keys in combos:
        rows = sum((idx[k] for k in keys), [])
        cons = LinearConstraint(A2[rows], b2[rows] - 1e-7, b2[rows] + 1e-7)
        res = milp(np.zeros(n), constraints=cons, bounds=Bounds(lo, hi), integrality=np.ones(n))
        out[lab] = res.status == 0
        if res.status == 0 and lab == "all":
            y = np.round(res.x); e = cvec - 2 * y
            Et = sum(int(ei) * mo for ei, (c, mo) in zip(e, mons))
            print(f"   解：R̃_adj = {sp.factor(Et)}", flush=True)
    print(f"T(2,{2*m+1})、s_W = {sW:+d}：R_adj の単項式 {n} 個（Σ|c| = {int(abs(cvec).sum())}）。|e| ≤ |c|、e ≡ c (mod 2) で満たせる条件：{out}", flush=True)
    return out


if __name__ == "__main__":
    for m in (1, 2):
        for sW in (1, -1):
            solve(m, sW)
