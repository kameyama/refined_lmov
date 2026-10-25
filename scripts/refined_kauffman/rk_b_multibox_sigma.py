r"""多箱の c = 1 部分の σ 非対称性の原因の切り分け（theory §11.16、progress.md 続き 55）。
  σ：r → 1/s、s → 1/r、ν 固定。c = 1：ν の偶数冪（分母に対する相対）。
  仮説：全 d の g に含まれる Adams 項 ½ψ²[g_1]（ψ²：r, s, ν → 2 乗）は、[1] の c = 2 部分（ν の奇数冪、σ 非対称）を ν の偶数冪に移すので、
  2 箱の「c = 1（ν の偶奇で定義）」に σ 非対称な寄与が混ざる。同様に積 g_1² も (c=2)² を偶数冪に移す。
  検査：p₁²・p₂ 方向、hatted、Adams 項・積を戻したもの、t = q の極限（q → 1/q、ν 固定）について c = 1 の σ 対称性。"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp

r_, s_, u_ = sp.symbols("r s u", positive=True)
nu = sp.Symbol("nu", positive=True)
loc = {"r": r_, "s": s_, "u": u_}
OUT = "data/refined_kauffman/daha"


def load(m):
    d = pickle.load(open(f"{OUT}/largeN_2box_m{m}.pkl", "rb"))
    conv = lambda v: sp.cancel(sp.sympify(v, locals=loc).subs(u_, nu * s_))
    return {k: conv(v) for k, v in d["g"].items()}, {k: conv(v) for k, v in d["f"].items()}


def sigma(e):
    return e.subs({r_: 1 / s_, s_: 1 / r_}, simultaneous=True)


def psi2(e):
    return e.subs({r_: r_**2, s_: s_**2, nu: nu**2}, simultaneous=True)


def c1(e):
    """ν の偶数冪の部分（分母の ν の冪に対する相対）。"""
    n_, d_ = sp.fraction(sp.cancel(e))
    P = sp.Poly(sp.expand(n_), nu)
    dn = sp.Poly(d_, nu).monoms()[0][0] if sp.Poly(d_, nu).is_monomial else 0
    keep = sum((c * nu**k[0] for k, c in zip(P.monoms(), P.coeffs()) if (k[0] - dn) % 2 == 0), sp.Integer(0))
    return sp.cancel(keep / d_)


def sym(e):
    return sp.cancel(e - sigma(e)) == 0


def nterms(e):
    return len(sp.Add.make_args(sp.expand(sp.fraction(sp.cancel(e))[0])))


def report(lab, e, tq=False):
    x = c1(e)
    if tq:   # t = q：s = r、σ は r → 1/r
        x = sp.cancel(x.subs(s_, r_))
        ok = sp.cancel(x - x.subs(r_, 1 / r_)) == 0
    else:
        ok = sym(x)
    anti = sp.cancel(x - (x.subs(r_, 1 / r_) if tq else sigma(x)))
    print(f"  {lab}: c=1 が σ 対称 {ok}" + ("" if ok else f"（非対称部分の分子 {nterms(anti)} 項、分母 {sp.factor(sp.fraction(anti)[1])}）"), flush=True)


if __name__ == "__main__":
    for m in (1, 2):
        g, f = load(m)
        g1, g2, g11 = g["(1,)"], g["(2,)"], g["(1, 1)"]
        fc1 = sp.cancel(f["(1,)"] - g1)                                  # ½f^{cov}_1
        odd1 = sp.cancel((g1 - g1.subs(nu, -nu)) / 2); ev1 = sp.cancel(g1 - odd1)
        zq, zqp, t = r_ - 1 / r_, r_ + 1 / r_, s_**2
        print(f"== T(2,{2*m+1})", flush=True)
        report("[1] g_1", g1)
        print(" p₁² 方向", flush=True)
        report("g_(2)+g_(1,1)", g2 + g11)
        report("(g_(2)+g_(1,1))/z_q（hatted）", (g2 + g11) / zq)
        report("(g_(2)+g_(1,1))/z_t", (g2 + g11) / (s_ - 1 / s_))
        report("g_(2)+g_(1,1) + g_1²", g2 + g11 + g1**2)
        report("g_(2)+g_(1,1) + odd(g_1)²", g2 + g11 + odd1**2)
        report("g_(2)+g_(1,1) − odd(g_1)²", g2 + g11 - odd1**2)
        report("g_(2)+g_(1,1) + 2 odd(g_1)·½f^cov_1", g2 + g11 + 2 * odd1 * fc1)
        print(" p₂ 方向", flush=True)
        report("g_(2)−g_(1,1)", g2 - g11)
        report("(1+t)(g_(2)−g_(1,1))", (1 + t) * (g2 - g11))
        report("(g_(2)−g_(1,1))/z_q⁺（hatted）", (g2 - g11) / zqp)
        for sg, lab in ((1, "+"), (-1, "−")):
            report(f"g_(2)−g_(1,1) {lab} ψ²[g_1]", g2 - g11 + sg * psi2(g1))
            report(f"g_(2)−g_(1,1) {lab} ψ²[odd g_1]", g2 - g11 + sg * psi2(odd1))
            report(f"g_(2)−g_(1,1) {lab} ψ²[even g_1]", g2 - g11 + sg * psi2(ev1))
            report(f"g_(2)−g_(1,1) {lab} ½ψ²[g_1]", g2 - g11 + sg * psi2(g1) / 2)
        print(" t = q（q → 1/q、ν 固定）", flush=True)
        report("[1] g_1", g1, tq=True)
        report("g_(2)+g_(1,1)", g2 + g11, tq=True)
        report("g_(2)−g_(1,1)", g2 - g11, tq=True)
        report("(1+t)(g_(2)−g_(1,1))", (1 + t) * (g2 - g11), tq=True)
