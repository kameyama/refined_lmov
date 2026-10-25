r"""タスク C (b)：大きな N の 2 箱のクロスキャップ部分 g（trefoil、全 d）を、自由エネルギーの右辺の Schur 基底から別の基底（Macdonald の P, Q、Hall–Littlewood、modified Macdonald）に替えたときの係数の分母（progress.md 続き 28）。"""
import sys, pickle; sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rcsknots.refined_kauffman.flint_rf3 import RF, from_sympy, conv as fconv
from rcsknots.macdonald import macdonald_P, macdonald_norm
from rcsknots.modified import modified_H
r, s, u = sp.symbols("r s u", positive=True)
q, t = r**2, s**2
conv = lambda f: sp.Rational(f.numerator, f.denominator)
d = pickle.load(open("data/refined_kauffman/daha/largeN_2box_m1.pkl", "rb"))
g = {eval(k): from_sympy(sp.sympify(v, locals={"r": r, "s": s, "u": u})) for k, v in d["g"].items()}
half = fconv(sp.Rational(1, 2)) if False else None
from fractions import Fraction
h = fconv(Fraction(1, 2))
c11, c2 = (g[(2,)] + g[(1, 1)]) * h, (g[(2,)] - g[(1, 1)]) * h
def in_basis(Pb):
    lams = [(2,), (1, 1)]
    Mx = sp.Matrix([[Pb[l].get(nu_, 0) for l in lams] for nu_ in ((1, 1), (2,))])
    Mi = Mx.inv().applyfunc(sp.cancel)                       # 小さな 2×2
    out = {}
    for i, l in enumerate(lams):
        a, b = from_sympy(sp.sympify(Mi[i, 0])), from_sympy(sp.sympify(Mi[i, 1]))
        out[l] = a * c11 + b * c2
    return out
def den(x):
    return sp.factor(RF(x.d).to_sympy() if False else sp.Poly(str(x.d).replace("^", "**"), r, s, u).as_expr())
cands = {}
for lab, (qq, tt) in {"P(q,t)": (q, t), "P(t,q)": (t, q), "P(q⁻¹,t⁻¹)": (1 / q, 1 / t), "P(t⁻¹,q⁻¹)": (1 / t, 1 / q)}.items():
    P = macdonald_P(2, qq, tt, conv)
    cands[lab] = P
    cands["Q" + lab[1:]] = {l: {k: v / macdonald_norm(l, qq, tt) for k, v in P[l].items()} for l in P}
cands["HL(t)"] = macdonald_P(2, sp.Integer(0), t, conv)
cands["HL(q)"] = macdonald_P(2, sp.Integer(0), q, conv)
cands["H̃(q,t)"] = modified_H(2, q, t, conv)
cands["H̃(t,q)"] = modified_H(2, t, q, conv)
cands["Schur"] = {(2,): {(1, 1): sp.Rational(1, 2), (2,): sp.Rational(1, 2)}, (1, 1): {(1, 1): sp.Rational(1, 2), (2,): -sp.Rational(1, 2)}}
for lab, Pb in cands.items():
    x = in_basis(Pb)
    print(f"{lab}: (2) 分母 {sp.factor(x[(2,)].to_sympy().as_numer_denom()[1])}、(1,1) 分母 {sp.factor(x[(1, 1)].to_sympy().as_numer_denom()[1])}", flush=True)
