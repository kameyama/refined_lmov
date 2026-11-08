"""brane/image を区別した記号ブロック（rk_f_n2_rows_lab、n2_rows_lab_m1_k14.pkl、theory §12.46、続き 116）の機能テスト：
ν^{−12} の X の 2 次の部分は −(s³/r)·X^b·X^i（brane の X ⊗ image の X）で、X^b = a₂²/((q−1)(t−1)) − w_b b₄/((q−1)(qt−1)) − w_c c₄/((t−1)(qt−1))。"""
import sys, pickle
import sympy as sp
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import rk_f_n2_rows_lab as L
import rk_f_n2_rows_sym as RS
from rk_f_n2_rows_sym import r, s, q, t


def _coeff(P, m):
    c = sp.sympify(P.coeff_monomial(sp.Poly.from_dict({m: 1}, *L.SYMS_LAB).as_expr()))
    return c.as_expr() if hasattr(c, "as_expr") else c


def test_square_is_brane_times_image():
    d = pickle.load(open("data/refined_kauffman/daha/n2_rows_lab_m1_k14.pkl", "rb"))
    loc = {str(x): x for x in L.SYMS_LAB}; loc.update({"r": r, "s": s})
    B, I = L.BR, L.IM
    for lab, wb, wc in (("2", 1, t), ("11", q, 1)):
        P = sp.Poly(sp.expand(sp.sympify(d["N"][(lab, 12)], locals=loc)), *L.SYMS_LAB, domain=RS.DOM)
        Xb = B["a2"]**2 / ((q - 1) * (t - 1)) - wb * B["b4"] / ((q - 1) * (q * t - 1)) - wc * B["c4"] / ((t - 1) * (q * t - 1))
        Xi = Xb.subs({B[k]: I[k] for k in B}, simultaneous=True)
        Q = sp.Poly(sp.expand(-(s**3 / r) * Xb * Xi), *L.SYMS_LAB, domain=RS.DOM)
        for m, c in Q.terms():
            assert sp.cancel(_coeff(P, m) - (sp.sympify(c).as_expr() if hasattr(c, "as_expr") else c)) == 0, (lab, m)
