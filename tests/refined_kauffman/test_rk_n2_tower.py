"""随伴の塔の母関数（rk_f_n2_tower、theory §12.51、続き 121）の機能テスト：閉じた形の κ_R(x)・κ_G(x)・背景の差（x = ν^{−2}）の級数が、
kmax = 20 の記号ブロックの σ の差 N_(2)∘σ − N_(1,1) の a₂²・gadj6・定数の係数（−(s³/r)(q−t) を外したもの）と ν^{−20} まで一致する。"""
import sys
import sympy as sp
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import rk_f_n2_grows as G
import rk_f_n2_rows_sym as RS
from rk_f_n2_rows_sym import r, s, q, t, a2, gadj, SYMS

x = sp.Symbol("x")
V, Vp = 1 + t + q * t, 1 + q + q * t
KAPPA_G = q**2 * x**3 * (q**2 * (Vp - t) - t**2 * x * (V - q)) / ((t * x - q) * (t**3 * x - q) * (t**2 * x - q**2) * (t * x - q**3))
BG = -q**2 * x**4 / ((x - q) * (t * x - q) * (t**2 * x - q) * (t * x - q**2) * (t * x - 1))
KAPPA_R_NUM = (r**12*s**2*x + 2*r**12*x - r**12 - r**10*s**6*x**2 - r**10*s**4*x - r**10*s**2*x**2 + r**10*x - r**8*s**6*x**3 + r**8*s**6*x**2 - r**8*s**4*x - r**8*s**2*x
               + r**6*s**10*x**3 - r**6*s**8*x**3 + r**6*s**8*x**2 - r**6*s**6*x**3 + r**6*s**6*x**2 - r**6*s**4*x**3 + r**6*s**4*x**2 - r**6*s**2*x**2 + r**4*s**10*x**4 + r**4*s**8*x**4
               - r**4*s**6*x**3 + r**4*s**6*x**2 - r**2*s**12*x**4 + r**2*s**10*x**3 + r**2*s**8*x**4 + r**2*s**6*x**3 + s**12*x**5 - 2*s**12*x**4 - s**10*x**4)
KAPPA_R = -r**4 * x**3 * KAPPA_R_NUM / ((x - q) * (t * x - q) * (t**2 * x - q) * (t**3 * x - q) * (t * x - q**2) * (t**2 * x - q**2) * (t * x - q**3) * (t * x - 1))


def _coeff(P, e):
    c = sp.sympify(P.coeff_monomial(e))
    return c.as_expr() if hasattr(c, "as_expr") else c


def test_tower_closed_forms_match_series():
    D = G.blocks("data/refined_kauffman/daha/n2_rows_sym_m1_k20.pkl")
    pref = -(s**3 / r) * (q - t)
    # sympy の series は遅いので、ν の有理式として rat_series（多項式の割り算）で展開する：x^n ↔ ν^{−2n}
    RS.SYMS = [sp.Symbol("_dummy")]; RS.ZERO = RS.P(0)
    nu = RS.nu
    ser = lambda K, n: {k // 2: v.as_expr() for k, v in RS.rat_series(K.subs(x, nu**-2), 2 * n).c.items()}
    sR, sG, sB = ser(KAPPA_R, 9), ser(KAPPA_G, 8), ser(BG, 11)
    class _S:
        def __init__(self, d): self.d = d
        def coeff(self, _x, i): return self.d.get(i, 0)
    serR, serG, serB = _S(sR), _S(sG), _S(sB)
    for k in range(10, 21, 2):
        Dk = G.sigma_difference(k, D)
        # a₂² は ν^{−k} ↔ x^{(k−10)/2+3}（κ_R ν⁴ の規格化：a₂²ν^{−4} × κ_R）、gadj6 ↔ x^{(k−12)/2+3}、定数 ↔ x^{(k−8)/2+4}
        cR = sp.cancel(_coeff(Dk, a2**2) / pref)
        assert sp.cancel(cR - serR.coeff(x, (k - 10) // 2 + 3)) == 0, ("κ_R", k)
        if k >= 12:
            cG = sp.cancel(_coeff(Dk, gadj[6]) / pref)
            assert sp.cancel(cG - serG.coeff(x, (k - 12) // 2 + 3)) == 0, ("κ_G", k)
        cB = sp.cancel(_coeff(Dk, sp.Integer(1)) / pref)
        assert sp.cancel(cB - serB.coeff(x, (k - 8) // 2 + 4)) == 0, ("bg", k)
