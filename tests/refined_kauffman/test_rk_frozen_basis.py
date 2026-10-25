"""凍結した状態の多重項の基底（scripts/refined_kauffman/rk_f_frozen_basis.py、theory §12）の機能テスト。
  (1) [1] T(2,3)：c=1 の ν ブロックは σ 対称で G 展開の符号が +, −, +；c=2 は (𝐚²−𝐪²)(𝐭²−1) × (−)。
  (2) trefoil の 2 箱 p₁² hatted の c=1（ν 奇）：σ 非対称なブロックで d_m = c_m − c_{−m} が (1−w²) × 同符号（全行同じ符号）、
      A⁻ の G 展開の符号は非対称部分と逆、β で交互。
  (3) Hopf：P_22|_{2β} = 単項式 × (C_11|_β/z_q)²（両端の β）。"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_f_frozen_basis as F


def test_one_box_trefoil_blocks():
    X = F.knot_one_box(1)
    res = F.test_c1("[1] T(2,3)", F.parity_split(X, True))
    assert all(v["sym"] for v in res.values())
    assert {b: v["Aminus"][0] for b, v in res.items()} == {-6: "+", -4: "−", -2: "+"}
    c2 = F.parity_split(X, False)
    Q = sp.cancel(c2 / ((F.nu**2 - F.t) * (F.r**2 - F.s**2)))
    assert sp.Poly(sp.fraction(Q)[1], F.r, F.s, F.nu).is_monomial
    assert F.sign_of(sp.Poly(sp.expand(sp.fraction(Q)[0]), F.r, F.s, F.nu).coeffs()) == "−"


def test_two_box_trefoil_asymmetric_rows():
    X = F.knot_two_box(1)
    res = F.test_c1("T(2,3) 2 箱", F.parity_split(X, False))
    asym = {b: v for b, v in res.items() if not v["sym"]}
    assert sorted(asym) == [-11, -9, -7, -5]
    assert all(v["uniform"] for v in asym.values())
    assert {b: v["asym_sign"] for b, v in asym.items()} == {-11: "−", -9: "+", -7: "−", -5: "+"}
    assert {b: v["Aminus"][0] for b, v in asym.items()} == {-11: "+", -9: "−", -7: "+", -5: "−"}
    assert res[-3]["sym"] and res[-3]["Aminus"][0] == "+"
    assert res[-1]["sym"] and res[-1]["Aminus"][1] == "+"


def test_hopf_sqrt_11_22():
    r, s, nu = sp.symbols("r s nu", positive=True)
    loc = {"r": r, "s": s, "nu": nu}
    g = {k: sp.sympify(v, locals=loc) for k, v in pickle.load(open("data/refined_kauffman/hopf_g_alld.pkl", "rb")).items()}
    f4 = {k: sp.sympify(v, locals=loc) for k, v in pickle.load(open("data/refined_kauffman/hopf_fcov_sector_A4.pkl", "rb")).items()}
    zq = r - 1 / r
    C = sp.cancel(g[(1, 1)] / zq)                           # z_q(ν~ − ν~⁻¹)/z_q
    P = sp.cancel((1 - s**2) * f4[(1, 1, 1, 1)] / zq**3)
    nut = sp.Symbol("nut", positive=True)
    Cn = sp.expand(sp.cancel(C.subs(nu, nut * s / r))); Pn = sp.expand(sp.cancel(P.subs(nu, nut * s / r)))
    assert Cn.coeff(nut, -1) == -1 and Cn.coeff(nut, 1) == 1
    assert sp.cancel(Pn.coeff(nut, -2) / Cn.coeff(nut, -1)**2) == r
    assert sp.cancel(Pn.coeff(nut, 2) / Cn.coeff(nut, 1)**2) == s**2 / r
