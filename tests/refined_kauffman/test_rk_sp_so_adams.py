"""SO と Sp の 2 箱の差（theory_orientifold_index.md §11.15）：
f^{Sp}_(2)(𝐚) − f^{SO}_(2)(−𝐚) = ψ²[δ]、f^{Sp}_(1,1)(𝐚) − f^{SO}_(1,1)(−𝐚) = −ψ²[δ]、δ = even_𝐚(g^{SO}_(1)) + 1（g でも同じ、被覆は共通）。"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rk_d_sp_so_adams import load, psi2, neg, r_, s_, nu, loc, OUT


def test_sp_so_two_box_difference_is_adams_of_even_part():
    g_sp, g_so = load()
    sp_f = pickle.load(open(f"{OUT}/largeN_2box_sp_m1.pkl", "rb"))["f"]
    so_f = pickle.load(open(f"{OUT}/largeN_2box_m1.pkl", "rb"))["f"]
    u_ = sp.Symbol("u", positive=True)
    f_sp = {k: sp.cancel(sp.sympify(v, locals=loc).subs(u_, nu / s_)) for k, v in sp_f.items()}
    f_so = {k: sp.cancel(sp.sympify(v, locals=loc).subs(u_, nu * s_)) for k, v in so_f.items()}
    g1 = g_so["(1,)"]
    delta = sp.cancel((g1 + neg(g1)) / 2 + 1)
    assert sp.cancel(delta - (nu**2 - s_**2) * (nu**2 * r_**2 - 1) * (nu**2 * r_**2 - s_**2) / (nu**6 * r_**4)) == 0
    D = psi2(delta)
    assert sp.cancel(f_sp["(2,)"] - neg(f_so["(2,)"]) - D) == 0
    assert sp.cancel(f_sp["(1, 1)"] - neg(f_so["(1, 1)"]) + D) == 0
    assert sp.cancel(g_sp["(2,)"] - neg(g_so["(2,)"]) - D) == 0
    # 冪和の成分：Z₁₁ = (f_(2)+f_(1,1))/2 + f₁²/2 は Sp(𝐚) = SO(−𝐚)、Z₂ = (f_(2)−f_(1,1))/2 + ψ²[f₁]/2 は Sp(𝐚) = SO(−𝐚) + 1
    Z2 = lambda f: sp.cancel((f["(2,)"] - f["(1, 1)"] + psi2(f["(1,)"])) / 2)
    Z11 = lambda f: sp.cancel((f["(2,)"] + f["(1, 1)"] + f["(1,)"]**2) / 2)
    assert sp.cancel(Z2(f_sp) - neg(Z2(f_so))) == 1
    assert sp.cancel(Z11(f_sp) - neg(Z11(f_so))) == 0
