"""refined Sp(2N) の trefoil の 2 箱（theory_orientifold_index.md §11.14、保存データ largeN_2box_sp_m1.pkl と largeN_2box_m1.pkl）。
(1) 整数性：g_(2)+g_(1,1) と (1+t)(g_(2)−g_(1,1)) が整数 Laurent、(2) E^{Sp}(t=−1) = E^{SO}(t=−1)、(3) hatted の p1² 成分で Sp(𝐚) = +SO(−𝐚)。"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp

r_, s_, u_ = sp.symbols("r s u", positive=True)
nu = sp.Symbol("nu", positive=True)
loc = {"r": r_, "s": s_, "u": u_}
OUT = "data/refined_kauffman/daha"


def _load():
    gsp = {k: sp.cancel(sp.sympify(v, locals=loc).subs(u_, nu / s_)) for k, v in pickle.load(open(f"{OUT}/largeN_2box_sp_m1.pkl", "rb"))["g"].items()}
    gso = {k: sp.cancel(sp.sympify(v, locals=loc).subs(u_, nu * s_)) for k, v in pickle.load(open(f"{OUT}/largeN_2box_m1.pkl", "rb"))["g"].items()}
    return gsp, gso


def _int_laurent(e):
    n_, d_ = sp.fraction(sp.cancel(e))
    return sp.Poly(d_, r_, s_, nu).is_monomial and all(c.is_integer for c in sp.Poly(sp.expand(n_ * (r_ * s_ * nu)**60), r_, s_, nu).coeffs())


def test_sp_two_box_integrality_and_E():
    gsp, gso = _load()
    t = s_**2
    assert _int_laurent(gsp["(2,)"] + gsp["(1, 1)"])
    assert _int_laurent((1 + t) * (gsp["(2,)"] - gsp["(1, 1)"]))
    Esp = sp.cancel((1 + t) * (gsp["(2,)"] - gsp["(1, 1)"])).subs(s_, sp.I)
    Eso = sp.cancel((1 + t) * (gso["(2,)"] - gso["(1, 1)"])).subs(s_, sp.I)
    assert sp.cancel(Esp - Eso) == 0


def test_sp_so_hatted_p1sq_relation():
    gsp, gso = _load()
    A_, Q_, T_ = sp.symbols("A Q T")
    bold = lambda e: sp.expand(sp.cancel(e.subs({r_: -T_ * Q_, nu: A_}).subs(s_, Q_)))
    zq = -T_ * Q_ + 1 / (T_ * Q_)
    hsp = sp.cancel(bold(gsp["(2,)"] + gsp["(1, 1)"]) / zq)
    hso = sp.cancel(bold(gso["(2,)"] + gso["(1, 1)"]) / zq)
    assert sp.simplify(hsp - hso.subs(A_, -A_)) == 0
