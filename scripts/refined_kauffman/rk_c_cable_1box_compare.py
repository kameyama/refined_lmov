r"""表現 [1] の大きな N の g の σ 対称性の比較（progress.md 続き 140）：T(2,3)（thin）、T(3,4)（thick なトーラス結び目）、ケーブル K1, K3。
t = q での c = 1 の q ↔ 1/q 対称性と refined の c = 1 の σ 対称性、T(3,4) の c = 1・c = 2 の判定。"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rcsknots.refined_kauffman.flint_rf3 import R, S, U, from_sympy, conv
from rcsknots.refined_kauffman.kauffman_hopf2 import K_rho
from rcsknots.refined_kauffman.kauffman_hopf import qs, ts, Us
import rk_c_2box_largeN as L
import rk_b_2box_structure as B
import rk_f_frozen_basis as FB
r_, s_, u_ = sp.symbols("r s u", positive=True)
def g_one(tag, knot, delta):
    c = pickle.load(open(f"data/refined_kauffman/daha/R_U_{tag}_1.pkl", "rb"))
    Rf = sum((from_sympy(sp.sympify(v, locals={"r": r_, "s": s_})) * U**int(j) for j, v in c.items()), conv(0))
    kr = K_rho((1,)); kr = kr.as_expr() if hasattr(kr, "as_expr") else kr
    return -Rf * from_sympy(sp.sympify(kr).subs({qs: r_**2, ts: s_**2, Us: u_})) + (S / R)**delta * L.Rdata_nu("1", knot) * L.Psym((1,))
for tag, knot, d in (("m1k3", "T3_4", 3), ("m1", "T2_3", 1), ("K1", "K1", 8), ("K3", "K3", 9)):
    g = g_one(tag, knot, d)
    gt = sp.cancel(-g.to_sympy().subs(u_, FB.nu * FB.s).subs({r_: FB.r, s_: FB.s}, simultaneous=True))
    c1 = FB.parity_split(gt, True)
    x = sp.cancel(c1.subs(FB.s, FB.r)); tq = sp.cancel(x - x.subs(FB.r, 1 / FB.r)) == 0
    print(f"{knot}：整数 {g.is_integer_laurent()}、t = q の c=1 が q↔1/q 対称 {tq}、refined の c=1 が σ 対称 {sp.cancel(c1 - c1.subs({FB.r: 1/FB.s, FB.s: 1/FB.r}, simultaneous=True)) == 0}", flush=True)
    if knot == "T3_4":
        FB.test_c1("T(3,4) の [1] の c=1", c1); FB.test_c2("T(3,4) の [1] の c=2", FB.parity_split(gt, False))
