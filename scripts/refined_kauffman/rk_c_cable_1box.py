r"""代数的なケーブルの表現 [1] の大きな N のクロスキャップ部分（progress.md 続き 140）。
g_(1) = f_(1) − ½f^{cov}_(1) = −W^{SO}_(1) + W^{U}_(1)（rk_c_2box_largeN.so_f・cover_f の [1] の部分）：
  W^{SO}_(1) = R_(1)(U) K_(1)(t^ρ)（R_U_{cable}_1.pkl、T_λ^{−F} の framing、F = p₂q₂）、
  W^{U}_(1) = (t/q)^{δ/2} R^{KN}_(1) P_(1)(t^ρ)（data/framed_knot/cables/{cable}_1.json、δ は rk_torus_tag.CABLE）。
判定（トーラス結び目の g_true = −g_(1) と同じ符号の規約、rk_b_uniform / rk_f_frozen_basis）：整数 Laurent、普遍因子、
c = 1（ν 偶）の G 展開と λ^β ごとの符号、c = 2（ν 奇）が −(𝐚²−𝐪²)(𝐭²−1) × 正 か。
使い方：python rk_c_cable_1box.py K1"""
import sys, io, pickle, contextlib
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rcsknots.refined_kauffman.flint_rf3 import R, S, U, from_sympy, conv
from rcsknots.refined_kauffman.kauffman_hopf2 import K_rho
from rcsknots.refined_kauffman.kauffman_hopf import qs, ts, Us
import rk_c_2box_largeN as L
import rk_torus_tag as TG

r_, s_, u_ = sp.symbols("r s u", positive=True)


def g_one(cable):
    c = pickle.load(open(f"data/refined_kauffman/daha/R_U_{cable}_1.pkl", "rb"))
    Rf = sum((from_sympy(sp.sympify(v, locals={"r": r_, "s": s_})) * U**int(j) for j, v in c.items()), conv(0))
    kr = K_rho((1,))
    kr = kr.as_expr() if hasattr(kr, "as_expr") else kr
    Wso = Rf * from_sympy(sp.sympify(kr).subs({qs: r_**2, ts: s_**2, Us: u_}))
    Wu = (S / R)**TG.delta(cable) * L.Rdata_nu("1", cable) * L.Psym((1,))
    return -Wso + Wu


if __name__ == "__main__":
    import rk_b_2box_structure as B
    import rk_f_frozen_basis as FB
    cable = sys.argv[1]
    g = g_one(cable)
    print(f"{cable}：g_(1) は整数 Laurent {g.is_integer_laurent()}", flush=True)
    gt = sp.cancel(-g.to_sympy().subs(u_, FB.nu * FB.s).subs({r_: FB.r, s_: FB.s}, simultaneous=True))     # g_true = −g_(1)、変数 (r, s, ν)
    pickle.dump({"g_true": str(gt)}, open(f"data/refined_kauffman/daha/largeN_1box_{cable}.pkl", "wb"))
    e = sp.cancel(gt.subs({FB.r: B.r, FB.s: B.s, FB.nu: B.nu}, simultaneous=True))
    print(f"  普遍因子 (q − t)(ν² − t) で refined だけの部分が割り切れる {B.divisible_universal(e)}", flush=True)
    sb = B.sign_blocks(e)
    print(f"  商の ā ブロックの符号：{sum(1 for n_, sg in sb.values() if sg != '混在')}/{len(sb)} {sb}", flush=True)
    FB.test_c1(f"{cable} の [1] の c=1（ν 偶）", FB.parity_split(gt, True), verbose=False)
    FB.test_c2(f"{cable} の [1] の c=2（ν 奇）", FB.parity_split(gt, False))
