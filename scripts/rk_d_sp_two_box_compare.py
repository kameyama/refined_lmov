r"""Sp の 2 箱（rk_d_sp_two_box、mono = 1）を確定し、SO（largeN_2box_m1.pkl）と比べる（progress.md 続き 52）。
  検査：g_(2)+g_(1,1) 整数 Laurent、(1+t)(g_(2)−g_(1,1)) 整数 Laurent、U = t², t³ の固定 N は無いので省略。
  比較：g^{Sp}_R(𝐚) = −g^{SO}_R(−𝐚)（R = (2), (1,1)）、E^{Sp}(t=−1) と E^{SO}(t=−1)、E を含まない方向 g_(2)+g_(1,1)。"""
import sys, os, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp
from fractions import Fraction
import rk_d_sp_two_box as T
from rcsknots.flint_rf3 import conv

r_, s_, u_ = T.r_, T.s_, T.u_
nu = T.nu
OUT = T.OUT


def main(m=1):
    avail = [N for N in range(2, 20) if os.path.exists(os.path.join(OUT, f"cn_twobox_N{N}_m{m}.pkl"))]
    usub = s_**4
    Rx = {}
    for lam in T.LAMS:
        c_, nodes, ok = T.interp_laurent({N: T.R_of_N(N, lam, usub, m) for N in avail}, avail)
        assert all(ok.values()); Rx[lam] = c_
    cf = T.coef_functions2(avail, usub, m)
    fc = T.cover_in_x(m)
    half = conv(Fraction(1, 2))
    f, W, co = T.assemble(m, Rx, cf, {lam: conv(1) for lam in T.LAMS})
    g = {k: f[k] - half * fc[k] for k in f}
    t = s_**2
    print("g_(2)+g_(1,1) 整数 Laurent:", (g[(2,)] + g[(1, 1)]).is_integer_laurent(), "；(1+t)(g_(2)−g_(1,1)) 整数 Laurent:", ((1 + t) * (g[(2,)] - g[(1, 1)])).is_integer_laurent() if hasattr(g[(2,)], "is_integer_laurent") else None, flush=True)
    gsp_str = {str(k): str(v.to_sympy()) for k, v in g.items()}
    pickle.dump({"f": {str(k): str(v.to_sympy()) for k, v in f.items()}, "g": gsp_str, "note": "U = x = t^N、ν = x s（Sp(2N)）"}, open(os.path.join(OUT, f"largeN_2box_sp_m{m}.pkl"), "wb"))
    gso = pickle.load(open(os.path.join(OUT, f"largeN_2box_m{m}.pkl"), "rb"))["g"]
    A_, Q_, T_ = sp.symbols("A Q T")
    gb_sp = {k: T.bold_nu(v)[0] for k, v in g.items()}
    gb_so = {k: T.bold_so(gso[str(k)]) for k in g}
    for k in g:
        print(f"R = {k}：項数 Sp {len(sp.Add.make_args(gb_sp[k]))}, SO {len(sp.Add.make_args(gb_so[k]))}；g^Sp(𝐚) = −g^SO(−𝐚)：{sp.expand(gb_sp[k] + gb_so[k].subs(A_, -A_)) == 0}；𝐭 → −𝐭 不変 Sp {sp.expand(gb_sp[k] - gb_sp[k].subs(T_, -T_)) == 0}", flush=True)
    ssum_sp = sp.expand(gb_sp[(2,)] + gb_sp[(1, 1)]); ssum_so = sp.expand(gb_so[(2,)] + gb_so[(1, 1)])
    print("E を含まない方向 g_(2)+g_(1,1)：Sp(𝐚) = −SO(−𝐚)：", sp.expand(ssum_sp + ssum_so.subs(A_, -A_)) == 0, flush=True)
    # E(t = −1)：(1+t)(g_(2) − g_(1,1)) を t = −1（𝐭 = −1 ではなく t = s² = −1、すなわち Q² = −1）で
    E_sp = sp.cancel(((1 + t) * (g[(2,)] - g[(1, 1)])).to_sympy().subs(u_, nu / s_))
    E_so = sp.cancel(((1 + t) * (sp.sympify(gso["(2,)"], locals=T.loc) - sp.sympify(gso["(1, 1)"], locals=T.loc))).subs(u_, nu * s_))
    Esp1, Eso1 = sp.factor(E_sp.subs(s_, sp.I)), sp.factor(E_so.subs(s_, sp.I))
    print("E^Sp(t=−1) =", Esp1, flush=True)
    print("E^SO(t=−1) =", Eso1, flush=True)
    for lab, sub in (("ν→−ν", {nu: -nu}), ("ν→ν", {}), ("ν→1/ν", {nu: 1 / nu}), ("ν→−1/ν", {nu: -1 / nu})):
        ratio = sp.cancel(Esp1 / Eso1.subs(sub))
        n_, d_ = sp.fraction(ratio)
        if sp.Poly(n_, r_, nu).is_monomial and sp.Poly(d_, r_, nu).is_monomial:
            print(f"   E^Sp(t=−1) = E^SO(t=−1)[{lab}] × {sp.factor(ratio)}", flush=True)
    return g


if __name__ == "__main__":
    main()
