r"""W_λ と冪和の係数 co の段階で SO と Sp を比べる（§11.15 の定数 1 の由来、progress.md 続き 54）。
  SO：rk_c_2box_largeN_m.so_parts(1)（変数 U_D = νs）。Sp：rk_d_sp_two_box の Rx・coef_functions2・assemble（変数 U = x = t^N、ν = xs）。
  検査：W^{Sp}_λ(𝐚) = ±W^{SO}_λ(−𝐚)、co の各係数の Sp と SO(−𝐚) の差、Z_2 の定数項 −½S2[()] の差。"""
import sys, os, pickle, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp
import rk_d_sp_two_box as T
import rk_c_2box_largeN_m as Lm
from rcsknots.flint_rf3 import conv

r_, s_, u_, nu = T.r_, T.s_, T.u_, T.nu
to_nu_sp = lambda e: sp.cancel(e.to_sympy().subs(u_, nu / s_))       # U = x、ν = xs
to_nu_so = lambda e: sp.cancel(e.to_sympy().subs(u_, nu * s_))       # U_D = νs
neg = lambda e: e.subs(nu, -nu)

if __name__ == "__main__":
    t0 = time.time()
    Wso, coso = Lm.so_parts(1)
    print(f"SO 側 {time.time()-t0:.0f}s", flush=True)
    avail = [N for N in range(2, 20) if os.path.exists(os.path.join(T.OUT, f"cn_twobox_N{N}_m1.pkl"))]
    usub = s_**4
    Rx = {}
    for lam in T.LAMS:                                               # R_x_sp_m1_*.pkl：{j: str}（変数 r, s）
        d = pickle.load(open(os.path.join(T.OUT, f"R_x_sp_m1_{''.join(map(str, lam))}.pkl"), "rb"))
        Rx[lam] = {j: T.from_sympy(sp.sympify(v, locals=T.loc)) for j, v in d.items()}
    cf = T.coef_functions2(avail, usub, 1)
    fsp, Wsp, cosp = T.assemble(1, Rx, cf, {lam: conv(1) for lam in T.LAMS})
    print(f"Sp 側 {time.time()-t0:.0f}s", flush=True)
    for lam in T.LAMS:
        a, b = to_nu_sp(Wsp[lam]), to_nu_so(Wso[lam])
        for eps in (1, -1):
            d = sp.cancel(a - eps * neg(b))
            print(f"W_{lam}: Sp(𝐚) − ({eps:+d})·SO(−𝐚) =", sp.factor(d) if len(sp.Add.make_args(sp.expand(sp.fraction(d)[0]))) <= 12 else f"(分子 {len(sp.Add.make_args(sp.expand(sp.fraction(d)[0])))} 項)", flush=True)
    for name in ("S1sq", "S2"):
        for k in cosp[name]:
            a = to_nu_sp(cosp[name][k]); b = to_nu_so(coso[name][k]) if k in coso[name] else sp.Integer(0)
            print(f"co[{name}][{k}]: Sp =", sp.factor(a), "; SO(−𝐚) =", sp.factor(neg(b)), "; 差 =", sp.factor(sp.cancel(a - neg(b))), flush=True)
