r"""タスク C：大きな N（U = t^N は形式変数）の refined SO(2N) trefoil の 2 箱のクロスキャップ部分 g = f − ½f^{cov}（一様な規格化）。
  SO：W_λ(U) = R_λ(U) K_λ(t^ρ)（R_λ は D_N の DAHA-Jones を U で補間、scripts/rk_c_interp.py）、核の係数は kauffman_hopf2.koornwinder_small。
  被覆 U(N_A)：A = ν² = U²/t、ν = U/s。普通の表現は W_μ = (t/q)^{|μ|/2} R_μ(KN の data、a^{1/2} = ν) P_μ(t^ρ)、
  随伴は W_adj = R_adj(ν) P_adj(t^ρ)（GL の DAHA-Jones を A^{1/2} で補間）、P_adj(t^ρ) = [N_A]² − κ。
変数は (r, s, u) = (q^{1/2}, t^{1/2}, U)（rcsknots.refined_kauffman.flint_rf3）。"""
import sys, json, pickle, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from fractions import Fraction
from rcsknots.refined_kauffman.flint_rf3 import RF, R, S, U, from_sympy, conv, ONE
from rcsknots.macdonald import principal_specialization

r_, s_, u_ = sp.symbols("r s u", positive=True)
q, t = R**2, S**2
nu = U / S
half = conv(Fraction(1, 2))


def so_parts():
    from rcsknots.refined_kauffman.kauffman_hopf2 import koornwinder_small, K_rho
    from rcsknots.refined_kauffman.kauffman_hopf import qs, ts, Us
    import rk_c_interp as I
    sub = {qs: r_**2, ts: s_**2, Us: u_}
    cv = lambda v: from_sympy(sp.sympify(v.as_expr() if hasattr(v, "as_expr") else v).subs(sub))
    ks = koornwinder_small()
    co = {name: {lam: cv(v) for lam, v in ks[name].items()} for name in ("S1sq", "S2")}
    assert cv(ks["kappa"]).is_zero()
    c1 = I.interpolate(1, (1,), [2, 3, 4, 5], -5, -2)
    Rl = {(1,): c1}
    saved = pickle.load(open("data/refined_kauffman/daha/R_U_m1_2box.pkl", "rb"))
    for lam in ((1, 1), (2,)):
        Rl[lam] = {j: from_sympy(sp.sympify(v, locals={"r": r_, "s": s_})) for j, v in saved[lam].items()}
    W = {}
    for lam, c in Rl.items():
        Rf = sum((v if isinstance(v, RF) else from_sympy(v.to_sympy())) * U**j for j, v in ((j, (from_sympy(v.to_sympy()) if not isinstance(v, RF) else v)) for j, v in c.items()))
        W[lam] = Rf * cv(K_rho(lam))
    return W, co


def so_f(W, co):
    Z1 = -W[(1,)]
    Z11 = half * (W[(2,)] + co["S1sq"][(1, 1)] * W[(1, 1)] + co["S1sq"][(1,)] * W[(1,)] + co["S1sq"][()])
    Z2 = -half * (W[(2,)] + co["S2"][(1, 1)] * W[(1, 1)] + co["S2"][(1,)] * W[(1,)] + co["S2"][()])
    f1 = Z1
    H11 = Z11 - Z1 * Z1 * half
    H2 = Z2 - f1.adams(2) * half
    return {(1,): f1, (2,): H11 + H2, (1, 1): H11 - H2}


def Rdata_nu(c, knot="T2_3"):
    import os
    p = f"data/refined_kauffman/torus/{knot}_{c}.json"
    d = json.load(open(p if os.path.exists(p) else f"data/framed_knot/torus/{knot}_{c}.json"))
    out = conv(0)
    for k, v in d["poly"].items():
        i, j, l = map(int, k.split(","))
        out = out + conv(v) * nu**(-i) * R**(-j) * S**(-l)      # 鏡像（rk_k5_trefoil_2box.Rdata と同じ）、a^{1/2} = ν
    return out


def Psym(lam):
    k = sum(lam)
    ps = principal_specialization(lam, nu * nu, q, t)            # a = A = ν²
    return S**k * nu**(-k) * ps


def cover_f():
    w = S / R
    Wc = {lam: w**sum(lam) * Rdata_nu("".join(map(str, lam))) * Psym(lam) for lam in ((1,), (1, 1), (2,))}
    radj = pickle.load(open("data/refined_kauffman/daha/Radj_A_m1.pkl", "rb"))
    Radj = sum((from_sympy(sp.sympify(v, locals={"r": r_, "s": s_})) * nu**j for j, v in radj.items()), conv(0))
    A = nu * nu
    qn = (nu - 1 / nu) / (S - 1 / S)
    kappa = (1 - A) * (1 - q) / ((1 - t) * (1 - q * A / t))
    Wadj = Radj * (qn * qn - kappa)
    c = (t - q) / (1 - q * t)
    A1 = -2 * Wc[(1,)]
    A2 = 2 * Wc[(1, 1)] + Wadj + kappa
    A11 = 2 * (Wc[(2,)] + c * Wc[(1, 1)]) + Wadj + kappa
    Z1, Z11, Z2 = A1, (A2 + A11) * half, (A2 - A11) * half
    f1 = Z1
    H11 = Z11 - Z1 * Z1 * half
    H2 = Z2 - f1.adams(2) * half
    return {(1,): f1, (2,): H11 + H2, (1, 1): H11 - H2}


if __name__ == "__main__":
    t0 = time.time()
    W, co = so_parts()
    f = so_f(W, co)
    fc = cover_f()
    g = {k: f[k] - half * fc[k] for k in f}
    print(f"組み立て {time.time()-t0:.0f}s", flush=True)
    pickle.dump({"f": {str(k): str(v.to_sympy()) for k, v in f.items()}, "fcov": {str(k): str(v.to_sympy()) for k, v in fc.items()},
                 "g": {str(k): str(v.to_sympy()) for k, v in g.items()}}, open("data/refined_kauffman/daha/largeN_2box_m1.pkl", "wb"))
    for k, v in g.items():
        print(f"g_{k}: Laurent {v.is_laurent()}、分母 {sp.factor(v.d and RF(v.d).to_sympy() if False else sp.fraction(sp.cancel(v.to_sympy()))[1])}", flush=True)
