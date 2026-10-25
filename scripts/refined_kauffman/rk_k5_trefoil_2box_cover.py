r"""固定 N の trefoil の 2 箱：被覆 U(N_A)（N_A = 2N − 1、A = ν² = t^{N_A}）とクロスキャップ部分 g = f − ½f^{cov}。
  Π_{α,i}(1 − y_α x_i)(1 − y_α/x_i) = Σ_{λ,μ} (−1)^{|λ|+|μ|} s_λ(Y) s_μ(Y) s_{λᵀ}(x) s_{μᵀ}(1/x)。|Y| ≤ 2 では
    s_1(Y)：−(e_1(x) + e_1(1/x))、s_2(Y)：e_2(x) + e_2(1/x) + e_1(x)e_1(1/x)、s_11(Y)：h_2(x) + h_2(1/x) + e_1(x)e_1(1/x)。
  U(N_A) の Macdonald で h_2 = P_(2) + c P_(1,1)（c = (t − q)/(1 − qt)）、e_1(x)e_1(1/x) = P_adj + κ（κ = (1−t^{N_A})(1−q)/((1−t)(1−qt^{N_A−1}))、
  P_adj = P_{(2,1^{N_A−2})}/e_{N_A}）。荷電共役 ⟨F(1/x)⟩ = ⟨F(x)⟩ で W_μ̄ = W_μ（data では R_μ̄/R_μ が単項式であることを SU(3), SU(4) で確認）。
  W_μ = m_μ R_μ P_μ(t^ρ)（R は data の鏡像、A_data = s^{N_A}；[1] は K3' の ½h^{ref} と比 1）。"""
import sys, time, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_k5_trefoil_2box as T
from rcsknots.refined_kauffman.kauffman_knot import r, s

q, t = r**2, s**2


def schur_to_f(A1, A2, A11):
    """Z = 1 + A1 s_1 + A2 s_2 + A11 s_11 → {(1): f, (2): f, (1,1): f}（全 d の Möbius、x^d：r → r^d、s → s^d）。"""
    Z1, Z11, Z2 = A1, (A2 + A11) / 2, (A2 - A11) / 2
    f1 = sp.cancel(Z1)
    H11 = Z11 - Z1**2 / 2
    H2 = Z2 - f1.subs({r: r**2, s: s**2}, simultaneous=True) / 2
    return {(1,): f1, (2,): sp.cancel(H11 + H2), (1, 1): sp.cancel(H11 - H2)}


ADJ = {3: "21", 5: "2111"}


def cover_W(N, m=None, unknot=False):
    NA = 2 * N - 1
    m = m or {}
    lam = {"1": (1,), "11": (1, 1), "2": (2,), "adj": (2,) + (1,) * (NA - 2)}
    W = {}
    for k, l in lam.items():
        R = sp.Integer(1) if unknot else T.Rdata(ADJ[NA] if k == "adj" else k, NA)
        W[k] = m.get(k, 1) * R * T.P_A(l, NA)
    return W


def cover_f(N, m=None, unknot=False):
    NA = 2 * N - 1
    W = cover_W(N, m, unknot)
    c = (t - q) / (1 - q * t)
    kappa = (1 - t**NA) * (1 - q) / ((1 - t) * (1 - q * t**(NA - 1)))
    A1 = -2 * W["1"]
    A2 = 2 * W["11"] + W["adj"] + kappa
    A11 = 2 * (W["2"] + c * W["11"]) + W["adj"] + kappa
    return schur_to_f(A1, A2, A11)


def so_f(N, unknot=False, m1=None):
    if unknot:
        return T.free_energy_2(N, {lam: T.unknot_K(N, lam) for lam in ((1,), (1, 1), (2,))})
    return T.free_energy_2(N, T.trefoil_W(N, T.M2[N], m1=m1))


def M_ADJ(N):
    """随伴の単項式：極 q t^{N_A−1} = 1 の上で決まる族 r s^{11}(qt²)^k（N_A=3）・r s^{29}(qt⁴)^k（N_A=5）のうち t = q で 1 のもの。"""
    return T.W_**(2 * N - 1)


def crosscap(N, m1=None):
    f = so_f(N, m1=m1)
    fc = cover_f(N, {"adj": M_ADJ(N)})
    return f, fc, {R: sp.cancel(f[R] - fc[R] / 2) for R in f}


def p2_divisible(e, x):
    """e の p_2 成分（(2) − (1,1)）の分子が x で割り切れるか。"""
    n = sp.fraction(sp.cancel(e[(2,)] - e[(1, 1)]))[0]
    return sp.rem(sp.Poly(n, r, s), sp.Poly(x, r, s)).is_zero


if __name__ == "__main__":
    from rcsknots.bps import M_matrix
    from rcsknots.refined_kauffman.kauffman_trefoil import g_ref, r as rT, s as sT, U as UT
    Minv = M_matrix(2, 1 / r).inv()
    zq = r - 1 / r
    out = {}
    for N in (2, 3):
        t0 = time.time()
        NA = 2 * N - 1
        fu = cover_f(N, unknot=True)
        gu = {R: sp.cancel(v - fu[R] / 2) for R, v in so_f(N, unknot=True).items()}
        print(f"N={N} unknot：f^cov_(1) = −2P_[1]: {sp.cancel(fu[(1,)] + 2 * T.P_A((1,), NA)) == 0}、f^cov_(2) = {fu[(2,)]}、f^cov_(1,1) = {fu[(1, 1)]}、g = {gu}", flush=True)
        print(f"N={N} 被覆の m_adj = 1 で Laurent: {all(T.laurent(v) for v in cover_f(N).values())}、m_adj = (t/q)^(N_A/2) で Laurent: {all(T.laurent(v) for v in cover_f(N, {'adj': M_ADJ(N)}).values())}", flush=True)
        out[N] = {}
        for lab, m1 in (("K3'", T.M1), ("hatted", T.M1_HAT)):
            f, fc, g = crosscap(N, m1)
            if lab == "K3'":
                print(f"N={N} g_(1) + g^ref_[1](U = t^N) = {sp.cancel(g[(1,)] + g_ref().subs({rT: r, sT: s, UT: s**(2 * N)}))}", flush=True)
            st = {R: (T.stats(g[R])["laurent"], T.stats(g[R])["integer"], T.stats(2 * g[R])["integer"]) for R in g}
            print(f"N={N} m_(1) = {m1}（{lab}）：g の (Laurent, 整数, 2g 整数) = {st}、g_(2)+g_(1,1) 整数 = {T.stats(g[(2,)] + g[(1, 1)])['integer']}", flush=True)
            for kind, e in (("f", f), ("f^cov", fc), ("g", g)):
                print(f"    {kind} の p_2 成分：1+q で割り切れる {p2_divisible(e, r**2 + 1)}、1+t で割り切れる {p2_divisible(e, s**2 + 1)}", flush=True)
            gh = Minv * sp.Matrix([g[(2,)], g[(1, 1)]])
            print(f"    ĝ = M^{{-1}}(q^{{-1}}) g の分母: {[sp.factor(sp.fraction(sp.cancel(gh[i]))[1]) for i in range(2)]}", flush=True)
            out[N][lab] = {"m1": str(m1), "g": {str(R): str(v) for R, v in g.items()}, "ghat": [str(sp.cancel(gh[i])) for i in range(2)]}
        out[N]["fcov"] = {str(R): str(v) for R, v in fc.items()}
        print(f"  N={N}: {time.time() - t0:.0f}s", flush=True)
    pickle.dump(out, open("data/refined_kauffman/trefoil_fixedN_2box_g.pkl", "wb"))
