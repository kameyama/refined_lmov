r"""(b)：大きな N の 2 箱で、KN 1703.05408 の refined LMOV の二つの枠（src/rcsknots/bps.py の f_formulas と同じ規約）
  q 枠：Σ_λ W_λ g_λ P_λ(x; q, t) = Exp(Σ_μ f^q_μ s_μ(x)/(q^{1/2} − q^{−1/2}))、g_λ = (q/t)^{|λ|}/⟨P_λ, P_λ⟩、f^q に (t/q)^{|μ|/2} を掛ける、
  t 枠：Σ_λ W_λ P_{λ^T}(−x; t, q) = Exp(Σ_μ f^t_μ s_μ(x)/(t^{−1/2} − t^{1/2}))、
  f^q = M(t) f̂、f^t = M(q^{−1}) f̂（同じ f̂）
を、SO（ベクトル表現の probe）と被覆（brane と image）とクロスキャップ部分 g = f − ½f^{cov} に当てはめる。
q 枠の核は Macdonald の Cauchy 核 exp(Σ_d (1/d) c_d (q/t)^d p_d(Y) P_d(X))、c_d = (1 − t^d)/(1 − q^d)、
P_d(X) = S_d（SO）または p_d(x) + p_d(1/x)（被覆）。t 枠は scripts/rk_c_2box_largeN.py の値に (t^{−1/2} − t^{1/2}) を掛けたもの。"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_c_2box_largeN as L
from rcsknots.refined_kauffman.flint_rf3 import R, S, U, conv, from_sympy

half = L.half
q, t = L.q, L.t


def plog2(Z1, Z11, Z2):
    f1 = Z1
    H11 = Z11 - Z1 * Z1 * half
    H2 = Z2 - f1.adams(2) * half
    return {(1,): f1, (2,): H11 + H2, (1, 1): H11 - H2}


def q_frame(P1, P1sq, P2):
    """⟨P_1⟩, ⟨P_1²⟩, ⟨P_2⟩ → q 枠の f^q（KN の規格化：(Q − 1/Q)(T/Q)^{|μ|}）。"""
    v = q / t
    c1 = (1 - t) / (1 - q)
    c2 = (1 - t * t) / (1 - q * q)
    Z1 = c1 * v * P1
    Z11 = half * c1 * c1 * v * v * P1sq
    Z2 = half * c2 * v * v * P2
    F = plog2(Z1, Z11, Z2)
    w = S / R
    return {k: (R - 1 / R) * w**sum(k) * x for k, x in F.items()}


def so_q(W, co):
    X = W[(2,)] + co["S1sq"][(1, 1)] * W[(1, 1)] + co["S1sq"][(1,)] * W[(1,)] + co["S1sq"][()]
    Xp = W[(2,)] + co["S2"][(1, 1)] * W[(1, 1)] + co["S2"][(1,)] * W[(1,)] + co["S2"][()]
    return q_frame(W[(1,)], X, Xp)


def cover_q():
    w = S / R
    Wc = {lam: w**sum(lam) * L.Rdata_nu("".join(map(str, lam))) * L.Psym(lam) for lam in ((1,), (1, 1), (2,))}
    radj = pickle.load(open("data/refined_kauffman/daha/Radj_A_m1.pkl", "rb"))
    Radj = sum((from_sympy(sp.sympify(v, locals={"r": L.r_, "s": L.s_})) * L.nu**j for j, v in radj.items()), conv(0))
    A = L.nu * L.nu
    qn = (L.nu - 1 / L.nu) / (S - 1 / S)
    kappa = (1 - A) * (1 - q) / ((1 - t) * (1 - q * A / t))
    Wadj = Radj * (qn * qn - kappa)
    c = (t - q) / (1 - q * t)
    P1 = 2 * Wc[(1,)]
    p1sq = Wc[(2,)] + (1 + c) * Wc[(1, 1)]                     # ⟨p_1(x)²⟩ = ⟨s_2 + s_11⟩
    P1sq = 2 * p1sq + 2 * (Wadj + kappa)
    P2 = 2 * (Wc[(2,)] + c * Wc[(1, 1)] - Wc[(1, 1)])          # ⟨p_2(x) + p_2(1/x)⟩
    return q_frame(P1, P1sq, P2)


def phiRF(rho, tt_sym):
    from rcsknots.characters import chi
    from rcsknots.partitions import partitions
    from rcsknots.bps import B
    return from_sympy(sum(chi(sg, rho) * B(sg, tt_sym) for sg in partitions(sum(rho))))


def hatted(g, tt_sym):
    p11 = (g[(2,)] + g[(1, 1)]) * half / phiRF((1, 1), tt_sym)
    p2 = (g[(2,)] - g[(1, 1)]) * half / phiRF((2,), tt_sym)
    return {(1,): g[(1,)], (2,): p11 + p2, (1, 1): p11 - p2}


if __name__ == "__main__":
    rr, ssy = sp.symbols("r s", positive=True)
    W, co = L.so_parts()
    zt_ = (1 / S - S)
    ft = {k: zt_ * v for k, v in L.so_f(W, co).items()}
    fct = {k: zt_ * v for k, v in L.cover_f().items()}
    fq = so_q(W, co)
    fcq = cover_q()
    for lab, A_, B_ in (("SO の f", fq, ft), ("被覆の f^cov", fcq, fct)):
        hq, ht = hatted(A_, ssy), hatted(B_, 1 / rr)
        print(f"{lab}：f̂^q = f̂^t {[(hq[k] - ht[k]).is_zero() for k in hq]}", flush=True)
    gq = {k: fq[k] - half * fcq[k] for k in fq}
    gt = {k: ft[k] - half * fct[k] for k in ft}
    hq, ht = hatted(gq, ssy), hatted(gt, 1 / rr)
    print(f"クロスキャップ部分：ĝ^q = ĝ^t {[(hq[k] - ht[k]).is_zero() for k in hq]}", flush=True)
    for k in hq:
        print(f"   ĝ_{k}：分母 {sp.factor(sp.fraction(sp.cancel(ht[k].to_sympy()))[1])}", flush=True)
    pickle.dump({"ghat": {str(k): str(v.to_sympy()) for k, v in ht.items()}, "ghat_q": {str(k): str(v.to_sympy()) for k, v in hq.items()}},
                open("data/refined_kauffman/daha/largeN_2box_m1_ghat_KN.pkl", "wb"))
