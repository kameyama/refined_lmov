r"""大きな N の trefoil の 3 箱のクロスキャップ部分 g = f − ½f^{cov}（progress.md 続き 33）。変数は (r, s, u) = (q^{1/2}, t^{1/2}, U)、rcsknots.refined_kauffman.flint_rf3。
SO：Z~(Y) = Σ_ν (−1)^{|ν|} s_ν(Y) ⟨s_{ν'}(x, x^{−1})⟩、s_{ν'}(z) を S_ρ = Π(p_d(x) + p_d(x^{−1})) で書き、D_N の K_λ（|λ| ≤ 3、
  data/refined_kauffman/koornwinder_DN、U は形式変数）の逆で K 基底に直し、⟨K_λ⟩ = W_λ = R_λ(U) K_λ(t^ρ)（R_λ は R_U_m1_{λ}.pkl）。
被覆 U(N_A)：Z~ = Σ_{λ,μ} (−1)^{|λ|+|μ|} ⟨P_λ(x) P_μ(x^{−1})⟩ P_{λ'}(Y; t, q) P_{μ'}(Y; t, q)。
  ⟨P_λ⟩ = (t/q)^{|λ|/2} R_λ(data) P_λ(t^ρ)、⟨P_1 P_1^{−1}⟩ = W_adj + κ、⟨P_α P_1^{−1}⟩ = W_{(α;1)} + d_α W_(1)、
  W_{(α;1)} = R_{(α;1)}(A^{1/2}) P_{(α;1)}(t^ρ)（R は scripts/rk_c_fast_composite.py の値を A^{1/2} で補間、P は kauffman_refined.dim_composite）、
  d_α = (P_α(t^ρ) P_1(t^ρ) − P_{(α;1)}(t^ρ)) / P_1(t^ρ)。
検算：U = t^N（N = 2, 3）で固定 N の SO・被覆の自由エネルギー（scripts/rk_k5_trefoil_3box.py）と全係数一致、g_R が固定 N の g と一致。"""
import sys, os, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from fractions import Fraction
from rcsknots.refined_kauffman.flint_rf3 import RF, R, S, U, from_sympy, conv
from rcsknots.characters import chi
from rcsknots.partitions import partitions, z as zlam
from rcsknots.macdonald import macdonald_P
import rk_c_2box_largeN as L
import rk_c_interp as I

r_, s_, u_ = sp.symbols("r s u", positive=True)
q, t = R * R, S * S
nu = U / S
OUT = "data/refined_kauffman/daha"
ZERO, ONE = conv(0), conv(1)
PARTS = [(), (1,), (2,), (1, 1), (3,), (2, 1), (1, 1, 1)]


# ---- 冪和基底の演算（キーは Y の冪和の分割） ----
def padd(A, B, c=ONE):
    out = dict(A)
    for k, v in B.items():
        out[k] = out.get(k, ZERO) + c * v
    return out


def pmul(A, B, kmax=3):
    out = {}
    for a, x in A.items():
        for b, y in B.items():
            if sum(a) + sum(b) <= kmax:
                k = tuple(sorted(a + b, reverse=True))
                out[k] = out.get(k, ZERO) + x * y
    return out


def schur_p(lam):
    n = sum(lam)
    return {rho: conv(Fraction(chi(lam, rho), zlam(rho))) for rho in partitions(n)} if n else {(): ONE}


def plog(Zp, kmax=3):
    X = {k: v for k, v in Zp.items() if k}
    X2 = pmul(X, X, kmax)
    X3 = pmul(X2, X, kmax)
    Lg = padd(padd(X, X2, conv(Fraction(-1, 2))), X3, conv(Fraction(1, 3)))
    f = dict(Lg)
    f1 = Lg.get((1,), ZERO)
    f[(2,)] = Lg.get((2,), ZERO) - conv(Fraction(1, 2)) * f1.adams(2)
    f[(3,)] = Lg.get((3,), ZERO) - conv(Fraction(1, 3)) * f1.adams(3)
    return f


def to_schur(F, k):
    return {lam: sum((F.get(rho, ZERO) * conv(chi(lam, rho)) for rho in partitions(k)), ZERO) for lam in partitions(k)}


# ---- SO 側 ----
def K_matrix():
    """K_λ = Σ_μ C[λ][μ] S_μ（|λ| ≤ 3）を読み、S_μ = Σ_λ Cinv[μ][λ] K_λ を返す。"""
    qs, ts, Us = sp.symbols("q t U")
    C = {(): {(): ONE}}
    for lam in PARTS[1:]:
        raw = pickle.load(open("data/refined_kauffman/koornwinder_DN/K_" + "_".join(map(str, lam)) + ".pkl", "rb"))
        C[lam] = {mu: from_sympy(sp.sympify(v, locals={"q": qs, "t": ts, "U": Us}).subs({qs: r_**2, ts: s_**2, Us: u_})) for mu, v in raw.items()}
    # 行列 C（行 λ、列 μ）を掃き出し法で逆にする：S_μ = Σ_λ Cinv[μ][λ] K_λ
    n = len(PARTS)
    M = [[C[lam].get(mu, ZERO) for mu in PARTS] + [ONE if i == j else ZERO for j in range(n)] for i, lam in enumerate(PARTS)]
    for col in range(n):
        piv = next(i for i in range(col, n) if not M[i][col].is_zero())
        M[col], M[piv] = M[piv], M[col]
        inv = 1 / M[col][col]
        M[col] = [x * inv for x in M[col]]
        for i in range(n):
            if i != col and not M[i][col].is_zero():
                f = M[i][col]
                M[i] = [x - f * y for x, y in zip(M[i], M[col])]
    # M の右半分 = C^{-1}（行 μ、列 λ）：K = C S ⇒ S = C^{-1} K
    Sinv = {mu: {lam: M[i][n + j] for j, lam in enumerate(PARTS) if not M[i][n + j].is_zero()} for i, mu in enumerate(PARTS)}
    return Sinv


def so_W():
    from rcsknots.refined_kauffman.kauffman_hopf2 import K_rho
    from rcsknots.refined_kauffman.kauffman_hopf import qs, ts, Us
    sub = {qs: r_**2, ts: s_**2, Us: u_}
    W = {(): ONE}
    for lam in PARTS[1:]:
        c = pickle.load(open(os.path.join(OUT, f"R_U_m1_{''.join(map(str, lam))}.pkl"), "rb"))
        Rf = sum((from_sympy(sp.sympify(v, locals={"r": r_, "s": s_})) * U**j for j, v in c.items()), ZERO)
        kr = K_rho(lam)
        W[lam] = Rf * from_sympy(sp.sympify(kr.as_expr() if hasattr(kr, "as_expr") else kr).subs(sub))
    return W


def so_Z(W, Sinv):
    """Z~(Y) を Y の冪和基底で（次数 ≤ 3）。⟨S_ρ⟩ = Σ_λ Sinv[ρ][λ] W_λ。"""
    Sexp = {rho: sum((c * W[lam] for lam, c in Sinv[rho].items()), ZERO) for rho in PARTS}
    Zp = {(): ONE}
    for n in (1, 2, 3):
        for lam in partitions(n):
            lamT = tuple(sum(1 for x in lam if x > j) for j in range(lam[0]))
            # ⟨s_{λ'}(z)⟩ = Σ_ρ χ_{λ'}(ρ)/z_ρ ⟨S_ρ⟩
            val = sum((conv(Fraction(chi(lamT, rho), zlam(rho))) * Sexp[rho] for rho in partitions(n)), ZERO)
            Zp = padd(Zp, schur_p(lam), conv((-1)**n) * val)
    return Zp


# ---- 被覆側 ----
def comp_interp(alpha):
    """R_{(α;1)}(N) を x = s^N（A^{1/2}）の Laurent 多項式として補間（傾きで窓、補間に使わない N で照合）。"""
    tag = "".join(map(str, alpha))
    avail = [N for N in range(2, 40) if os.path.exists(os.path.join(OUT, f"Rcomp_N{N}_m1_{tag}.pkl"))]
    vals = {N: sp.sympify(pickle.load(open(os.path.join(OUT, f"Rcomp_N{N}_m1_{tag}.pkl"), "rb"))["R"], locals={"r": r_, "s": s_}) for N in avail}
    a, b, c3 = avail[-3:]
    rng = {N: I.s_range(vals[N]) for N in (a, b, c3)}
    jmin, jmax = rng[c3][0] - rng[b][0], rng[c3][1] - rng[b][1]
    assert rng[b][0] - rng[a][0] == jmin and rng[b][1] - rng[a][1] == jmax, rng
    n = jmax - jmin + 1
    nodes = avail[-n:]
    c = I.interpolate_generic({N: vals[N] for N in nodes}, 1, jmin, jmax)
    from rcsknots.refined_kauffman.flint_rf import from_sympy as fs2
    ok = {N: (I.predict_generic(c, 1, N) - fs2(vals[N])).is_zero() for N in avail if N not in nodes}
    print(f"composite ({tag};1)：A^{{1/2}} の冪 {jmin}〜{jmax}（ノード N_A = {nodes}）、照合 {ok}", flush=True)
    pickle.dump({j: str(v.to_sympy()) for j, v in c.items()}, open(os.path.join(OUT, f"Rcomp_A_m1_{tag}.pkl"), "wb"))
    return sum((from_sympy(v.to_sympy()) * nu**j for j, v in c.items()), ZERO)


def cover_parts():
    from rcsknots import kauffman_refined as KR
    w = S / R
    Wo = {lam: w**sum(lam) * L.Rdata_nu("".join(map(str, lam))) * L.Psym(lam) for lam in PARTS[1:]}
    radj = pickle.load(open(os.path.join(OUT, "Radj_A_m1.pkl"), "rb"))
    Radj = sum((from_sympy(sp.sympify(v, locals={"r": r_, "s": s_})) * nu**j for j, v in radj.items()), ZERO)
    A = nu * nu
    qn = (nu - 1 / nu) / (S - 1 / S)
    kappa = (1 - A) * (1 - q) / ((1 - t) * (1 - q * A / t))
    Wadj = Radj * (qn * qn - kappa)
    pair = {((1,), (1,)): Wadj + kappa}
    for alpha in ((2,), (1, 1)):
        Rc = comp_interp(alpha)
        dimc = sp.sympify(KR.dim_composite(alpha, (1,)))
        dimc = from_sympy(sp.cancel(dimc.subs({KR.q: r_**2, KR.t: s_**2, KR.A: (u_ / s_)**2})))
        Wcomp = Rc * dimc
        d = (L.Psym(alpha) * L.Psym((1,)) - dimc) / L.Psym((1,))
        pair[(alpha, (1,))] = Wcomp + d * Wo[(1,)]
        pair[((1,), alpha)] = pair[(alpha, (1,))]
    return Wo, pair


def cover_Z(Wo, pair):
    PY = {}
    for n in (1, 2, 3):
        P = macdonald_P(n, t, q, lambda f: conv(f))                       # P_λ(Y; t, q) の冪和表示
        for lam, v in P.items():
            PY[lam] = v
    PY[()] = {(): ONE}
    T_ = lambda lam: tuple(sum(1 for x in lam if x > j) for j in range(lam[0])) if lam else ()
    Zp = {(): ONE}
    for lam in PARTS:
        for mu in PARTS:
            n = sum(lam) + sum(mu)
            if n == 0 or n > 3:
                continue
            if not lam:
                val = Wo[mu]
            elif not mu:
                val = Wo[lam]
            else:
                val = pair[(lam, mu)]
            term = pmul(PY[T_(lam)], PY[T_(mu)])
            Zp = padd(Zp, term, conv((-1)**n) * val)
    return Zp


def main():
    Sinv = K_matrix()
    W = so_W()
    Zso = so_Z(W, Sinv)
    Wo, pair = cover_parts()
    Zc = cover_Z(Wo, pair)
    F, Fc = plog(Zso), plog(Zc)
    G = {k: F.get(k, ZERO) - conv(Fraction(1, 2)) * Fc.get(k, ZERO) for k in set(F) | set(Fc)}
    # 固定 N との照合
    import rk_k5_trefoil_3box as X
    fx = pickle.load(open("data/refined_kauffman/T2_3_fixedN_3box_g.pkl", "rb"))
    for N in (2, 3):
        Fx = X.plog(X.so_Z(N, X.so_uniform(N)))
        Fcx = X.plog(X.cover_Z(N, {}, uniform=True))
        sub = lambda e: sp.cancel(e.to_sympy().subs(u_, s_**(2 * N)))
        okso = all(sp.cancel(sub(F.get(k, ZERO)) - Fx.get(k, X.ZERO).to_sympy()) == 0 for k in Fx if sum(k) <= 3)
        okc = all(sp.cancel(sub(Fc.get(k, ZERO)) - Fcx.get(k, X.ZERO).to_sympy()) == 0 for k in Fcx if sum(k) <= 3)
        g3 = to_schur(G, 3)
        okg = all(sp.cancel(sub(g3[R_]) - sp.sympify(fx[N]["g_schur"][str(R_)], locals={"r": r_, "s": s_})) == 0 for R_ in g3)
        print(f"N = {N}：SO の f 一致 {okso}、被覆の f 一致 {okc}、g_R（|R| = 3）一致 {okg}", flush=True)
    g3 = to_schur(G, 3)
    for R_, v in g3.items():
        print(f"g_{R_}：分母 {sp.factor(sp.fraction(sp.cancel(v.to_sympy()))[1])}、整数 Laurent {v.is_integer_laurent()}", flush=True)
    pickle.dump({"f": {str(k): str(v.to_sympy()) for k, v in F.items()}, "fcov": {str(k): str(v.to_sympy()) for k, v in Fc.items()},
                 "g": {str(k): str(v.to_sympy()) for k, v in G.items()}}, open(os.path.join(OUT, "largeN_3box_m1.pkl"), "wb"))
    return F, Fc, G


if __name__ == "__main__":
    main()
