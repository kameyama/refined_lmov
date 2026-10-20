r"""案 D：C_N の DAHA-Jones JD(ω_1)（N = 2..8、u は記号）から Cherednik の hyper-polynomial ℋ^C_{r,2}(ω_1; q,t,u,a) を復元し、
  u = 1（SO(2N)）で既知の g^{true}_[1] を再現する規格化を固定したうえで、u = t²（AM の refined Sp(2N)）の g^{Sp}_[1] を作る（progress.md 続き 50）。
  (1) JD(N)·q^α t^{β(N−1)} u^γ = ℋ(a = −t^{N−1})、ℋ は a の多項式（次数 ≤ r、係数は Z[q,t,u]、定数項 1）。(α,β,γ) を走査し、N = 2..7 で当てはめ、N = 8 で照合。
  (2) SO：W = ℋ(u=1, a = −ν/t^{1/2})(鏡像の有無、単項式) × dim_V(ν)、g = −W + W^{cov}_1(ν)（被覆 U(2N−1)、A = ν²）。既知の g^{true}_[1] と比べて規約を固定。
  (3) Sp：同じ規約で u = t²、dim_V^{Sp}、ν_Sp = t^{N+1/2}（A = t^{2N+1}）。g^{Sp}_[1] と g^{SO}_[1] を比べる。"""
import sys, pickle, itertools
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp

q, t, u, a = sp.symbols("q t u a")
r_, s_, u_ = sp.symbols("r s u", positive=True)
OUT = "data/refined_kauffman/daha"


def jd(N, m):
    d = pickle.load(open(f"{OUT}/cn_omega1_N{N}_m{m}.pkl", "rb"))
    loc = {"r": r_, "s": s_, "u": u_}
    Z = sp.sympify(d["Z"], locals=loc); X1 = sp.sympify(d["X1"], locals=loc)
    return sp.cancel((Z / X1).subs({r_: sp.sqrt(q), s_: sp.sqrt(t), u_: sp.sqrt(u)}))


def fit_hyper(m, Ns=range(2, 8), check=8):
    """ℋ(a) = Σ_{j≤r} h_j(q,t,u) a^j を、a_N = −t^{N−1} での値 JD(N)·M(N) から解く（Vandermonde、r + 1 個の N）。"""
    r = 2 * m + 1
    vals = {N: jd(N, m) for N in list(Ns) + [check]}
    for al, be, ga in itertools.product(range(0, 13), range(0, 13), range(0, 7)):
        M = lambda N: q**al * t**(be * (N - 1)) * u**ga
        nodes = list(Ns)[: r + 1]
        aN = [-t**(N - 1) for N in nodes]
        rhs = [sp.cancel(vals[N] * M(N)) for N in nodes]
        V = sp.Matrix([[aN[i]**j for j in range(r + 1)] for i in range(r + 1)])
        try:
            h = V.LUsolve(sp.Matrix(rhs))
        except Exception:
            continue
        h = [sp.cancel(x) for x in h]
        if not all(sp.Poly(sp.fraction(x)[1], q, t, u).is_monomial for x in h):
            continue
        if sp.cancel(h[0] - 1) != 0:
            continue
        H = sum(h[j] * a**j for j in range(r + 1))
        if sp.cancel(H.subs(a, -t**(check - 1)) - vals[check] * M(check)) != 0:
            continue
        print(f"T({r},2)：M(N) = q^{al} t^{{{be}(N−1)}} u^{ga}、ℋ^C = {sp.expand(H)}", flush=True)
        return H, (al, be, ga)
    raise RuntimeError("当てはまらない")


def dim_V(kind, nu):
    """大きな N の fundamental の次元（ν の関数、変数 r = q^{1/2}, s = t^{1/2}）。
    SO(2N)：Σ_{i=0}^{N−1}(t^i + t^{−i})、U = t^N = ν s。Sp(2N)：Σ_{j=1}^{N}(t^j + t^{−j})、U = t^N = ν/s（A = t^{2N+1} = ν²）。"""
    t2 = s_**2
    if kind == "SO":
        Uv = nu * s_
        return sp.cancel((Uv - 1) / (t2 - 1) + (1 / Uv - 1) / (1 / t2 - 1))
    Uv = nu / s_
    return sp.cancel((t2 * Uv - t2) / (t2 - 1) + (1 / (t2 * Uv) - 1 / t2) / (1 / t2 - 1))


def cover_W1(m, nu):
    """被覆 U(N_A) の [1]：(t/q)^{δ/2} R₁(A) P₁(t^ρ)、A = ν²（rk_c_2box_largeN.Rdata_nu・Psym と同じ）。変数 r, s, ν。"""
    import rk_c_2box_largeN as L
    from rcsknots.flint_rf3 import S, R
    w = (S / R)**m
    W = (w * L.Rdata_nu("1", f"T2_{2*m+1}") * L.Psym((1,))).to_sympy()
    uu = sp.Symbol("u", positive=True)
    return sp.cancel(W.subs(uu, nu * s_))


def g_known(m):
    """既知の g^{true}_[1]（一様な規格化、rk_b_uniform）。変数 r, s, ν。"""
    import rk_b_uniform as Bu
    from rcsknots.kauffman_trefoil import r as r0, s as s0, U as U0
    nu = sp.Symbol("nu", positive=True)
    return sp.cancel(Bu.g_true(m).subs(U0, nu * s0).subs({r0: r_, s0: s_})), nu


def H_rs(H, uval):
    """ℋ(q,t,u,a) を変数 r, s に（u は uval：1 なら SO、s_**4 なら Sp）。"""
    return sp.expand(H.subs({q: r_**2, t: s_**2, u: uval}))


def is_monomial_rs(e, nu):
    n_, d_ = sp.fraction(sp.cancel(e))
    try:
        return sp.Poly(n_, r_, s_, nu).is_monomial and sp.Poly(d_, r_, s_, nu).is_monomial
    except sp.PolynomialError:
        return False


def calibrate_SO(m, H):
    """g = −W + W^{cov}、W = mono · ℋ(u=1, a ↦ a_map(ν)) (鏡像あり/なし) · dim_V^{SO}(ν)。既知の g^{true} と比べ、mono が単項式になる規約を探す。"""
    gk, nu = g_known(m)
    Wc = cover_W1(m, nu)
    HD = H_rs(H, 1)
    hits = []
    for mirror in (False, True):
        for sg, k in itertools.product((1, -1), (1, -1, 3, -3)):
            for inv in (1, -1):
                Hs = HD.subs(a, sg * nu**inv * s_**k)
                if mirror:
                    Hs = Hs.subs({r_: 1 / r_, s_: 1 / s_, nu: 1 / nu}, simultaneous=True)
                W0 = sp.cancel(Hs * dim_V("SO", nu))
                for gsign in (1, -1):                      # g_(1)（全 d の定義）= −g^{true}_[1]（report §9）なので W_SO = W^{cov} + g^{true}
                    mono = sp.cancel((Wc + gsign * gk) / W0)
                    if is_monomial_rs(mono, nu):
                        hits.append((mirror, sg, k, inv, gsign, sp.factor(mono)))
    print(f"T(2,{2*m+1})、SO の規約の候補 (鏡像, a の符号, s の冪, ν の冪, g の符号, 単項式) = {hits}", flush=True)
    return hits


if __name__ == "__main__":
    for m in (1, 2):
        H, M = fit_hyper(m)
        pickle.dump({"H": str(H), "M": M}, open(f"{OUT}/cn_hyper_omega1_m{m}.pkl", "wb"))
        calibrate_SO(m, H)
