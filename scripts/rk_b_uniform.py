r"""タスク B の判定式を一様な規格化（data の R_μ × (t/q)^{δ|μ|/2}、conventions.md「一様な規格化」）で作り直す。
[1]：g^{true}_[1] = (t/q)^{δ/2} g^{ref}_[1]（T(2,2m+1)、δ = m、m = 1..5）について
  (N) SO(2) の点 𝐚² = 𝐪²（U = t）での値、(P) 𝐭 → −𝐭 の偶奇、(c=1) G 展開の ν^β ごとの分類、(c=2) 普遍的な因子と符号。
固定 N の 2 箱・3 箱（ν = t^{N−1/2} が q, t と混ざるので G 展開はできない）：g_R の r → −r（q^{1/2} → −q^{1/2}）の偶奇を、規格化を変えて比べる。"""
import sys, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
from fractions import Fraction
import sympy as sp
import rk_k3_torus_family as Tf
from rk_k4_crosscap import to_VW, split_parity, in_G, V, W, nu, G
from rcsknots.kauffman_trefoil import r, s, U

Q, Tb, A = sp.symbols("Q T A")


def g_true(m):
    return sp.cancel(Tf.g_family(m, 7 * m + 4) * (s / r)**m)


def bold(e):
    return sp.expand(sp.cancel(e.subs({r: -Tb * Q, U: A * Q}).subs(s, Q)))


def classify_c1(m):
    """c = 1 部分を G 展開し、ν^β ごとに 'coherent±' か 'tfactor±'（W⁴ − 1 = 𝐭² − 1 で割って同符号）か 'mixed'。"""
    E, _ = split_parity(to_VW(g_true(m)))
    P1 = sp.expand(in_G(E))
    out = {}
    for beta in sorted({t.as_coeff_exponent(nu)[1] for t in sp.Add.make_args(P1)}):
        blk = sp.expand(P1.coeff(nu, beta))
        cs = sp.Poly(sp.expand(blk * W**200), G, W).coeffs()
        if all(c > 0 for c in cs) or all(c < 0 for c in cs):
            out[beta] = "coherent" + ("+" if cs[0] > 0 else "−"); continue
        quo = sp.cancel(blk / (W**4 - 1))
        if sp.Poly(sp.fraction(quo)[1], G, W).is_monomial:
            qc = sp.Poly(sp.expand(quo * W**200), G, W).coeffs()
            if all(c > 0 for c in qc) or all(c < 0 for c in qc):
                out[beta] = "tfactor" + ("+" if qc[0] > 0 else "−"); continue
        out[beta] = "mixed"
    return out


def check_one(m):
    g = g_true(m)
    so2 = {sp.simplify(g.subs(U, sg * s * s)) for sg in (1, -1)}
    gb = bold(g)
    even_t = sp.expand(gb - gb.subs(Tb, -Tb)) == 0
    odd = sp.expand((gb - gb.subs(A, -A)) / 2)
    Pi, (Q0, T0, A0) = Tf.refined_only_part(m)
    pre = -(A**2 - Q**2) * (Tb**2 - 1) * A**(-(4 * m + 1)) * Q**(-(2 * m - 1)) * Tb**(-4 * m)
    c2 = sp.expand(odd - pre * Pi.subs({Q0: Q, T0: Tb, A0: A})) == 0
    return so2, even_t, c2, Pi.subs({Q0: Q, T0: Tb, A0: A})


def fixedN_parity(N):
    import rk_k5_trefoil_3box as X
    W_ = X.W_
    rr = sp.Symbol("r", positive=True)
    cases = {"一様": (X.so_uniform(N), {}, True),
             "被覆 m=1・SO 一様": (X.so_uniform(N), X.cover_framing(N), False),
             "SO K3'・被覆 m=1": (X.so_framing(N, W_, (W_**5, W_**5 if N == 3 else W_**3)), X.cover_framing(N), False),
             "SO K3'・被覆一様": (X.so_framing(N, W_, (W_**5, W_**5 if N == 3 else W_**3)), {}, True)}
    out = {}
    for lab, (fr, mc, uni) in cases.items():
        F = X.plog(X.so_Z(N, fr)); Fc = X.plog(X.cover_Z(N, mc, uniform=uni))
        Gd = {k: F.get(k, X.ZERO) - X.conv(Fraction(1, 2)) * Fc.get(k, X.ZERO) for k in set(F) | set(Fc)}
        par = []
        for k in (1, 2, 3):
            for R_, v in X.to_schur(Gd, k).items():
                e = v.to_sympy()
                par.append(sp.cancel(e - e.subs(rr, -rr)) == 0)
        out[lab] = (sum(par), len(par))
        if lab == "一様":
            hat = {}
            for k in (2, 3):
                for R_, v in X.hatted(Gd, k).items():
                    e = v.to_sympy()
                    hat[R_] = sp.cancel(e - (-1)**(k - 1) * e.subs(rr, -rr)) == 0
            out["ĝ_R(−r) = (−1)^(|R|−1) ĝ_R(r)"] = hat
    return out


if __name__ == "__main__":
    for m in range(1, 6):
        t0 = time.time()
        so2, even_t, c2, Pi = check_one(m)
        print(f"m={m}: SO(2) の点の値 {so2}、𝐭 について偶 {even_t}、c=2 = −(𝐚²−𝐪²)(𝐭²−1)𝐚^-(4m+1)𝐪^-(2m-1)𝐭^-4m Π_m {c2}（Π_m の係数 {set(sp.Poly(Pi, Q, Tb, A).coeffs())}）", flush=True)
        print(f"   c=1 の G 展開の分類: {classify_c1(m)}（{time.time()-t0:.0f}s）", flush=True)
    for N in (2, 3):
        t0 = time.time()
        print(f"N={N} 固定 N の g_R（|R| ≤ 3）の r → −r での不変性（成り立つ数/全数）: {fixedN_parity(N)}（{time.time()-t0:.0f}s）", flush=True)
