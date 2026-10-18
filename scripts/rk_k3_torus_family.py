r"""A：K3' の trefoil [1] の h/g 分解を T(2,2m+1)（m = 1..5）に広げる。
  W_SO = Z_V(T(2,2m+1))（AS §6 の一般式、unknot で正規化）× mono × unknot_SO、
  ½h^{ref} = 𝒫_m(𝐚_U, 𝐪, 𝐭) の鏡像 × (ν − ν^{−1})/(𝐪 − 𝐪^{−1})、𝐚_U² = A (t/q)^{1/2}、A = ν²（K3' と同じ）、
  𝒫_m（DGR の T_{2,2m+1}）= 𝐚^{2m} Σ_{i=0}^{m} 𝐪^{−2m+4i} 𝐭^{2i} + 𝐚^{2m+2} Σ_{i=0}^{m−1} 𝐪^{−2m+2+4i} 𝐭^{2i+3}（m = 1 で K3' の P_DGR）。
  g^{ref} = W_SO − ½h^{ref} が Laurent になる単項式 mono = ±r^c s^d U^e を走査し、仮説 T_[1]^{−n}(q/t)^{1/2}（n = 2(2m+1)、T_[1] = q^{1/2}t^{N−1}）と比べる。
  条件：N := g (s − 1/s) が s = ±1 で 0（g の分母の候補は s − 1/s だけ）。"""
import sys
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp
from rcsknots.kauffman_trefoil import Z_V_T2, unknot_SO, P_DGR, r, s, U, bq, bt, ba, is_laurent

nu = ba


def P_T2(m, A2, Q, T):
    return A2**m * sum(Q**(-2 * m + 4 * i) * T**(2 * i) for i in range(m + 1)) + \
        A2**(m + 1) * sum(Q**(-2 * m + 2 + 4 * i) * T**(2 * i + 3) for i in range(m))


def h_half(m):
    A2U = nu**2 * s / r
    return P_T2(m, 1 / A2U, 1 / bq, 1 / bt) * (nu - 1 / nu) / (bq - 1 / bq)


def scan(m, crange=range(-30, 4), erange=range(-24, 1), drange=range(-4, 50)):
    Z = sp.expand(Z_V_T2(m))
    P = sp.expand(P_T2(m, s / (nu**2 * s / r) * 1 / s, 1 / bq, 1 / bt))      # 鏡像の 𝒫（𝐚_U² = A s/r の逆数）
    found = []
    for sg in (1, -1):
        for c in crange:
            for e in erange:
                N1 = sg * Z * r**c * U**e * ((nu - 1 / nu) + (s - 1 / s)) - P * (nu - 1 / nu)
                if sp.expand(N1.subs(s, 1)) != 0:
                    continue
                for d in drange:
                    N = sg * Z * r**c * s**d * U**e * ((nu - 1 / nu) + (s - 1 / s)) - P * (nu - 1 / nu)
                    if sp.expand(N.subs(s, -1)) == 0 and is_laurent(N / (s - 1 / s)):
                        found.append((sg, c, d, e))
    return found


if __name__ == "__main__":
    assert sp.expand(P_T2(1, *sp.symbols("a q t")) - P_DGR(*sp.symbols("a q t"))) == 0
    for m in range(1, 6):
        n = 2 * (2 * m + 1)
        pred = (1, 1 - n, 2 * n - 1, -n)
        f = scan(m)
        print(f"T(2,{2 * m + 1}): Laurent になる (符号, c, d, e) = {f}   仮説 T^{{-{n}}}(q/t)^{{1/2}} = {pred} {'含む' if pred in f else '含まない'}", flush=True)


def kauffman_q_power(m):
    """t = q で W = Z_V·unknot が r^x 倍で r → −1/r 不変（z = r − 1/r の多項式、Kauffman 多項式の形）になる x。"""
    nuS = sp.Symbol("nuS", positive=True)
    Wq = sp.cancel((Z_V_T2(m) * unknot_SO).subs(U, nuS * s).subs(s, r))
    for x in range(-60, 61):
        e = Wq * r**x
        if sp.cancel(e - e.subs(r, -1 / r)) == 0:
            return x
    return None


def g_family(m, d):
    n = 2 * (2 * m + 1)
    mono = r**(-(3 * m + 2)) * s**d * U**(-n)
    return sp.cancel(Z_V_T2(m) * mono * unknot_SO - h_half(m))


def refined_only_part(m):
    """g^{ref}（太字 𝐪 = Q, 𝐭 = T, 𝐚 = A）の 𝐚 の奇数冪の部分を (𝐚²−𝐪²)(𝐭²−1) · (−1)^{m+1} 𝐚^{−(4m+1)} 𝐪^{−(2m−1)} 𝐭^{−3m} · Π_m に分け、Π_m を返す。"""
    Q, Tb, A = sp.symbols("Q T A")
    gb = sp.expand(sp.cancel(g_family(m, 7 * m + 4).subs({r: -Tb * Q, U: A * Q}).subs(s, Q)))
    odd = sp.expand((gb - gb.subs(A, -A)) / 2)
    pre = (A**2 - Q**2) * (Tb**2 - 1) * (-1)**(m + 1) * A**(-(4 * m + 1)) * Q**(-(2 * m - 1)) * Tb**(-3 * m)
    return sp.expand(sp.cancel(odd / pre)), (Q, Tb, A)


def c1_blocks(m):
    """c = 1 部分（G, W, ν；G = z_q z_t、W = (q/t)^{1/4}）の ν^β ごとのブロックを分類する：
    'coherent±'（全係数が同符号）か 'tfactor±'（(W⁴ − 1) = (𝐭² − 1) で割り切れ、商の係数が同符号）か 'mixed'。"""
    from rk_k4_crosscap import to_VW, split_parity, in_G, W, nu, G
    E, _ = split_parity(to_VW(g_family(m, 7 * m + 4)))
    P1 = sp.expand(in_G(E))
    out = {}
    for beta in sorted({t.as_coeff_exponent(nu)[1] for t in sp.Add.make_args(P1)}):
        blk = sp.expand(P1.coeff(nu, beta))
        cs = sp.Poly(sp.expand(blk * W**80), G, W).coeffs()
        if all(c > 0 for c in cs) or all(c < 0 for c in cs):
            out[beta] = "coherent" + ("+" if cs[0] > 0 else "−")
            continue
        quo = sp.cancel(blk / (W**4 - 1))
        if sp.fraction(quo)[1].as_poly(G, W).is_monomial:
            qc = sp.Poly(sp.expand(quo * W**80), G, W).coeffs()
            if all(c > 0 for c in qc) or all(c < 0 for c in qc):
                out[beta] = "tfactor" + ("+" if qc[0] > 0 else "−")
                continue
        out[beta] = "mixed"
    return out


def bold_blocks(m):
    """G 展開の前：g^{ref}（太字 𝐪 = Q, 𝐭 = T, 𝐚 = A）の 𝐚 の偶数冪（c = 1）を 𝐚^β ごとに、c1_blocks と同じ基準で分類する（因子は 𝐭² − 1）。"""
    Q, Tb, A = sp.symbols("Q T A")
    gb = sp.expand(sp.cancel(g_family(m, 7 * m + 4).subs({r: -Tb * Q, U: A * Q}).subs(s, Q)))
    even = sp.expand((gb + gb.subs(A, -A)) / 2)
    out = {}
    for beta in sorted({t.as_coeff_exponent(A)[1] for t in sp.Add.make_args(even)}):
        blk = sp.expand(even.coeff(A, beta))
        cs = sp.Poly(sp.expand(blk * Q**80 * Tb**80), Q, Tb).coeffs()
        if all(c > 0 for c in cs) or all(c < 0 for c in cs):
            out[beta] = "coherent" + ("+" if cs[0] > 0 else "−"); continue
        quo = sp.cancel(blk / (Tb**2 - 1))
        if sp.Poly(sp.fraction(quo)[1], Q, Tb).is_monomial:
            qc = sp.Poly(sp.expand(quo * Q**80 * Tb**80), Q, Tb).coeffs()
            if all(c > 0 for c in qc) or all(c < 0 for c in qc):
                out[beta] = "tfactor" + ("+" if qc[0] > 0 else "−"); continue
        out[beta] = "mixed"
    return out


def special_points(m):
    """(SO(2) での値 g^{ref}(𝐚 = ±𝐪)、𝐭 → −𝐭 の偶奇の符号、すなわち g(−𝐭)/g(𝐭))。予想：((−𝐭)^m, (−1)^m)。"""
    Q, Tb, A = sp.symbols("Q T A")
    gb = sp.expand(sp.cancel(g_family(m, 7 * m + 4).subs({r: -Tb * Q, U: A * Q}).subs(s, Q)))
    so2 = {sp.factor(sp.expand(gb.subs(A, sg * Q))) for sg in (1, -1)}
    par = sp.cancel(gb.subs(Tb, -Tb) / gb)
    return so2, par, (Q, Tb, A)
