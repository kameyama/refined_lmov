r"""固定 N（N = 2, 3）の refined SO(2N) trefoil の 3 箱（probe 多数、Y の Schur 基底 |R| ≤ 3）と、被覆 U(N_A)（N_A = 2N − 1）。
SO 側は型 A の経路で計算する：D_3 = A_3、D_2 = A_1 × A_1 なので refined SO(6) = refined SU(4)、SO(4) = SU(2)²（Macdonald 多項式と測度が一致）。
  Π_{α,v}(1 − y_α v)（v はベクトル表現の重み）= Σ_λ (−1)^{|λ|} s_λ(Y) s_{λᵀ}[V]、V = e_2(z)（SU(4)）、V = p_1(Z)p_1(W)（SU(2)²）。
  s_{λᵀ}[V] を Macdonald P で展開し（Schur → P は安定、三角形の解法）、ℓ > n の項を落として全列を外す。⟨P_μ⟩ = W_μ = m_μ R_μ P_μ(t^ρ)。
被覆：Π_{α,i}(1 − y_α x_i)(1 − y_α/x_i) = Σ_{λ,μ} (−1)^{|λ|+|μ|} s_λ(Y)s_μ(Y) s_{λᵀ}(x)s_{μᵀ}(1/x)。
  ⟨P_α(x)P_β(1/x)⟩：β = ∅ は W_α、(1;1) は W_adj + κ、(α;1)（|α| = 2）は W_{(α;1)} + d_α W_(1)（d_α は主特殊化から、構造は数値で確認）。
  荷電共役で W_μ̄ = W_μ。composite (α;1) の SU(N_A) の分割は α + 1、(1^{N_A−2}) を付けたもの。
自由エネルギーは全 d の Möbius 反転 F = Σ_d μ(d)/d ψ_d(log Z)（ψ_d：r → r^d、s → s^d、p_k → p_{dk}）。
有理関数は rcsknots.flint_rf.RF（sympy の gcd は次数 6 で失敗する）。"""
import sys, time, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
from fractions import Fraction
from functools import lru_cache
import sympy as sp
from rcsknots.flint_rf import RF, R, S, ONE, conv, from_sympy
from rcsknots.macdonald import macdonald_P, _p_in_m
from rcsknots.characters import chi
from rcsknots.partitions import partitions, z, conjugate
from rcsknots.bps import B
import rk_k5_trefoil_2box as T

Q, Tt = R**2, S**2
W_ = S / R                      # (t/q)^{1/2}
ZERO = conv(0)


# ---- 冪和の多項式（{分割: 係数}） ----
def padd(acc, b, c=1):
    for k, v in b.items():
        acc[k] = acc.get(k, ZERO) + c * v
    return acc


def pmul(a, b, dmax=None):
    out = {}
    for x, u in a.items():
        for y, v in b.items():
            k = tuple(sorted(x + y, reverse=True))
            if dmax is not None and sum(k) > dmax:
                continue
            out[k] = out.get(k, ZERO) + u * v
    return out


def schur_p(lam):
    n = sum(lam)
    return {nu: conv(Fraction(chi(lam, nu), z(nu))) for nu in partitions(n)} if n else {(): ONE}


@lru_cache(maxsize=None)
def P_m(n):
    """{λ: {μ: P_λ の m_μ の係数}}（RF）。"""
    P = macdonald_P(n, Q, Tt, conv)
    Rm = _p_in_m(n)
    return {lam: {mu: sum((P[lam].get(nu, ZERO) * Rm[nu][mu] for nu in partitions(n)), ZERO) for mu in partitions(n)} for lam in partitions(n)}


def to_P(F, n):
    """冪和の次数 n の F → {λ: c}（F = Σ c_λ P_λ）。P_λ = m_λ + 辞書式で下の項、の三角形の解法。"""
    if n == 0:
        return {(): F.get((), ZERO)}
    Rm = _p_in_m(n)
    Fm = {mu: sum((F.get(nu, ZERO) * Rm[nu][mu] for nu in partitions(n)), ZERO) for mu in partitions(n)}
    PM = P_m(n)
    out = {}
    for lam in partitions(n):          # 辞書式の降順
        c = Fm[lam]
        if c.is_zero():
            continue
        out[lam] = c
        for mu, v in PM[lam].items():
            Fm[mu] = Fm[mu] - c * v
    return out


# ---- SO 側（型 A の経路） ----
def _e2_pleth(m):
    return {(m, m): conv(Fraction(1, 2)), (2 * m,): conv(Fraction(-1, 2))}


@lru_cache(maxsize=None)
def so_kernel(N, k):
    """{λ ⊢ k: {表現: 係数}}。N = 3：表現は SU(4) の分割、N = 2：SU(2)² の (a, b)。"""
    out = {}
    for lam in partitions(k):
        lt = conjugate(lam)
        acc = {}
        if N == 3:
            G = {}
            for rho in partitions(k):
                c = conv(Fraction(chi(lt, rho), z(rho)))
                term = {(): ONE}
                for part in rho:
                    term = pmul(term, _e2_pleth(part))
                padd(G, term, c)
            for mu, c in to_P(G, 2 * k).items():
                if len(mu) > 4:
                    continue
                red = tuple(x - (mu[3] if len(mu) == 4 else 0) for x in mu)
                red = tuple(x for x in red if x)
                acc[red] = acc.get(red, ZERO) + c
        else:
            for rho in partitions(k):
                c = conv(Fraction(chi(lt, rho), z(rho)))
                b = to_P({rho: ONE}, k)
                for ka, u in b.items():
                    for kb, v in b.items():
                        if len(ka) > 2 or len(kb) > 2:
                            continue
                        a1 = ka[0] - (ka[1] if len(ka) == 2 else 0)
                        b1 = kb[0] - (kb[1] if len(kb) == 2 else 0)
                        acc[(a1, b1)] = acc.get((a1, b1), ZERO) + c * u * v
        out[lam] = {key: v for key, v in acc.items() if not v.is_zero()}
    return out


@lru_cache(maxsize=None)
def RP(lab, n):
    """R_μ P_μ(t^ρ)（SU(n)、data の鏡像）。μ = () は 1。"""
    if not lab:
        return ONE
    return from_sympy(T.Rdata("".join(map(str, lab)), n)) * from_sympy(T.P_A(lab, n))


@lru_cache(maxsize=None)
def P_only(lab, n):
    return from_sympy(T.P_A(lab, n)) if lab else ONE


def so_W(N, rep, fr, unknot=False):
    if N == 3:
        base = P_only(rep, 4) if unknot else RP(rep, 4)
    else:
        a, b = rep
        f = (lambda x: P_only((x,) if x else (), 2)) if unknot else (lambda x: RP((x,) if x else (), 2))
        base = f(a) * f(b)
    return base if unknot else fr[rep] * base


def so_Z(N, fr, kmax=3, unknot=False):
    """Y の冪和基底の Z~（次数 ≤ kmax）。"""
    Zp = {(): ONE}
    for k in range(1, kmax + 1):
        for lam, coeffs in so_kernel(N, k).items():
            Zl = sum((c * so_W(N, rep, fr, unknot) for rep, c in coeffs.items()), ZERO)
            padd(Zp, schur_p(lam), (-1)**k * Zl)
    return Zp


# ---- 被覆 ----
@lru_cache(maxsize=None)
def schur_to_P(kappa):
    return to_P(schur_p(kappa), sum(kappa)) if kappa else {(): ONE}


def comp_shift(alpha, NA):
    """composite (α;1) の SU(N_A) の分割。"""
    return tuple(sorted([a + 1 for a in alpha] + [1] * (NA - 1 - len(alpha)), reverse=True))


def cover_parts(N, mc, unknot=False):
    NA = 2 * N - 1
    one = (1,)

    def Word(a):
        # 被覆は U(N_A)：ℓ(a) = N_A（N_A = 3 の (1,1,1) = det）も data の U(N_A) の値をそのまま使う
        if not a:
            return ONE
        return P_only(a, NA) if unknot else RP(a, NA)

    def Wcomp(alpha):
        sh = comp_shift(alpha, NA)
        if unknot:
            return P_only(sh, NA)
        return mc[alpha] * RP(sh, NA)

    kappa = (1 - Tt**NA) * (1 - Q) / ((1 - Tt) * (1 - Q * Tt**(NA - 1)))
    d = {al: (P_only(al, NA) * P_only(one, NA) - P_only(comp_shift(al, NA), NA)) / P_only(one, NA) for al in ((2,), (1, 1))}

    def Wc(a, b):
        if not a and not b:
            return ONE
        if not b:
            return Word(a)
        if not a:
            return Word(b)
        if a == one and b == one:
            return Wcomp(one) + kappa
        if b == one:
            return Wcomp(a) + d[a] * Word(one)
        if a == one:
            return Wcomp(b) + d[b] * Word(one)
        raise ValueError((a, b))
    return Wc, kappa, d


def cover_Z(N, mc, kmax=3, unknot=False):
    Wc, _, _ = cover_parts(N, mc, unknot)
    Zp = {(): ONE}
    for n in range(1, kmax + 1):
        for i in range(0, n + 1):
            for lam in (partitions(i) if i else [()]):
                for mu in (partitions(n - i) if n - i else [()]):
                    coefY = pmul(schur_p(lam), schur_p(mu))
                    E = ZERO
                    for al, u in schur_to_P(conjugate(lam) if lam else ()).items():
                        for be, v in schur_to_P(conjugate(mu) if mu else ()).items():
                            E = E + u * v * Wc(al, be)
                    padd(Zp, coefY, (-1)**n * E)
    return Zp


# ---- plethystic log（全 d）と Schur 基底 ----
def plog(Zp, kmax=3):
    X = {k: v for k, v in Zp.items() if k}
    X2 = pmul(X, X, kmax)
    X3 = pmul(X2, X, kmax)
    L = padd(padd(dict(X), X2, conv(Fraction(-1, 2))), X3, conv(Fraction(1, 3)))
    F = dict(L)
    L1 = L.get((1,), ZERO)
    for d_, mu_ in ((2, -1), (3, -1)):
        if d_ <= kmax:
            key = (d_,)
            F[key] = F.get(key, ZERO) + conv(Fraction(mu_, d_)) * L1.adams(d_)
    return F


def to_schur(F, k):
    return {lam: sum((F.get(rho, ZERO) * chi(lam, rho) for rho in partitions(k)), ZERO) for lam in partitions(k)}


def phi(rho):
    """M(q^{-1}) の冪和基底での固有値 φ_ρ = Σ_σ χ_σ(ρ) B_σ(q^{-1/2})（f_ρ = φ_ρ f̂_ρ）。"""
    rr = sp.Symbol("r", positive=True)
    return from_sympy(sum(chi(sg, rho) * B(sg, 1 / rr) for sg in partitions(sum(rho))))


def p_components(fR, k):
    """Schur 基底の {R: f_R} → 冪和の係数 {ρ: f_ρ}（Σ f_R s_R = Σ f_ρ p_ρ）。"""
    return {rho: sum((fR[lam] * conv(Fraction(chi(lam, rho), z(rho))) for lam in partitions(k)), ZERO) for rho in partitions(k)}


# ---- framing と 3 箱の判定 ----
def so_framing(N, m1, new):
    """既知の単項式（2 箱・1 probe の規則）と 3 箱の新しい表現の単項式 new = (m_a, m_b)。W_ = (t/q)^{1/2}。"""
    W = W_
    if N == 3:
        return {(): ONE, (1, 1): m1, (2, 1, 1): W**4, (2, 2): W**4, (2,): m1, (2, 2, 2): W**4 * m1, (3, 3): new[0], (3, 2, 1): new[1]}
    return {(0, 0): ONE, (1, 1): m1, (2, 0): W**2, (0, 2): W**2, (2, 2): W**4, (3, 3): new[0], (3, 1): new[1], (1, 3): new[1]}


def cover_framing(N):
    NA = 2 * N - 1
    return {(1,): W_**NA, (2,): W_**NA, (1, 1): W_**NA}


def Z3(Zp):
    return to_schur({k: v for k, v in Zp.items() if sum(k) == 3}, 3)


def scan_pairs(build, rng=range(-8, 17)):
    """build(m_a, m_b) → Zp。Z の 3 次の全成分を Laurent にする (β_a, β_b)（m = W_^β、t = q で 1 の族）。"""
    import itertools
    A0 = Z3(build(ZERO, ZERO))
    Aa = {k: v - A0[k] for k, v in Z3(build(ONE, ZERO)).items()}
    Ab = {k: v - A0[k] for k, v in Z3(build(ZERO, ONE)).items()}
    return [(a, b) for a, b in itertools.product(rng, rng)
            if all((A0[k] + W_**a * Aa[k] + W_**b * Ab[k]).is_laurent() for k in A0)]


def so_new_framing(N, m1):
    hits = scan_pairs(lambda ma, mb: so_Z(N, so_framing(N, m1, (ma, mb))))
    return hits


def crosscap3(N, m1, new):
    F = plog(so_Z(N, so_framing(N, m1, new)))
    Fc = plog(cover_Z(N, cover_framing(N)))
    G = {k: F.get(k, ZERO) - conv(Fraction(1, 2)) * Fc.get(k, ZERO) for k in set(F) | set(Fc)}
    return F, Fc, G


FACT = {"z_q": R**2 - 1, "1+q": R**2 + 1, "1+q+q²": R**4 + R**2 + 1, "z_t": S**2 - 1, "1+t": S**2 + 1, "1+t+t²": S**4 + S**2 + 1}


def hatted(E, k):
    """Schur 基底の f̂_R = Σ_S M^{-1}_{RS}(q^{-1}) f_S（冪和の成分を φ_ρ で割る）。"""
    pc = p_components(to_schur(E, k), k)
    return {R_: sum((pc[rho] / phi(rho) * chi(R_, rho) for rho in partitions(k)), ZERO) for R_ in partitions(k)}


def hatted_by_rho(E, k):
    pc = p_components(to_schur(E, k), k)
    return {rho: (pc[rho] / phi(rho)).is_laurent() for rho in partitions(k)}


def at_t_eq_q(e):
    rr, ss = sp.symbols("r s", positive=True)
    return from_sympy(e.to_sympy().subs(ss, rr))


if __name__ == "__main__":
    import math, random
    from rcsknots.macdonald import macdonald_P as _mP
    import rk_k5_trefoil_2box_cover as C
    out = {}
    W = W_
    for N in (2, 3):
        t0 = time.time()
        NA = 2 * N - 1
        # 検証：unknot、2 箱の再現、composite の分解
        fu, fcu = plog(so_Z(N, None, unknot=True)), plog(cover_Z(N, None, unknot=True))
        print(f"N={N} unknot：|R| = 2, 3 で SO・被覆とも f = 0: {all(v.is_zero() for k in (2, 3) for v in list(to_schur(fu, k).values()) + list(to_schur(fcu, k).values()))}", flush=True)
        fr2 = {k: v for k, v in so_framing(N, from_sympy(T.M1), (ZERO, ZERO)).items()}
        F2 = plog(so_Z(N, fr2, kmax=2), kmax=2)
        old = C.so_f(N)
        Fc2 = plog(cover_Z(N, cover_framing(N), kmax=2), kmax=2)
        oldc = C.cover_f(N, {"adj": C.M_ADJ(N)})
        print(f"N={N} 2 箱：SO 側（型 A の経路）が D_N の経路と一致 {all((to_schur(F2, sum(R_))[R_] - from_sympy(old[R_])).is_zero() for R_ in old)}、"
              f"被覆が 2 箱の実装と一致 {all((to_schur(Fc2, sum(R_))[R_] - from_sympy(oldc[R_])).is_zero() for R_ in oldc)}", flush=True)
        random.seed(N)
        q0, t0_ = Fraction(9, 49), Fraction(4, 25)
        Pn = {n: _mP(n, q0, t0_, Fraction) for n in range(1, 2 * NA + 1)}
        ev = lambda lam, x: sum(c * math.prod(sum(xi**k for xi in x) for k in nu) for nu, c in Pn[sum(lam)][lam].items())
        _, _, dd = cover_parts(N, {}, unknot=True)
        ok = True
        for al in ((2,), (1, 1)):
            x = [Fraction(random.randint(2, 9), random.randint(2, 9)) for _ in range(NA)]
            dv = sp.nsimplify(dd[al].to_sympy().subs({sp.Symbol("r", positive=True): sp.Rational(3, 7), sp.Symbol("s", positive=True): sp.Rational(2, 5)}))
            ok &= ev(al, x) * ev((1,), [1 / xi for xi in x]) == ev(comp_shift(al, NA), x) / math.prod(x) + Fraction(int(dv.p), int(dv.q)) * ev((1,), x)
        print(f"N={N} composite の分解 P_α(x)P_1(1/x) = P_(α;1) + d_α P_1（数値）: {ok}", flush=True)
        # 被覆の新しい単項式
        hc = scan_pairs(lambda ma, mb: cover_Z(N, {(1,): W**NA, (2,): ma, (1, 1): mb}))
        print(f"N={N} 被覆：(2;1), (11;1) の単項式 (t/q)^(β/2) の β = {hc}（(t/q)^(N_A/2) は β = {NA}）", flush=True)
        # m_(1) = (t/q)^(β1/2) ごとに：新しい SO の単項式、g の整数性、ĝ の Laurent 性
        out[N] = {}
        for b1 in range(-2, 7):
            m1 = W**b1
            hits = so_new_framing(N, m1)
            if len(hits) != 1:
                print(f"N={N} β1 = {b1}: SO の新しい単項式が一意でない {hits}", flush=True)
                continue
            new = tuple(W**b for b in hits[0])
            F, Fc, G = crosscap3(N, m1, new)
            g3 = to_schur(G, 3)
            gp = p_components(g3, 3)
            hat = {rho: (gp[rho] / phi(rho)).is_laurent() for rho in partitions(3)}
            hat_q = {rho: (at_t_eq_q(gp[rho]) / at_t_eq_q(phi(rho))).is_laurent() for rho in partitions(3)}
            so_hat = all(v.is_integer_laurent() for k in (2, 3) for v in hatted(F, k).values())
            print(f"N={N} m_(1) = (t/q)^({b1}/2)：新しい SO の β = {hits[0]}、SO の f̂（|R| ≤ 3）は整数 Laurent {so_hat}、g_R（|R|=3）は整数 Laurent {all(v.is_integer_laurent() for v in g3.values())}、"
                  f"ĝ_ρ の Laurent 性 {hat}、t = q で {hat_q}", flush=True)
            if b1 == 2:
                for k in (2, 3):
                    print(f"    被覆の f̂^cov_ρ の Laurent 性（|R|={k}）: {hatted_by_rho(Fc, k)}、SO の f̂_ρ: {hatted_by_rho(F, k)}", flush=True)
            if b1 in (1, 2):
                out[N][b1] = {"new": hits[0], "g3": {str(k): str(v.to_sympy()) for k, v in g3.items()},
                              "p_div": {str(rho): {nm: gp[rho].num_divisible(x) for nm, x in FACT.items()} for rho in partitions(3)}}
                print(f"    p 成分の割り切れ: {out[N][b1]['p_div']}", flush=True)
        print(f"  N={N}: {time.time() - t0:.0f}s", flush=True)
    pickle.dump(out, open("data/refined_kauffman/trefoil_fixedN_3box_g.pkl", "wb"))
