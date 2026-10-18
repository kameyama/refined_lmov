r"""最小版（2026-10-18）：Hopf 絡み目で両成分とも多 probe（成分 2 は次数 ≤ 2）。SO 側（D_N、U = t^N を形式変数）。

成分 2 の核（virtual 正規化の fermionic な核）：Π_{α,i}(1 − z_α x_i)(1 − z_α/x_i) = exp(−Σ_d p_d(Z) S_d(x)/d)、S_d(x) = Σ_i (x_i^d + x_i^{−d})。
次数 ≤ 2：1 − p_1 S_1 + ½ p_1² S_1² − ½ p_2 S_2。これを D_N の Koornwinder 多項式 K_λ（(±1, ±q^{1/2})）で展開し、
  ⟨K_λ(U_2) F(U_1)⟩ = K_λ(t^ρ) F(t^ρ q^λ)（評価の対称性、K2'(b) と同じ）
で成分 1（多 probe、exp(−Σ_d p_d(Y) S_d(x)/d)）と組み合わせる。x = t^ρ q^λ の x_i = t^{N−i} q^{λ_i}。
K_(1) = S_1 + κ、K_(1,1) = ½(S_1² − S_2) + b_1 S_1 + b_0、K_(2) = (1 − c/2) S_2 + (c/2) S_1² + c_1 S_1 + c_0 の未定係数を、
ずらした点 t^ρ q^ν（ν = ∅, (1), (1,1)）での値（評価の対称性）から決める。N が単独で現れない形なので U を形式変数にできる。"""
import sympy as sp

from .kauffman_hopf import K, q, t, U, qs, ts, Us, _sym_mul, _exp_powersums
from .kauffman_refined import koornwinder_eval
from .kauffman_refined import q as _qs, t as _ts, U as _Us

_TBAR = (1, -1, sp.sqrt(_qs), -sp.sqrt(_qs))


def K_rho(lam):
    """K_λ(t^ρ)（U の有理式、K の元）。"""
    e = sp.cancel(koornwinder_eval(tuple(x for x in lam if x) if lam else (), _TBAR)).subs({_qs: qs, _ts: ts, _Us: Us})
    return K.from_expr(e)


def S(d, lam):
    """S_d(t^ρ q^λ) = Σ_{i=1}^N (x_i^d + x_i^{−d})、x_i = t^{N−i} q^{λ_i}（λ は上の数行だけ非零）。"""
    tot = (K.one - U**d) / (K.one - t**d) + (K.one - U**-d) / (K.one - t**-d)
    for i, li in enumerate(lam, start=1):
        if li:
            base = U * t**(-i)                      # t^{N−i}
            tot += (base * q**li)**d + (base * q**li)**-d - base**d - base**-d
    return tot


def _solve(M, rhs):
    from .kauffman_hopf import _solve_field
    return _solve_field([list(row) for row in M], rhs)


def koornwinder_small():
    """K_(1), K_(1,1), K_(2) の (S_1, S_2) 表示の係数と、S_1², S_2 の K 展開。"""
    Kr = {lam: K_rho(lam) for lam in ((), (1,), (1, 1), (2,))}
    pts = [(), (1,), (1, 1), (2,)]
    S1 = {nu: S(1, nu) for nu in pts}
    S2 = {nu: S(2, nu) for nu in pts}
    kappa = Kr[(1,)] - S1[()]
    K1 = lambda nu: S1[nu] + kappa                                   # K_(1)(t^ρ q^ν)
    # K_(1,1)：ν = ∅, (1) の二点
    rhs11 = [Kr[(1, 1)], Kr[(1, 1)] * K1((1, 1)) / Kr[(1,)]]
    M11 = [[S1[nu], K.one] for nu in ((), (1,))]
    b1, b0 = _solve(M11, [rhs11[i] - (S1[nu]**2 - S2[nu]) / 2 for i, nu in enumerate(((), (1,)))])
    K11 = lambda nu: (S1[nu]**2 - S2[nu]) / 2 + b1 * S1[nu] + b0
    # 検算：K_(1,1)(t^ρ q^{(2)}) = K_(1,1)(t^ρ) K_(2)(t^ρ q^{(1,1)}) / K_(2)(t^ρ) は K_(2) を求めた後に行う
    # K_(2)：ν = ∅, (1), (1,1)
    nus = ((), (1,), (1, 1))
    rhs2 = [Kr[(2,)], Kr[(2,)] * K1((2,)) / Kr[(1,)], Kr[(2,)] * K11((2,)) / Kr[(1, 1)]]
    M2 = [[(S1[nu]**2 - S2[nu]) / 2, S1[nu], K.one] for nu in nus]
    c, c1, c0 = _solve(M2, [rhs2[i] - S2[nu] for i, nu in enumerate(nus)])
    K2 = lambda nu: S2[nu] + c * (S1[nu]**2 - S2[nu]) / 2 + c1 * S1[nu] + c0
    check = K11((2,)) - Kr[(1, 1)] * K2((1, 1)) / Kr[(2,)]
    # 逆：S_1² = K_(2) + A K_(1,1) + B K_(1) + C、S_2 = K_(2) + A' K_(1,1) + B' K_(1) + C'
    #   K_(2) = (1 − c/2) S_2 + (c/2) S_1² + c1 S_1 + c0、K_(1,1) = ½ S_1² − ½ S_2 + b1 S_1 + b0、K_(1) = S_1 + κ
    # 基底 (S_1², S_2, S_1, 1) の係数で解く
    basis_K = {(2,): [c / 2, 1 - c / 2, c1, c0], (1, 1): [K.one / 2, -K.one / 2, b1, b0], (1,): [K.zero, K.zero, K.one, kappa], (): [K.zero, K.zero, K.zero, K.one]}
    order = [(2,), (1, 1), (1,), ()]
    Mt = [[basis_K[lam][j] for lam in order] for j in range(4)]       # 列が K_λ、行が基底の成分
    sq = _solve(Mt, [K.one, K.zero, K.zero, K.zero])                  # S_1²
    s2 = _solve(Mt, [K.zero, K.one, K.zero, K.zero])                  # S_2
    return dict(Kr=Kr, kappa=kappa, b=(b1, b0), c=(c, c1, c0), check=check,
                S1sq={lam: v for lam, v in zip(order, sq)}, S2={lam: v for lam, v in zip(order, s2)})


def so_Z2(AY):
    """Z~(Y; Z)：{(μ_Y, μ_Z): K の元}、|μ_Y| ≤ AY、|μ_Z| ≤ 2（冪和の分割）。"""
    ks = koornwinder_small()
    kap = ks["kappa"]
    # 成分 2 の核の K_λ 展開の係数（Z の冪和の多項式）：{λ: {μ_Z: 係数}}
    coeff = {(): {(): K.one, (1,): kap, (1, 1): ks["S1sq"][()] / 2, (2,): -ks["S2"][()] / 2},
             (1,): {(1,): -K.one, (1, 1): ks["S1sq"][(1,)] / 2, (2,): -ks["S2"][(1,)] / 2},
             (1, 1): {(1, 1): ks["S1sq"][(1, 1)] / 2, (2,): -ks["S2"][(1, 1)] / 2},
             (2,): {(1, 1): ks["S1sq"][(2,)] / 2, (2,): -ks["S2"][(2,)] / 2}}
    # 1 + p_1 κ：S_1 = K_(1) − κ なので −p_1 S_1 = −p_1 K_(1) + p_1 κ
    Z = {}
    for lam, cz in coeff.items():
        ey = _exp_powersums({d: -S(d, lam) / d for d in range(1, AY + 1)}, AY)
        for (muY, _), vy in ey.items():
            for muZ, c in cz.items():
                key = (muY, muZ)
                Z[key] = Z.get(key, K.zero) + ks["Kr"][lam] * c * vy
    return Z, ks


def _mul2(P, Q, AY, AZ):
    R = {}
    for (a1, b1), v1 in P.items():
        for (a2, b2), v2 in Q.items():
            if sum(a1) + sum(a2) <= AY and sum(b1) + sum(b2) <= AZ:
                key = (tuple(sorted(a1 + a2, reverse=True)), tuple(sorted(b1 + b2, reverse=True)))
                R[key] = R.get(key, K.zero) + v1 * v2
    return R


def log2(Z, AY, AZ, zero=None, one=None):
    zero = K.zero if zero is None else zero
    one = K.one if one is None else one
    X = {k: v for k, v in Z.items() if k != ((), ())}
    L, P = {}, dict(X)
    for m in range(1, AY + AZ + 1):
        for k, v in P.items():
            L[k] = L.get(k, zero) + v * sp.Rational((-1)**(m + 1), m)
        P = _mul2(P, X, AY, AZ)
    return L


def plethystic2(L, AY, AZ, adams):
    """log Z = Σ_d (1/d) Σ_{R,S} f_{R,S}(x^d) s_R(Y^d) s_S(Z^d) の f（sympy の式）。adams(e, d) は変数の d 乗。"""
    from .characters import chi
    from .partitions import partitions, z as z_lam
    Le = {k: v.as_expr() for k, v in L.items()}
    f = {}
    for n1 in range(AY + 1):
        for n2 in range(AZ + 1):
            if n1 + n2 == 0:
                continue
            H = {(m1, m2): Le.get((m1, m2), 0) for m1 in partitions(n1) for m2 in partitions(n2)}
            for d in sp.divisors(int(sp.gcd(n1, n2))):
                if d == 1:
                    continue
                for R in partitions(n1 // d):
                    for Sx in partitions(n2 // d):
                        fv = f.get((R, Sx))
                        if not fv:
                            continue
                        fa = adams(fv, d) / d
                        for m1 in partitions(n1 // d):
                            for m2 in partitions(n2 // d):
                                c = sp.Rational(chi(R, m1), z_lam(m1)) * sp.Rational(chi(Sx, m2), z_lam(m2))
                                H[(tuple(d * x for x in m1), tuple(d * x for x in m2))] -= fa * c
            for R in partitions(n1):
                for Sx in partitions(n2):
                    f[(R, Sx)] = sp.cancel(sum(chi(R, m1) * chi(Sx, m2) * H[(m1, m2)] for (m1, m2) in H))
    return f


def so_f2(AY):
    Z, ks = so_Z2(AY)
    L = log2(Z, AY, 2)
    ad = lambda e, d: e.subs({qs: qs**d, ts: ts**d, Us: Us**d}, simultaneous=True)
    return plethystic2(L, AY, 2, ad)


# ---------------------------------------------------------------------------
# 被覆側（U(N_A)、A = ν² を形式変数、体 KC = Q(q, s, ν)）。成分 2 は brane 側 Z と image 側 Z' の多 probe で全次数 ≤ 2。
#   核 = exp(−Σ_d [p_d(Z) p_d(x) + p_d(Z') p_d(1/x)]/d) = Σ_μ coef_μ(Z, Z') P_μ(x)（μ は composite）
#   p_1(x)² = P_(2) + (2 − c) P_(1,1)、p_2(x) = P_(2) − c P_(1,1)（c = (1+q)(1−t)/(1−qt)）、p_1(x) p_1(1/x) = ψ_adj P_adj + ψ_0
#   ⟨P_μ(U_2) F(U_1)⟩ = P_μ(t^ρ) F(t^ρ q^μ)、x_i = ν s^{1−2i} q^{μ_i}。
from .kauffman_cover import KC, qc, sc, nuc, _psi_sym, _principal_sym, _exp_two, _mul4


def _pp(d, top, bottom):
    """(Σ_i x_i^d, Σ_i x_i^{−d})、x_i = ν s^{1−2i} q^{μ_i}；top = 上の行の μ、bottom = 下の行の μ（下から）。"""
    out = []
    for e in (d, -d):
        base = nuc**e * sc**(-e) * (1 - nuc**(-2 * e)) / (1 - sc**(-2 * e))
        for i, m in enumerate(top, start=1):
            if m:
                base += (nuc * sc**(1 - 2 * i))**e * (qc**(m * e) - 1)
        for k, m in enumerate(bottom):
            if m:
                base += (nuc**-1 * sc**(1 + 2 * k))**e * (qc**(m * e) - 1)
        out.append(base)
    return out


def _P_rho(lam):
    """対称な正規化の主特殊化 (t/A)^{|λ|/2} Π (t^{l'} − A q^{a'})/(1 − q^a t^{l+1})（A = ν²）。"""
    from .partitions import conjugate
    t = sc**2; A = nuc**2
    out = (sc / nuc)**sum(lam)
    conj = conjugate(lam)
    for i in range(len(lam)):
        for j in range(lam[i]):
            out *= (t**i - A * qc**j) / (1 - qc**(lam[i] - j - 1) * t**(conj[j] - i))
    return out


def cover_Z2(AY):
    """{(μY, μY', μZ, μZ'): KC}、|μY| + |μY'| ≤ AY、|μZ| + |μZ'| ≤ 2。"""
    t = sc**2
    c = (1 + qc) * (1 - t) / (1 - qc * t)
    psi_adj, psi_0 = _psi_sym(1, 1, 1), _psi_sym(1, 0, 0)
    one = KC.one
    # 成分 2 の composite：(top, bottom, P(t^ρ), coef{(μZ, μZ'): 係数})
    reps = [
        ((), (), one, {((), ()): one, ((1,), (1,)): psi_0}),
        ((1,), (), _P_rho((1,)), {((1,), ()): -one}),
        ((), (-1,), _P_rho((1,)), {((), (1,)): -one}),
        ((2,), (), _P_rho((2,)), {((1, 1), ()): one / 2, ((2,), ()): -one / 2}),
        ((1, 1), (), _P_rho((1, 1)), {((1, 1), ()): (2 - c) / 2, ((2,), ()): c / 2}),
        ((), (-2,), _P_rho((2,)), {((), (1, 1)): one / 2, ((), (2,)): -one / 2}),
        ((), (-1, -1), _P_rho((1, 1)), {((), (1, 1)): (2 - c) / 2, ((), (2,)): c / 2}),
        ((1,), (-1,), _principal_sym(1, 1), {((1,), (1,)): psi_adj}),
    ]
    Z = {}
    for top, bottom, Prho, coef in reps:
        pp = {d: _pp(d, top, bottom) for d in range(1, AY + 1)}
        ey = _exp_two({d: -pp[d][0] / d for d in pp}, {d: -pp[d][1] / d for d in pp}, AY)
        for (m1, m2, _, _), v in ey.items():
            for (mz, mz2), cf in coef.items():
                key = (m1, m2, mz, mz2)
                Z[key] = Z.get(key, KC.zero) + Prho * cf * v
    return Z


def _mul4b(P, Q, AY, AZ):
    R = {}
    for k1, v1 in P.items():
        for k2, v2 in Q.items():
            if sum(k1[0]) + sum(k1[1]) + sum(k2[0]) + sum(k2[1]) <= AY and sum(k1[2]) + sum(k1[3]) + sum(k2[2]) + sum(k2[3]) <= AZ:
                key = tuple(tuple(sorted(a + b, reverse=True)) for a, b in zip(k1, k2))
                R[key] = R.get(key, KC.zero) + v1 * v2
    return R


def cover_f2(AY):
    """f^{cov}_{ρ,ρ',σ,σ'}（4 方向の plethystic な反転、sympy の式、記号 q, s, nu）。"""
    from .characters import chi
    from .partitions import partitions, z as z_lam
    from .kauffman_cover import qs_, ss_, nus_
    Z = cover_Z2(AY)
    one = ((), (), (), ())
    assert Z[one] == KC.one
    X = {k: v for k, v in Z.items() if k != one}
    L, P = {}, dict(X)
    for m in range(1, AY + 3):
        for k, v in P.items():
            L[k] = L.get(k, KC.zero) + KC(sp.Rational((-1)**(m + 1), m)) * v
        P = _mul4b(P, X, AY, 2)
    Le = {k: v.as_expr() for k, v in L.items()}
    ad = lambda e, d: e.subs({qs_: qs_**d, ss_: ss_**d, nus_: nus_**d}, simultaneous=True)
    f = {}
    import itertools
    for degs in itertools.product(range(AY + 1), range(AY + 1), range(3), range(3)):
        if degs[0] + degs[1] > AY or degs[2] + degs[3] > 2 or sum(degs) == 0:
            continue
        H = {ms: Le.get(ms, 0) for ms in itertools.product(*[partitions(n) for n in degs])}
        g = 0
        for x in degs:
            g = sp.gcd(g, x)
        for d in sp.divisors(int(g)):
            if d == 1:
                continue
            for Rs in itertools.product(*[partitions(n // d) for n in degs]):
                fv = f.get(Rs)
                if not fv:
                    continue
                fa = ad(fv, d) / d
                for ms in itertools.product(*[partitions(n // d) for n in degs]):
                    cc = 1
                    for Rk, mk in zip(Rs, ms):
                        cc *= sp.Rational(chi(Rk, mk), z_lam(mk))
                    H[tuple(tuple(d * x for x in mk) for mk in ms)] -= fa * cc
        for Rs in itertools.product(*[partitions(n) for n in degs]):
            tot = 0
            for ms, h in H.items():
                cc = 1
                for Rk, mk in zip(Rs, ms):
                    cc *= chi(Rk, mk)
                tot += cc * h
            f[Rs] = sp.cancel(tot)
    return f
