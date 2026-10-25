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
    from ..characters import chi
    from ..partitions import partitions, z as z_lam
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
    from ..partitions import conjugate
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
    from ..characters import chi
    from ..partitions import partitions, z as z_lam
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


# ---------------------------------------------------------------------------
# 拡張版（|λ| ≤ nmax）：K_λ を S_μ = Π_i S_{μ_i}（|μ| ≤ |λ|）の多項式として評価の対称性で決める。
#   K_λ = m^{GL}_λ(y) + Σ_{ν ⊢ |λ|, ν <lex λ} c_ν m^{GL}_ν(y) + Σ_{|μ| < |λ|} a_μ S_μ、y = (x_i, x_i^{−1})（m^{GL} は 2N 変数の単項対称関数、p_ρ(y) = S_ρ）。
#   未知数の数 = それまでに決めた K の数（すべての小さい形と、同じ大きさで lex で小さい形）で、その点 t^ρ q^ν で
#   K_λ(t^ρ q^ν) = K_λ(t^ρ) K_ν(t^ρ q^λ) / K_ν(t^ρ) を課す。
def _order(nmax):
    from ..partitions import partitions
    out = []
    for n in range(nmax + 1):
        out += list(reversed(partitions(n)))      # lex の昇順
    return out


def _eval_S(poly, nu, cache):
    """poly = {μ: 係数}（S_μ の線形結合）を t^ρ q^ν で評価。"""
    tot = K.zero
    for mu, c in poly.items():
        v = K.one
        for d in mu:
            key = (d, nu)
            if key not in cache:
                cache[key] = S(d, nu)
            v *= cache[key]
        tot += c * v
    return tot


def koornwinder_general(nmax, log=None, store=None, workers=1):
    """{λ: {μ: 係数}}：K_λ の S_μ 表示（|λ| ≤ nmax）、と K_λ(t^ρ)。
    store：保存先のディレクトリ（λ ごとに保存し、あれば読み込んで飛ばす）。workers > 1 なら行列式を並列に計算する。"""
    import os, pickle
    from ..macdonald import m_in_p
    order = _order(nmax)
    Kr = {lam: K_rho(lam) for lam in order}
    Kfun, cache = {(): {(): K.one}}, {}
    for idx, lam in enumerate(order):
        if lam == ():
            continue
        path = os.path.join(store, "K_" + "_".join(map(str, lam)) + ".pkl") if store else None
        if path and os.path.exists(path):
            Kfun[lam] = {mu: K.from_expr(sp.sympify(v, locals={"q": qs, "t": ts, "U": Us})) for mu, v in pickle.load(open(path, "rb")).items()}
            if log:
                log(f"{lam} 読み込み")
            continue
        n = sum(lam)
        top = m_in_p(n)
        # D_N のパラメータは x → −x で不変なので K_λ(−x) = (−1)^{|λ|} K_λ(x)：|μ| ≡ |λ| (mod 2) の S_μ だけが現れる
        prev = [nu for nu in order[:idx] if (sum(nu) - n) % 2 == 0]   # 評価点 = 未知数の数
        same = [nu for nu in prev if sum(nu) == n]
        lower = [mu for mu in order if sum(mu) < n and (sum(mu) - n) % 2 == 0]
        base = {rho: K(c) for rho, c in top[lam].items()}
        unknown_funcs = [{rho: K(c) for rho, c in top[nu].items()} for nu in same] + [{mu: K.one} for mu in lower]
        assert len(unknown_funcs) == len(prev), (lam, len(unknown_funcs), len(prev))
        M, rhs = [], []
        for nu in prev:
            target = Kr[lam] * _eval_S(Kfun[nu], lam, cache) / Kr[nu]
            rhs.append(target - _eval_S(base, nu, cache))
            M.append([_eval_S(f, nu, cache) for f in unknown_funcs])
        from .flint_solve import solve as _fsolve, solve_parallel as _psolve
        if len(M) <= 2:
            sol = _solve(M, rhs)
        else:
            sol = _psolve(M, rhs, K, workers) if workers > 1 else _fsolve(M, rhs, K)
        poly = dict(base)
        for f, c in zip(unknown_funcs, sol):
            for mu, v in f.items():
                poly[mu] = poly.get(mu, K.zero) + c * v
        Kfun[lam] = {mu: v for mu, v in poly.items() if v != 0}
        if path:
            pickle.dump({mu: str(v.as_expr()) for mu, v in Kfun[lam].items()}, open(path, "wb"))
        if log:
            log(f"{lam} 済")
    return Kfun, Kr, cache


# ---------------------------------------------------------------------------
# 拡張版の S2：被覆側（U(N_A)、A = ν² 形式）。成分 2 の brane・image を Z' = Z とおいた次数 ≤ AZ、成分 1 も Y' = Y。
#   composite の P_{[α,β]} = P_α(x)P_β(1/x) + Σ_{k≥1} Σ c · T_{μ'ν'}（T_{μν} = p_μ(x)p_ν(1/x)、|μ'| = |α|−k、|ν'| = |β|−k）
#   の c を評価の対称性 P_λ(t^ρ q^μ)P_μ(t^ρ) = P_μ(t^ρ q^λ)P_λ(t^ρ)（composite の点）で決める。P_{[α,β]}(t^ρ) は dim_composite。
def _comp_point_pp(d, gam, dlt):
    """x = t^ρ q^{[γ,δ]} での (Σ x_i^d, Σ x_i^{−d})。"""
    return _pp(d, tuple(gam), tuple(-x for x in dlt))


def _T_eval(mu, nu, gam, dlt, cache):
    v = KC.one
    for d in mu:
        key = (d, gam, dlt)
        if key not in cache:
            cache[key] = _comp_point_pp(d, gam, dlt)
        v *= cache[key][0]
    for d in nu:
        key = (d, gam, dlt)
        if key not in cache:
            cache[key] = _comp_point_pp(d, gam, dlt)
        v *= cache[key][1]
    return v


def _dim_comp(alpha, beta):
    from .kauffman_refined import dim_composite, q as qR, t as tR, A as AR
    from .kauffman_cover import qs_, ss_, nus_
    e = dim_composite(list(alpha), list(beta))
    e = e.subs({qR: qs_, tR: ss_**2, AR: nus_**2})
    return KC.from_expr(sp.powdenest(sp.expand_power_base(e, force=True), force=True))


def composite_P(nmax):
    """{(α, β): {(μ, ν): 係数}}（|α| + |β| ≤ nmax）と P_{[α,β]}(t^ρ)。"""
    from ..partitions import partitions
    from ..macdonald import macdonald_P
    from .kauffman_cover import qs_ as _q, ss_ as _s, nus_ as _n
    conv = lambda c: KC(c)
    Pm = {n: macdonald_P(n, qc, sc**2, conv) for n in range(nmax + 1)}
    comps = [(a, b) for n in range(nmax + 1) for na in range(n + 1) for a in partitions(na) for b in partitions(n - na)]
    Pr = {c: (KC.one if c == ((), ()) else _dim_comp(*c)) for c in comps}
    P, cache = {}, {}
    for (al, be) in comps:
        top = {}
        for mu, x in Pm[sum(al)][al].items() if al else [((), KC.one)]:
            for nu, y in Pm[sum(be)][be].items() if be else [((), KC.one)]:
                top[(mu, nu)] = top.get((mu, nu), KC.zero) + x * y
        lowers = [(g, d) for k in range(1, min(sum(al), sum(be)) + 1)
                  for g in partitions(sum(al) - k) for d in partitions(sum(be) - k)]
        if not lowers:
            P[(al, be)] = top
            continue
        M, rhs = [], []
        for (g, d) in lowers:                       # 評価点 = 下位の composite
            target = Pr[(al, be)] * sum((c * _T_eval(mu, nu, al, be, cache) for (mu, nu), c in P[(g, d)].items()), KC.zero) / Pr[(g, d)]
            rhs.append(target - sum((c * _T_eval(mu, nu, g, d, cache) for (mu, nu), c in top.items()), KC.zero))
            M.append([_T_eval(mu, nu, g, d, cache) for (mu, nu) in lowers])
        from .flint_solve import solve as _fsolve
        sol = _fsolve(M, rhs, KC)
        poly = dict(top)
        for (mu, nu), c in zip(lowers, sol):
            poly[(mu, nu)] = poly.get((mu, nu), KC.zero) + c
        P[(al, be)] = poly
    return P, Pr, cache


def _mul2f(P, Q, AY, AZ, zero):
    R = {}
    for (a1, b1), v1 in P.items():
        for (a2, b2), v2 in Q.items():
            if sum(a1) + sum(a2) <= AY and sum(b1) + sum(b2) <= AZ:
                key = (tuple(sorted(a1 + a2, reverse=True)), tuple(sorted(b1 + b2, reverse=True)))
                R[key] = R.get(key, zero) + v1 * v2
    return R


def log2f(Z, AY, AZ, field):
    X = {k: v for k, v in Z.items() if k != ((), ())}
    L, P = {}, dict(X)
    for m in range(1, AY + AZ + 1):
        for k, v in P.items():
            L[k] = L.get(k, field.zero) + v * field(sp.Rational((-1)**(m + 1), m))
        P = _mul2f(P, X, AY, AZ, field.zero)
    return L


def _exp_ps(coef, A, field):
    """exp(Σ_d coef[d] p_d) の {μ: 係数}（|μ| ≤ A）。"""
    out = {(): field.one}
    for d, c in coef.items():
        term = {(): field.one}
        acc = field.one
        for m in range(1, A // d + 1):
            acc = acc * c / m
            term[(d,) * m] = acc
        new = {}
        for a, va in out.items():
            for b, vb in term.items():
                if sum(a) + sum(b) <= A:
                    k = tuple(sorted(a + b, reverse=True))
                    new[k] = new.get(k, field.zero) + va * vb
        out = new
    return out


def cover_Z_LR(AY, AZ):
    """被覆の Z^{cov}(Y, Y'=Y; Z, Z'=Z)：{(μY, μZ): KC}（Σ N N f^{cov} を直接与える母関数）。"""
    from ..partitions import partitions, z as z_lam
    P, Pr, cache = composite_P(AZ)
    comps = list(P)
    # T_{μν} を composite の P で表す（ブロックごとに上位の行列を逆に解き、下位を再帰的に代入）
    Tin = {}
    blocks = sorted({(sum(a), sum(b)) for a, b in comps})
    for (na, nb) in blocks:
        cs = [(a, b) for a, b in comps if sum(a) == na and sum(b) == nb]
        ts = [(mu, nu) for mu in (partitions(na) if na else [()]) for nu in (partitions(nb) if nb else [()])]
        X = [[P[c].get(tk, KC.zero) for tk in ts] for c in cs]          # 行：composite、列：上位の T
        # P_c − (下位の T の部分) = Σ_tk X[c][tk] T_tk  ⇒ T_tk = Σ_c Xinv[tk][c] (P_c − 下位_c)
        from .flint_solve import solve as _fsolve
        for j, tk in enumerate(ts):
            e = [KC.one if i == j else KC.zero for i in range(len(ts))]
            col = _fsolve([[X[c_i][t_i] for c_i in range(len(cs))] for t_i in range(len(ts))], e, KC) if len(ts) > 1 else [KC.one / X[0][0]]
            expr = {}
            for c_i, coeff in enumerate(col):
                if coeff == 0:
                    continue
                c = cs[c_i]
                expr[c] = expr.get(c, KC.zero) + coeff
                for (mu, nu), v in P[c].items():
                    if (sum(mu), sum(nu)) != (na, nb):
                        for c2, w in Tin[(mu, nu)].items():
                            expr[c2] = expr.get(c2, KC.zero) - coeff * v * w
            Tin[tk] = {c: v for c, v in expr.items() if v != 0}
    # 成分 2 の核：Σ_ρ ε_ρ p_ρ(Z)/z_ρ Π_i (p_{ρ_i}(x) + p_{ρ_i}(1/x))
    import itertools
    kern = {}                                                            # {composite: {ρ_Z: 係数}}
    for n in range(AZ + 1):
        for rho in (partitions(n) if n else [()]):
            pref = KC((-1)**len(rho)) / KC(z_lam(rho)) if rho else KC.one
            for choice in itertools.product((0, 1), repeat=len(rho)):
                mu = tuple(sorted((d for d, c in zip(rho, choice) if c == 0), reverse=True))
                nu = tuple(sorted((d for d, c in zip(rho, choice) if c == 1), reverse=True))
                for comp, w in Tin[(mu, nu)].items():
                    kern.setdefault(comp, {})
                    kern[comp][rho] = kern[comp].get(rho, KC.zero) + pref * w
    Z = {}
    for comp, cz in kern.items():
        al, be = comp
        pp = {d: _comp_point_pp(d, al, be) for d in range(1, AY + 1)}
        ey = _exp_ps({d: -(pp[d][0] + pp[d][1]) / d for d in pp}, AY, KC)
        for muY, vy in ey.items():
            for rhoZ, c in cz.items():
                if c == 0:
                    continue
                key = (muY, rhoZ)
                Z[key] = Z.get(key, KC.zero) + Pr[comp] * c * vy
    return Z


def cover_f_LR(AY, AZ):
    """Σ N^R_{ρρ'} N^S_{σσ'} f^{cov}_{ρρ'σσ'}（sympy の式、記号 q, s, nu）。"""
    from .kauffman_cover import qs_, ss_, nus_
    Z = cover_Z_LR(AY, AZ)
    assert Z[((), ())] == KC.one
    L = log2f(Z, AY, AZ, KC)
    ad = lambda e, d: e.subs({qs_: qs_**d, ss_: ss_**d, nus_: nus_**d}, simultaneous=True)
    return plethystic2(L, AY, AZ, ad)


def so_Z_general(AY, AZ, store="data/refined_kauffman/koornwinder_DN"):
    """SO 側 Z~(Y; Z)、|μ_Y| ≤ AY、|μ_Z| ≤ AZ（K_λ は koornwinder_general の保存を使う）。"""
    from ..partitions import partitions, z as z_lam
    Kf, Kr, cache = koornwinder_general(AZ, store=store, workers=8)
    order = [lam for lam in _order(AZ)]
    # S_ρ を K_λ で表す：K_λ = Σ_μ a_{λμ} S_μ を大きさの順に逆に解く
    Sin = {(): {(): K.one}}
    for n in range(1, AZ + 1):
        lams = [lam for lam in order if sum(lam) == n]
        rhos = list(reversed(partitions(n)))
        X = [[Kf[lam].get(rho, K.zero) for rho in rhos] for lam in lams]
        from .flint_solve import solve as _fsolve
        for j, rho in enumerate(rhos):
            e = [K.one if i == j else K.zero for i in range(len(rhos))]
            col = _fsolve([[X[l_i][r_i] for l_i in range(len(lams))] for r_i in range(len(rhos))], e, K) if len(rhos) > 1 else [K.one / X[0][0]]
            expr = {}
            for l_i, coeff in enumerate(col):
                if coeff == 0:
                    continue
                lam = lams[l_i]
                expr[lam] = expr.get(lam, K.zero) + coeff
                for mu, v in Kf[lam].items():
                    if sum(mu) < n:
                        for lam2, w in Sin[mu].items():
                            expr[lam2] = expr.get(lam2, K.zero) - coeff * v * w
            Sin[rho] = {lam: v for lam, v in expr.items() if v != 0}
    kern = {}
    for n in range(AZ + 1):
        for rho in (partitions(n) if n else [()]):
            pref = K((-1)**len(rho)) / K(z_lam(rho)) if rho else K.one
            for lam, w in Sin[rho].items():
                kern.setdefault(lam, {})
                kern[lam][rho] = kern[lam].get(rho, K.zero) + pref * w
    Z = {}
    for lam, cz in kern.items():
        ey = _exp_powersums({d: -S(d, lam) / d for d in range(1, AY + 1)}, AY)
        for (muY, _), vy in ey.items():
            for rhoZ, c in cz.items():
                if c == 0:
                    continue
                key = (muY, rhoZ)
                Z[key] = Z.get(key, K.zero) + Kr[lam] * c * vy
    return Z


def so_f_general(AY, AZ):
    Z = so_Z_general(AY, AZ)
    assert Z[((), ())] == K.one
    L = log2f(Z, AY, AZ, K)
    ad = lambda e, d: e.subs({qs: qs**d, ts: ts**d, Us: Us**d}, simultaneous=True)
    return plethystic2(L, AY, AZ, ad)
