r"""K4'：被覆（U(N) refined CS、composite 表現）の Hopf 絡み目の自由エネルギー（有限 N、厳密な有理数）。

Z^{cov}(y, z) = ⟨ Π_i (1 − y x_i)(1 − y/x_i) · Π_i (1 − z x'_i)(1 − z/x'_i) ⟩_{Hopf}
  Π_i (1 − y x_i)(1 − y/x_i) = Σ_{k,l} (−y)^{k+l} e_k(x) e_l(1/x)、e_k e_l(1/x) = e_N^{−1} e_k e_{N−l}
  = Σ_j ψ'_{λ/1^k} P_{λ−1^N}(x)、λ = 2^j 1^{k+N−l−2j}（Pieri、縦帯）。
  ⟨P_μ(U_1) F(U_2)⟩_{Hopf} / S_00 = P_μ(t^ρ) F(t^ρ q^μ)、μ = λ − 1^N（composite）、P_μ(t^ρ) = P_λ(t^ρ)（e_N(t^ρ) = 1）。
変数：s = t^{1/2}（有理数）、q（有理数）、x_i = s^{N+1−2i} q^{μ_i}、ν = A^{1/2} = s^N。
Pieri：ψ'_{λ/μ} = Π_{s ∈ C_{λ/μ}∖R_{λ/μ}} b_λ(s)/b_μ(s)（数値で特定、tests/test_rk_cover.py）。
"""
from fractions import Fraction as F

from .partitions import conjugate


def _b(lam, i, j, q, t):
    a = lam[i] - j - 1
    l = conjugate(lam)[j] - i - 1
    return (1 - q**a * t**(l + 1)) / (1 - q**(a + 1) * t**l)


def pieri_vertical(lam, mu, q, t):
    """P_μ e_r = Σ_λ ψ'_{λ/μ} P_λ（λ/μ は縦帯）の ψ'。"""
    mu2 = tuple(mu) + (0,) * (len(lam) - len(mu))
    strip = [(i, j) for i in range(len(lam)) for j in range(mu2[i], lam[i])]
    rows = {i for i, _ in strip}
    cols = {j for _, j in strip}
    out = F(1)
    for i in range(len(mu)):
        for j in range(mu[i]):
            if j in cols and i not in rows:
                out *= _b(lam, i, j, q, t) / _b(mu, i, j, q, t)
    return out


def principal(lam, N, q, t, s):
    """P_λ(t^ρ)、t^ρ = (s^{N−1}, s^{N−3}, …, s^{1−N})：(t/A)^{|λ|/2} Π_s (t^{l'} − A q^{a'})/(1 − q^a t^{l+1})、A = s^{2N}。"""
    A = s**(2 * N)
    out = F(1)
    conj = conjugate(lam)
    for i in range(len(lam)):
        for j in range(lam[i]):
            out *= (t**i - A * q**j) / (1 - q**(lam[i] - j - 1) * t**(conj[j] - i))
    n = sum(lam)
    return out * s**n / s**(N * n)


def zprod(mu, N, q, s, B):
    """Π_i (1 − z x_i)(1 − z/x_i)、x_i = s^{N+1−2i} q^{μ_i}、z^B まで。"""
    p = {0: F(1)}
    for i in range(1, N + 1):
        x = s**(N + 1 - 2 * i) * q**mu[i - 1]
        X = x + 1 / x
        new = {}
        for b, v in p.items():
            for db, f in ((0, F(1)), (1, -X), (2, F(1))):
                if b + db <= B:
                    new[b + db] = new.get(b + db, F(0)) + v * f
        p = new
    return p


def cover_Z(N, q, s, A_, B):
    """Z^{cov}(y,z) mod (y^{A_+1}, z^{B+1})：{(a,b): Fraction}。"""
    t = s**2
    Z = {}
    for k in range(0, A_ + 1):
        for l in range(0, A_ + 1 - k):
            for j in range(0, k + 1):
                m = k + N - l - 2 * j
                z0 = l - k + j
                if m < 0 or z0 < 0 or N - l - j < 0:
                    continue
                lam = (2,) * j + (1,) * m + (0,) * z0
                lam_p = tuple(x for x in lam if x > 0)
                mu_k = (1,) * k
                psi = pieri_vertical(lam_p, mu_k, q, t)
                w = psi * principal(lam_p, N, q, t, s)
                mu = tuple(x - 1 for x in lam)
                pz = zprod(mu, N, q, s, B)
                for b, vb in pz.items():
                    Z[(k + l, b)] = Z.get((k + l, b), F(0)) + (-1)**(k + l) * w * vb
    return Z


def log_coeffs(Z, A_, B):
    assert Z[(0, 0)] == 1
    X = {ab: v for ab, v in Z.items() if ab != (0, 0)}

    def mul(P, Q):
        R = {}
        for (a1, b1), v1 in P.items():
            for (a2, b2), v2 in Q.items():
                if a1 + a2 <= A_ and b1 + b2 <= B:
                    R[(a1 + a2, b1 + b2)] = R.get((a1 + a2, b1 + b2), F(0)) + v1 * v2
        return R

    L = {}
    P = dict(X)
    for m in range(1, A_ + B + 1):
        for ab, v in P.items():
            L[ab] = L.get(ab, F(0)) + F((-1)**(m + 1), m) * v
        P = mul(P, X)
    return L


# ---------------------------------------------------------------------------
# A = ν² を形式変数とした記号計算（Q(q, s, ν)、s = t^{1/2}）。
# λ = 2^j 1^m 0^{z0}（c1 = j + m = N − z0 行）について
#   ψ'_{λ/1^k} = Π_{i=j}^{k−1} (1 − A t^{−z0−i})(1 − q t^{k−i−1}) / ((1 − q A t^{−z0−i−1})(1 − t^{k−i}))
#   P_λ(t^ρ) = s^{z0²} ν^{−z0} s^{j−j²} s^{j} ν^{−j} Π_{i<j}(t^i − A)/(1 − qA t^{−z0−i}) Π_{i<j}(t^i − Aq)/(1 − t^{j−i}) Π_{i=1}^{z0}(1 − A t^{i−z0−j})/(1 − t^i)
#   p_d = Σ_i (x_i^d + x_i^{−d})、x_i = ν s^{1−2i} q^{μ_i}（μ = λ − 1^N）：三つの等比級数の和（中央のブロックの長さは ν^{±2d} で閉じる）。
# いずれも有限 N の cover_Z と比較して検証する（tests/test_rk_cover.py）。
import sympy as sp
from sympy import QQ, mobius
from sympy.polys.fields import field

KC, qc, sc, nuc = field("q,s,nu", QQ)
qs_, ss_, nus_ = KC.symbols


def _psi_sym(k, j, z0):
    t = sc**2; A = nuc**2
    out = KC.one
    for i in range(j, k):
        out *= (1 - A * t**(-z0 - i)) * (1 - qc * t**(k - i - 1)) / ((1 - qc * A * t**(-z0 - i - 1)) * (1 - t**(k - i)))
    return out


def _principal_sym(j, z0):
    t = sc**2; A = nuc**2
    out = sc**(z0 * z0) * nuc**(-z0) * sc**(j - j * j) * sc**j * nuc**(-j)
    for i in range(j):
        out *= (t**i - A) / (1 - qc * A * t**(-z0 - i))
        out *= (t**i - A * qc) / (1 - t**(j - i))
    for i in range(1, z0 + 1):
        out *= (1 - A * t**(i - z0 - j)) / (1 - t**i)
    return out


def _p_d_sym(d, j, z0):
    """Σ_i (x_i^d + x_i^{−d})、x_i = ν s^{1−2i} q^{μ_i}：行 1..j は q、中央は q^0、最後の z0 行は q^{−1}。"""
    tot = KC.zero
    for sign in (1, -1):
        e = sign * d
        nu_e, q_e = nuc**e, qc**e
        # 行 1..j
        tot += sum((nu_e * q_e * sc**(e * (1 - 2 * i)) for i in range(1, j + 1)), KC.zero)
        # 中央 i = j+1 .. N−z0：Σ s^{e(1−2i)} = s^{e(1−2(j+1))} (1 − s^{−2e(N−z0−j)})/(1 − s^{−2e})、s^{−2eN} = ν^{−2e}
        ratio = sc**(-2 * e)
        first = sc**(e * (1 - 2 * (j + 1)))
        tot += nu_e * first * (1 - nuc**(-2 * e) * sc**(2 * e * (z0 + j))) / (1 - ratio)
        # 最後の z0 行 i = N−z0+1..N：s^{e(1−2i)} = s^{e(1−2N)} s^{2e(N−i)} = ν^{−2e} s^{e} s^{2e(N−i)}
        tot += sum((nu_e / q_e * nuc**(-2 * e) * sc**e * sc**(2 * e * m) for m in range(0, z0)), KC.zero)
    return tot


def _zprod_sym(j, z0, B):
    g = {d: -_p_d_sym(d, j, z0) / d for d in range(1, B + 1)}
    f = [KC.one]
    for m in range(1, B + 1):
        f.append(sum((d * g[d] * f[m - d] for d in range(1, m + 1)), KC.zero) / m)
    return f


def cover_Z_sym(A_, B):
    """Z^{cov}(y,z) mod (y^{A_+1}, z^{B+1})、係数は Q(q, s, ν)。N → ∞ の形（z0 = l − k + j ≥ 0、j ≤ k）。"""
    Z = {}
    for k in range(0, A_ + 1):
        for l in range(0, A_ + 1 - k):
            for j in range(0, k + 1):
                z0 = l - k + j
                if z0 < 0:
                    continue
                w = _psi_sym(k, j, z0) * _principal_sym(j, z0)
                pz = _zprod_sym(j, z0, B)
                for b in range(B + 1):
                    Z[(k + l, b)] = Z.get((k + l, b), KC.zero) + (-1)**(k + l) * w * pz[b]
    return Z


def log_coeffs_sym(Z, A_, B):
    assert Z[(0, 0)] == KC.one
    X = {ab: v for ab, v in Z.items() if ab != (0, 0)}

    def mul(P, Q):
        R = {}
        for (a1, b1), v1 in P.items():
            for (a2, b2), v2 in Q.items():
                if a1 + a2 <= A_ and b1 + b2 <= B:
                    R[(a1 + a2, b1 + b2)] = R.get((a1 + a2, b1 + b2), KC.zero) + v1 * v2
        return R

    L = {}
    P = dict(X)
    for m in range(1, A_ + B + 1):
        for ab, v in P.items():
            L[ab] = L.get(ab, KC.zero) + KC(sp.Rational((-1)**(m + 1), m)) * v
        P = mul(P, X)
    return L


def cover_f_sym(A_, B):
    """f^{cov}_{ab}（Möbius 反転済み、sympy の式、記号 q, s, nu）。"""
    L = log_coeffs_sym(cover_Z_sym(A_, B), A_, B)
    Le = {ab: v.as_expr() for ab, v in L.items()}
    out = {}
    for (a, b) in Le:
        g = sp.gcd(a, b) if a and b else (a or b)
        tot = 0
        for d in sp.divisors(int(g)):
            tot += sp.Rational(mobius(d), d) * Le[(a // d, b // d)].subs({qs_: qs_**d, ss_: ss_**d, nus_: nus_**d}, simultaneous=True)
        out[(a, b)] = sp.cancel(tot)
    return out


# ---------------------------------------------------------------------------
# K4'(a)：セクター別の被覆。probe を x 側（y1, z1）と 1/x 側（y2, z2）に分け、
#   Z^{cov}(y1,y2;z1,z2) = ⟨ Π_i(1 − y1 x_i)(1 − y2/x_i) · Π_i(1 − z1 x'_i)(1 − z2/x'_i) ⟩
# の係数を (k,l,m,n)（= 成分 1 の (|R1|,|R2|)、成分 2 の (|S1|,|S2|)）で持つ。
# 軌道体の作用 Ĩ は (k,l,m,n) → (l,k,n,m) なので、Ĩ で固定される 1 粒子状態は k = l, m = n のセクターにしかない。

def _p_pm_sym(d, j, z0):
    """(Σ_i x_i^d, Σ_i x_i^{−d})、x_i = ν s^{1−2i} q^{μ_i}（_p_d_sym の符号ごとの部分）。"""
    out = []
    for sign in (1, -1):
        e = sign * d
        nu_e, q_e = nuc**e, qc**e
        tot = sum((nu_e * q_e * sc**(e * (1 - 2 * i)) for i in range(1, j + 1)), KC.zero)
        ratio = sc**(-2 * e)
        first = sc**(e * (1 - 2 * (j + 1)))
        tot += nu_e * first * (1 - nuc**(-2 * e) * sc**(2 * e * (z0 + j))) / (1 - ratio)
        tot += sum((nu_e / q_e * nuc**(-2 * e) * sc**e * sc**(2 * e * m) for m in range(0, z0)), KC.zero)
        out.append(tot)
    return out


def _exp_series(g, B):
    """exp(Σ_d g_d z^d) mod z^{B+1}。"""
    f = [KC.one]
    for m in range(1, B + 1):
        f.append(sum((d * g[d] * f[m - d] for d in range(1, m + 1)), KC.zero) / m)
    return f


def cover_Z_sym4(A_, B):
    """Z^{cov}(y1,y2;z1,z2) mod (y^{A_+1}, z^{B+1})：{(k,l,m,n): Q(q,s,ν) の元}（k+l ≤ A_、m+n ≤ B）。"""
    Z = {}
    for k in range(0, A_ + 1):
        for l in range(0, A_ + 1 - k):
            for j in range(0, k + 1):
                z0 = l - k + j
                if z0 < 0:
                    continue
                w = _psi_sym(k, j, z0) * _principal_sym(j, z0)
                pp = {d: _p_pm_sym(d, j, z0) for d in range(1, B + 1)}
                f1 = _exp_series({d: -pp[d][0] / d for d in pp}, B)
                f2 = _exp_series({d: -pp[d][1] / d for d in pp}, B)
                for m in range(B + 1):
                    for n in range(B + 1 - m):
                        Z[(k, l, m, n)] = Z.get((k, l, m, n), KC.zero) + (-1)**(k + l) * w * f1[m] * f2[n]
    return Z


def log_coeffs_sym4(Z, A_, B):
    assert Z[(0, 0, 0, 0)] == KC.one
    X = {key: v for key, v in Z.items() if key != (0, 0, 0, 0)}

    def mul(P, Q):
        R = {}
        for k1, v1 in P.items():
            for k2, v2 in Q.items():
                key = tuple(a + b for a, b in zip(k1, k2))
                if key[0] + key[1] <= A_ and key[2] + key[3] <= B:
                    R[key] = R.get(key, KC.zero) + v1 * v2
        return R

    L = {}
    P = dict(X)
    for m in range(1, A_ + B + 1):
        for key, v in P.items():
            L[key] = L.get(key, KC.zero) + KC(sp.Rational((-1)**(m + 1), m)) * v
        P = mul(P, X)
    return L


def cover_f_sym4(A_, B):
    """f^{cov}_{(k,l),(m,n)}（全 d の Möbius 反転、sympy の式、記号 q, s, nu）。"""
    L = log_coeffs_sym4(cover_Z_sym4(A_, B), A_, B)
    Le = {key: v.as_expr() for key, v in L.items()}
    out = {}
    for key in Le:
        g = 0
        for x in key:
            g = sp.gcd(g, x)
        tot = 0
        for d in sp.divisors(int(g)):
            sub = tuple(x // d for x in key)
            tot += sp.Rational(mobius(d), d) * Le[sub].subs({qs_: qs_**d, ss_: ss_**d, nus_: nus_**d}, simultaneous=True)
        out[key] = sp.cancel(tot)
    return out


# ---------------------------------------------------------------------------
# K4'(a)：成分 1 に probe を無限個（brane 側の冪和 p_d(Y)、image 側の冪和 p_d(Y')）、成分 2 は 1 変数（z1 = x 側、z2 = 1/x 側）。
# Hopf の対称性で、composite の和を成分 2 の核 Π(1 − z1 x_i)(1 − z2/x_i) から取り、成分 1 は exp(−Σ_d [p_d(Y) p_d(x) + p_d(Y') p_d(1/x)]/d)。
# 鍵は (μ, μ', k, l)：p_μ(Y) p_μ'(Y') z1^k z2^l。

def _mul4(P, Q, A_, B):
    R = {}
    for (m1, m2, k1, l1), v1 in P.items():
        for (n1, n2, k2, l2), v2 in Q.items():
            if sum(m1) + sum(m2) + sum(n1) + sum(n2) <= A_ and k1 + k2 + l1 + l2 <= B:
                key = (tuple(sorted(m1 + n1, reverse=True)), tuple(sorted(m2 + n2, reverse=True)), k1 + k2, l1 + l2)
                R[key] = R.get(key, KC.zero) + v1 * v2
    return R


def _exp_two(c1, c2, A_):
    """exp(Σ_d c1[d] p_d(Y) + c2[d] p_d(Y'))、全次数 ≤ A_。"""
    out = {((), (), 0, 0): KC.one}
    for side, coef in ((0, c1), (1, c2)):
        for d, c in coef.items():
            term = {((), (), 0, 0): KC.one}
            acc = KC.one
            for m in range(1, A_ // d + 1):
                acc = acc * c / m
                mu = (d,) * m
                term[(mu, (), 0, 0) if side == 0 else ((), mu, 0, 0)] = acc
            out = _mul4(out, term, A_, 0)
    return out


def cover_Z_multi(A_, B):
    """Z^{cov}(Y, Y'; z1, z2)：{(μ, μ', k, l): Q(q,s,ν)}、|μ| + |μ'| ≤ A_、k + l ≤ B。"""
    Z = {}
    for k in range(0, B + 1):
        for l in range(0, B + 1 - k):
            for j in range(0, k + 1):
                z0 = l - k + j
                if z0 < 0:
                    continue
                w = _psi_sym(k, j, z0) * _principal_sym(j, z0)
                pp = {d: _p_pm_sym(d, j, z0) for d in range(1, A_ + 1)}
                ey = _exp_two({d: -pp[d][0] / d for d in pp}, {d: -pp[d][1] / d for d in pp}, A_)
                for (m1, m2, _, _), v in ey.items():
                    key = (m1, m2, k, l)
                    Z[key] = Z.get(key, KC.zero) + (-1)**(k + l) * w * v
    return Z


def log_coeffs_multi(Z, A_, B):
    one = ((), (), 0, 0)
    assert Z[one] == KC.one
    X = {key: v for key, v in Z.items() if key != one}
    L, P = {}, dict(X)
    for m in range(1, A_ + B + 1):
        for key, v in P.items():
            L[key] = L.get(key, KC.zero) + KC(sp.Rational((-1)**(m + 1), m)) * v
        P = _mul4(P, X, A_, B)
    return L


def cover_f_multi(A_, B):
    """log Z^{cov} = Σ_d (1/d) Σ f_{ρ,ρ',k,l}(x^d) s_ρ(Y^d) s_ρ'(Y'^d) z1^{kd} z2^{ld} を満たす f（sympy の式、記号 q, s, nu）。"""
    from .characters import chi
    from .partitions import partitions, z as z_lam
    L = log_coeffs_multi(cover_Z_multi(A_, B), A_, B)
    Le = {key: v.as_expr() for key, v in L.items()}

    def ad(e, d):
        return e.subs({qs_: qs_**d, ss_: ss_**d, nus_: nus_**d}, simultaneous=True)

    f = {}
    for n1 in range(A_ + 1):
        for n2 in range(A_ + 1 - n1):
            for k in range(B + 1):
                for l in range(B + 1 - k):
                    if n1 + n2 + k + l == 0:
                        continue
                    H = {(m1, m2): Le.get((m1, m2, k, l), 0) for m1 in partitions(n1) for m2 in partitions(n2)}
                    g = 0
                    for x in (n1, n2, k, l):
                        g = sp.gcd(g, x)
                    for d in sp.divisors(int(g)):
                        if d == 1:
                            continue
                        for r1 in partitions(n1 // d):
                            for r2 in partitions(n2 // d):
                                fv = f.get((r1, r2, k // d, l // d))
                                if not fv:
                                    continue
                                fa = ad(fv, d) / d
                                for m1 in partitions(n1 // d):
                                    for m2 in partitions(n2 // d):
                                        c = sp.Rational(chi(r1, m1), z_lam(m1)) * sp.Rational(chi(r2, m2), z_lam(m2))
                                        key = (tuple(d * x for x in m1), tuple(d * x for x in m2))
                                        H[key] = H[key] - fa * c
                    for r1 in partitions(n1):
                        for r2 in partitions(n2):
                            f[(r1, r2, k, l)] = sp.cancel(sum(chi(r1, m1) * chi(r2, m2) * H[(m1, m2)] for (m1, m2) in H))
    return f
