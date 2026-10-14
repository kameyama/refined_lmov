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
