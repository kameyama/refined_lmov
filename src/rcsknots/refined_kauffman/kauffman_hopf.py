r"""K2'(a)：双対 Cauchy（t-brane）型の refined orientifold OV 核での Hopf 絡み目（probe 各 1 個）を、
U = t^N を形式変数として Q(q,t,U) の上で計算する。

virtual 正規化 Z~(y,z) := y^N z^N Z(y,z)
  = Σ_{k=0}^{N} (−1)^k K_{1^k}(t^ρ; q,t; t̄) · [y^N K_{(N−k)}(y; t,q; t̄)] · Π_{i=1}^{N} (1 − z x_i)(1 − z/x_i)|_{x_i = t^{N−i} q^{[i≤k]}}
は y, z の多項式で Z~(0,0) = 1。各因子の U 依存性：
  - K_{1^k}(t^ρ)：van Diejen–Sahi の評価式（kauffman_refined.koornwinder_eval）。
  - y^N K_{(N−k)}(y) = y^k Σ_j c_{2j}(N−k) y^{2j}（対称性 c_{2j}(n) = c_{2(n−j)}(n)）、c_{2j} は askey_wilson.top_coeff_closed。
  - z 側の積 = exp(−Σ_d p_d z^d/d)、p_d = Σ_i (x_i^d + x_i^{−d}) は等比級数の和で U^{±d} の有理式。
y^a の次数 a に寄与するのは k ≤ a だけなので、和は N によらず有限。
log Z~ = Σ L_{ab} y^a z^b、f_{ab} = Σ_{d|(a,b)} (μ(d)/d) L_{a/d,b/d}(q^d,t^d,U^d)（Möbius 反転）。
"""
import sympy as sp
from sympy import QQ, mobius
from sympy.polys.fields import field

from .askey_wilson import top_coeff_closed
from .kauffman_refined import koornwinder_eval
from .kauffman_refined import q as _qs, t as _ts, U as _Us

K, q, t, U = field("q,t,U", QQ)
qs, ts, Us = K.symbols
TBAR = (1, -1, sp.sqrt(_qs), -sp.sqrt(_qs))


def _W(k):
    e = sp.cancel(koornwinder_eval((1,) * k if k else (), TBAR)).subs({_qs: qs, _ts: ts, _Us: Us})
    return K.from_expr(e)


def _p_z(d, k):
    td, Ud, qd = t**d, U**d, q**d
    s1 = qd * Ud * sum((td**(-i) for i in range(1, k + 1)), K.zero)
    s2 = (K.one - Ud * td**(-k)) / (K.one - td)
    s3 = Ud**-1 / qd * sum((td**i for i in range(1, k + 1)), K.zero)
    s4 = (K.one - Ud**-1 * td**k) / (K.one - td**-1)
    return s1 + s2 + s3 + s4


def _zprod(k, B):
    """exp(Σ_d g_d z^d) mod z^{B+1}、g_d = −p_d/d。f_m = (1/m) Σ_{d≤m} d g_d f_{m−d}。"""
    g = {d: -_p_z(d, k) / d for d in range(1, B + 1)}
    f = [K.one]
    for m in range(1, B + 1):
        f.append(sum((d * g[d] * f[m - d] for d in range(1, m + 1)), K.zero) / m)
    return f


def hopf_Ztilde(A, B):
    """Z~(y,z) mod (y^{A+1}, z^{B+1})：{(a,b): Q(q,t,U) の元}。"""
    Z = {}
    for k in range(0, A + 1):
        Wk = _W(k)
        Un = U * t**(-k)
        ky = {k + 2 * j: top_coeff_closed(j, Un, q, t) for j in range((A - k) // 2 + 1)}
        pz = _zprod(k, B)
        for a, va in ky.items():
            for b in range(B + 1):
                Z[(a, b)] = Z.get((a, b), K.zero) + (-1)**k * Wk * va * pz[b]
    return Z


def log_coeffs(Z, A, B):
    """log Z~ = Σ L_{ab} y^a z^b（Z~(0,0) = 1 を仮定）。"""
    assert Z[(0, 0)] == K.one
    X = {ab: v for ab, v in Z.items() if ab != (0, 0)}

    def mul(P, Q):
        R = {}
        for (a1, b1), v1 in P.items():
            for (a2, b2), v2 in Q.items():
                if a1 + a2 <= A and b1 + b2 <= B:
                    R[(a1 + a2, b1 + b2)] = R.get((a1 + a2, b1 + b2), K.zero) + v1 * v2
        return R

    L = {}
    P = dict(X)
    for m in range(1, A + B + 1):
        for ab, v in P.items():
            L[ab] = L.get(ab, K.zero) + K(sp.Rational((-1)**(m + 1), m)) * v
        P = mul(P, X)
    return L


def hopf_log_coeffs(A, B):
    """L_{ab}（a ≤ A, b ≤ B）を sympy の式（記号 q, t, U）で返す。"""
    L = log_coeffs(hopf_Ztilde(A, B), A, B)
    return {ab: v.as_expr() for ab, v in L.items()}


def mobius_f(L, a, b):
    """f_{ab} = Σ_{d|gcd(a,b)} (μ(d)/d) L_{a/d,b/d}(q^d, t^d, U^d)。"""
    g = sp.gcd(a, b)
    tot = 0
    for d in sp.divisors(int(g)):
        tot += sp.Rational(mobius(d), d) * L[(a // d, b // d)].subs({qs: qs**d, ts: ts**d, Us: Us**d}, simultaneous=True)
    return sp.cancel(tot)


# ---------------------------------------------------------------------------
# K2'(b)：成分 1 に probe を無限個（冪和 P_d = p_d(y) を形式変数）、成分 2 に probe 1 個（z）。
# S 行列の対称性 S_{λ'μ} = S_{μλ'} を使って成分 1 の和を双対 Cauchy 恒等式で縮約すると
#   Z~(P; z) = Σ_{k=0}^{B} (−1)^k K_{1^k}(t^ρ) \hat K_{(k)}(z;U) · exp(−Σ_d p_d^{(k)} P_d / d)
# となり、多変数の Koornwinder 多項式は不要になる（p_d^{(k)} は _p_z と同じ冪和）。
# 自由エネルギーは log Z~ = Σ_d (1/d) Σ_{R,b} f_{R,b}(q^d,t^d,U^d) B_R(y^d) z^{bd}（B_R は Schur か Macdonald P_R(y;t,q)）。

from fractions import Fraction

from ..characters import chi
from ..macdonald import macdonald_P
from ..partitions import partitions, z as z_lam


def _sym_mul(P, Q, A):
    R = {}
    for (mu1, b1), v1 in P.items():
        for (mu2, b2), v2 in Q.items():
            if sum(mu1) + sum(mu2) <= A:
                key = (tuple(sorted(mu1 + mu2, reverse=True)), b1 + b2)
                R[key] = R.get(key, K.zero) + v1 * v2
    return R


def _exp_powersums(coef, A):
    """exp(Σ_d coef[d] P_d) を |μ| ≤ A で打ち切った {(μ, 0): 係数}。"""
    out = {((), 0): K.one}
    for d, c in coef.items():
        term = {((), 0): K.one}
        acc = K.one
        for m in range(1, A // d + 1):
            acc = acc * c / m
            term[((d,) * m, 0)] = acc
        out = _sym_mul(out, term, A)
    return out


def hopf_Ztilde_multi(A, B):
    """Z~(P; z) mod (|μ| ≤ A, z^{B+1})：{(μ, b): Q(q,t,U) の元}。"""
    Z = {}
    for k in range(0, B + 1):
        Wk = _W(k)
        Un = U * t**(-k)
        kz = {k + 2 * j: top_coeff_closed(j, Un, q, t) for j in range((B - k) // 2 + 1)}
        ey = _exp_powersums({d: -_p_z(d, k) / d for d in range(1, A + 1)}, A)
        for (mu, _), vy in ey.items():
            for b, vb in kz.items():
                Z[(mu, b)] = Z.get((mu, b), K.zero) + (-1)**k * Wk * vb * vy
    return Z


def log_coeffs_multi(Z, A, B):
    assert Z[((), 0)] == K.one
    X = {key: v for key, v in Z.items() if key != ((), 0)}
    L = {}
    P = dict(X)
    for m in range(1, A + B + 1):
        for key, v in P.items():
            L[key] = L.get(key, K.zero) + K(sp.Rational((-1)**(m + 1), m)) * v
        P = {key: v for key, v in _sym_mul(P, X, A).items() if key[1] <= B}
    return L


def _basis_in_p(n, basis, d):
    """{R: {ρ: 係数}}：次数 n の基底 B_R を冪和 p_ρ で表す。d 乗のプレシズム後のパラメータ (t^d, q^d)。"""
    if basis == "schur":
        return {R: {rho: K(Fraction(chi(R, rho), z_lam(rho))) for rho in partitions(n)} for R in partitions(n)}
    if basis == "macdonald":  # probe 側の極限基底 P_R(y; t, q)：Macdonald の (q,t) に (t^d, q^d) を入れる
        return macdonald_P(n, t**d, q**d, lambda c: K(c))
    raise ValueError(basis)


def _adams(v, d):
    return K.from_expr(v.as_expr().subs({qs: qs**d, ts: ts**d, Us: Us**d}, simultaneous=True))


def _solve_field(M, rhs):
    """M x = rhs（K の元の正方行列）。"""
    n = len(M)
    Aug = [row[:] + [r] for row, r in zip(M, rhs)]
    for c in range(n):
        p = next(i for i in range(c, n) if Aug[i][c] != 0)
        Aug[c], Aug[p] = Aug[p], Aug[c]
        inv = K.one / Aug[c][c]
        Aug[c] = [v * inv for v in Aug[c]]
        for i in range(n):
            if i != c and Aug[i][c] != 0:
                f = Aug[i][c]
                Aug[i] = [vi - f * vc for vi, vc in zip(Aug[i], Aug[c])]
    return [Aug[i][n] for i in range(n)]


def plethystic_f(L, A, B, basis="schur"):
    """log Z~ = Σ_d (1/d) Σ_{R,b} f_{R,b}(q^d,t^d,U^d) B_R(y^d) z^{bd} を満たす f_{R,b}（|R| ≤ A, b ≤ B）。"""
    f = {}
    for n in range(0, A + 1):
        for b in range(0, B + 1):
            if n == 0 and b == 0:
                continue
            H = {rho: L.get((rho, b), K.zero) for rho in partitions(n)}
            g = sp.gcd(n, b) if n and b else (n or b)
            for d in sp.divisors(int(g)):
                if d == 1:
                    continue
                Bd = _basis_in_p(n // d, basis, d)
                for R, vec in Bd.items():
                    fv = f.get((R, b // d))
                    if fv is None or fv == 0:
                        continue
                    fa = _adams(fv, d) / d
                    for rho, c in vec.items():
                        drho = tuple(d * x for x in rho)
                        H[drho] = H[drho] - fa * c
            B1 = _basis_in_p(n, basis, 1)
            parts = partitions(n)
            M = [[B1[R].get(rho, K.zero) for R in parts] for rho in parts]
            sol = _solve_field(M, [H[rho] for rho in parts])
            for R, v in zip(parts, sol):
                f[(R, b)] = v
    return f
