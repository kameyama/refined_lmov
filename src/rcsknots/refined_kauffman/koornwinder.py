r"""Koornwinder 多項式 K_λ(x; q, t; t_0, t_1, t_2, t_3)（BC_n 対称 Laurent 多項式、小さい n 用）。

Koornwinder (1992) の q-差分作用素
  D = Σ_i [Φ_i(x) (T_{q,i} − 1) + Φ_i(x^{-1}) (T_{q,i}^{-1} − 1)],
  Φ_i(x) = Π_r (1 − t_r x_i) / ((1 − x_i²)(1 − q x_i²)) × Π_{j≠i} (1 − t x_i x_j)(1 − t x_i/x_j) / ((1 − x_i x_j)(1 − x_i/x_j))
の固有関数で、固有値は
  E_λ = Σ_i [ q^{-1} t_0 t_1 t_2 t_3 t^{2n−i−1} (q^{λ_i} − 1) + t^{i−1} (q^{−λ_i} − 1) ]
（Rains–Warnaar 1506.02755 §2.4、Koornwinder の原論文の規約）。K_λ = m^W_λ + Σ_{μ<λ} c_{λμ} m^W_μ を線形方程式で解く。
パラメータ q, t, t_r は数値（有理数）でも記号でもよいが、記号だと遅い。
"""
import itertools
from functools import lru_cache

import sympy as sp

from ..partitions import partitions, dominates


def _symbols(n):
    return sp.symbols(f"x1:{n+1}")


def monomial_W(lam, xs):
    """BC_n の軌道の単項対称関数 m^W_λ = Σ_{α ∈ W·λ} x^α。"""
    n = len(xs)
    lam = tuple(lam) + (0,) * (n - len(lam))
    terms = set()
    for perm in itertools.permutations(lam):
        for signs in itertools.product((1, -1), repeat=n):
            terms.add(tuple(s * p for s, p in zip(signs, perm)))
    return sum(sp.prod([x**a for x, a in zip(xs, alpha)]) for alpha in terms)


def _phi(i, xs, q, t, tbar):
    xi = xs[i]
    num = sp.prod([1 - tr * xi for tr in tbar])
    den = (1 - xi**2) * (1 - q * xi**2)
    for j, xj in enumerate(xs):
        if j != i:
            num *= (1 - t * xi * xj) * (1 - t * xi / xj)
            den *= (1 - xi * xj) * (1 - xi / xj)
    return num / den


def apply_D(f, xs, q, t, tbar):
    out = 0
    for i, xi in enumerate(xs):
        phi = _phi(i, xs, q, t, tbar)
        phi_inv = phi.subs(xi, 1 / xi)
        out += phi * (f.subs(xi, q * xi) - f) + phi_inv * (f.subs(xi, xi / q) - f)
    return sp.cancel(sp.together(out))


def _coeffs_in_mW(expr, xs, parts):
    """BC_n 対称な Laurent 多項式を m^W_μ（μ ∈ parts）で展開した係数。"""
    expr = sp.expand(expr)
    P = sp.Poly(expr * sp.prod(xs) ** 40, *xs)  # 負冪を避けるためのずらし
    d = P.as_dict()
    out = {}
    for mu in parts:
        n = len(xs)
        alpha = tuple(mu) + (0,) * (n - len(mu))
        key = tuple(a + 40 for a in alpha)
        out[mu] = d.get(key, 0)
    return out


def eigenvalue(lam, n, q, t, tbar):
    t0, t1, t2, t3 = tbar
    lam = tuple(lam) + (0,) * (n - len(lam))
    return sum(t0 * t1 * t2 * t3 / q * t ** (2 * n - i - 1) * (q ** lam[i - 1] - 1) + t ** (i - 1) * (q ** (-lam[i - 1]) - 1)
               for i in range(1, n + 1))


def koornwinder(lam, n, q, t, tbar, xs=None):
    """K_λ(x_1..x_n; q, t; t̄)。λ の長さ ≤ n。"""
    xs = xs or _symbols(n)
    lam = tuple(x for x in lam if x > 0)
    assert len(lam) <= n
    lower = [mu for mu in _parts_below(lam, n) if mu != lam]
    if not lower:
        return sp.expand(monomial_W(lam, xs))
    allp = [lam] + lower
    Dm = {mu: _coeffs_in_mW(apply_D(monomial_W(mu, xs), xs, q, t, tbar), xs, allp) for mu in allp}
    E = eigenvalue(lam, n, q, t, tbar)
    cs = sp.symbols(f"c0:{len(lower)}")
    eqs = []
    for nu in lower:
        eq = Dm[lam][nu] + sum(c * (Dm[mu][nu] - (E if mu == nu else 0)) for c, mu in zip(cs, lower))
        eqs.append(eq)
    sol = sp.solve(eqs, cs, dict=True)
    assert len(sol) == 1, sol
    K = monomial_W(lam, xs) + sum(sol[0][c] * monomial_W(mu, xs) for c, mu in zip(cs, lower))
    return sp.expand(K)


def _parts_below(lam, n):
    """|μ| ≤ |λ|、長さ ≤ n、μ ≤ λ（BC_n の dominance：μ_1+…+μ_k ≤ λ_1+…+λ_k、|μ| ≡ |λ| mod 2 を含めて全部）。"""
    out = []
    lam_full = tuple(lam) + (0,) * (n - len(lam))
    for m in range(sum(lam), -1, -1):
        for mu in partitions(m):
            if len(mu) > n:
                continue
            mu_full = tuple(mu) + (0,) * (n - len(mu))
            if all(sum(mu_full[:k]) <= sum(lam_full[:k]) for k in range(1, n + 1)):
                out.append(tuple(mu))
    return out
