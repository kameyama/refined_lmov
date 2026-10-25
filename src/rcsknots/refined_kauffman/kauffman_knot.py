r"""(b)：N を固定した refined SO(2N) の結び目の自由エネルギー（probe 1 個、双対 Cauchy 型、virtual 正規化）。

  Z~(y) = Σ_{k=0}^{N} (−1)^k W_{1^k} · [y^N K_{(N−k)}(y; t, q; 1, −1, q^{1/2}, −q^{1/2})]
W_{1^k} は結び目の unnormalized な不変量（unknot では K_{1^k}(t^ρ)）。log Z~ = Σ_a L_a y^a、f_a = Σ_{d|a} (μ(d)/d) L_{a/d}(q^d, t^d)。
変数は r = q^{1/2}、s = t^{1/2}（sympy の記号）。"""
import sympy as sp
from sympy import mobius

from .askey_wilson import askey_wilson_monic
from .kauffman_refined import koornwinder_eval, q as _q, t as _t, U as _U

r, s = sp.symbols("r s", positive=True)


def unknot_W(N, k):
    """K_{1^k}(t^ρ; q, t; 1, −1, q^{1/2}, −q^{1/2})（D_N、U = t^N）。"""
    e = koornwinder_eval((1,) * k if k else (), (1, -1, sp.sqrt(_q), -sp.sqrt(_q)))
    return sp.cancel(e.subs({_q: r**2, _t: s**2, _U: s**(2 * N)}))


def Ztilde(N, W):
    """W = [W_0, …, W_N] から Z~(y) の係数 {a: 式}（a = 0..2N）。"""
    aw = askey_wilson_monic(N, s**2, sp.Integer(1), sp.Integer(-1), r, -r)
    Z = {}
    for k in range(N + 1):
        n = N - k
        for j in range(n + 1):
            c = aw[n].get(n - 2 * j, 0)
            if c:
                Z[k + 2 * j] = Z.get(k + 2 * j, 0) + (-1)**k * W[k] * c
    return {a: sp.cancel(v) for a, v in Z.items()}


def _mul(P, Q, amax):
    R = {}
    for a1, v1 in P.items():
        for a2, v2 in Q.items():
            if a1 + a2 <= amax:
                R[a1 + a2] = R.get(a1 + a2, 0) + v1 * v2
    return {a: sp.cancel(v) for a, v in R.items()}


def free_energy(N, W, amax):
    """f_a（a = 1..amax）。"""
    Z = Ztilde(N, W)
    assert sp.cancel(Z[0] - 1) == 0, Z[0]
    X = {a: v for a, v in Z.items() if 0 < a <= amax}
    L, P = {}, dict(X)
    for m in range(1, amax + 1):
        for a, v in P.items():
            L[a] = L.get(a, 0) + sp.Rational((-1)**(m + 1), m) * v
        P = _mul(P, X, amax)
    f = {}
    for a in range(1, amax + 1):
        tot = 0
        for d in sp.divisors(a):
            tot += sp.Rational(mobius(d), d) * L.get(a // d, 0).subs({r: r**d, s: s**d}, simultaneous=True)
        f[a] = sp.cancel(tot)
    return f
