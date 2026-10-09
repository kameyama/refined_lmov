"""ダブルツイスト結び目 K_{m,n} の refined 不変量の「チャネル分解」（別研究 [2,1] の X2）。

KN の変数 (A, Q, T) = (a, q, t)^{1/2} で書く。X1 で確かめた読み替え（conventions.md §16）：
  rCS_{[r^s]}(K_{m,n}) = N(m,n) Σ_{λ ⊆ [r^s]} Z_λ τ_λ(m) τ_λ(n) / τ_λ(1)，  N = (T/Q)^{rs (sgn m + sgn n)/2}
  Z_λ = 𝐃^{(s)}_λ(Q, T) 𝐃^{(r)}_{λ^T}(1/T, 1/Q) Π_{□∈λ} {A Q^{r+a'} T^{−l'}}{A Q^{a'} T^{−s−l'}}（4_1 の項）
  τ_λ(m) は KNTZ の 𝐓𝐰^{(m)}_λ（(t/q)^{rs} を除く、𝐅 の中で A → A(T/Q)^{1/2}）。τ_λ(−1) = 1, τ_λ(0) = 0。
τ_λ(m) = Σ_{μ⊆λ} U_{λμ} e_μ^m は全ての m で一つの指数和（refined の 𝓑 の行）で、固有値は e_μ = (A²Q/T)^{|μ|} T_μ²（T_μ = Q^{2n(μ')} T^{−2n(μ)}）。
"""
from functools import lru_cache

import sympy as sp

from .double_twist_refined import A, Bref, Dref, Fref, br, q as Q, subdiagrams, t as T
from .partitions import arm, boxes, conjugate, leg, n_stat

S = sp.sqrt(T / Q)


def eig(mu):
    """e_μ = (A²Q/T)^{|μ|} Q^{4n(μ')} T^{−4n(μ)}。"""
    return (A**2 * Q / T) ** sum(mu) * Q ** (4 * n_stat(conjugate(mu))) * T ** (-4 * n_stat(mu))


def tre0(lam):
    e = (-A**2 * Q / T) ** sum(lam)
    for i, j in boxes(lam):
        e *= Q ** (2 * arm(lam, i, j)) * T ** (-2 * leg(lam, i, j))
    return e


def tau(lam, m):
    """τ_λ(m)（KNTZ の RefinedTwFactor から (t/q)^{rs} を除き、KN の変数に読み替えたもの）。"""
    if m == -1:
        return sp.Integer(1)
    if m == 0:
        return sp.Integer(0) if lam else sp.Integer(1)
    if m == 1:
        return tre0(lam)
    if m < 0:
        return Fref(lam, -m).subs({A: 1 / (A * S), Q: 1 / Q, T: 1 / T}, simultaneous=True)
    return Fref(lam, m).subs(A, A * S) * tre0(lam)


@lru_cache(maxsize=None)
def U_row(lam):
    """τ_λ(m) = Σ_μ U[μ] e_μ^m の係数 {μ: U}（m = −1, 0, 1, …, k − 2 で解く。指数和の形が全ての m で成り立つことはテストで確認）。"""
    if not lam:
        return {(): sp.Integer(1)}
    mus = subdiagrams(lam)
    cs = sp.symbols(f"c0:{len(mus)}")
    eqs = [sum(c * eig(mu) ** m for c, mu in zip(cs, mus)) - tau(lam, m) for m in range(-1, len(mus) - 1)]
    sol = sp.solve(eqs, cs, dict=True)[0]
    return {mu: sp.factor(sol[c]) for c, mu in zip(cs, mus)}


def Z_rect(lam, r, s):
    """長方形 [r^s] の Z_λ（KNTZ の 4_1 の項、KN の変数）。"""
    z = Dref(s, lam, Q, T) * Dref(r, conjugate(lam), 1 / T, 1 / Q)
    for i, j in boxes(lam):
        z *= br(A * Q ** (r + j) * T ** (-i)) * br(A * Q ** j * T ** (-s - i))
    return sp.factor(z)


def M_rect(r, s):
    """長方形 [r^s] の係数 M[(μ, ν)]：rCS = N Σ M[μ,ν] e_μ^m e_ν^n。"""
    M = {}
    for lam in subdiagrams((r,) * s):
        U = U_row(lam)
        for mu, u in U.items():
            for nu, v in U.items():
                M[(mu, nu)] = M.get((mu, nu), 0) + Z_rect(lam, r, s) * u * v / tre0(lam)
    return M


def value(M, m, n, rs):
    """N(m, n) Σ M e^m e^n（rs は N の指数の長方形の箱の数）。"""
    sg = (sp.sign(m) + sp.sign(n)) / 2
    return (T / Q) ** (rs * sg) * sum(c * eig(mu) ** m * eig(nu) ** n for (mu, nu), c in M.items())
