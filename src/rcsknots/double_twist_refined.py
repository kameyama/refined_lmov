"""KNTZ 1902.02275 の refined 版：ダブルツイスト結び目 K_{m,n} の長方形 [r^s] の Poincaré 多項式（別研究 [2,1] の X1）。

基準は本文（Cyclotomic-Rectangular.tex :706-761、付録 A :779-）ではなく、著者の補助ファイル
references/1902.02275/src/anc/HOMFLYPTinvariants.m（以下「パッケージ」。行番号は CR を改行に直した HOMFLYPTinvariants.lf.m の行：
補間 Macdonald :443-460、BinomialFormulaqt :487、reprefactor :511、RefinedFFactor :565、reFtrefoil :569、RefinedTwFactor :573、Poincaré の DoubleBraidPolynomial :665）。本文との違い：
- 置換 (cov) q → −t_c, t → 1/t_r, A → 𝐚√(−t_r t_c), σ → t_r^{−s}𝐐^{−1} は最後にまとめて掛ける。Tw の中の 𝐅 だけは先に A → 𝐚（m < 0 では A → 1/𝐚, q → 1/q, t → 1/t）にする。
  本文の Tw の 𝐅(𝐚, 𝐭_c, 𝐭_r) とは q, t の符号と逆数が違う。
- 三葉結び目の因子 reFtrefoil に (t/q)^{rs} が付く（本文の前因子 (−𝐭_r𝐭_c)^{−rs(sgn m + sgn n)/2} の代わり）。
- 付録 A の P*_λ の定義の積 Π_{□∈μ} は Π_{□∈λ} の誤植（パッケージの xTableShift）。
変数：計算の途中は (A, q, t, σ)、結果は (𝐚, 𝐐, 𝐭_r, 𝐭_c)（記号 a, Q, t_r, t_c）。照合は scripts/x1_check.py とテスト。
"""
from functools import lru_cache
from itertools import product

import sympy as sp

from .partitions import arm, boxes, conjugate, leg

A, q, t, sig = sp.symbols("A q t sigma")
a, Qb, tr, tc = sp.symbols("a Q t_r t_c")


def br(x):
    return x - 1 / x


def _contains(lam, mu):
    return len(mu) <= len(lam) and all(m <= l for m, l in zip(mu, lam))


def subdiagrams(lam):
    """λ に含まれるヤング図形（∅ と λ を含む）。"""
    out = []
    for vec in product(*[range(x + 1) for x in lam]):
        if all(vec[i] >= vec[i + 1] for i in range(len(vec) - 1)):
            out.append(tuple(x for x in vec if x))
    return sorted(set(out), key=lambda m: (sum(m), m))


# ---------------------------------------------------------------- 補間 Macdonald 多項式（付録 A、パッケージ :443-460）
def _b(lam, i, j, qq, tt):
    return (1 - qq ** arm(lam, i, j) * tt ** (leg(lam, i, j) + 1)) / (1 - qq ** (arm(lam, i, j) + 1) * tt ** leg(lam, i, j))


def _psi(lam, mu, qq, tt):
    """ψ_{λ/μ}：λ/μ の行に入り列に入らない箱で b_μ/b_λ の積（λ/μ は水平帯）。"""
    rows = {i for i, j in boxes(lam) if not (i < len(mu) and j < mu[i])}
    cols = {j for i, j in boxes(lam) if not (i < len(mu) and j < mu[i])}
    out = sp.Integer(1)
    for i, j in boxes(lam):
        if i in rows and j not in cols:
            out *= _b(mu, i, j, qq, tt) / _b(lam, i, j, qq, tt)
    return out


def _ssyt(shape, n):
    """形 shape、成分 1..n の半標準盤（行は広義増加、列は狭義増加）。{(i, j): 値}。"""
    cells = boxes(shape)
    out = []

    def rec(k, T):
        if k == len(cells):
            out.append(dict(T))
            return
        i, j = cells[k]
        lo = 1
        if j > 0:
            lo = max(lo, T[(i, j - 1)])
        if i > 0:
            lo = max(lo, T[(i - 1, j)] + 1)
        for v in range(lo, n + 1):
            T[(i, j)] = v
            rec(k + 1, T)
            del T[(i, j)]
    rec(0, {})
    return out


def interp_P(mu, xs, qq=q, tt=t):
    """P*_μ(x_1, …, x_n; q, t)（n = len(xs)）。半標準盤 T について ψ_T と、逆盤 R = n + 1 − T の Π t^{1−R}(x_R − q^{a'} t^{−l'})。"""
    n = len(xs)
    if not mu:
        return sp.Integer(1)
    total = 0
    for T in _ssyt(mu, n):
        shapes = []
        for k in range(n, -1, -1):
            shapes.append(tuple(x for x in (sum(1 for j in range(r) if T[(i, j)] <= k) for i, r in enumerate(mu)) if x))
        w = sp.Integer(1)
        for big, small in zip(shapes, shapes[1:]):
            w *= _psi(big, small, qq, tt)
        for (i, j), v in T.items():
            R = n + 1 - v
            w *= tt ** (1 - R) * (xs[R - 1] - qq ** j * tt ** (-i))
        total += w
    return total


# ---------------------------------------------------------------- 𝐁, 𝐃, 𝐅（パッケージ :487-565）
def Dref(s, lam, qq=q, tt=t):
    """𝐃_λ^{(s)}(q, t) = Π {t^{s − l'} q^{a'}}/{t^{l+1} q^{a}}。"""
    out = sp.Integer(1)
    for i, j in boxes(lam):
        out *= br(tt ** (s - i) * qq ** j) / br(tt ** (leg(lam, i, j) + 1) * qq ** arm(lam, i, j))
    return out


def binom_qt(lam, mu, s):
    """パッケージの BinomialFormulaqt[λ, μ, s]。"""
    num = interp_P(mu, [q ** x for x in lam])
    den = interp_P(mu, [q ** x for x in mu])
    first = (num / den).subs({q: q**2, t: t**2}, simultaneous=True)
    second = (interp_P(lam, [0] * s) / interp_P(mu, [0] * s)).subs({q: q**-2, t: t**-2}, simultaneous=True)
    return first * second


@lru_cache(maxsize=None)
def Bref(lam, mu, s):
    """𝐁^{(s)}_{λ,μ}(A, q, t)（本文 :726-733、パッケージの reblock）。"""
    if not _contains(lam, mu):
        return sp.Integer(0)
    pre = (-t ** (1 - s)) ** (sum(lam) - sum(mu)) * (A * q / t) ** (2 * sum(mu))
    for i, j in boxes(lam):
        pre *= q ** (2 * j) * t ** (-2 * i)
    for i, j in boxes(mu):
        pre *= q ** (2 * j) * t ** (-2 * i)
    return sp.cancel(pre * binom_qt(lam, mu, s))


@lru_cache(maxsize=None)
def Fref(lam, m, s=None):
    """𝐅^{(m)}_λ(A, q, t)（m ≥ 1）= (𝐁^{m−1}𝐃)_λ / 𝐃_λ（入れ子の和、パッケージの RefinedFFactor）。"""
    s = s or max(len(lam), 1)
    subs = subdiagrams(lam)
    v = {mu: Dref(s, mu) for mu in subs}
    for _ in range(m - 1):
        v = {nu: sp.cancel(sum(Bref(nu, mu, s) * v[mu] for mu in subs if _contains(nu, mu))) for nu in subs}
    return sp.factor(sp.cancel(v[lam] / Dref(s, lam)))


def Tw_ref(lam, m, r, s):
    """パッケージの RefinedTwFactor[λ, m, r, s]（:569-573）。"""
    tre = (t / q) ** (r * s) * (-A**2 * q / t) ** sum(lam)
    for i, j in boxes(lam):
        tre *= q ** (2 * arm(lam, i, j)) * t ** (-2 * leg(lam, i, j))
    if m == -1:
        return sp.Integer(1)
    if m == 1:
        return tre
    if m < 0:
        return Fref(lam, -m, s).subs({A: 1 / a, t: 1 / t, q: 1 / q}, simultaneous=True)
    return Fref(lam, m, s).subs(A, a) * tre


def summand_pre(r, s, lam):
    """和の各項の Tw 以外の部分（(A, q, t, σ)、置換の前）。"""
    term = Dref(s, lam) * Dref(r, conjugate(lam)).subs({q: 1 / t, t: 1 / q}, simultaneous=True)
    for i, j in boxes(lam):
        term *= br(A * q ** (r + j) * t ** (-i) / sig) * br(A * q ** j * t ** (-s - i) * sig)
    return term


def poincare_pre(r, s, m, n):
    """パッケージの DoubleBraidPolynomial[r, s, m, n] の置換の前の式（(A, q, t, σ, a)）。"""
    total = 0
    for lam in subdiagrams((r,) * s):
        tw = Tw_ref(lam, m, r, s) / Tw_ref(lam, -1, r, s) * Tw_ref(lam, n, r, s) / Tw_ref(lam, 1, r, s)
        total += summand_pre(r, s, lam) * tw
    return total


def cov(expr, s):
    """(cov)：q → −t_c, t → 1/t_r, A → 𝐚√(−t_r t_c), σ → t_r^{−s}𝐐^{−1}（パッケージ :665 の末尾、本文 :721-722）。"""
    return expr.subs({q: -tc, t: 1 / tr, A: a * sp.sqrt(-tr * tc), sig: tr ** (-s) / Qb}, simultaneous=True)


@lru_cache(maxsize=None)
def poincare(r, s, m, n):
    """𝒫_{[r^s]}(K_{m,n}; 𝐚, 𝐐, 𝐭_r, 𝐭_c)（パッケージの DoubleBraidPolynomial、KNTZ (final) :746-758）。"""
    return sp.expand(sp.cancel(sp.together(cov(sp.cancel(sp.together(poincare_pre(r, s, m, n))), s))))
