"""ダブルツイスト結び目 K_{m,n} の彩色 HOMFLY-PT（t = q）の文献の式（別研究 [2,1] の X0）。

変数は文献の (A, q)：skein A H(+) − A^{-1} H(−) = (q − q^{-1}) H(0)、自明な結び目で 1（reduced）、framing 0（KNTZ :537-552）。
本研究の (a, Q) へは A² = a, q = Q（KN の t = q、Q = q_KN^{1/2}。conventions.md §14）。
{x} = x − x^{-1}、D_n = {A q^n}、[n] = {q^n}/{q}。行番号は references/<id>/src の TeX（1903.00259, 1612.00422 は CR を改行に直した *.lf.tex）。

- KNTZ 1902.02275（Cyclotomic-Rectangular.tex）：長方形 [r^s] の式。B :612-626、F :631-640、Tw :644-648、H :649-656。
- Morozov 1903.00259（main.lf.tex）：[2,1] の 6×6 の三角行列 𝓑 :322-339、F_X^{(m)} = Σ_Y (𝓑^{m+1})_{XY} :154-157、Z 因子 :385-390、H = Σ_X Z^X F_X^{(m)} F_X^{(n)} :139-151（Z を F^{(1)}F^{(−1)} で割って読む、H21_1903）。
- Morozov 1612.00422（1612.00422v4.lf.tex）：[2,1] の微分展開 (H21de) :1249-1271、Z 因子 (Zf21) :1275-1289、F 因子 (Ftwist21) :1309-1343、固有値 :1345-1349。
"""
from functools import lru_cache

import sympy as sp

from .partitions import boxes, arm, leg, conjugate

A, q = sp.symbols("A q", positive=True)


def br(x):
    return x - 1 / x


def Dn(n):
    return br(A * q**n)


def qnum(n):
    return br(q**n) / br(q)


def _contains(lam, mu):
    return len(mu) <= len(lam) and all(m <= l for m, l in zip(mu, lam))


def subdiagrams(lam):
    """λ に含まれる全てのヤング図形（∅ を含む、大きさの昇順）。"""
    out = [()]
    for k in range(1, sum(lam) + 1):
        from .partitions import partitions
        out += [mu for mu in partitions(k) if _contains(lam, mu)]
    return out


# ---------------------------------------------------------------- KNTZ（長方形）
def Dprinc(lam, s, qq=q):
    """D_λ^{(s)}(q) = s_λ(q^{s−1}, …, q^{1−s}) = Π {q^{s − l' + a'}}/{q^{l + a + 1}}（KNTZ :572-579）。"""
    out = sp.Integer(1)
    for i, j in boxes(lam):
        out *= br(qq ** (s - i + j)) / br(qq ** (leg(lam, i, j) + arm(lam, i, j) + 1))
    return out


def gauss(n, k, v):
    """ガウスの二項係数 [n choose k]_v（v = q²。B_{λ,∅} = D_λ になるのはこの規約、テストで確認）。"""
    if k < 0 or k > n:
        return sp.Integer(0)
    num = den = sp.Integer(1)
    for i in range(k):
        num *= 1 - v ** (n - i)
        den *= 1 - v ** (i + 1)
    return num / den


def Bkntz(lam, mu, s):
    """B^{(s)}_{λ,μ}(A, q)（KNTZ :612-626）。"""
    if not _contains(lam, mu):
        return sp.Integer(0)
    L = list(lam) + [0] * (s - len(lam))
    M = list(mu) + [0] * (s - len(mu))
    mat = sp.Matrix(s, s, lambda i, j: gauss(L[i] + s - (i + 1), M[j] + s - (j + 1), q**2))
    pre = A ** (2 * sum(mu)) * q ** ((1 - s) * (sum(lam) - sum(mu)))
    for i, j in boxes(mu):
        pre *= q ** (4 * (j - i))
    return pre * mat.det()


@lru_cache(maxsize=None)
def F_kntz(lam, m, s=None):
    """F^{(m)}_λ(A, q)（m ≥ 1、KNTZ :631-640）。入れ子の和 = (B^{m−1} D)_λ / D_λ。"""
    s = s or max(len(lam), 1)
    subs = subdiagrams(lam)
    v = {mu: Dprinc(mu, s) for mu in subs}
    for _ in range(m - 1):
        v = {nu: sum(Bkntz(nu, mu, s) * v[mu] for mu in subs if _contains(nu, mu)) for nu in subs}
    return sp.factor(sp.cancel(v[lam] / Dprinc(lam, s)))


def Tw(lam, m):
    """Tw^{(m)}_λ（KNTZ :644-648）。"""
    if m > 0:
        pre = (-A**2) ** sum(lam)
        for i, j in boxes(lam):
            pre *= q ** (2 * arm(lam, i, j) - 2 * leg(lam, i, j))
        return pre * F_kntz(lam, m)
    return F_kntz(lam, -m).subs({A: 1 / A, q: 1 / q}, simultaneous=True)


def H_rect(r, s, m, n):
    """[r^s] 彩色の reduced HOMFLY-PT of K_{m,n}（KNTZ :649-656。Fang–Zhou 2610.00303 で定理）。"""
    R = (r,) * s
    total = 0
    for lam in subdiagrams(R):
        term = Tw(lam, m) * Tw(lam, n) / (Tw(lam, 1) * Tw(lam, -1))
        term *= Dprinc(lam, s) * Dprinc(conjugate(lam), r, 1 / q)
        for i, j in boxes(lam):
            term *= br(A * q ** (r + j - i)) * br(A * q ** (j - i - s))
        total += term
    return sp.expand(sp.cancel(sp.together(total)))


# ---------------------------------------------------------------- Morozov 1903.00259（[2,1]）
CH21 = ["0", "1", "11", "2", "21", "X2"]   # ∅, [1], [1,1], [2], [2,1], X̃₂ の順（main.lf.tex :322-339）


def B21():
    """𝓑^{[2,1]}（1903.00259 main.lf.tex :322-339 を書き写したもの）。"""
    q2 = qnum(2)
    q3 = qnum(3)
    return sp.Matrix([
        [1, 0, 0, 0, 0, 0],
        [-A**2, A**2, 0, 0, 0, 0],
        [A**4 / q**2, -q2 * A**4 / q**3, A**4 / q**4, 0, 0, 0],
        [q**2 * A**4, -q2 * q**3 * A**4, 0, q**4 * A**4, 0, 0],
        [-A**6, q3 * A**6, -q3 * A**6 / (q2 * q), -q3 * q * A**6 / q2, A**6, 0],
        [-A**6, q3 * A**6, (A**2 - q**6) * A**4 / (q**2 * (q**4 - 1)), -(A**2 * q**6 - 1) * A**4 / (q**4 - 1), 0, A**4],
    ])


def Z21_1903():
    """Z^X_{[2,1]}（1903.00259 main.lf.tex :385-390）。"""
    q2, q3 = qnum(2), qnum(3)
    return [sp.Integer(1),
            (q3 * Dn(0)**2 + q3**2 * Dn(2) * Dn(-2)) / q2**2,
            q3 / q2 * Dn(2) * Dn(0) * Dn(-2) * Dn(-3),
            q3 / q2 * Dn(3) * Dn(2) * Dn(0) * Dn(-2),
            Dn(3) * Dn(2) * Dn(1) * Dn(-1) * Dn(-2) * Dn(-3),
            -q3**2 * br(q)**4 * Dn(2) * Dn(-2)]


@lru_cache(maxsize=None)
def F21_1903(m):
    """F_X^{(m)} = Σ_Y (𝓑^{m+1})_{XY}（X は CH21 の順。m < −1 は逆行列の冪）。"""
    P = B21() ** (m + 1)
    return [sp.factor(sp.cancel(sum(P[i, j] for j in range(6)))) for i in range(6)]


def H21_1903(m, n):
    """1903.00259 の (difexpan) を、Z^X を F^{(1)}_X F^{(−1)}_X で割る正規化で読んだもの（K_{m,n} は KNTZ と同じ付け方）。

    箱の式（main.lf.tex :139-151）には分母 F^{(1)}F^{(−1)} がコメントアウトされて残っている（:148）。
    書かれたとおり Σ Z F^{(m)} F^{(n)} と読むと、一致するのは H^{(m,−1)} = K_{m,1}（F^{(−1)} = 1 のため）だけで、
    K_{m,2} では合わない（X0 の計算、conventions.md §14）。F^{(1)}_X は単項式なので、これは Z の正規化の違い。"""
    Z = Z21_1903()
    Fm, Fn, F1, Fm1 = F21_1903(m), F21_1903(n), F21_1903(1), F21_1903(-1)
    return sp.expand(sp.cancel(sp.together(sum(Z[i] * Fm[i] * Fn[i] / (F1[i] * Fm1[i]) for i in range(6)))))


# ---------------------------------------------------------------- Morozov 1612.00422（[2,1]）
def _lam2():
    """Λ² の値（1612.00422v4.lf.tex :1345-1349）。"""
    return {"0": sp.Integer(1), "1": A**2, "X2": A**4, "2": q**4 * A**4, "11": A**4 / q**4, "21": A**6}


@lru_cache(maxsize=None)
def F21_1612(m):
    """(Ftwist21)（1612.00422v4.lf.tex :1309-1343）。順は ∅, [1], [1,1], [2], [2,1], X̃₂（X̃₂ は \\tF_{[2,1]}）。"""
    L = {k: v**m for k, v in _lam2().items()}
    q2, q3, q4 = qnum(2), qnum(3), qnum(4)
    D = Dn
    F1 = A * (L["0"] / D(0) - L["1"] / D(0))
    F2 = q * A**2 * (L["0"] / (D(1) * D(0)) - q2 * L["1"] / (D(2) * D(0)) + L["2"] / (D(2) * D(1)))
    F11 = A**2 / q * (L["0"] / (D(0) * D(-1)) - q2 * L["1"] / (D(0) * D(-2)) + L["11"] / (D(-1) * D(-2)))
    F21 = A**3 * (L["0"] / (D(1) * D(0) * D(-1)) - q3 * L["1"] / (D(2) * D(0) * D(-2))
                  + q3 / q2 * L["2"] / (D(2) * D(1) * D(-2)) + q3 / q2 * L["11"] / (D(2) * D(-1) * D(-2))
                  - L["21"] / (D(2) * D(0) * D(-2)))
    tF21 = A**3 * (L["0"] / (D(1) * D(0) * D(-1)) - q4 / q2 * L["1"] / (D(2) * D(0) * D(-2))
                   - 1 / (q2**2 * br(q)**2) * (D(3) * L["2"] / (D(2) * D(1)) - 2 * L["X2"] / D(0)
                                               + D(-3) * L["11"] / (D(-1) * D(-2))))
    return [sp.Integer(1)] + [sp.factor(sp.cancel(sp.together(f))) for f in (F1, F11, F2, F21, tF21)]


def Z21_1612():
    """(Zf21)（1612.00422v4.lf.tex :1275-1289。下線の微分は defect 0 の結び目では残す）と、(H21de) の係数 [3]/[2]（[2], [1,1] の項）。"""
    q2, q3 = qnum(2), qnum(3)
    D = Dn
    return [sp.Integer(1),
            D(3) * D(-3) + D(2) * D(0) + D(0) * D(-2),
            q3 / q2 * D(2) * D(0) * D(-2) * D(-3),
            q3 / q2 * D(3) * D(2) * D(0) * D(-2),
            D(3) * D(2) * D(1) * D(-1) * D(-2) * D(-3),
            -q3**2 * br(q)**4 * D(2) * D(-2)]


def H21_1612(m, n):
    """(H21de)。[2] の項の分母は本文では F^{(2)}_{[2]} F^{(−1)}_{[1]} と印刷されているが、他の項と同じ F^{(1)} F^{(−1)} と読む（誤植と判断、テストで確認）。"""
    Z = Z21_1612()
    Fm, Fn, F1, Fm1 = F21_1612(m), F21_1612(n), F21_1612(1), F21_1612(-1)
    tot = sum(Z[i] * Fm[i] * Fn[i] / (F1[i] * Fm1[i]) for i in range(6))
    return sp.expand(sp.cancel(sp.together(tot)))


def to_kn(expr):
    """文献の (A, q) → 本研究の (a, Q)：A² = a, q = Q（conventions.md §14）。"""
    from .homfly import a, Q
    return sp.expand(expr.subs({A: sp.sqrt(a), q: Q}, simultaneous=True))
