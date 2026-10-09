"""open_problems D1：Mellit 1704.07630 の Cor (daha) の文面 (−1)^n (a(qt)^{−1/2})^δ (1/(1−t)) P_{m,n}(1)[a−1]/(a−1) について、
P_{m,n}(1) を標準の正規化の EHA の元（BGLX の Q_{M,N}·1, Q_{N,M}·1、GN の \tilde P·1）と読み、変数の取り方 8 通り × 代入 5 通りで
Thm (dyck paths) と単項式倍で一致するものを総当たりで探す（T(2,3), T(3,2), T(2,5), T(3,4), T(4,3)）。
結果（logs/d1_mellit_cor.log）：一致するのは代入 [a/q−1] で Q（どちらの添字の順でも、(q,t) でも (t,q) でも）だけ、比はちょうど −1。[a−1] では一致する読み方が無い。
使い方：PYTHONPATH=src:scripts python scripts/d1_mellit_cor.py > logs/d1_mellit_cor.log
"""
from fractions import Fraction as Fr
import math
from d1_bglx import Qop, gn_vector
from d1_ops2 import pleth, dyck

def is_power(r, base):
    """r = ±base^i（i 整数, |i| ≤ 40）なら i を返す。"""
    for i in range(-40, 41):
        if abs(r) == base ** i:
            return i
    return None

def ratio(F, M, N, a, q, t, sq, apl):
    d = (M - 1) * (N - 1) // 2
    cor = (a / sq) ** d / (1 - t) * pleth(F(q, t), apl(a, q, t)) / (apl(a, q, t) - 1)
    return dyck(M, N, a, q, t, sq) / cor

def monomial(F, M, N, apl):
    # q = s², t = u² で √(qt) = s·u を有理数に
    pts = [(Fr(2, 7), Fr(2, 3), Fr(3, 5))]
    a, s, u = pts[0]
    r0 = ratio(F, M, N, a, s * s, u * u, s * u, apl)
    ra = ratio(F, M, N, 2 * a, s * s, u * u, s * u, apl)
    rs = ratio(F, M, N, a, (2 * s) ** 2, u * u, 2 * s * u, apl)
    ru = ratio(F, M, N, a, s * s, (2 * u) ** 2, s * 2 * u, apl)
    i, j, k = is_power(ra / r0, 2), is_power(rs / r0, 2), is_power(ru / r0, 2)
    if None in (i, j, k):
        return None
    # 第 2 の点で単項式の予測を確かめる
    a2, s2, u2 = Fr(5, 11), Fr(3, 7), Fr(7, 4)
    c = r0 / (a ** i * s ** j * u ** k)
    r2 = ratio(F, M, N, a2, s2 * s2, u2 * u2, s2 * u2, apl)
    return (c, i, j, k) if r2 == c * a2 ** i * s2 ** j * u2 ** k else None

cands = {}
for name, fn in {"Q_{M,N}": lambda M, N: (lambda q, t: Qop(M, N, {(): Fr(1)}, q, t)),
                 "Q_{N,M}": lambda M, N: (lambda q, t: Qop(N, M, {(): Fr(1)}, q, t)),
                 "GN(M,N)": lambda M, N: (lambda q, t: gn_vector(M, N, q, t)),
                 "GN(N,M)": lambda M, N: (lambda q, t: gn_vector(N, M, q, t))}.items():
    for vname, var in {"(q,t)": lambda q, t: (q, t), "(t,q)": lambda q, t: (t, q), "(1/q,1/t)": lambda q, t: (1 / q, 1 / t),
                       "(1/q,t)": lambda q, t: (1 / q, t), "(q,1/t)": lambda q, t: (q, 1 / t), "(1/t,1/q)": lambda q, t: (1 / t, 1 / q),
                       "(t,1/q)": lambda q, t: (t, 1 / q), "(1/t,q)": lambda q, t: (1 / t, q)}.items():
        cands[(name, vname)] = (fn, var)

for apname, apl in {"[a−1]": lambda a, q, t: a, "[a/q−1]": lambda a, q, t: a / q, "[aq−1]": lambda a, q, t: a * q,
                    "[a/t−1]": lambda a, q, t: a / t, "[at−1]": lambda a, q, t: a * t}.items():
    for (name, vname), (fn, var) in cands.items():
        res = []
        for M, N in [(2, 3), (3, 2), (2, 5), (3, 4), (4, 3)]:
            try:
                F = (lambda f, v: (lambda q, t: f(*v(q, t))))(fn(M, N), var)
                res.append(monomial(F, M, N, apl))
            except Exception as ex:
                res.append(f"err {ex}")
        if all(isinstance(r, tuple) for r in res):
            print(apname, name, vname, res, flush=True)
print("done")
