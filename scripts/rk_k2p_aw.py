r"""K2'(a)：probe 側の Askey–Wilson 多項式 K_{(n)}(y; t, q; 1, −1, q^{1/2}, −q^{1/2}) の上位係数 c_{2j}(n)（y^{n−2j} の係数）。
(1) 三項漸化式（askey_wilson_monic）が q-差分作用素の固有関数（koornwinder）と一致するか（n ≤ 4）。
(2) c_{2j}(n) が U_n = t^n の有理式か：有理数 (q,t) で n ≤ 40 の値を有理補間し、次数と分母を見る。
(3) 閉じた形 top_coeff_closed（予想）を複数の有理点で j ≤ 5、n ≤ 30 で照合。
"""
import sys
sys.path.insert(0, "src")
from fractions import Fraction as F
import sympy as sp
from rcsknots.askey_wilson import askey_wilson_monic, top_coeff_closed
from rcsknots.koornwinder import koornwinder


def solve(M, b):
    n = len(M); m = len(M[0])
    A = [row[:] + [bb] for row, bb in zip(M, b)]
    piv = []; r = 0
    for c in range(m):
        p = next((i for i in range(r, n) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        inv = 1 / A[r][c]; A[r] = [v * inv for v in A[r]]
        for i in range(n):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [vi - f * vr for vi, vr in zip(A[i], A[r])]
        piv.append(c); r += 1
        if r == n: break
    if any(A[i][m] != 0 for i in range(r, n)): return None
    x = [F(0)] * m
    for i, c in enumerate(piv): x[c] = A[i][m]
    return x


def ratfit(pts, dp, dq):
    """v (1 + Σ_{k≥1} q_k U^k) = Σ_j p_j U^j を dp+dq+1 点で解き、全点で検証する。"""
    M = [[U**j for j in range(dp + 1)] + [-v * U**k for k in range(1, dq + 1)] for U, v in pts[: dp + dq + 1]]
    x = solve(M, [v for _, v in pts[: dp + dq + 1]])
    if x is None: return None
    P, Q = x[: dp + 1], [F(1)] + x[dp + 1:]
    ok = all(sum(P[j] * U**j for j in range(dp + 1)) == v * sum(Q[k] * U**k for k in range(dq + 1)) for U, v in pts)
    return (P, Q) if ok else None


def find_rational(pts, maxdeg):
    for tot in range(0, 2 * maxdeg + 1):
        for dq in range(0, min(tot, maxdeg) + 1):
            dp = tot - dq
            if dp > maxdeg or dp + dq + 1 > len(pts) - 2: continue
            r = ratfit(pts, dp, dq)
            if r: return (dp, dq) + r
    return None


if __name__ == "__main__":
    q = sp.Rational(9, 49); sq = sp.Rational(3, 7); t = sp.Rational(5, 11); y = sp.Symbol("y")
    aw = askey_wilson_monic(4, t, sp.Integer(1), sp.Integer(-1), sq, -sq)
    for n in range(1, 5):
        K = sp.expand(koornwinder((n,), 1, t, q, (1, -1, sq, -sq), (y,)))
        print(f"(1) n={n}: 漸化式 = 固有関数:", sp.expand(K - sum(v * y**j for j, v in aw[n].items())) == 0)
    qf, sqf, tf = F(9, 49), F(3, 7), F(5, 11)
    aw = askey_wilson_monic(40, tf, F(1), F(-1), sqf, -sqf)
    Us = sp.Symbol("U")
    for j in (1, 2, 3, 4):
        pts = [(tf**n, aw[n].get(n - 2*j, F(0))) for n in range(2*j, 41)]
        dp, dq, P, Q = find_rational(pts, 8)
        Qd = sp.factor(sum(sp.Rational(c.numerator, c.denominator) * Us**i for i, c in enumerate(Q)))
        print(f"(2) j={j}: 有理式、次数 (分子,分母) = ({dp},{dq})、分母 = {Qd}")
    for q_, sq_, t_ in [(F(9,49), F(3,7), F(5,11)), (F(4,9), F(2,3), F(7,13)), (F(25,169), F(5,13), F(3,17)), (F(49,4), F(7,2), F(8,5))]:
        aw = askey_wilson_monic(30, t_, F(1), F(-1), sq_, -sq_)
        bad = [(j, n) for j in range(6) for n in range(2*j, 31) if aw[n].get(n - 2*j, F(0)) != top_coeff_closed(j, t_**n, q_, t_)]
        print(f"(3) q={q_} t={t_}: 閉じた形との不一致 {bad}")
