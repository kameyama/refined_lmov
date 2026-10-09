"""別研究 [2,1] の X2：正則なチャネルの Z̃_X（X = [1], [1,1], [2], [2,1]）を、特殊化での値と対称性から Laurent 多項式として解き、解の空間の次元を数える。

Z̃_X = Σ c_{ijk} A^i Q^j T^k（i は Morozov の Z_X と同じ A の次数の範囲、|j|, |k| ≤ K）について、次の線形条件を課す：
  (a) t = q（Q = T）：Morozov 1612.00422 の Z_X（KN の変数で A² = a, q_Mor = Q）
  (b) EGP（Q = 1）と (d) U(2)（A = T²）：scripts/x2_specializations.py で一意に決まった値
  (c), (e)：(b), (d) を転置したもの（T = 1、A = Q^{−2}）。転置 Z̃_X(A, 1/T, 1/Q) = Z̃_{X^T}(A, Q, T)
  (f) ミラー Z̃_X(1/A, 1/Q, 1/T) = Z̃_X(A, Q, T)
  (g) 転置で自分に戻るチャネル（[1], [2,1]）では c_{ijk} = c_{i,−k,−j}
K を増やして、解があるか・解の空間の次元を出す。[1,1] は [2] から転置で決まるので解かない。
使い方：python scripts/x2_window.py KMAX [チャネル ...]（例：1 21。[2] は Morozov の Z に [3]/[2] の分母があり Laurent の仮定に合わない）
"""
import sys
import time

import flint
import sympy as sp

from rcsknots import double_twist as U
from rcsknots.double_twist_family import A, Q, T
from x2_specializations import solve_at

Z_MOR = dict(zip([(), (1,), (1, 1), (2,), (2, 1)], U.Z21_1612()[:5]))
TRANS = {(1,): (1,), (2, 1): (2, 1), (2,): (1, 1), (1, 1): (2,)}


def terms(expr, gens):
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    P, D = sp.Poly(num, *gens), sp.Poly(den, *gens)
    assert len(D.terms()) == 1, f"分母が単項式でない：{den}"
    (dexp, dc), = D.terms()
    return {tuple(x - y for x, y in zip(k, dexp)): sp.Rational(c, dc) for k, c in P.terms()}


def specialization_data():
    egp = solve_at({Q: 1}, "EGP")[0]
    u2 = solve_at({A: T**2}, "U(2)")[0]
    z = lambda X: sp.Symbol(f"z_{''.join(map(str, X)) or '0'}")  # noqa: E731
    return {X: egp[z(X)] for X in TRANS}, {X: u2[z(X)] for X in TRANS}


DEN = {(2,): Q**2 * T**2 - 1}    # [2] の分母（t = q で [2] = {q²}/{q} の分母を与える refined の候補）。N = Z̃·DEN を Laurent として解く


def solve_channel(X, K, egp, u2):
    den = DEN.get(X, sp.Integer(1))
    denT = den.subs({Q: 1 / T, T: 1 / Q}, simultaneous=True)       # 転置した分母（[2] ↔ [1,1] で同じ形）
    mor = Z_MOR[X].subs({U.A: A, U.q: Q}) * den.subs(T, Q)
    Ad = max(abs(k[0]) for k in terms(mor, (A, Q)))
    mons = [(i, j, k) for i in range(-Ad, Ad + 1) for j in range(-K, K + 1) for k in range(-K, K + 1)]
    idx = {m: n for n, m in enumerate(mons)}
    rows, rhs = [], []

    def add(target_map, data):
        """target_map(i, j, k) → 像の単項式。data：{像の単項式: 係数}。"""
        groups = {}
        for m in mons:
            groups.setdefault(target_map(*m), []).append(idx[m])
        for key in set(groups) | set(data):
            if key not in groups:
                rows.append({}), rhs.append(None)      # 窓に入らない：解なし
                continue
            rows.append({n: 1 for n in groups[key]})
            rhs.append(data.get(key, 0))

    XT = TRANS[X]
    add(lambda i, j, k: (i, j + k), terms(mor, (A, Q)))                                            # (a) Q = T
    add(lambda i, j, k: (i, k), terms(egp[X] * den.subs(Q, 1), (A, T)))                           # (b) Q = 1
    add(lambda i, j, k: (i, j), terms(egp[XT].subs(T, 1 / Q) * den.subs(T, 1), (A, Q)))             # (c) T = 1
    add(lambda i, j, k: (j, 2 * i + k), terms(u2[X] * den, (Q, T)))                                 # (d) A = T²
    add(lambda i, j, k: (j - 2 * i, k), terms(u2[XT].subs({Q: 1 / T, T: 1 / Q}, simultaneous=True) * den, (Q, T)))  # (e) A = 1/Q²
    # (f) ミラー：Z̃(1/x) = Z̃(x)。N = Z̃·den で den(1/x) = sg·Q^a T^b·den(x) なら c_{ijk} = sg·c_{−i, −j−a, −k−b}
    ratio = sp.factor(sp.cancel(den.subs({A: 1 / A, Q: 1 / Q, T: 1 / T}, simultaneous=True) / den))
    ((_, ex, ey), sg), = terms(ratio, (A, Q, T)).items()
    for m in mons:
        mm = (-m[0], -m[1] - ex, -m[2] - ey)
        if mm in idx and (m < mm or (m == mm and sg == -1)):
            rows.append({idx[m]: 1, idx[mm]: -sg} if m != mm else {idx[m]: 1}), rhs.append(0)
    if XT == X:                                                                                     # (g) 転置
        for (i, j, k) in mons:
            t = (i, -k, -j)
            if (i, j, k) < t and t in idx:
                rows.append({idx[(i, j, k)]: 1, idx[t]: -1}), rhs.append(0)
    if any(r is None for r in rhs):
        return None
    M = flint.fmpq_mat(len(rows), len(mons) + 1)
    for r, (row, b) in enumerate(zip(rows, rhs)):
        for c, v in row.items():
            M[r, c] = int(v)
        M[r, len(mons)] = flint.fmpq(int(sp.numer(b)), int(sp.denom(b)))
    R, rank = M.rref()
    pivots = []
    for r in range(rank):
        c = next(c for c in range(len(mons) + 1) if R[r, c] != 0)
        pivots.append(c)
    if len(mons) in pivots:
        return None                                                                                 # 矛盾
    dim = len(mons) - rank
    sol = {mons[c]: R[r, len(mons)] for r, c in enumerate(pivots) if R[r, len(mons)] != 0}
    return dim, sol


if __name__ == "__main__":
    KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    chans = [tuple(int(c) for c in x) for x in sys.argv[2:]] or [(1,), (2, 1)]
    egp, u2 = specialization_data()
    for X in chans:
        for K in range(1, KMAX + 1):
            t0 = time.time()
            res = solve_channel(X, K, egp, u2)
            if res is None:
                print(f"{X} K={K}：解なし（{time.time() - t0:.1f}s）", flush=True)
                continue
            dim, sol = res
            print(f"{X} K={K}：解の空間の次元 {dim}（{time.time() - t0:.1f}s）", flush=True)
            if dim == 0:
                expr = sum(sp.Rational(int(v.p), int(v.q)) * A**i * Q**j * T**k for (i, j, k), v in sol.items())
                print(f"   Z̃_{X} = {sp.factor(expr)}", flush=True)
