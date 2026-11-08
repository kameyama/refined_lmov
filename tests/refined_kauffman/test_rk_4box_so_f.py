"""SO の自由エネルギーの次数 4 までの Möbius 反転（scripts/refined_kauffman/rk_c_4box_so_f.plog、progress 続き 72）の機能テスト。
  (1) 次数 ≤ 3 の成分が 3 箱の plog（rk_c_3box_largeN.plog）と一致すること。
  (2) 一粒子の f から作った Z = exp(Σ_d ψ_d f / d) に plog を当てると f に戻ること（次数 4 まで）。"""
import sys
from fractions import Fraction
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import rk_c_4box_so_f as F4
import rk_c_3box_largeN as T3
from rcsknots.refined_kauffman.flint_rf3 import R, S, U, conv
from rcsknots.partitions import partitions

ONE = conv(1)


def _toy_f():
    q, t = R * R, S * S
    return {(1,): R + U / S, (2,): q - 2 * t, (1, 1): U * R, (3,): S, (2, 1): R / U, (1, 1, 1): 3 * q,
            (4,): t * U, (3, 1): R * S, (2, 2): 1 / (1 - t), (2, 1, 1): U, (1, 1, 1, 1): R - S}


def _exp_pleth(f):
    """Z = exp(Σ_d ψ_d(f)/d)（次数 ≤ 4）。"""
    L = {}
    for d in (1, 2, 3, 4):
        for k, v in f.items():
            if d * sum(k) <= 4:
                kk = tuple(sorted((d * x for x in k), reverse=True))
                L[kk] = L.get(kk, conv(0)) + conv(Fraction(1, d)) * v.adams(d)
    Z, term = {(): ONE}, {(): ONE}
    for n in range(1, 5):
        term = T3.pmul(term, L, 4)
        Z = T3.padd(Z, term, conv(Fraction(1, 1)) / conv(_fact(n)))
    return Z


def _fact(n):
    return 1 if n <= 1 else n * _fact(n - 1)


def test_plog_inverts_exp():
    f = _toy_f()
    got = F4.plog(_exp_pleth(f))
    assert all((got.get(k, conv(0)) - v).is_zero() for k, v in f.items())


def test_plog_matches_3box_up_to_degree_3():
    Z = _exp_pleth(_toy_f())
    Z3 = {k: v for k, v in Z.items() if sum(k) <= 3}
    a, b = F4.plog(Z), T3.plog(Z3)
    assert all((a[k] - b[k]).is_zero() for n in (1, 2, 3) for k in partitions(n))
