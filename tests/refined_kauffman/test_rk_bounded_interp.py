"""大きな N の補間に要る N を減らす一次方程式（scripts/refined_kauffman/rk_c_bounded_interp.py、progress 続き 71）の機能テスト。
  (1) 人工の c_j(r, s) から作った R(N) で、使う N の個数が α ごとの U の冪の項数以上なら一意に元へ戻り、足りなければ「一意でない」になること。
  (2) 保存値（SO の 3 箱 (3)、trefoil）で N = 2〜6 から一意に決まり、使わない N = 7 で全係数一致すること。"""
import sys
from fractions import Fraction
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import rk_c_bounded_interp as B


def _synthetic():
    # r^0：U^{-3}〜U^{-1}（3 項）、r^2：U^{-2}〜U^{-1}（2 項）
    sol = {-3: {(0, 4): Fraction(1), (0, 6): Fraction(-2)}, -2: {(0, 2): Fraction(3), (2, 2): Fraction(1)},
           -1: {(0, 0): Fraction(1), (2, 0): Fraction(5), (2, 2): Fraction(-1)}}
    T = {N: B.predict(sol, N) for N in range(2, 9)}
    return sol, T


def test_synthetic_recovery():
    sol, T = _synthetic()
    got, st, _ = B.solve_bounded(T, [3, 4, 5], -3, -1, "alpha")
    assert st == "一意"
    assert {j: c for j, c in got.items()} == sol
    assert all(B.predict(got, N) == T[N] for N in T)


def test_synthetic_too_few_N():
    _, T = _synthetic()
    _, st, _ = B.solve_bounded(T, [3, 4], -3, -1, "alpha")
    assert st == "一意でない"


def test_so_three_box_sym():
    T = {N: B.terms_true(N, 1, (3,)) for N in range(2, 8)}
    rng, d = B.slopes(T)
    assert d[-1] == (-15, -6)
    sol, st, _ = B.solve_bounded(T, [2, 3, 4, 5, 6], -15, -6, "alpha")
    assert st == "一意"
    assert B.predict(sol, 7) == T[7]
