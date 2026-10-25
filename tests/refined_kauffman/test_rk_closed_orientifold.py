"""閉弦の refined orientifold（AgSc §6、AM 2205.12832 §3）の分解と、開弦の判定式のずれの向き（progress.md 続き 38）。"""
import sys
sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_b_closed_orientifold as C
import rk_b_shift_direction as S


def test_agsc_large_n_matches_product_formula():
    assert C.check_A(range(2, 6))


def test_am_klein_bottle_is_frozen_plus_invariant():
    assert C.check_B()


def test_shift_direction_unique_for_trefoil_one_box():
    import rk_b_uniform as Bu
    from rcsknots.refined_kauffman.kauffman_trefoil import r as r0, s as s0, U as U0
    g = sp.cancel(Bu.g_true(1).subs({r0: S.r, s0: S.s}).subs(U0, S.nu * S.s))
    ok = [k for k in range(-2, 3) if all(S.divisible(g, k).values())]
    assert ok == [1]
