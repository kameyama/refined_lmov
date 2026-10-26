"""Laurent 多項式の係数（flint_lp）と評価の汎関数による高速版（daha_fast_lp）の機能テスト。
  (1) LP の演算が RF と同じ値を与えること。
  (2) composite (2;2), (2;11), (11;11) の生の値が、旧版（rk_c_fast_composite22.comp_raw2、RF）と全係数で一致すること（N = 3, 4、final = poly と functional）。"""
import sys
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
from rcsknots.refined_kauffman import flint_lp as LPm, flint_rf as RFm
from rcsknots.refined_kauffman import daha_fast_lp as L


def test_lp_arithmetic():
    R, S = LPm.R, LPm.S
    e = (S - S ** -1) * (R ** 3 + 2) - R ** -2 * S + 1 / (R * S ** 2)
    r, s = RFm.R, RFm.S
    ref = (s - 1 / s) * (r ** 3 + 2) - r ** -2 * s + 1 / (r * s ** 2)
    assert (e.to_rf() - ref).is_zero()
    assert (e - e).is_zero() and not e.is_zero()


def test_composite22_matches_rf_version():
    import rk_c_fast_composite22 as C
    P2 = [(2,), (1, 1)]
    for N in (3, 4):
        olds = {ab: C.comp_raw2(N, *ab) for ab in (((2,), (2,)), ((2,), (1, 1)), ((1, 1), (1, 1)))}
        for final in ("poly", "functional"):
            V = L.basis_values(N, P2, P2, final=final)
            for ab, old in olds.items():
                assert (L.composite_raw(N, *ab, V) - old).is_zero(), (N, final, ab)
