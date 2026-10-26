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


def test_composite22_decomposition_and_normalization():
    """二つの評価から解いた分解の係数が多項式の分解と一致し（N = 3）、規格化した値が KN の data と比 1（U(4) の [3,1], [2,2]）。"""
    import sympy as sp
    import rk_c_composite22_run as C
    import rk_k5_trefoil_2box as T
    from rcsknots.refined_kauffman.daha_dn import GLN, gl_P_poly, gl_adj_poly
    from rcsknots.refined_kauffman.daha_fast import solve_symmetric, poly_mul, gl_adj_fast
    N = 3
    D = GLN(N)
    inv = lambda f: {tuple(-x for x in a): v for a, v in f.items()}
    for alpha, beta in (((2,), (2,)), ((2,), (1, 1))):
        b = C.bvec(N, alpha, beta); k = -min(b)
        Pb = {tuple(x - k for x in a): v for a, v in gl_P_poly(D, tuple(x + k for x in b if x + k > 0)).items()}
        lhs = poly_mul(gl_P_poly(D, alpha), inv(gl_P_poly(D, beta)))
        sol = solve_symmetric(D, lhs, {"b": Pb, "adj": gl_adj_poly(D), "0": {tuple([0] * N): RFm.ONE}})
        c1, c0, Kb, _ = C.decomposition(N, alpha, beta)
        assert sol["b"] == RFm.conv(1) and (sol["adj"] - c1).is_zero() and (sol["0"] - c0).is_zero()
        assert (D.evaluate(Pb) - Kb).is_zero()
    N = 4
    V = L.basis_values(N, C.P2, C.P2)
    Zadj = gl_adj_fast(N, 1)
    for ab, sh in ((((2,), (1, 1)), "31"), (((1, 1), (1, 1)), "22")):
        Rb, _ = C.reduced(N, V, *ab, Zadj)
        W = (C.ss / C.rr)**(sum(int(c) for c in sh)) * T.Rdata(sh, N)
        assert sp.cancel(Rb / W) == 1
