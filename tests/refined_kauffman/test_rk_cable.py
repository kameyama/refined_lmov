"""代数的なケーブル（K1 = trefoil の (2,13)-ケーブル、K3 = (2,15)-ケーブル）の iterated の DAHA-Jones（progress.md 続き 140）。
機能テスト：GL_N の値が framed_knot の GN の値と単項式倍で一致、D_N の高速版が旧版と一致、
D_N の値が SU(2)²・SU(4) の一様な規格化の固定 N の値と framing の単項式 T_[1]^{−F}（F = p₂q₂）だけ違うこと。
companion の段の約分なしの版（続き 141）が旧版と一致すること、GL の規格化が T_b^{−F} × (|b| だけによる単項式) の形であること。
results：大きな N の表現 [1] の g が整数 Laurent で普遍因子で割り切れること。"""
import json
import sys

import pytest
import sympy as sp

sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
from rcsknots.refined_kauffman.daha_iter import (iter_jones_gl, iter_jones_dn, iter_jones_dn_fast, iter_jones_dn_lp,
                                                 iter_adj_gl, iter_adj_gl_lp)

r, s = sp.symbols("r s", positive=True)
mir = lambda e: e.subs({r: 1 / r, s: 1 / s}, simultaneous=True)


def _kn(cable, c, n):
    d = json.load(open(f"data/framed_knot/cables/{cable}_{c}.json"))
    A = s**n
    return sp.expand(sum(v * A**(-int(k.split(",")[0])) * r**(-int(k.split(",")[1])) * s**(-int(k.split(",")[2])) for k, v in d["poly"].items()))


def _is_monomial(e):
    n_, d_ = sp.fraction(sp.factor(sp.cancel(e)))
    return sp.Poly(n_, r, s).is_monomial and sp.Poly(d_, r, s).is_monomial


@pytest.mark.parametrize("cable", ["K1", "K3"])
def test_gl_matches_gn(cable):
    for N, lam in ((2, (1,)), (3, (1,)), (2, (2,)), (3, (1, 1))):
        Z, K = iter_jones_gl(N, lam, cable)
        red = mir(sp.cancel((Z / K).to_sympy()))
        assert _is_monomial(red / _kn(cable, "".join(map(str, lam)), N)), (cable, N, lam)


def test_dn_fast_matches_slow():
    for cable in ("K1", "K3"):
        for N in (2, 3, 4):
            Zs, Ks = iter_jones_dn(N, (1,), cable)
            Zf, Kf = iter_jones_dn_fast(N, (1,), cable)
            assert (Zs - Zf).is_zero() and (Ks - Kf).is_zero()


def test_dn_matches_fixed_N_uniform():
    from rk_k5_framing_scan import P_sym
    for cable, F, delta in (("K1", 26, 8), ("K3", 30, 9)):
        w = (s / r)**delta
        Wso = {2: w**2 * (_kn(cable, "1", 2) * P_sym((1,), 2))**2, 3: w**2 * _kn(cable, "11", 4) * P_sym((1, 1), 4)}
        for N in (2, 3):
            Z, _ = iter_jones_dn_fast(N, (1,), cable)
            rat = sp.cancel(mir(Z.to_sympy()) / Wso[N])
            assert rat == r**F * s**(2 * F * (N - 1)), (cable, N, rat)        # 鏡像(T_[1]^{−F})、T_[1] = q^{1/2} t^{N−1}


@pytest.mark.parametrize("cable", ["K1", "K3"])
def test_lp_matches_old(cable):
    for N in (2, 3):
        for lam in ((2,), (1, 1)):
            Zs, Ks = iter_jones_dn_fast(N, lam, cable)
            Zl, Kl = iter_jones_dn_lp(N, lam, cable)
            assert (Zs - Zl).is_zero() and (Ks - Kl).is_zero(), (N, lam)
        assert (iter_adj_gl(N, cable) - iter_adj_gl_lp(N, cable)).is_zero(), N


def test_gl_rule_small():
    """M_b = Z_b / 鏡像(一様な規格化の W_b) = T_b^{−F} t^{c|b|(N−1)}（rk_c_cable_gl_rule.py の小さい場合）。"""
    import rk_c_cable_gl_rule as G
    import rk_c_2box_largeN as L2
    u = sp.Symbol("u", positive=True)
    for cable, F, delta, c in (("K1", 26, 8, 10), ("K3", 30, 9, 12)):
        for N, lam in ((2, (1,)), (2, (2,)), (3, (1, 1))):
            Z, _ = iter_jones_gl(N, lam, cable)
            W = (L2.S / L2.R)**(delta * sum(lam)) * L2.Rdata_nu("".join(map(str, lam)), cable) * L2.Psym(lam)
            W = sp.cancel(W.to_sympy().subs(u, s**(N + 1)))
            b = list(lam) + [0] * (N - len(lam))
            bb = sum(x * x for x in b); br = sum(x * (N - 1 - i) for i, x in enumerate(b))
            assert sp.cancel(Z.to_sympy() / mir(W) - r**(-F * bb) * s**(-2 * F * br + 2 * c * sum(lam) * (N - 1))) == 0, (cable, N, lam)


@pytest.mark.results
@pytest.mark.parametrize("cable", ["K1", "K3"])
def test_large_N_one_box(cable):
    import rk_c_cable_1box as C
    import rk_b_2box_structure as B
    g = C.g_one(cable)
    assert g.is_integer_laurent()
    u = sp.Symbol("u", positive=True)
    e = sp.cancel(g.to_sympy().subs(u, B.nu * B.s).subs({r: B.r, s: B.s}, simultaneous=True))
    assert B.divisible_universal(e)


def test_structure_flint_matches_sympy():
    """flint 版の 2 箱の判定（rk_c_2box_structure_flint.py、続き 141）が sympy 版（rk_b_2box_structure）と同じ答えを返す（trefoil）。"""
    import rk_c_2box_structure_flint as FL
    import rk_b_2box_structure as B
    import rk_b_2box_structure_m as BM
    import rk_c_2box_largeN as L
    import rk_c_2box_largeN_m as M
    W, co = M.so_parts(1)
    f = L.so_f(W, co)
    fc = M.cover_f(1, "T2_3")
    G = FL.hat({k: FL.to_nu(f[k] - FL.half * fc[k]) for k in f})
    Gs = BM.hat(BM.load(1)[2])
    for k in ((1,), "sum"):
        assert FL.c1_symmetric(G[k]) == B.c1_symmetric(Gs[k])
        assert FL.divisible_universal(G[k]) == B.divisible_universal(Gs[k])
        assert FL.sign_blocks(G[k]) == B.sign_blocks(Gs[k])


@pytest.mark.results
def test_large_N_two_box_K1():
    """K1 の大きな N の 2 箱：ĝ_(2) + ĝ_(1,1) が整数 Laurent、ĝ の refined だけの部分が普遍因子で割り切れる。"""
    import rk_c_2box_structure_flint as FL
    import rk_c_2box_largeN as L
    import rk_c_2box_largeN_m as M
    W, co = M.so_parts("K1")
    f = L.so_f(W, co)
    fc = M.cover_f("K1", "K1")
    G = FL.hat({k: FL.to_nu(f[k] - FL.half * fc[k]) for k in f})
    assert G["sum"].is_integer_laurent()
    assert all(FL.divisible_universal(G[k]) for k in ((1,), (2,), (1, 1), "sum"))
