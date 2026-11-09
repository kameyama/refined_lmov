"""代数的なケーブル（K1 = trefoil の (2,13)-ケーブル、K3 = (2,15)-ケーブル）の iterated の DAHA-Jones（progress.md 続き 140）。
機能テスト：GL_N の値が framed_knot の GN の値と単項式倍で一致、D_N の高速版が旧版と一致、
D_N の値が SU(2)²・SU(4) の一様な規格化の固定 N の値と framing の単項式 T_[1]^{−F}（F = p₂q₂）だけ違うこと。
results：大きな N の表現 [1] の g が整数 Laurent で普遍因子で割り切れること。"""
import json
import sys

import pytest
import sympy as sp

sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
from rcsknots.refined_kauffman.daha_iter import iter_jones_gl, iter_jones_dn, iter_jones_dn_fast

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
