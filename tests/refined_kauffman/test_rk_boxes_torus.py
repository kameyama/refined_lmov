"""箱を増やした検証（progress.md 続き 139）のテスト。
機能テスト：composite (α;1̄)（|α| = 3）の分解（Pieri の ψ'、psi_A）と KN の data との一致、|α| = 2 の既存の値の再現、
  U(N_A) での固定 N の分解の重みの絞り込み。
results：T(2,5) の 3 箱と trefoil の 4 箱の大きな N の g の整数性・(1 + t) の極の位置・次数 ≤ 3 の一致（保存済みの結果の再判定）。"""
import os
import pickle
import sys

import pytest
import sympy as sp

sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman"); sys.path.insert(1, "scripts/framed_knot")
import rk_c_fast_composite3 as C3

rr, ss = C3.rr, C3.ss


def test_lower_and_psi_A():
    assert C3.lower((2, 1)) == [(0, (1, 1)), (1, (2,))]
    assert C3.lower((3,)) == [(0, (2,))]
    assert C3.lower((1, 1, 1)) == [(2, (1, 1))]
    A = sp.Symbol("A")
    for alpha in ((3,), (2, 1), (1, 1, 1), (2,), (1, 1)):
        for i, mu in C3.lower(alpha):
            for N in (4, 5, 6):
                cA = sp.sympify(C3.psi_A(alpha, i, rr**2, ss**2, A)).subs(A, ss**(2 * N))
                assert sp.cancel(cA - C3.c_pieri(alpha, i, N)) == 0, (alpha, mu, N)


def test_composite3_matches_polynomial_and_old_values():
    # N = 4 で多項式の分解と照合（comp_reduced の check）、|α| = 2 は保存済みの値を再現
    for alpha in ((2, 1), (1, 1, 1)):
        C3.comp_reduced(4, alpha, 1, check=True)
    for alpha in ((2,), (1, 1)):
        e = C3.comp_reduced(4, alpha, 1, check=True)
        old = sp.sympify(pickle.load(open(f"data/refined_kauffman/daha/Rcomp_N4_m1_{''.join(map(str, alpha))}.pkl", "rb"))["R"],
                         locals={"r": rr, "s": ss})
        assert sp.expand(e - old) == 0


def test_composite3_matches_kn_data():
    import rk_k5_trefoil_2box as T
    for N, alpha, sh in ((2, (3,), "4"), (4, (2, 1), "321"), (4, (1, 1, 1), "222")):
        e = C3.comp_reduced(N, alpha, 1)
        W = (ss / rr)**sum(int(c) for c in sh) * T.Rdata(sh, N)
        assert sp.cancel(e / W) == 1, (N, alpha, sh)


def test_fixedN_weights_dominance():
    import rk_c_4box_cover_fixedN as X
    ws = X.weights(3, 0, 2, 2, (2, 0, -2))
    assert set(ws) == {(2, 0, -2), (1, 1, -2), (2, -1, -1), (1, 0, -1), (0, 0, 0)}


def _load3(m):
    from rcsknots.refined_kauffman.flint_rf3 import from_sympy
    d = pickle.load(open(f"data/refined_kauffman/daha/largeN_3box_m{m}.pkl", "rb"))
    r_, s_, u_ = sp.symbols("r s u", positive=True)
    return {eval(k): from_sympy(sp.sympify(v, locals={"r": r_, "s": s_, "u": u_})) for k, v in d["g"].items()}


@pytest.mark.results
def test_T25_3box_integer():
    from rcsknots.characters import chi
    from rcsknots.partitions import partitions
    from rcsknots.refined_kauffman.flint_rf3 import conv
    g = _load3(2)
    for R_ in partitions(3):
        gR = sum((g[rho] * conv(chi(R_, rho)) for rho in partitions(3)), conv(0))
        assert gR.is_integer_laurent(), R_


@pytest.mark.results
def test_trefoil_4box_poles_and_lower_degrees():
    import rk_c_4box_largeN_g as G4
    from rcsknots.characters import chi
    from rcsknots.partitions import partitions
    from rcsknots.refined_kauffman.flint_rf3 import conv, S
    d = pickle.load(open("data/refined_kauffman/daha/largeN_4box_m1.pkl", "rb"))["g"]
    g = {eval(k): G4.deser(v) for k, v in d.items()}
    g3 = _load3(1)
    for k, v in g3.items():
        assert (g[k] - v).is_zero(), k
    t = S * S
    for rho in partitions(4):
        assert g[rho].is_laurent() == (rho not in ((4,), (2, 2))), rho      # (1 + t) は冪和がすべて偶数の方向だけ
    for R_ in partitions(4):
        gR = sum((g[rho] * conv(chi(R_, rho)) for rho in partitions(4)), conv(0))
        assert (2 * (1 + t) * gR).is_integer_laurent(), R_          # g = ½E/(1+t) + C（2 箱と同じく分母は 2 × 単項式 × (1 + t)）
