"""別研究 [2,1] の X2：ダブルツイスト結び目の refined 不変量のチャネル分解と、[2,1] の未知数・条件の数え上げ（機能テスト）。"""
import itertools
import sys
from pathlib import Path

import pytest
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from rcsknots.double_twist_family import A, Q, T, M_rect, U_row, eig, tau, value  # noqa: E402
from rcsknots.double_twist_refined import a, poincare_pre, sig  # noqa: E402


def zero(x):
    return sp.expand(sp.cancel(sp.together(x))) == 0


@pytest.mark.parametrize("rs", [(1, 1), (2, 1), (1, 2)])
def test_channel_decomposition_reproduces_kntz(rs):
    """KNTZ の置換の前の式で σ = 1, 𝐚 = A(T/Q)^{1/2} とした値 = KN の rCS（X1 で確認）が、チャネル分解 N Σ M e^m e^n と一致。"""
    M = M_rect(*rs)
    for mn in ((1, 1), (-1, 1), (2, 1), (-2, -1), (3, -2)):
        direct = poincare_pre(*rs, *mn).subs({sig: 1, a: A * sp.sqrt(T / Q)}, simultaneous=True)
        assert zero(value(M, *mn, rs[0] * rs[1]) - direct), mn


@pytest.mark.parametrize("lam", [(1,), (2,), (1, 1), (2, 1)])
def test_tau_is_single_exponential_sum(lam):
    """τ_λ(m) = Σ_μ U_{λμ} e_μ^m が負の m も含めて成り立つ（refined の 𝓑 の正則な行）。当てはめに使わない m で確かめる。"""
    sp.core.cache.clear_cache()      # 前のテストで膨らんだ sympy のキャッシュで数十倍遅くなるため
    Urow = U_row(lam)
    for m in (-3, -2, len(Urow) - 1):
        assert zero(sum(u * eig(mu) ** m for mu, u in Urow.items()) - tau(lam, m)), m


@pytest.mark.parametrize("X,XT", [((2,), (1, 1)), ((1,), (1,)), ((2, 1), (2, 1))])
def test_tau_transposition_and_mirror(X, XT):
    """τ_X(A, 1/T, 1/Q) = τ_{X^T}(A, Q, T)、τ_X(−m)(x) = τ_X(m)(1/x)/τ_X(1)(1/x)。"""
    inv = {A: 1 / A, Q: 1 / Q, T: 1 / T}
    for m in (-2, 1, 2):
        assert zero(tau(X, m).subs({Q: 1 / T, T: 1 / Q}, simultaneous=True) - tau(XT, m))
    for m in (1, 2):
        assert zero(tau(X, -m) - tau(X, m).subs(inv, simultaneous=True) / tau(X, 1).subs(inv, simultaneous=True))


def test_egp_and_u2_are_consistent_and_fix_channels():
    """EGP（Q = 1）と U(2)（A = T²）は、どちらも 25 本の式で 5 つの z_X を一意に決める（過剰決定だが矛盾しない）。"""
    from x2_specializations import solve_at
    for sub in ({Q: 1}, {A: T**2}):
        sol = solve_at(sub, "test")
        assert len(sol) == 1 and len(sol[0]) == 5


def test_channel_21_is_not_a_product_of_refined_differentials():
    """Z̃_{[2,1]} を {A Q^α T^β} の 6 個の積（t = q で D₃D₂D₁D₋₁D₋₂D₋₃、q = 1 で EGP の値）とすると、転置で閉じるものは無い。"""
    betas, ns = [0, 0, -1, -1, -2, -3], [3, 2, 1, -1, -2, -3]
    cands = {tuple(sorted((n - b, b) for n, b in zip(ns, p))) for p in itertools.permutations(betas)}
    assert len(cands) == 180
    assert not [p for p in cands if tuple(sorted((-b, -a_) for a_, b in p)) in cands]


def test_channel_1_is_not_unique():
    """Z̃_{[1]} は窓 K = 3 で初めて解をもち、解の空間の次元は 8（一意でない）。"""
    from x2_window import solve_channel, specialization_data
    egp, u2 = specialization_data()
    assert solve_channel((1,), 2, egp, u2) is None
    dim, _ = solve_channel((1,), 3, egp, u2)
    assert dim == 8


def test_explicit_freedom_preserves_all_conditions():
    """f = (Q−T)(Q−1)(T−1)(A²−T⁴)(A²Q⁴−1)/(A²Q³T³) と g = Q/T − T/Q で δ = f·g は、5 つの特殊化で消え、転置で不変、反転で不変。"""
    f = (Q - T) * (Q - 1) * (T - 1) * (A**2 - T**4) * (A**2 * Q**4 - 1) / (A**2 * Q**3 * T**3)
    d = f * (Q / T - T / Q)
    for sub in ({Q: T}, {Q: 1}, {T: 1}, {A: T**2}, {A: 1 / Q**2}):
        assert zero(d.subs(sub))
    assert zero(d.subs({Q: 1 / T, T: 1 / Q}, simultaneous=True) - d)
    assert zero(d.subs({A: 1 / A, Q: 1 / Q, T: 1 / T}, simultaneous=True) - d)
    assert not zero(d)
