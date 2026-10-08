"""別研究 [2,1] の X1：KNTZ 1902.02275 の refined 𝐁（補間 Macdonald）の実装 rcsknots.double_twist_refined の機能テスト。

照合：KNTZ の補助ファイル（Wolfram で実行した値）、GNSSS 1512.07883 の式（第 2 期の入力）、t = q で X0 の式、三葉結び目の EHA、
第 2 期 T3 のツイスト結び目の入力、Morozov 1903.00259 の 𝓑^{[2,1]} の正則な行（t = q）。
"""
import sys
from pathlib import Path

import pytest
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from rcsknots import double_twist as U  # noqa: E402
from rcsknots.double_twist_refined import (A, Bref, Dref, Fref, interp_P, poincare, q, t, tc, tr)  # noqa: E402


def zero(x):
    return sp.expand(sp.cancel(sp.together(x))) == 0


def test_interpolation_vanishing():
    """補間 Macdonald 多項式：P*_μ(q^λ) = 0 ⇔ μ ⊄ λ（|μ| ≤ |λ|、Okounkov の性質）。"""
    shapes = [(1,), (2,), (1, 1), (2, 1), (3,)]
    for lam in shapes:
        for mu in shapes:
            if sum(mu) > sum(lam):
                continue
            n = max(len(lam), len(mu))
            v = sp.cancel(interp_P(mu, [q ** x for x in list(lam) + [0] * (n - len(lam))]))
            inside = len(mu) <= len(lam) and all(m <= l for m, l in zip(mu, lam))
            assert (v != 0) == inside, (lam, mu)


@pytest.mark.parametrize("lam", [(1,), (2,), (1, 1), (2, 1), (2, 2)])
def test_B_empty_F_laurent_and_t_eq_q(lam):
    """KNTZ :734-744：𝐁_{λ,∅} = 𝐃_λ、𝐅 は s によらない Laurent 多項式。t = q で非 refined の F に戻る。"""
    for s in (len(lam), len(lam) + 1):
        assert zero(Bref(lam, (), s) - Dref(s, lam))
    for m in (2, 3):
        f = Fref(lam, m)
        assert zero(f - Fref(lam, m, len(lam) + 1))
        assert len(sp.Poly(sp.fraction(sp.cancel(f))[1], A, q, t).terms()) == 1
        assert zero(f.subs(t, q) - U.F_kntz(lam, m).subs({U.A: A, U.q: q}))


def test_kntz_package_values():
    """KNTZ のパッケージ（scripts/x1_kntz_package.wls）の 𝐅 と、小さな 𝒫 が全係数一致。全件は results。"""
    from x1_check import load_package
    n = 0
    keep = lambda kind, args: kind == "F" or (kind == "P" and args[0] * args[1] <= 2 and abs(args[2]) <= 2 and args[3] == 1)  # noqa: E731
    for kind, args, val in load_package(keep):
        ours = Fref(tuple(args[0]), args[1], args[2]) if kind == "F" else poincare(*args)
        assert zero(ours - val), (kind, args)
        n += 1
    assert n >= 20


@pytest.mark.parametrize("mn", [(1, 1), (-1, 1), (2, 1), (2, 2)])
def test_gnsss_symmetric_and_antisymmetric(mn):
    """[r]：GNSSS の四重次数の式（t3_pretzel.quadruple）と一致。[1^r]：GGS のミラー対称性（t_r ↔ t_c）を介して一致。"""
    from t3_pretzel import quadruple
    for r in (1, 2):
        g = sp.expand(sp.cancel(quadruple(*mn, r)))
        assert sp.expand(poincare(r, 1, *mn) - g) == 0
        assert sp.expand(poincare(1, r, *mn) - g.subs({tr: tc, tc: tr}, simultaneous=True)) == 0


@pytest.mark.parametrize("rs", [(2, 1), (1, 2), (2, 2)])
def test_t_eq_q_is_x0_formula(rs):
    """置換の前の式で t = q, σ = 1, a = A とすると、KNTZ の非 refined の式（X0 でケーブルと一致したもの）に戻る。"""
    from x1_check import tq_limit
    for mn in ((1, 1), (-1, 1), (2, 1)):
        assert tq_limit(*rs, *mn)


def test_trefoil_matches_eha():
    """K_{1,1} = T(2,3) の [1], [2], [1,1], [2,2]：t_c 次数の三重次数（GGS の 𝐪 = R𝐐 − 𝐭_r + 𝐭_c）を KN の rCS に直すと EHA の値（(N2') で正規化）に一致。"""
    from x1_check import y6_trefoil
    assert y6_trefoil(lambda s: None) == 0


def test_twist_knots_match_phase2_inputs():
    """K_{±1,1}, K_{±2,1} の [1], [2], [1,1] が第 2 期 T3 の入力（NRZS の式＋EGP で決めた δ'）と一致。"""
    from x1_check import y6_twist
    assert y6_twist([-2, -1, 1, 2], lambda s: None) == 0


@pytest.mark.parametrize("m", [-3, -2, -1, 1, 2, 3])
def test_morozov_regular_rows_are_kntz_twist(m):
    """t = q：1903.00259 の 𝓑^{[2,1]} の正則な行の F^{(m)}_λ（λ = [1], [1,1], [2], [2,1]）は KNTZ の Tw^{(m)}_λ に等しい。"""
    F = U.F21_1903(m)
    for i, lam in {1: (1,), 2: (1, 1), 3: (2,), 4: (2, 1)}.items():
        assert zero(F[i] - U.Tw(lam, m)), lam


@pytest.mark.results
def test_kntz_package_all():
    from x1_check import check_package
    assert check_package(lambda s: None) == 0


@pytest.mark.results
def test_twist_knots_match_phase2_inputs_wide():
    from x1_check import y6_twist
    assert y6_twist([p for p in range(-6, 7) if p], lambda s: None) == 0
