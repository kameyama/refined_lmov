"""第 4 期 V4 の続き：pretzel K_{s,t} の \\hat f_[2], \\hat f_[1,1] の 2 変数の閉じた形と正値性の証明に使う事実。

入力は GNSSS の予想式（scripts/framed_knot/t3_pretzel.py、浮動小数点の混入を修正した版）。閉じた形の c_kl は各領域の |s|, |t| = 1..5 の 25 点から作り、
照合は使っていない点で全係数を厳密に比較する。
"""
import sys
from pathlib import Path

import pytest
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/framed_knot"))

import t3_pretzel  # noqa: E402
from t3_pretzel import tc_graded  # noqa: E402
from v4_closed import RHOS, hat_from_inputs  # noqa: E402
from v4_pretzel import ORDER, PROOF_ORDERS, REGIONS, closed_eval, closed_form, euler_q1, sigma  # noqa: E402
from rcsknots.bps import g_expand  # noqa: E402

_hat = {}


def hat(s, t):
    k = (min(s, t), max(s, t))
    if k not in _hat:
        _hat[k] = hat_from_inputs(tc_graded(*k, 1), tc_graded(*k, 2))
    return _hat[k]


def coeffs(region, orders=None):
    return closed_form(region, orders, hat=hat)


def test_no_float_and_special_polynomial():
    """修正後の式は厳密（Float を含まない）。旧版で「不備」とされた点でも P_r(𝐪=1, 𝐭=−1) = (1 − (1−𝐚^{2s})(1−𝐚^{2t}))^r、次元は二乗。"""
    a = sp.Symbol("a")
    assert not t3_pretzel.quadruple(2, 5, 2).atoms(sp.Float)
    for s, t in [(2, 5), (-3, 3), (-4, 4), (-5, -2)]:
        P1, P2 = tc_graded(s, t, 1), tc_graded(s, t, 2)
        assert sum(P2.values()) == sum(P1.values()) ** 2
        for r, P in ((1, P1), (2, P2)):
            e = sp.Poly(sp.expand((1 - (1 - a ** (2 * s)) * (1 - a ** (2 * t))) ** r * a ** 100), a)
            assert euler_q1(P) == {k - 100: int(c) for (k,), c in e.terms()}


def test_mirror_line():
    """K_{−1,t} の入力は K_{1,−t} の入力の鏡像（全変数の反転）。"""
    inv = lambda P: {tuple(-e for e in k): c for k, c in P.items()}  # noqa: E731
    for t in (-3, 2, 4):
        for r in (1, 2):
            assert tc_graded(-1, t, r) == inv(tc_graded(1, -t, r))


@pytest.mark.parametrize("region,pt", [("pp", (6, 7)), ("nn", (-7, -6)), ("np", (-6, 7)), ("np", (-7, 2))])
def test_closed_form_matches_direct(region, pt):
    """c_kl の作成に使っていない点で、閉じた形が直接計算と全係数一致（証明に使う並びでも、既定の並びでも）。"""
    for orders in (None, PROOF_ORDERS if region != "np" else None):
        CS = coeffs(region, orders)
        for rho in RHOS:
            pair = (ORDER, ORDER) if orders is None else orders[(region, rho)]
            assert closed_eval(CS[rho], region, *pt, order_pair=pair) == hat(*pt)[rho]


@pytest.mark.parametrize("region", ["pp", "nn"])
def test_newton_coefficients_uniform_sign(region):
    """PROOF_ORDERS で c_kl の係数は全て < 0、𝐪 の指数は偶数。h は非負係数なので領域内の全ての (s, t) で \\hat f ≤ 0。"""
    CS = coeffs(region, PROOF_ORDERS)
    for rho in RHOS:
        assert all(v < 0 for row in CS[rho] for c in row for v in c.values())
        assert all(e[1] % 2 == 0 for row in CS[rho] for c in row for e in c)


def test_np_not_uniform_with_default_order():
    """np では既定の並びの c_kl に正の係数が現れる（2 変数のニュートン形では示せない。直線ごとに示す）。"""
    CS = coeffs("np")
    assert any(v > 0 for row in CS[(2,)] for c in row for v in c.values())


@pytest.mark.parametrize("region", list(REGIONS))
def test_sigma_invariance_on_grid(region):
    """25 点で σ 不変・G 展開が閉じる。根の集合は σ で閉じるので、各変数 5 点で全ての (s, t) に従う。"""
    ss, st = REGIONS[region]
    for i in range(1, 6):
        for j in range(1, 6):
            for rho in RHOS:
                P = hat(ss * i, st * j)[rho]
                assert P == sigma(P)
                g_expand(P)


def test_np_lines():
    """np の直線：K_{−3,t}（t ≥ 1）は [2] で基点 t = 3、[1,1] で基点 t = 2 から c_l < 0（並び (1, X, X², Y, Ȳ)）。基点より前は個別に負。"""
    from v4_pretzel_np import line_proof
    CS = coeffs("np")
    for rho, base in (((2,), 3), ((1, 1), 2)):
        b, pm, vals = line_proof(lambda n: closed_eval(CS[rho], "np", -3, n), 1, 6)
        assert (b, pm) == (base, ("1", "X", "X2", "Y", "Yb"))
        assert all(v <= 0 for n in range(1, b) for v in vals[n].values())
