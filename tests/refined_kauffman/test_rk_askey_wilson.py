"""refined Kauffman K2'(a)：probe 側の Askey–Wilson 多項式（基底 t、パラメータ (1, −1, q^{1/2}, −q^{1/2})）。"""
import sys
from fractions import Fraction as F
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from rcsknots.refined_kauffman.askey_wilson import askey_wilson_monic, top_coeff_closed  # noqa: E402
from rcsknots.refined_kauffman.koornwinder import koornwinder  # noqa: E402

q = sp.Rational(9, 49); sq = sp.Rational(3, 7); t = sp.Rational(5, 11)


def test_recurrence_matches_koornwinder_solver():
    """三項漸化式の monic 化が、q-差分作用素の固有関数として解いた 1 変数 Koornwinder と一致する（n ≤ 3）。"""
    y = sp.Symbol("y")
    aw = askey_wilson_monic(3, t, sp.Integer(1), sp.Integer(-1), sq, -sq)
    for n in range(1, 4):
        K = sp.expand(koornwinder((n,), 1, t, q, (1, -1, sq, -sq), (y,)))
        K2 = sum(v * y**j for j, v in aw[n].items())
        assert sp.expand(K - K2) == 0


def test_top_coefficients_closed_form():
    """y^{n−2j} の係数の閉じた形（U_n = t^n）が漸化式の値と一致する（j ≤ 3、n ≤ 12）。"""
    qf, sqf, tf = F(9, 49), F(3, 7), F(5, 11)
    aw = askey_wilson_monic(12, tf, F(1), F(-1), sqf, -sqf)
    for j in range(0, 4):
        for n in range(2 * j, 13):
            assert aw[n].get(n - 2 * j, F(0)) == top_coeff_closed(j, tf**n, qf, tf)


def test_top_coefficients_limit():
    """U_n → 0 で (q y²; t²)_∞/(t y²; t²)_∞ の係数（q-二項定理）に一致する。"""
    qf, tf = F(9, 49), F(5, 11)
    for j in range(0, 4):
        lim = tf**j
        for i in range(j):
            lim *= (1 - qf / tf * tf**(2 * i)) / (1 - tf**(2 * i + 2))
        assert top_coeff_closed(j, F(0), qf, tf) == lim
