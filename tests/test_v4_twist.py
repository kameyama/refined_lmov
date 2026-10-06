"""第 4 期 V4：ツイスト結び目 K_p の \\hat f_[2], \\hat f_[1,1] の閉じた形（ニュートン形）と正値性の証明に使う事実。

入力は GNSSS の予想式（scripts/t3_pretzel.py の tc_graded）。閉じた形の c_k は |p| = 1..5 の直接計算から作り、
照合は c_k の作成に使っていない |p| = 6, 7, 9 で全係数を厳密に比較する。
"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from t3_pretzel import tc_graded  # noqa: E402
from v4_closed import RHOS, L, closed_eval, hat_from_inputs, mono, newton_coeffs, roots_for  # noqa: E402
from rcsknots.bps import g_expand  # noqa: E402


@pytest.fixture(scope="module")
def direct():
    out = {}
    for sign in (1, -1):
        for n in (1, 2, 3, 4, 5, 6, 7, 9):
            p = sign * n
            out[p] = hat_from_inputs(tc_graded(1, p, 1), tc_graded(1, p, 2))
    return out


def coeffs(direct, sign):
    return {rho: newton_coeffs([direct[sign * n][rho] for n in range(1, 6)], roots_for(sign)) for rho in RHOS}


@pytest.mark.parametrize("sign", [1, -1])
def test_closed_form_matches_direct(direct, sign):
    CS = coeffs(direct, sign)
    for n in (6, 7, 9):                                    # c_k の作成に使っていない p
        for rho in RHOS:
            assert closed_eval(CS[rho], sign, sign * n) == direct[sign * n][rho]


@pytest.mark.parametrize("sign", [1, -1])
def test_newton_coefficients_nonpositive(direct, sign):
    """c_k の係数はすべて ≤ 0（全体の符号 −1 で非負）、𝐪 の指数は偶数。h_n は非負係数なので全ての p で \\hat f ≤ 0 が従う。"""
    for cs in coeffs(direct, sign).values():
        assert all(v < 0 for c in cs for v in c.values())
        assert all(e[1] % 2 == 0 for c in cs for e in c)


@pytest.mark.parametrize("sign", [1, -1])
def test_mirror_symmetry(direct, sign):
    """σ：(e_𝐚, e_𝐪, e_𝐭) ↦ (e_𝐚, −e_𝐪, e_𝐭 − e_𝐪)（u = 𝐪𝐭^{1/2} ↦ u^{-1}）で不変。根の集合が σ で閉じるので |p| = 1..5 で全ての p に従う。"""
    for n in range(1, 6):
        for rho in RHOS:
            P = direct[sign * n][rho]
            assert P == {(a, -q, t - q): c for (a, q, t), c in P.items()}
            g_expand(P)                                    # G 展開が閉じる


def test_negative_p_is_t2_times_continuation(direct):
    """p < 0 の \\hat f は、p > 0 の指数多項式を負の p に延長したものの 𝐭² 倍（p = 0 での延長は 0）。"""
    ms = roots_for(1)
    char = [L.from_dict({(0, 0, 0): 1})]                   # Π (x − m_j) の係数（低次から）
    for m in ms:
        new = [L.from_dict({}) for _ in range(len(char) + 1)]
        for i, c in enumerate(char):
            new[i + 1] = new[i + 1] + c
            new[i] = new[i] - mono(m) * c
        char = new
    for rho in RHOS:
        s = {n: L.from_dict(direct[n][rho]) for n in range(1, 6)}
        for n in range(0, -4, -1):                         # Σ_i κ_i s_{n+i} = 0 を後ろ向きに解く
            acc = L.from_dict({})
            for i in range(1, 6):
                acc = acc + char[i] * s[n + i]
            s[n] = (L.from_dict({}) - acc) / char[0]
        assert s[0].is_zero()
        for n in (-1, -2, -3):
            assert (mono((0, 0, 2)) * s[n]).to_dict() == direct[n][rho]
