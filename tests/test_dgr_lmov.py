"""DGR / GGS の微分 d_{n|m} の類似（src/rcsknots/dgr.py、docs/dgr_lmov.md）の機能テスト。小さい結び目だけで数秒。"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from framing_check import load_reduced  # noqa: E402
from framing_twist import hat_framed  # noqa: E402
from v4_closed import mir  # noqa: E402
from rcsknots.bps import g_expand, to_bold  # noqa: E402
from rcsknots.dgr import (divides, d11_top_genus, hook_rule, judge_quotient, quotient, sigma,  # noqa: E402
                          specialize, u)
from rcsknots.extract import A_SHIFT  # noqa: E402
from rcsknots.lmov import RHOS, load_lmov  # noqa: E402

UNKNOT = {(0, 0, 0): 1}
TREFOIL_P1 = {(2, -2, 0): 1, (2, 2, 2): 1, (4, 0, 3): 1}       # DGR の三葉結び目


def hat(K):
    F = load_lmov(K)
    return {rho: to_bold(F[("hat", rho)], A_SHIFT["nontorus"][sum(rho)]) for rho in RHOS}


def trefoil_p2():
    from t3_pretzel import tc_graded
    return tc_graded(1, 1, 2)                                   # ツイスト結び目 K_1 = 三葉結び目の [2]（GNSSS、t_c 次数）


def test_sigma_swaps_d_nm_and_d_mn():
    for n in range(3):
        for m in range(3):
            assert sigma({u(n, m): 1}) == {u(m, n): 1}


@pytest.mark.parametrize("K", ["T2_3", "4_1"])
def test_universal_factor_and_positive_quotient(K):
    """|ρ| ≥ 2 は d_{0|0}, d_{1|0}, d_{0|1} の因子で割り切れ、商は spin・genus 表示とも非負（[1] は d_{0|0} だけ）。"""
    for rho, b in hat(K).items():
        n = sum(rho)
        assert divides(b, 0, 0)
        assert divides(b, 1, 0) == (n > 1) and divides(b, 0, 1) == (n > 1)
        sn, _, gn, _, same = judge_quotient(b, n)
        assert (sn, gn, same) == (0, 0, True)


def test_quotient_of_trefoil_rho2_explicit():
    """三葉結び目：R_{[2]} = −(𝐚²𝐪⁻² + 𝐚²𝐪²𝐭²)(1 + 𝐚²𝐭³)。"""
    R = quotient(hat("T2_3")[(2,)], 2)
    assert R == {(2, -2, 0): -1, (2, 2, 2): -1, (4, -2, 3): -1, (4, 2, 5): -1}


def test_framed_congruence_with_unknot():
    """framing x では \\hat f^{(x)}(K) ≡ c^{|ρ|} \\hat f^{(x)}(○) mod (1 + u_{1|0}), (1 + u_{0|1})、c = [1] の特殊化。"""
    P1, P2 = TREFOIL_P1, trefoil_p2()
    (c_key, c), = specialize(P1, 1, 0).items()
    assert c_key[0] == 0
    for x in (-2, -1, 1, 3):
        H, U = hat_framed(P1, P2, x), hat_framed(UNKNOT, UNKNOT, x)
        for rho in H:
            for nm in ((1, 0), (0, 1)):
                lhs = specialize(H[rho], *nm)
                rhs = {(k[0], k[1] + 2 * c_key[1]): c ** 2 * v for k, v in specialize(U[rho], *nm).items()}
                assert lhs == rhs


def test_framed_unknot_rho2_is_a_single_unknot_state():
    """framing −1 の自明な結び目の \\hat f_{[2]} は −(1 + 𝐚²𝐭) × 単項式（d_{1|0}, d_{0|1} の因子を持たない）。"""
    b = hat_framed(UNKNOT, UNKNOT, -1)[(2,)]
    assert len(b) == 2 and divides(b, 0, 0) and not divides(b, 1, 0)


@pytest.mark.parametrize("K", ["T2_3", "T2_5", "4_1"])
def test_hook_rule_at_d11(K):
    """hook λ で rCS_λ|_{d11} = c · A_K(𝐪^{|λ|}𝐭^{arm}, 𝐭)。"""
    R = load_reduced(K, 3)
    bold = {lam: to_bold(P, 0 if all(k[0] % 2 == 0 for k in P) else 1) for lam, P in R.items()}
    A1 = specialize(bold[(1,)], 1, 1)
    for lam, b in bold.items():
        assert hook_rule(A1, specialize(b, 1, 1), lam) is not None, lam


def test_hook_rule_fails_for_wrong_shift():
    R = load_reduced("T2_5", 2)
    A1 = specialize(to_bold(R[(1,)], 0), 1, 1)
    assert hook_rule(A1, specialize(to_bold(R[(2,)], 0), 1, 1), (1, 1)) is None   # [2] を arm 0 として読むと合わない


def test_d11_top_genus_linear_in_seifert_genus():
    """T(2,5)（g = 2）：d_{1|1} で残る最大 genus は |ρ|g − c_ρ（c = 3, 2, 4, 3, 3）。特殊化の前は |ρ|² 程度で大きい。"""
    c = {(2,): 3, (1, 1): 2, (3,): 4, (2, 1): 3, (1, 1, 1): 3}
    for rho, b in hat("T2_5").items():
        R = quotient(b, sum(rho))
        if rho != (1,):
            assert d11_top_genus(R) == 2 * sum(rho) - c[rho]
            assert max(k[0] for k in g_expand(R)) > d11_top_genus(R)


def test_sigma_orbit_factors_positive_in_genus_rep():
    """(1 + u_{n|m})(1 + u_{m|n}) は genus 表示で非負係数、最大 genus は |n − m|。"""
    for n in range(3):
        for m in range(n, 4):
            e = u(n, m)
            s = u(m, n)
            prod = {(0, 0, 0): 1, e: 1} if e == s else {(0, 0, 0): 1, e: 1, s: 1, tuple(e[i] + s[i] for i in range(3)): 1}
            N = g_expand(prod)
            assert min(N.values()) > 0 and max(k[0] for k in N) == m - n


def test_unrefined_double_zero_at_a1():
    """t = q（𝐭 = −1）では d_{0|0} と d_{1|1} が重なり、|ρ| ≥ 2 の \\hat f は (1 − 𝐚²)² で割り切れる。"""
    for rho, b in hat("T2_3").items():
        un = {}
        for (ea, eq, et), cc in b.items():
            un[(ea, eq)] = un.get((ea, eq), 0) + cc * (-1) ** (et % 2)
        # (1 − 𝐚²) で割れる回数 = 𝐚² = 1 での零点の位数：導関数で数える
        import sympy as sp
        x, y = sp.symbols("x y")
        P = sp.expand(sum(v * x ** (k[0] // 2) * y ** k[1] for k, v in un.items()) * x ** 20 * y ** 40)
        order = 0
        while sp.expand(P.subs(x, 1)) == 0:
            P = sp.diff(P, x)
            order += 1
        assert order == (1 if rho == (1,) else 2)
