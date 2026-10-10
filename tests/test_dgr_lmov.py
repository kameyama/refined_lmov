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


# ---- 段階 D2：d_{1|1} の普遍式（scripts/dgr_d11_universal.py） ----

def test_d11_universal_formula_matches_twist_knots():
    """hook 性質から作る \\hat f_ρ|_{d11}（A_K だけで決まる）が実際の値と全係数で一致する。"""
    from dgr_d11_universal import universal
    from v4_closed import hat_from_inputs, inputs
    for p in (-2, 2):
        P1, P2 = inputs(p)
        Sigma = -next(iter(specialize(P1, 1, 0)))[1]
        real = {rho: specialize(b, 1, 1) for rho, b in hat_from_inputs(P1, P2).items()}
        assert universal(specialize(P1, 1, 1), Sigma) == real


def test_d11_top_genus_drop_rule():
    """genus 2、Σ = 2：最高 genus の係数 α_2 が 𝐭^{−4} なら [2] だけ、1 なら [1,1] だけが一般の 2g−1 = 3 から 1 落ちる。"""
    from dgr_d11_universal import from_genus, universal_qdeg
    low = [{-1: 1, 0: 2}, {1: -1, 0: 1}]
    top = lambda al: {r: v[0] for r, v in universal_qdeg(from_genus(low + [al]), 2).items()}  # noqa: E731
    assert top({-4: 1}) == {(2,): 2, (1, 1): 3}
    assert top({0: 1}) == {(2,): 3, (1, 1): 2}
    assert top({-4: 2}) == {(2,): 3, (1, 1): 3}


# ---- H1：(S_{1|0}), (S_{0|1}) の導出の補題（scripts/dgr_h1_lemma.py） ----

def _gn():
    from fractions import Fraction as Fr
    from rcsknots.gn import GN
    q, t = Fr(3, 7), Fr(5, 11)
    return q, t, GN(q, t, t, lambda x: x), GN(q, t, 1 / q, lambda x: x)


def test_h1_kernels_of_one_variable_evaluations():
    """H̃_μ[p_k ↦ 1 − u^k] = 0 ⇔ ℓ(μ) ≥ 2（u = t）、μ₁ ≥ 2（u = 1/q）、μ ⊇ (2,2)（u = t/q）。"""
    from dgr_h1_lemma import check_kernel
    q, t, Gt, _ = _gn()
    assert check_kernel(Gt, t, lambda mu: len(mu) >= 2, 4)
    assert check_kernel(Gt, 1 / q, lambda mu: mu[0] >= 2, 4)
    assert check_kernel(Gt, t / q, lambda mu: len(mu) >= 2 and mu[1] >= 2, 4)


def test_h1_quotient_is_quantum_torus():
    """P̃_{kn,km} は核を保ち、商に x^j ↦ γ^k q^{(k²−k)nm/2 + kmj} x^{j+kn} で作用する（γ = (−t)^{n−1}q^{(n−1)(m−1)/2}）。"""
    from fractions import Fraction as Fr
    from dgr_h1_lemma import scalars
    q, t, Gt, _ = _gn()
    for n, m in ((1, 1), (2, 1)):
        gam = (-t) ** (n - 1) * q ** Fr((n - 1) * (m - 1), 2)
        assert scalars(Gt, t, n, m, 1, 2) == [gam * q ** (m * j) for j in range(3)]


def test_h1_row_and_column_values_trefoil():
    """reduced rCS_{[r]}|_{a=t} = γ^r q^{nm r(r−1)/2}、二次の部分は (N2') の T_{[r]}^{nm}。列も a = 1/q で同様に c^r 型。"""
    from fractions import Fraction as Fr
    q, t, Gt, Gq = _gn()
    n, m = 2, 3
    gam = (-t) * q
    for r in (1, 2):
        assert Gt.torus((r,), n, m) / Gt.torus((r,), 1, 0) == gam ** r * q ** Fr(n * m * r * (r - 1), 2)
    v = [Gq.torus((1,) * r, n, m) / Gq.torus((1,) * r, 1, 0) / t ** Fr(-n * m * r * (r - 1), 2) for r in (1, 2)]
    assert v[1] == v[0] ** 2


# ---- D8・D9：GT1/GT2 の逆解き（src/rcsknots/large_n_inverse.py） ----

def test_inverse_round_trip_trefoil():
    """\\hat f（data/lmov）から GT1・GT2 を逆に解いた reduced の rCS が、既存の rCS に点で一致する（|λ| ≤ 3）。"""
    from fractions import Fraction as Fr
    from rcsknots.large_n_inverse import ev_poly, rcs_from_hat, reduced
    F = load_lmov("T2_3")
    hat_ = {rho: F[("hat", rho)] for rho in RHOS}
    R = load_reduced("T2_3", 3)
    P = (Fr(2, 3), Fr(5, 7), Fr(3, 11))
    for which in ("q", "t"):
        red = reduced(rcs_from_hat(hat_, *P, which), *P)
        assert all(red[lam] == ev_poly(R[lam], *P) for lam in red)


def test_gt1_gt2_agree_for_any_hat_f():
    """\\hat f を摂動しても GT1 と GT2 の逆解きは一致する（2 つの式は同値で、\\hat f への制約にならない）。"""
    from fractions import Fraction as Fr
    from rcsknots.large_n_inverse import rcs_from_hat
    F = load_lmov("T2_3")
    hat_ = {rho: dict(F[("hat", rho)]) for rho in RHOS}
    hat_[(2,)][(2, 1, 1)] = hat_[(2,)].get((2, 1, 1), 0) + 1
    hat_[(2, 1)][(3, 1, 1)] = hat_[(2, 1)].get((3, 1, 1), 0) + 1
    P = (Fr(2, 3), Fr(5, 7), Fr(3, 11))
    assert rcs_from_hat(hat_, *P, "q") == rcs_from_hat(hat_, *P, "t")


def test_d7_true_answer_is_feasible_and_not_unique():
    """D7 の整数計画：三葉結び目の正解の rCS_{[2]} は全ての条件を満たすが、総数が同じ別解もある。"""
    import numpy as np
    from dgr_d7_construct import build, parity_class, solve, solution, unrefined
    from v4_closed import inputs
    P1, P2 = inputs(1)
    true = {k: int(v) for k, v in P2.items() if v}
    ts = [k[2] for k in P2]
    lin, mons = build(P1, unrefined(P2), (min(ts) - 2, max(ts) + 3), "b", (-1, -1), genus=True, egp=True, parity=parity_class(P1))
    obj = np.array([1.0 if n[0] == "h" else 0.0 for n in lin.col])
    res = solve(lin, objective=obj)
    assert round(res.fun) == 96
    assert solution(lin, res, mons) != true                     # 最小解は正解と限らない（𝐭 のずらしなど）
    for m in mons:
        lin.eq({("x", m): 1}, true.get(m, 0))
    assert solve(lin).x is not None                              # 正解は実行可能
