"""open_problems D1：正規化 (N1)(N2) で決めた単項式 M_□ の閉じた式。

M_□ = (−q/t)^{b−1}·(a (q/t)^{1/2})^{δ}。b はブレイドの本数（Π p_i）、δ は反復ごとに δ ↦ p δ + (p−1)(q−1)/2（トーラス結び目 T(m,n) は δ = (m−1)(n−1)/2）。
正のブレイド表示では δ = (w − b + 1)/2（w は writhe）で、Markov 跡による正規化の refined 版と読める。保存済みの全トーラス結び目と全ケーブルで照合する。
"""
import json
from fractions import Fraction
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def steps_of(d):
    if "torus" in d:
        return [tuple(d["torus"])]
    if "steps" in d:
        return [tuple(s) for s in d["steps"]]
    return [tuple(d["cable"]["companion"]), tuple(d["cable"]["pattern"])]


def predicted(steps):
    b, delta = 1, Fraction(0)
    for p, q in steps:
        b, delta = b * p, p * delta + Fraction((p - 1) * (q - 1), 2)
    # M = sign · a^δ · (q/t)^{δ/2 + b − 1}
    return (-1) ** (b - 1), delta, delta / 2 + b - 1


FILES = sorted(ROOT.glob("data/torus/*_1.json")) + sorted(ROOT.glob("data/cables/K*_1.json"))


@pytest.mark.parametrize("path", FILES, ids=[f.stem for f in FILES])
def test_box_monomial_formula(path):
    d = json.loads(path.read_text())
    m = d["monomial_applied"]
    sign, a_exp, qt_exp = predicted(steps_of(d))
    assert m["sign"] == sign
    assert Fraction(m["a"]) == a_exp
    assert Fraction(m["q"]) == qt_exp and Fraction(m["t"]) == -qt_exp


@pytest.mark.parametrize("lam,p,n", [((1,), 2, -3), ((2,), 2, -1), ((1, 1), 3, -2), ((2, 1), 2, 1), ((2,), 3, -1)])
def test_framing_operator_shifts_slope(lam, p, n):
    """framing 演算子 τ: P_μ ↦ T_μ P_μ（T_μ = q^{n(μ')} t^{−n(μ)}）について τ^k·P_λ[p_k ↦ P_{k(p,n)}]·1 = P_λ[p_k ↦ P_{k(p,n+kp)}]·1（比はちょうど 1）。

    ケーブルの正規化 M(C_{p,q}K) = M(K)^p·M(T(p,q)) の導出に使う（conventions.md §13）。有理数の点で厳密に比較する。
    """
    from rcsknots.gn import GN
    from rcsknots.modified import solve_in_basis
    from rcsknots.partitions import n_stat, conjugate
    q, t, a = Fraction(3, 7), Fraction(5, 11), Fraction(2, 13)
    G = GN(q, t, a, lambda x: x)
    N = p * sum(lam)
    F0 = solve_in_basis(G.pattern(lam, p, n), G.P(N), Fraction(0))
    for k in (1, -1, 2):
        F1 = solve_in_basis(G.pattern(lam, p, n + k * p), G.P(N), Fraction(0))
        tau = {mu: c * (q ** n_stat(conjugate(mu)) * t ** (-n_stat(mu))) ** k for mu, c in F0.items()}
        assert {mu: c for mu, c in tau.items() if c} == {mu: c for mu, c in F1.items() if c}


@pytest.mark.parametrize("m,n", [(2, 3), (3, 2), (3, 4), (4, 3), (2, 5)])
def test_bottom_coefficient_is_rational_catalan(m, n):
    """GN の正規化前の reduced 値の a⁰ の係数は (−t/q)^{m−1} C_{m,n}(1/q, t)（C は Dyck 路の和、GN の Conjecture catalan）。"""
    import sys
    import sympy as sp
    sys.path.insert(0, str(ROOT / "scripts"))
    from d1_catalan import catalan, q, t
    from rcsknots.exact import reconstruct
    raw = reconstruct(lambda G: G.torus_reduced_apoly((1,), m, -n), log=lambda s: None)
    R = sum(c * q ** eq * t ** et for (ea, eq, et), c in raw.items() if ea == 0)
    assert sp.expand(R - (-t / q) ** (m - 1) * catalan(m, n).subs(q, 1 / q)) == 0
