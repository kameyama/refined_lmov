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


@pytest.mark.parametrize("M,N", [(2, 3), (3, 2), (2, 5)])
def test_mellit_cor_daha_normalization(M, N):
    """Mellit 1704.07630 の Thm (dyck paths) の 𝒫_{M,N} は −(a(qt)^{−1/2})^δ (1/(1−t)) (Q_{M,N}·1)[a/q−1]/(a/q−1)（BGLX の Q）。
    Cor (daha) の文面どおりの代入 [a−1] では、どの符号でも一致しない（scripts/d1_mellit_cor.py の総当たりの一部）。"""
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    from d1_bglx import Qop
    from d1_ops2 import dyck, pleth
    a, s, u = Fraction(2, 7), Fraction(2, 3), Fraction(3, 5)
    q, t, sq = s * s, u * u, s * u
    d = (M - 1) * (N - 1) // 2
    V = Qop(M, N, {(): Fraction(1)}, q, t)
    D = dyck(M, N, a, q, t, sq)
    pref = (a / sq) ** d / (1 - t)
    assert D == -pref * pleth(V, a / q) / (a / q - 1)
    assert D not in (pref * pleth(V, a) / (a - 1), -pref * pleth(V, a) / (a - 1))


def test_n1n2_is_cd14_tilde_times_corner():
    """(N1)(N2) の [1] と CD14（1408.4348）の tilde 正規化の対応：係数が非負の全ての保存値で、太字の P は
    𝐚^{2δ}𝐪^{−2δ} × (最低 𝐚 の行の最低の項が 1 になる tilde 正規化) で、DGR の対称性 (e_𝐚, e_𝐪, e_𝐭) ↔ (e_𝐚, −e_𝐪, e_𝐭 − e_𝐪)
    と 𝐪 の指数の偶数性（命題 1 の条件）をもつ。CD14 の明示例との比も 𝐚^{2δ}𝐪^{−2δ}（scripts/check_cd14_tilde.py）。"""
    from rcsknots.bps import to_bold
    files = sorted(ROOT.glob("data/torus/*_1.json")) + sorted(ROOT.glob("data/cables/*_1.json"))
    n = 0
    for f in files:
        d = json.loads(f.read_text())
        P = to_bold({tuple(map(int, k.split(","))): v for k, v in d["poly"].items()}, 0)
        if any(v < 0 for v in P.values()):
            continue
        _, delta, _ = predicted(steps_of(d))
        amin = min(k[0] for k in P)
        low = {k: v for k, v in P.items() if k[0] == amin}
        assert amin == 2 * delta and low.get((2 * delta, -2 * delta, 0)) == 1, d["knot"]
        assert all(k[1] >= -2 * delta and k[2] >= 0 for k in low), d["knot"]
        assert all(P.get((ea, -eq, et - eq)) == c and eq % 2 == 0 for (ea, eq, et), c in P.items()), d["knot"]
        n += 1
    assert n >= 17
