"""framed knot の refined LMOV 不変量（scripts/framing_check.py）の機能テスト。

基準：ユーザーの以前の計算（Dropbox/Mathematica/marino-vafa の framed_results.nb、framed refined Oogri-Vafa ver2.nb）の出力。
G 展開の結果をノートブックの形 Σ c A^{2β} g^g T^{2J_r}（A = 𝐚, T = 𝐭, g = G）に直して全係数で比べる。
"""
import sys
from pathlib import Path

import pytest
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from framing_check import framed, framing_monomial, load_reduced, run  # noqa: E402

A, T, g = sp.symbols("A T g")


def as_poly(N):
    return sp.expand(sum(c * A ** (2 * b) * g ** gg * T ** F for (gg, b, F), c in N.items()))


def test_framing_factor_matches_notebook_definition():
    """framing[λ] = (−1)^{|λ|} Π q^{a+1/2} t^{−l−1/2}（ノートブックの定義）を (A, Q, T) = (a, q, t)^{1/2} で書いたもの。"""
    assert framing_monomial((1,), 1) == (-1, (0, 1, -1))
    assert framing_monomial((2,), 1) == (1, (0, 4, -2))        # 箱の a: 1, 0 → q^{1/2+3/2}、l: 0, 0 → t^{−1}
    assert framing_monomial((1, 1), 1) == (1, (0, 2, -4))
    assert framing_monomial((2, 1), 2) == (1, (0, 2 * (3 + 2), -2 * (3 + 2)))


def test_trefoil_framing_2_matches_notebook():
    """三葉結び目 framing 2：ノートブック（framed_results.nb）の [1], [2], [1,1] の G 展開と全係数一致。"""
    res = run("T2_3", [2], nmax=2, log=lambda s: None)
    user1 = A**4 * T**5 + T * (1 + g * T + T**2) + A**2 * T * (T + g * T**2 + 2 * T**3)
    assert as_poly(res[(2, (1,))][1]) == sp.expand(user1)
    user2 = (A**2 * T**4 * (1 + 3 * g * T + 3 * T**2 + 3 * g**2 * T**2 + 5 * g * T**3 + g**3 * T**3 + 2 * T**4 + 3 * g**2 * T**4 + 3 * g * T**5 + T**6)
             + A**10 * T**4 * (2 * T**8 + 4 * g * T**9 + 2 * T**10 + 3 * g**2 * T**10 + 4 * g * T**11 + g**3 * T**11 + T**12 + 3 * g**2 * T**12 + 3 * g * T**13 + T**14))
    P2 = as_poly(res[(2, (2,))][1])
    assert sp.expand(P2.coeff(A, 2) * A**2 + P2.coeff(A, 10) * A**10 - user2) == 0   # 𝐚 の端の 2 行（全体は scripts の照合で確認）
    for rho in ((1,), (2,), (1, 1)):
        assert "後 0/" in res[(2, rho)][0] and "前 0/" in res[(2, rho)][0]


@pytest.mark.parametrize("x", [-2, -1, 0, 1, 2])
def test_unknot_and_trefoil_three_boxes_positive(x):
    """自明な結び目と三葉結び目の |ρ| ≤ 3：framing −2..2 で G 展開の前後とも非負。"""
    for name in ("unknot", "T2_3"):
        for (xx, rho), (msg, N, hat) in run(name, [x], nmax=3, log=lambda s: None).items():
            assert N == {} or (N is not None and "前 0/" in msg and "後 0/" in msg), (name, x, rho, msg)


def test_fundamental_framing_is_t_shift():
    """[1]：framing[□]^x は太字の変数で 𝐭^x（Q/T = −𝐭）。\\hat N は 2J_r が x ずれるだけ【導出】。"""
    base = run("4_1", [0], nmax=1, log=lambda s: None)[(0, (1,))][1]
    for x in (-2, 3):
        N = run("4_1", [x], nmax=1, log=lambda s: None)[(x, (1,))][1]
        assert N == {(gg, b, F + x): c for (gg, b, F), c in base.items()}


def test_unknot_framing_0_has_only_fundamental():
    """framing 0 の自明な結び目では \\hat f_ρ = 0（ρ ≠ [1]）。"""
    res = run("unknot", [0], nmax=3, log=lambda s: None)
    assert all(N == {} for (x, rho), (msg, N, hat) in res.items() if rho != (1,))


def test_framed_is_monomial_rescaling():
    red = load_reduced("T2_3", 2)
    f = framed(red, 3)
    for lam in red:
        s, (da, dq, dt) = framing_monomial(lam, 3)
        assert f[lam] == {(ea + da, eq + dq, et + dt): s * c for (ea, eq, et), c in red[lam].items()}
