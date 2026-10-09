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


# ---------------------------------------------------------------- framing を記号として扱う（構造の導出）
def test_two_box_denominator_loci_framing_factors_agree():
    """2 箱の公式の分母の各既約因子の零点で、残る 2 つの入力の framing の因子が一致する ⇒ Laurent 性と命題 2 は framing によらない。"""
    from v4_closed import formula
    Asym, Q, T = sp.symbols("A Q T")
    mult = [Q**4 * T**-2, Q**2 * T**-4, Q**2 * T**-2, -Q**2 * T**-2]          # r2, r11, r1², r1(Adams) の framing[λ]（x = 1）
    for rho, (ns, den) in formula().items():
        polys = [sum(c * Asym**i * Q**j * T**k for (i, j, k), c in n.to_dict().items()) for n in ns]
        for locus in ({Q: 1}, {Q: -1}, {Q: sp.I}, {T: 1}, {T: -1}, {T: sp.I}, {Q: 1 / T}, {Q: -1 / T}):
            alive = [i for i, P in enumerate(polys) if sp.simplify(P.subs(locus)) != 0]
            assert len(alive) == 2
            assert len({sp.simplify(mult[i].subs(locus)) for i in alive}) == 1, (rho, locus)


def test_framing_factor_transposition():
    """framing[λ](1/t, 1/q) = framing[λ^T](q, t) ⇒ framing を付けても転置対称性（G 展開が閉じる条件の σ）は保たれる。"""
    from rcsknots.partitions import conjugate
    for lam in [(1,), (2,), (1, 1), (3,), (2, 1), (1, 1, 1)]:
        s, (_, eq, et) = framing_monomial(lam, 1)
        sT, (_, eqT, etT) = framing_monomial(conjugate(lam), 1)
        assert s == sT and (-et, -eq) == (eqT, etT)          # (Q, T) → (1/T, 1/Q)


def test_two_box_framing_is_exponential_sum():
    """𝐭^{−2x}\\hat f^{(x)} は x の指数和（根 𝐪²𝐭², 𝐪⁻², 1, −1）：x = 0..3 で係数を決め、x = 4, −1, −2 を予測（三葉結び目、有理数の点で）。"""
    from fractions import Fraction as Fr
    from framing_newton import S_values
    red = load_reduced("T2_3", 2)
    pt = (Fr(2), Fr(3), Fr(5))
    roots = [pt[1] ** 2 * pt[2] ** 2, pt[1] ** -2, Fr(1), Fr(-1)]
    for rho in [(2,), (1, 1)]:
        v = {x: sum(Fr(c) * pt[0] ** i * pt[1] ** j * pt[2] ** k for (i, j, k), c in S_values("T2_3", 2, x, red)[rho].items())
             for x in (-2, -1, 0, 1, 2, 3, 4)}
        c = sp.Matrix([[r ** x for r in roots] for x in range(4)]).LUsolve(sp.Matrix([v[x] for x in range(4)]))
        assert all(sum(c[i] * roots[i] ** x for i in range(4)) == v[x] for x in (4, -1, -2))


def test_trefoil_two_box_all_framings_newton():
    """三葉結び目の 2 箱：4 つの (偶奇, 向き) でニュートン形の係数が一様な符号になる起点があり、起点と 0 の間も直接判定で正。"""
    from framing_newton import run as nrun
    res = nrun("T2_3", 2, 4, log=lambda s: None)
    for key, found in res.items():
        assert found and found[1] is not None and found[3], key
        assert all(g[0] == found[2] and g[1] is True for x, g in found[4]), key


# ---------------------------------------------------------------- t = q での Mariño–Vafa との整合性
def test_t_eq_q_limit_is_marino_vafa_framing():
    """t = q で framing[λ]^x = (−1)^{ℓx} q^{xκ_λ/2}：Mariño–Vafa（hep-th/0108064、式 \\framefin と \\kapar）の
    W^{(p)}_R = (−1)^{ℓp} q^{pκ_R/2} W_R、κ_R = ℓ + Σ_i (l_i² − 2i l_i)（i は 1 から）と p = x で一致する。
    (A, Q, T) = (a, q, t)^{1/2} なので t = q は T = Q、q^{κ/2} = Q^κ。"""
    for lam in [(1,), (2,), (1, 1), (3,), (2, 1), (1, 1, 1), (3, 2, 2, 1)]:
        ell = sum(lam)
        kappa = ell + sum(l * l - 2 * i * l for i, l in enumerate(lam, start=1))
        for x in (-3, -1, 1, 2):
            s, (ea, eq, et) = framing_monomial(lam, x)
            assert (s, ea, eq + et) == ((-1) ** (ell * x), 0, x * kappa)


def test_t_eq_q_unknot_matches_marino_vafa_variables():
    """変数の向き：t = q の本研究の unknot（主特殊化）は MV の量子次元 Π (λ^{1/2}q^{c/2} − λ^{−1/2}q^{−c/2})/(q^{h/2} − q^{−h/2})
    と a = λ、q = q_MV で単項式倍を除いて一致する（q ↔ 1/q では [2] と [1,1] が入れ替わり一致しない）。
    したがって framing の向き x = p も含めて MV と同じ規約。a = A²、q = Q² として整数冪で比べる。"""
    from rcsknots.macdonald import principal_specialization
    from rcsknots.partitions import conjugate
    Asym, Q = sp.symbols("A Q")

    def is_monomial(e):
        num, den = sp.fraction(sp.cancel(sp.together(e)))
        return all(len(sp.Poly(sp.expand(z), Asym, Q).terms()) == 1 for z in (num, den))

    for lam in [(2,), (1, 1), (2, 1), (3,), (3, 1)]:
        lamT = conjugate(lam)
        mv = sp.Integer(1)
        for i, row in enumerate(lam):
            for j in range(row):
                h = (row - j - 1) + (lamT[j] - i - 1) + 1
                mv *= (Asym * Q ** (j - i) - Q ** (i - j) / Asym) / (Q ** h - Q ** -h)
        assert is_monomial(principal_specialization(lam, Asym**2, Q**2, Q**2) / mv)
        if lamT != lam:                                    # 自己共役な [2,1] は q ↔ 1/q で区別できない
            assert not is_monomial(principal_specialization(lam, Asym**2, Q**-2, Q**-2) / mv)


def test_marino_vafa_sign_is_needed_for_laurent():
    """MV の符号 (−1)^{|λ|x} を外し T_λ^x だけを掛けると、x が奇数のとき 2 箱の公式の分母の零点 Q = i, T = i で
    残る 2 入力の因子が食い違う（Adams の写しの項だけ符号が変わる）⇒ 一般の入力で \\hat f_{[2]}, \\hat f_{[1,1]} が Laurent でなくなる。
    箱あたりの (q/t)^{1/2} は外しても食い違わない（Adams 作用と両立）。"""
    from v4_closed import formula
    Asym, Q, T = sp.symbols("A Q T")
    one = sp.Integer(1)

    def mismatched(mult):
        bad = set()
        for rho, (ns, den) in formula().items():
            polys = [sum(c * Asym**i * Q**j * T**k for (i, j, k), c in n.to_dict().items()) for n in ns]
            for locus in ({Q: 1}, {Q: -1}, {Q: sp.I}, {T: 1}, {T: -1}, {T: sp.I}, {Q: 1 / T}, {Q: -1 / T}):
                alive = [i for i, P in enumerate(polys) if sp.simplify(P.subs(locus)) != 0]
                if len({sp.simplify(mult[i].subs(locus)) for i in alive}) != 1:
                    bad.add((rho, tuple(locus.items())[0]))
        return bad

    for x in (1, 3, -1):                                   # r2, r11, r1², r1(Adams) の因子
        assert not mismatched([Q ** (2 * x), T ** (-2 * x), one, (-one) ** x])
        assert mismatched([Q ** (2 * x), T ** (-2 * x), one, one]) == {
            ((2,), (Q, sp.I)), ((2,), (T, sp.I)), ((1, 1), (Q, sp.I)), ((1, 1), (T, sp.I))}
    assert not mismatched([Q**4, T**-4, one, one])         # x = 2：符号は効かない


@pytest.mark.parametrize("x", [-3, 2])
def test_framing_big_matches_framing_check(x):
    """scripts/framing_big.py（\\hat f だけを 1 回の補間で）が framing_check.py（ρ ごとの補間）と全係数で一致する（三葉結び目、3 箱）。"""
    from framing_big import hat_f_all
    R = framed(load_reduced("T2_3", 3), x)
    hats = hat_f_all(R, 3, log=lambda s: None, workers=1, cache=None)
    ref = run("T2_3", [x], nmax=3, log=lambda s: None)
    for rho, hat in hats.items():
        assert hat == ref[(x, rho)][2]
