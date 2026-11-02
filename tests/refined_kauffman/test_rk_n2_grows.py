"""n = 2 の行の母関数の記号版ブロック（rk_f_n2_rows_sym、n2_rows_sym_m1.pkl、theory §12.40–12.41）の機能テスト。
(1) ν^{−12} の行の最上部の 2 次形式は −(s³/r)X_λ²。(2) ν^{−14} の最上行×次の行の交差項は −(s³/r)·2X₀X₁（X₁ は次の行の同じ線形形式）。
(3) λ = (1,1) の ν^{−12} の結び目の部分は −(s³/r)[(X+Z)² + ((q−t)/t)Z²]、Z = −t(1+qt)/q³、X + Z = −(1+q)(1+t)(1+qt)/q³（§12.27 の Y）。"""
import sys
import sympy as sp
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import rk_f_n2_grows as G
import rk_f_n2_rows_sym as RS
from rk_f_n2_rows_sym import r, s, q, t, a2, a4, b4, b6, c4, c6, SYMS


def _poly(lab, k):
    return sp.Poly(sp.expand(G.blocks()["N"][(lab, k)]), *SYMS, domain=RS.DOM)


def _coeff(P, e):
    c = sp.sympify(P.coeff_monomial(e))
    return c.as_expr() if hasattr(c, "as_expr") else c


def test_nu12_top_square():
    tops = [a2**2, b4, c4]
    for lab in ("2", "11"):
        N12 = _poly(lab, 12)
        Q = sp.Poly(sp.expand(-(s**3 / r) * G.X_of(lab)**2), *SYMS, domain=RS.DOM)
        for i in range(3):
            for j in range(i, 3):
                e = tops[i] * tops[j]
                assert sp.cancel(_coeff(N12, e) - _coeff(Q, e)) == 0, (lab, e)


def test_nu14_cross_terms():
    tops, nexts = [a2**2, b4, c4], [a2 * a4, b6, c6]
    for lab in ("2", "11"):
        N14 = _poly(lab, 14)
        X0 = G.X_of(lab)
        X1 = X0.subs({a2: sp.sqrt(2 * a2 * a4), b4: b6, c4: c6}, simultaneous=True)   # a₂² → 2a₂a₄（(a₂ + a₄ν⁻²)² の交差項）
        Q = sp.Poly(sp.expand(-(s**3 / r) * 2 * X0 * X1), *SYMS, domain=RS.DOM)
        for x in tops:
            for y in nexts:
                assert sp.cancel(_coeff(N14, x * y) - _coeff(Q, x * y)) == 0, (lab, x * y)


def test_nu12_knot_part_11_square_completion():
    """λ = (1,1)：K₁₂ = −(s³/r)[Y² + (q−t)(t/q⁴)a₂²]、Y = X + Z = −(1+q)(1+t)(1+qt)/q³。
    G の項（G の記号を含む単項式）＋結び目の Adams が −X² を −Y² に補完し、(q−t) の補正は記号ブロックの a₂² の項そのもの。"""
    val = RS.actual_rows()
    valrf = G.rows_rf(val)
    X = sp.cancel(G.X_of("11").subs(val))
    K = G.to_sym(G.knot_rf("11", 12, valrf))
    Z = -t * (1 + q * t) / q**3
    Y = -(1 + q) * (1 + t) * (1 + q * t) / q**3
    assert sp.cancel(X + Z - Y) == 0
    assert sp.cancel(K + (s**3 / r) * (Y**2 + (q - t) * t * val[a2]**2 / q**4)) == 0
    # 記号ブロックの a₂² の項の係数は −(s³/r)(q−t)t/q⁴
    c = _coeff(_poly("11", 12), a2**2)
    assert sp.cancel(c + (s**3 / r) * (q - t) * t / q**4) == 0


def test_sigma_covariance_nu12():
    """N_(2)∘σ − N_(1,1)（σ：2 箱の行と G の重みの付け替え）は ν^{−12} で随伴セクターの単項式（a₂², a₂a₄, gadj6, 定数）だけで、全係数が (q−t) で割れる。"""
    Dk = G.sigma_difference(12)
    mons = {sp.Poly.from_dict({mon: 1}, *SYMS).as_expr() for mon, _ in Dk.terms()}
    assert mons == {a2**2, a2 * a4, G.gadj[6], sp.Integer(1)}
    for _, c in Dk.terms():
        cq = sp.cancel(_coeff_expr(c) / (q - t))
        assert sp.fraction(cq)[1].is_polynomial(r, s) and sp.fraction(cq)[1].subs(s, r) != 0      # (q−t) を分母に残さない


def _coeff_expr(c):
    c = sp.sympify(c)
    return c.as_expr() if hasattr(c, "as_expr") else c


def test_nu12_lambda2_square_by_sigma_split():
    """λ = (2) の ν^{−12}：σ 共変な部分 W が完全平方 Y_(2)²、Y_(2) = [(1+q)(1+t)(1+qt)² + qt(q−t)]/q⁴。補正 C は随伴セクターの (q−t) の項。"""
    val = RS.actual_rows()
    valrf = G.rows_rf(val)
    for lab in ("2", "11"):
        W, C = G.split_sigma(lab, 12, valrf)
        assert sp.cancel(G.to_sym(W) - G.Y_LAM[lab]**2) == 0
        K = G.to_sym(G.knot_rf(lab, 12, valrf))
        assert sp.cancel(K + (s**3 / r) * (G.to_sym(W) + G.to_sym(C))) == 0
        assert sp.cancel(G.to_sym(C) / (q - t)).is_polynomial(s)                                   # C ∝ (q−t)
