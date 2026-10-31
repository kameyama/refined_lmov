"""n = 2（trefoil）の N_λ の結び目のデータごとの分解（rk_f_n2_rterms、theory §12.26、progress.md 続き 92）の機能テスト。
  (1) 基準 R = 0 の背景 T(∅) の閉じた形（6 セルの Fock 空間）。
  (2) 包除：[α;1] だけの項 T(c) は 0、λ = (2) の上端 2 ブロックは N = 背景 + 結び目の部分 −(q − t)(1 + qt)²。"""
import sys
import sympy as sp
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import pytest
import rk_f_n2_rterms as RT
import rk_f_cover_fixed22 as CF
import rk_c_2box_largeN as L
from rk_f_n2_univ import hatted_N

r_, s_ = RT.r_, RT.s_


@pytest.fixture(autouse=True)
def _restore_knot_data():
    """set_groups は CF.load_poly_A と L.Rdata_nu を差し替えるので、各テストの後に元に戻す（他のテストファイルに影響させない）。"""
    yield
    CF.load_poly_A, L.Rdata_nu = RT._ORIG


def _N(on, base):
    RT.set_groups(on, base=base)
    Ns, _ = hatted_N()
    return Ns


def test_background_closed_form():
    B = _N((), 0)
    q, t, nu, A = CF.q, CF.t, CF.nu, CF.A
    den = (A - t) * (A * q - 1) * (A * q - t) * (A * q - t**2) * (A * q**2 - t) * (A**2 * q**2 - t**2)
    assert (B["2"] + CF.R**7 * CF.S**3 * (q - t) * nu**6 / den).is_zero()
    assert (B["11"] + CF.R**3 * CF.S**7 * (q - t) * nu**2 / den).is_zero()
    blocks = {k: RT.nu_block(v, k) for k, v in RT.series_u(B["2"], 14).items()}
    assert [len(blocks[k].n.coeffs()) for k in (8, 10, 12, 14)] == [2, 6, 12, 20]       # k(k+1)
    assert all(blocks[k].is_laurent() for k in blocks)


def test_terms_and_top_blocks():
    B, Nc, N = _N((), 0), _N(("c",), 0), _N(RT.GROUPS, 0)
    assert (Nc["2"] - B["2"]).is_zero() and (Nc["11"] - B["11"]).is_zero()            # T(c) = 0
    assert N["2"].is_laurent() and N["11"].is_laurent()
    q, t = r_**2, s_**2
    K = {k: RT.nu_block(v, k) for k, v in RT.series_u(N["2"] - B["2"], 10).items()}
    assert 8 not in K or K[8].is_zero()
    assert sp.expand(RT.to_sym(K[10]) + r_**-11 * s_**3 * (q - t) * (1 + q * t)**2) == 0
    Nb = {k: RT.nu_block(v, k) for k, v in RT.series_u(N["2"], 10).items()}
    V, Vp = 1 + t + q * t, 1 + q + q * t
    assert sp.expand(RT.to_sym(Nb[10]) + r_**-11 * s_**3 * (q - t) * V * Vp) == 0


def test_nu12_square_structure():
    """ν^{−12} の結び目の部分（続き 93、theory §12.27）：λ = (1,1) は −Y⊗Y − (q − t)·t·(1 + qt)⊗(1 + qt)、Y = (1 + q)(1 + t)(1 + qt)。
    λ = (2) は t = q で −(1 + q)⁴(1 + q²)⁴、Y₂ = Y(1 + qt) で割った補正 R は 17 状態・Σ|R| = 82。"""
    N, B = _N(RT.GROUPS, 0), _N((), 0)
    q, t = r_**2, s_**2
    K11 = RT.to_sym(RT.nu_block(RT.series_u(N["11"] - B["11"], 12)[12], 12))
    Y = (1 + q) * (1 + t) * (1 + q * t)
    assert sp.expand(K11 + r_**-13 * s_**3 * (Y**2 + (q - t) * t * (1 + q * t)**2)) == 0
    K2 = sp.expand(RT.to_sym(RT.nu_block(RT.series_u(N["2"] - B["2"], 12)[12], 12)) * r_**17 / s_**3)
    assert sp.expand(K2.subs(s_, r_) + (1 + r_**2)**4 * (1 + r_**4)**4) == 0
    Q, R = sp.div(RT._qt_poly(sp.expand(-K2 - (Y * (1 + q * t))**2)), sp.Poly(RT.QQ - RT.TT, RT.QQ, RT.TT))
    assert R.is_zero and len(Q.as_dict()) == 17 and sum(abs(v) for v in Q.as_dict().values()) == 82


def test_nu12_hypothesis_rows():
    """上端の E と λ = (1,1) の ν^{−12} の E（全因子に入れ替え、σ = −1、真空 (−,−)）を固定すると (i)+(ii) の行に λ = (2) の解があり、
    σ = +1・真空 (+,+) の λ = (1,1) では解が無い（hyp12 の 187 通り中 1 通り）。"""
    mult, keys, A, b, idx, sel, rows_, _ = RT._n2_rows(True, -12)
    unknown = [k for k in sel if k[4] == -12 and k[0] == "2"]
    assert RT._ranges(mult, keys, A, b, rows_, idx, unknown, {**RT._top_E(), **RT._E11_12()}) is not None
    q, t = r_**2, s_**2
    wrong = RT._e_of(r_**-13 * s_**3 * (-(1 + q**2) * (1 + t**2) * (1 + q**2 * t**2) - (q - t) * t * (1 + q**2 * t**2) - q * t**2 * (q - t)), "11", -12)
    assert RT._ranges(mult, keys, A, b, rows_, idx, unknown, {**RT._top_E(), **wrong}) is None
