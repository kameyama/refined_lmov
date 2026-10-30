"""n = 2（trefoil）の N_λ の結び目のデータごとの分解（rk_f_n2_rterms、theory §12.26、progress.md 続き 92）の機能テスト。
  (1) 基準 R = 0 の背景 T(∅) の閉じた形（6 セルの Fock 空間）。
  (2) 包除：[α;1] だけの項 T(c) は 0、λ = (2) の上端 2 ブロックは N = 背景 + 結び目の部分 −(q − t)(1 + qt)²。"""
import sys
import sympy as sp
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import rk_f_n2_rterms as RT
import rk_f_cover_fixed22 as CF
from rk_f_n2_univ import hatted_N

r_, s_ = RT.r_, RT.s_


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
