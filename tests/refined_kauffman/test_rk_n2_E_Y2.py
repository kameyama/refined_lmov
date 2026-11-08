"""E への適用（rk_f_n2_E_Y2、theory §12.43、続き 113）の機能テスト：λ = (2) の ν^{−12} の結び目の部分が
K₁₂ = −r^{−17}s³[Ŷ² + (q−t)R]（Ŷ = (1+q)(1+t)(1+qt)² + qt(q−t)）と書け、R は 15 状態・Σ|R| = 52 で、行の 4 つの部分の和に一致する。"""
import sys
import sympy as sp
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import rk_f_n2_E_Y2 as EY


def test_setup_R_and_parts():
    K2, EB, Yh, R, parts = EY.setup()          # setup の中で (q−t) で割り切れることと parts の和 = R を assert している
    Rm = {m: int(v) for m, v in EY._qt_poly(R).as_dict().items()}
    assert len(Rm) == 15 and sum(abs(v) for v in Rm.values()) == 52
    q, t = EY.q, EY.t
    assert sp.expand(Yh.subs(t, q) - (1 + q)**2 * (1 + q**2)**2) == 0
    assert sp.expand(K2 + Yh**2 + (q - t) * R) == 0
