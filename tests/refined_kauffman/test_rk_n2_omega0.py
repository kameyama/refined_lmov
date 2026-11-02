"""単一の brane の λ 射影キュムラント Ω^{(0)}_λ（rk_f_n2_omega0、theory §12.49、続き 119）の機能テスト：
ν^{−4} の係数が両 λ で Ŷ_λ に比例し、比例定数は fac = −r s²/(t−1)（λ によらない）。ν^{−2} の係数は (1,1) で −a₂/q、(2) で −a₂²/q。"""
import sys
import sympy as sp
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import rk_f_n2_omega0 as O
import rk_f_n2_rows_sym as RS
import rk_f_n2_grows as G
from rk_f_n2_rows_sym import r, s, q, t


def test_second_level_is_Yhat():
    RS.SYMS = [sp.Symbol("_dummy")]; RS.ZERO = RS.P(0)
    Om, f, Zp = O.omega0()
    val = RS.actual_rows()
    fac = -r * s**2 / (t - 1)
    for lab, lvl2 in (("11", -val[RS.a2] / q), ("2", -val[RS.a2]**2 / q)):
        ser = O.top(Om[lab], 8)
        assert sp.cancel(ser[4] - fac * G.Y_LAM[lab]) == 0
        assert sp.cancel(ser[2] - fac * lvl2) == 0
