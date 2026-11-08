"""λ = (1,1) の ν^{−12} の E の余分な (q−t)qt² は背景の E（theory §12.45、続き 115）：
B₁₂(1,1) = −r^{−13}s³(q−t)qt²、§12.27 の E₁₂(1,1) = r^{−13}s³[ψ²Ŷ + (q−t)t ψ²(1+qt)] − B₁₂(1,1)（結び目の部分は全部入れ替え）。"""
import sys
import sympy as sp
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import rk_f_n2_grows as G
import rk_f_n2_rterms as NR
from rk_f_n2_rows_sym import r, s, q, t


def test_extra_term_is_background():
    P = r**-13 * s**3
    B = G.background("11", 12)
    assert sp.cancel(B + P * (q - t) * q * t**2) == 0
    psi2 = lambda F: sp.expand(F.subs({r: r**2, s: s**2}, simultaneous=True))
    Yh = (1 + q) * (1 + t) * (1 + q * t)
    E_saved = NR._poly_of(NR._E11_12(), "11", -12)
    assert sp.expand(E_saved - P * (psi2(Yh) + (q - t) * t * psi2(1 + q * t)) + B) == 0
