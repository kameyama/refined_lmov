"""案 B：ORS のセル分解（Γ = ⟨2,3⟩）が T(2,3) の DGR の unreduced superpolynomial を再現する（theory §11.18）。"""
import sys
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rk_b_planB_cells import P_alg, a, q, t


def test_cells_reproduce_trefoil():
    lmax = 5
    Pa = P_alg(2, 3, lmax)
    target = sp.expand(sp.series((a * t + 1 / a) * q / (1 - q**2) * (a**2 * q**-2 + a**4 * t**3 + a**2 * q**2 * t**2), q, 0, 2 * lmax).removeO())
    diff = sp.expand((Pa - target) * q)
    low = sum(term for term in sp.Add.make_args(diff) if sp.Poly(term, q).degree() <= 2 * lmax - 2)
    assert sp.expand(low) == 0
