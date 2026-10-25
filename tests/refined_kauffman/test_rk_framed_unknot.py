"""framed unknot（framing 2）の固定セクターで E が一意に決まること（theory_orientifold_index.md §11.8、Young の loop quiver の予言）。"""
import sys, io, contextlib
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import pytest
import rk_b_framed_unknot as F


@pytest.mark.results
def test_framed_unknot_x2_E_unique():
    with contextlib.redirect_stdout(io.StringIO()):
        res = F.solve_E(2, "A", quiet=True)
    assert res is not None and res["feasible"] and res["free"] == []
    r, s, nu = sp.symbols("r s nu", positive=True)
    E = res["E"].subs({F.r_: r, F.s_: s})
    expected = (nu**2 + 1) * (r**2 + 1) * (nu**2 * r**2 - s**2) / s**2
    assert sp.cancel(sp.sympify(str(E), locals={"r": r, "s": s, "nu": nu}) - expected) == 0


def test_framed_unknot_x1_E_unique_by_enumeration():
    """framing 1：(ii) の形を仮定せず、(i) と「t = q で C が整数」だけで E が一意（16 候補の全数）。"""
    with contextlib.redirect_stdout(io.StringIO()):
        ok = F.enumerate_E(1)
    r, s, nu = sp.symbols("r s nu", positive=True)
    assert len(ok) == 1
    E = sp.sympify(str(ok[0]), locals={"r": r, "s": s, "nu": nu})
    assert sp.cancel(E + (nu**2 + 1) * (r**2 + 1)) == 0
