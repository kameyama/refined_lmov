"""固定セクターの捻ったトレースによる E の整数計画（progress.md 続き 39）。"""
import sys, io, contextlib
sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import pytest
import rk_b_twisted_multiplet as M


@pytest.mark.results
def test_trefoil_E_unique_with_sum_boundary_factor():
    """trefoil の 2 箱：前の因子 z_q^+ と Pk の単項式ごとの符号で (i)+(ii) の解が一意、z_q では解なし。"""
    with contextlib.redirect_stdout(io.StringIO()):
        res = M.knot(1, twists=("mono",))
    assert ("mono", "z_q^+") in res and res[("mono", "z_q^+")]["free"] == []
    assert ("mono", "z_q") not in res
