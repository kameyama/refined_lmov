"""refined Kauffman K2：refined unknot の閉じた式の検査。

- D_N refined unknot（AS）の安定な式 dim_DN が N = 5, 6 の一般 t の評価式と一致。
- Koornwinder の主特殊化（van Diejen–Sahi）がパラメータ (1, −1, q^{1/2}, −q^{1/2}) で dim_DN と一致。
- U(M) composite の Macdonald 次元の安定な式が M = 8 の直接計算と一致。
- 極：dim_DN([2]) の分母に (Uq − t)、dim_DN([1,1]) に (U²q − t³)、composite ([1],[1]) に (Aq − t)。
"""
import importlib.util
import pathlib
import sys

import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))  # worktree の src を優先（.venv は本体の src を指す）
from rcsknots.refined_kauffman.kauffman_refined import (A, U, composite_diagram, dim_A_sym, dim_DN, dim_composite,
                                       koornwinder_eval, metric_DN, q, t)


def _eval_generic():
    spec = importlib.util.spec_from_file_location("m", ROOT / "scripts" / "refined_kauffman" / "rk_so_unknot_check.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_dim_DN_matches_generic_evaluation():
    m = _eval_generic()
    for lam in ((1,), (2,), (1, 1), (2, 1)):
        for N in (5, 6):
            d = sp.simplify(dim_DN(lam).subs(U, t**N) - m.eval_generic(lam, N).subs({m.q: q, m.t: t}))
            assert d == 0


def test_koornwinder_evaluation_is_DN():
    sq = sp.sqrt(q)
    for lam in ((1,), (2,), (1, 1), (3,), (2, 1)):
        assert sp.simplify(koornwinder_eval(lam, (1, -1, sq, -sq)) - dim_DN(lam)) == 0


def test_composite_dimension_stable_form():
    for lam, mu in (((1,), (1,)), ((2,), (1,)), ((1,), (1, 1)), ((2, 1), (1,))):
        e = dim_composite(lam, mu)
        d = sp.simplify(e.subs(A, t**8) - dim_A_sym(composite_diagram(lam, mu, 8), 8))
        assert d == 0


def _den_factors(expr, var):
    num, den = sp.fraction(sp.factor(expr))
    return {sp.factor(f) for f, m in sp.factor_list(den)[1] if sp.Poly(f, var).degree() > 0}


def test_poles():
    assert sp.factor(U*q - t) in _den_factors(dim_DN((2,)), U)
    assert sp.factor(U**2*q - t**3) in _den_factors(dim_DN((1, 1)), U)
    assert sp.factor(A*q - t) in _den_factors(dim_composite((1,), (1,)), A)
    # 計量は同じ因子を 2 乗でもつ
    num, den = sp.fraction(sp.factor(metric_DN((2,))))
    assert dict(sp.factor_list(den)[1]).get(sp.factor(U*q - t)) == 2
