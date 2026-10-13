"""Stage 0: 変数の対応（1703.05408 ↔ DGR ↔ DAHA）を三葉結び目で確認する。"""
import sympy as sp

from rcsknots.notebook_ref import torusknotsym

a, q, t = sp.symbols("a q t", positive=True)
A, Q, T = sp.symbols("A Q T", positive=True)  # DGR (bold) 変数 a, q, t
DGR_TREFOIL = A**2 * Q**-2 + A**2 * Q**2 * T**2 + A**4 * T**3  # math/0505662, Example T23


def kn_to_dgr(expr):
    """1703.05408 (cov): a = -A² T, q^{1/2} = -Q T, t^{1/2} = Q."""
    hq, ht = sp.symbols("hq ht", positive=True)
    e = expr.subs({q: hq**2, t: ht**2})
    e = sp.powsimp(sp.expand_power_base(e, force=True), force=True)
    return sp.simplify(e.subs({a: -A**2 * T, hq: -Q * T, ht: Q}))


def test_notebook_trefoil_is_dgr():
    """ノートブックの reduced 三葉結び目（T(2,3) = torusknotsym[1,1] / torusknotsym[1,0]）は DGR の式に一致する。"""
    reduced = torusknotsym(1, 1, a, q, t, sp.sqrt(a), sp.sqrt(q), sp.sqrt(t)) / \
        torusknotsym(1, 0, a, q, t, sp.sqrt(a), sp.sqrt(q), sp.sqrt(t))
    reduced = sp.factor(sp.simplify(reduced))
    assert sp.simplify(reduced - a * (1 + q * t - a * q) / sp.sqrt(q * t)) == 0
    assert sp.expand(kn_to_dgr(reduced) - DGR_TREFOIL) == 0


def test_daha_trefoil_matches_up_to_monomial():
    """DAHA (CD14 qtareli): t = Q², q = Q² T², a_DAHA = A² T。すなわち a_DAHA = −a_KN。"""
    daha = 1 + q * t + a * q
    in_dgr = sp.expand(daha.subs({a: A**2 * T, q: Q**2 * T**2, t: Q**2}))
    assert sp.expand(in_dgr * A**2 * Q**-2 - DGR_TREFOIL) == 0
    kn = a * (1 + q * t - a * q) / sp.sqrt(q * t)
    assert sp.simplify(kn / daha.subs(a, -a) - a / sp.sqrt(q * t)) == 0
