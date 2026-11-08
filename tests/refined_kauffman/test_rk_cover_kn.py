"""被覆の brane–brane セクター（rk_b_cover_split.cover_sectors）が [KN17] の U(N_A) の refined LMOV（framed_knot の再現データ）に一致すること（続き 131）。
(t^{-1/2} − t^{1/2}) f^{bb}_R × (q/t)^{δ|R|/2} = −f^{t̄}_R（鏡像：(A,Q,T) → (ν,r,s)^{-1}、a = ν² = t^{N_A}）、R = [1], [2], [1,1]。
hatted（p₁² 成分を z_q、p₂ 成分を q^{1/2}+q^{-1/2} で割る M⊗M）は、鏡像のため R ↔ R^T で KN の f̂ に一致する。trefoil と T(2,5)。"""
import json, sys
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import pytest
import sympy as sp

r_, s_, u_ = sp.symbols("r s u", positive=True)
nu = sp.Symbol("nu", positive=True)


def _kn(d):
    return sum(sp.Integer(c) * nu**(-int(k.split(",")[0])) * r_**(-int(k.split(",")[1])) * s_**(-int(k.split(",")[2])) for k, c in d.items())


def _tosym(x):
    return sp.cancel(sp.cancel(x.to_sympy()).subs(u_, nu * s_))


@pytest.mark.parametrize("m,knot", [(1, "T2_3"), (2, "T2_5")])
def test_brane_brane_equals_kn(m, knot):
    import rk_b_cover_split as CS
    from rcsknots.refined_kauffman.flint_rf3 import R
    sec = CS.cover_sectors(m)
    w = (s_ / r_)**m                                   # 一様な規格化 (t/q)^{δ/2}、δ = m
    pref = 1 / s_ - s_                                 # KN の t 側の因子 (t^{-1/2} − t^{1/2})
    l1 = json.load(open(f"data/framed_knot/lmov/{knot}_1.json"))["f"]
    l2 = json.load(open(f"data/framed_knot/lmov/{knot}_2.json"))["f"]
    f1 = pref * _tosym(sec["f1"]) / w
    f2 = pref * _tosym(sec["p11"] + sec["p2"]) / w**2
    f11 = pref * _tosym(sec["p11"] - sec["p2"]) / w**2
    assert sp.cancel(f1 + _kn(l1["ft_1"])) == 0
    assert sp.cancel(f2 + _kn(l2["ft_2"])) == 0
    assert sp.cancel(f11 + _kn(l2["ft_11"])) == 0
    zq = R - 1 / R
    h2 = pref * _tosym(sec["p11"] / zq + sec["p2"] / (R + 1 / R)) / w**2
    h11 = pref * _tosym(sec["p11"] / zq - sec["p2"] / (R + 1 / R)) / w**2
    assert sp.cancel(h2 - _kn(l2["hat_11"])) == 0          # 鏡像で R ↔ R^T
    assert sp.cancel(h11 - _kn(l2["hat_2"])) == 0
    assert sp.cancel(f1 + _kn(l1["hat_1"])) == 0
