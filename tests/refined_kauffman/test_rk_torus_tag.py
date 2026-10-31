"""rk_torus_tag（続き 102）：T(mk+1, k) のタグと δ。δ は data（KN の規格化）の U(2) の行列式表現 R_(1,1)(a = t) = (t/q)^δ から読む（conventions.md「一様な規格化」）。"""
import sys, json
import sympy as sp
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import rk_torus_tag as TG

a, q, t = sp.symbols("a q t")


def _R11(knot):
    d = json.load(open(f"data/framed_knot/torus/{knot}_11.json"))
    return sum(v * a**i * q**sp.Rational(j, 2) * t**sp.Rational(l, 2) for k, v in d["poly"].items() for i, j, l in [map(int, k.split(","))])


def test_tags():
    assert TG.tag(3) == "m3" and TG.tag(1, 3) == "m1k3"
    assert TG.knot_name(7) == "T2_15" and TG.knot_name(1, 3) == "T3_4" and TG.knot_name(2, 3) == "T3_7"
    assert TG.rs(7) == 30 and TG.rs(1, 3) == 12


def test_delta_from_data():
    for m, k in ((1, 2), (2, 2), (1, 3), (1, 4), (2, 3)):
        assert sp.simplify(_R11(TG.knot_name(m, k)).subs(a, t) - (t / q)**TG.delta(m, k)) == 0
