"""ステージ T4：2 箱（[2], [1,1]）の入力の検証と、9_42 の破れの原因。"""
import sys
from collections import Counter
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/framed_knot"))
sys.path.insert(0, str(ROOT / "data/framed_knot/reference/gs_1112.0030"))

import s2_thick as gs  # noqa: E402
from colored_from_cables import colored, load  # noqa: E402
from t4_colored import tc_colored  # noqa: E402
from t4_egp import square_at_q1  # noqa: E402
from rcsknots.bps import to_bold  # noqa: E402
from rcsknots.homfly import Q, a  # noqa: E402
from rcsknots.knot_table import load_homology, superpoly  # noqa: E402
from rcsknots.laurent import load_table  # noqa: E402

H = load_homology()
R = load_table(ROOT / "data/framed_knot/reference/notebook_torus.txt")
LINE2 = gs.terms((1 + gs.q) * (1 + gs.a**-1 * gs.q * gs.t**-1) * (1 + gs.a**-1 * gs.q**-2 * gs.t**-3)
                 * (1 + gs.a**-1 * gs.q**-3 * gs.t**-3) * (gs.a**2 * gs.q**2 * gs.t**6 + gs.a * gs.q**3 * gs.t**5))
PRINTED = gs.terms(gs.P2_942)
FIXED = {k: PRINTED[k] - LINE2.get(k, 0) for k in PRINTED if PRINTED[k] - LINE2.get(k, 0)}


def kn_t_equals_q(knot, color):
    A_, Q_, T_ = sp.symbols("A Q T")
    P = R[(knot, color)]
    e = sp.expand(sum(c * A_**i * Q_**j * T_**k for (i, j, k), c in P.items()).subs(T_, Q_))
    return e.subs({A_: sp.sqrt(a), Q_: Q})


def test_cable_inversion_reproduces_known_colored_homfly():
    """2 本ケーブル（s = 1, 3）の HOMFLY-PT からケーブル公式を逆に解くと、三葉結び目と 4_1 の [2], [1,1] が一致する。"""
    data = load(ROOT / "data/framed_knot/reference/two_cable_homfly.txt")
    for name, knot in [("3_1", "T2_3"), ("4_1", "4_1")]:
        C = colored(name, data)
        assert sp.simplify(C[(2,)] - kn_t_equals_q(knot, "2")) == 0
        assert sp.simplify(C[(1, 1)] - kn_t_equals_q(knot, "11")) == 0


def euler_bold(gens):
    x, y = sp.symbols("x y")
    return sp.expand(sum(c * (-1) ** (k % 2) * x ** (2 * i) * y ** (2 * j) for (i, j, k), c in gens.items())), x, y


def test_gs_942_s2_misprint():
    """GS 1112.0030 の 9_42 の S² は 2 行目を除くと Euler 標数が正しい [2] 彩色 HOMFLY-PT になる（印刷どおりでは違う）。"""
    C = colored("9_42", load(ROOT / "data/framed_knot/reference/two_cable_homfly.txt"))
    for gens, ok in [(PRINTED, False), (FIXED, True)]:
        e, x, y = euler_bold(gens)
        assert (sp.simplify(e - C[(2,)].subs({a: x**2, Q: y})) == 0) == ok
    assert sum(PRINTED.values()) == 433 and sum(FIXED.values()) == 401


def test_tc_conversion_reproduces_torus_knots():
    """t_r → t_c の変換で三葉結び目（δ' = −1 一定）が KN の rCS_[2], rCS_[1,1] に一致。8_19 は (𝐚, 𝐪) の重複集合が一致。"""
    two, one_one = tc_colored(gs.terms(gs.P2_31), lambda x: -1)
    assert two == to_bold(R[("T2_3", "2")], 0) and one_one == to_bold(R[("T2_3", "11")], 0)
    gs819 = Counter()
    for (i, j, k), c in gs.terms(gs.P2_819).items():
        gs819[(2 * i, 2 * j)] += c
    ours = Counter()
    for (e_a, e_q, _), c in to_bold(R[("T3_4", "2")], 0).items():
        ours[(e_a, e_q)] += c
    assert gs819 == ours


def test_egp_forces_exponential_growth_of_rank():
    """命題 2：非負の入力で \\hat f_[2] が Laurent なら dim H^[2] = (dim H^[1])²。右辺の打ち消しの無い大きさを確認。"""
    for name, gens in [("3_1", gs.terms(gs.P2_31)), ("8_19", gs.terms(gs.P2_819)), ("9_42", FIXED)]:
        P1 = superpoly(H[name])
        need = sum(abs(v) for v in square_at_q1(P1).values())
        assert need == sum(P1.values()) ** 2
        assert (sum(gens.values()) == need) == (name != "9_42")
