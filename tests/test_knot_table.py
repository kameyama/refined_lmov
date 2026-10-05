"""ステージ T0：結び目表の HOMFLY-PT ホモロジー（yui-kr）の規約の確定。"""
import json
from pathlib import Path

import pytest
import sympy as sp

from rcsknots.bps import to_bold
from rcsknots.knot_table import DATA, homfly, judge1, load_homology, mirror, superpoly
from rcsknots.laurent import load_table
from rcsknots.tables import normalize_sign

ROOT = Path(__file__).resolve().parent.parent
H = load_homology()
P = {name: superpoly(h) for name, h in H.items()}


def test_all_knots_present():
    assert len(H) == 801


def test_trefoil_is_dgr():
    """正のブレイドの三葉結び目 = DGR の 𝐚²𝐪⁻² + 𝐚²𝐪²𝐭² + 𝐚⁴𝐭³（第 1 期の T(2,3)）。"""
    assert P["3_1"] == {(2, -2, 0): 1, (2, 2, 2): 1, (4, 0, 3): 1}
    sol = json.loads((ROOT / "data/superpoly/T2_3.json").read_text())["solution"]
    assert P["3_1"] == {tuple(map(int, k.split(","))): v for k, v in sol.items()}


# ---- Euler 標数 = ブレイドから Sage で計算した HOMFLY-PT -----------------------------------------
L, M, A, Qs = sp.symbols("L M A Q")


def sage_to_dgr(text):
    """Sage (L P+ + L^{-1} P- + M P0 = 0) の値を DGR の HOMFLY-PT {(e_𝐚, e_𝐪): c} にする。

    tests/test_stage4_cable_homfly.py の sage_to_kn と同じ変換（L = i𝐚, M = −i(𝐪 − 𝐪^{-1}) の後、鏡像 (𝐚, 𝐪) → (𝐚^{-1}, 𝐪^{-1})）。
    """
    e = sp.sympify(text, locals={"L": L, "M": M})
    e = e.subs({L: sp.I / A, M: -sp.I * (1 / Qs - Qs)}, simultaneous=True)
    poly = sp.Poly(sp.expand(e * A ** 60 * Qs ** 60), A, Qs)
    return {(i - 60, j - 60): int(c) for (i, j), c in poly.terms()}


SAGE = dict(line.rstrip("\n").split("\t") for line in (DATA / "sage_homfly.txt").open())


def test_euler_characteristic_is_braid_homfly():
    bad = [name for name in P if homfly(P[name]) != sage_to_dgr(SAGE[name])]
    assert bad == []


# ---- DGR の対称性（Gorsky–Hogancamp–Mellit arXiv:2103.01212）-------------------------------------
def test_dgr_symmetry():
    """(e_𝐚, e_𝐪, e_𝐭) ↔ (e_𝐚, −e_𝐪, e_𝐭 − e_𝐪)（第 1 期の superpolynomial の構成と同じ形）。"""
    bad = [name for name, p in P.items() if {(a, -q, t - q): c for (a, q, t), c in p.items()} != p]
    assert bad == []


# ---- トーラス結び目 = 第 1 期の rCS_□（KN のノートブックと一致済み）------------------------------
TORUS = {"3_1": "T2_3", "5_1": "T2_5", "7_1": "T2_7", "9_1": "T2_9", "8_19": "T3_4"}


@pytest.mark.parametrize("name", list(TORUS))
def test_torus_knots_match_rcs(name):
    rcs = load_table(ROOT / "data/reference/notebook_torus.txt")[(TORUS[name], "1")]
    b = to_bold(rcs, next(iter(rcs))[0] % 2)
    assert P[name] in (b, {k: -v for k, v in b.items()})


# ---- KN 付録 C の [1] の表 ------------------------------------------------------------------------
PAPER = json.load(open(ROOT / "data/reference/kn_tables.json"))


def kn_table(knot):
    d = next(x for x in PAPER if x["knot"] == knot and x["rho"] == ["1"])["data"]
    return {tuple(map(int, k.split(","))): v for k, v in d.items()}


def bps1(p):
    return normalize_sign(judge1(p)["after"])[1]


@pytest.mark.parametrize("name,knot,chirality", [
    ("5_1", "T2_5", "K"), ("4_1", "4_1", "K"), ("5_2", "5_2", "K"), ("6_3", "6_3", "K"),
    ("6_1", "6_1", "K*"),     # KN の 6_1 は KnotInfo の鏡像
])
def test_kn_tables_reproduced(name, knot, chirality):
    p = P[name] if chirality == "K" else mirror(P[name])
    assert bps1(p) == kn_table(knot)


def test_kn_trefoil_table_uses_the_other_a_shift():
    """KN の T(2,3) の表だけ a 次数シフトが逆（第 1 期で記録済み）。β と 2J_r が 1 ずつずれる。"""
    assert {(g, b + 1, F + 1): v for (g, b, F), v in bps1(P["3_1"]).items()} == kn_table("T2_3")


# ---- 鏡像で判定が変わらない（G は 𝐪 → 𝐪⁻¹, 𝐭 → 𝐭⁻¹ で不変）------------------------------------
def test_mirror_invariance():
    bad = []
    for name, p in P.items():
        r, s = judge1(p), judge1(mirror(p))
        if r["closes"] != s["closes"] or (r["closes"] and (
                sorted(r["after"].values()) != sorted(s["after"].values())
                or r["parity_aligned"] != s["parity_aligned"])):
            bad.append(name)
    assert bad == []
