"""ステージ T0：結び目表の HOMFLY-PT ホモロジー（yui-kr）の規約の確定。"""
import json
from pathlib import Path

import pytest
import sympy as sp

from rcsknots.bps import to_bold
from rcsknots.knot_table import DATA, homfly, judge1, load_homology, mirror, parity_uniform, superpoly
from rcsknots.laurent import load_table
from rcsknots.tables import normalize_sign

ROOT = Path(__file__).resolve().parents[2]
H = load_homology()
P = {name: superpoly(h) for name, h in H.items()}


def test_all_knots_present():
    assert len(H) == 801


def test_trefoil_is_dgr():
    """正のブレイドの三葉結び目 = DGR の 𝐚²𝐪⁻² + 𝐚²𝐪²𝐭² + 𝐚⁴𝐭³（第 1 期の T(2,3)）。"""
    assert P["3_1"] == {(2, -2, 0): 1, (2, 2, 2): 1, (4, 0, 3): 1}
    sol = json.loads((ROOT / "data/framed_knot/superpoly/T2_3.json").read_text())["solution"]
    assert P["3_1"] == {tuple(map(int, k.split(","))): v for k, v in sol.items()}


# ---- Euler 標数 = ブレイドから Sage で計算した HOMFLY-PT -----------------------------------------
L, M, A, Qs = sp.symbols("L M A Q")


def sage_to_dgr(text):
    """Sage (L P+ + L^{-1} P- + M P0 = 0) の値を DGR の HOMFLY-PT {(e_𝐚, e_𝐪): c} にする。

    tests/framed_knot/test_stage4_cable_homfly.py の sage_to_kn と同じ変換（L = i𝐚, M = −i(𝐪 − 𝐪^{-1}) の後、鏡像 (𝐚, 𝐪) → (𝐚^{-1}, 𝐪^{-1})）。
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
    rcs = load_table(ROOT / "data/framed_knot/reference/notebook_torus.txt")[(TORUS[name], "1")]
    b = to_bold(rcs, next(iter(rcs))[0] % 2)
    assert P[name] in (b, {k: -v for k, v in b.items()})


# ---- KN 付録 C の [1] の表 ------------------------------------------------------------------------
PAPER = json.load(open(ROOT / "data/framed_knot/reference/kn_tables.json"))


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


# ---- ステージ T1：[1] の正値性は「P ≥ 0 と DGR の対称性」から従う ---------------------------------
def test_lemma_chebyshev_expansion_is_positive():
    """y^k + y^{−k} を G = y + y^{−1} − s − s^{−1} の多項式に展開すると、係数は s の非負係数ローラン多項式（k ≤ 20）。"""
    y, s, G, w = sp.symbols("y s G w")
    for k in range(1, 21):
        e = sp.expand(sp.expand(2 * sp.chebyshevt(k, w / 2)).subs(w, G + s + 1 / s))
        for g in range(k + 1):
            c = sp.Poly(sp.expand(e.coeff(G, g) * s ** k), s)
            assert all(x >= 0 for x in c.coeffs())
            assert all((m[0] - k - (k - g)) % 2 == 0 for m in c.monoms())   # 𝐭 の指数 ≡ k − g（偶奇の条件の根拠）


def test_positive_symmetric_polynomials_give_positive_bps():
    """非負係数・偶数の 𝐪 指数・DGR の対称性をもつ任意の P で、G 展開は閉じて \\hat N ≥ 0（乱数で 300 例）。"""
    import random
    rng = random.Random(0)
    for _ in range(300):
        P = {}
        for _ in range(rng.randint(1, 8)):
            a, q, t = 2 * rng.randint(-3, 3), 2 * rng.randint(-4, 4), rng.randint(-6, 6)
            c = rng.randint(1, 3)
            for e in {(a, q, t), (a, -q, t - q)}:
                P[e] = P.get(e, 0) + c
        r = judge1(P)
        assert r["closes"] and r["after_neg"] == 0 and r["before_neg"] == 0
        assert r["parity_aligned"] == parity_uniform(P)


def test_knot_table_all_positive_and_parity_criterion():
    for name, p in P.items():
        r = judge1(p)
        assert r["closes"] and r["before_neg"] == 0 and r["after_neg"] == 0, name
        assert r["parity_aligned"] == parity_uniform(p), name


def test_parity_criterion_equals_chandler_gorsky():
    """Chandler–Gorsky（2209.13058 main.tex:334）の parity knot：HHH 全体で t/2（Nakagane–Sano の次数）の偶奇が一定。
    𝐚 次数の台に隙間が無ければ、本研究の条件（parity_uniform）と同じ。801 個すべてで台に隙間が無く、判定が一致する。"""
    for name, h in load_homology().items():
        p = superpoly(h)
        a_degs = sorted({a for a, q, t in p})
        assert all(y - x == 2 for x, y in zip(a_degs, a_degs[1:])), name
        cg = len({k // 2 % 2 for (i, j, k) in h}) == 1
        assert cg == (len({(t - a // 2) % 2 for a, q, t in p}) == 1) == parity_uniform(p), name
