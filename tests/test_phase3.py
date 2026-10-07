"""第 3 期：3 本ケーブルによる |λ| = 3 の彩色 HOMFLY-PT、[2,1] の構造的な性質、U4 の 3 箱の結果。"""
import json
import sys
from pathlib import Path

import pytest
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from colored_from_cables import colored_p  # noqa: E402
from rcsknots.bps import g_expand, to_bold  # noqa: E402
from rcsknots.extract import A_SHIFT  # noqa: E402
from rcsknots.homfly import Q as HQ, a as Ha  # noqa: E402
from rcsknots.laurent import load_table  # noqa: E402

R = load_table(ROOT / "data/reference/notebook_torus.txt")
A, Q, T = sp.symbols("A Q T")


def ex(P):
    return sum(c * A ** i * Q ** j * T ** k for (i, j, k), c in P.items())


@pytest.mark.parametrize("name,knot", [("3_1", "T2_3"), pytest.param("4_1", "4_1", marks=pytest.mark.results)])
def test_three_cable_inversion(name, knot):
    """(3, s)-ケーブル（s = 1, 2, 4）の HOMFLY-PT から [3], [2,1], [1,1,1] が KN の rCS の t = q 極限と一致。"""
    C = colored_p(name, 3, ROOT / "data/reference/cable_homfly_p3.txt")
    for lam, c in [((3,), "3"), ((2, 1), "21"), ((1, 1, 1), "111")]:
        ref = sp.expand(ex(R[(knot, c)]).subs(T, Q)).subs({A: sp.sqrt(Ha), Q: HQ})
        assert sp.simplify(C[lam] - ref) == 0


@pytest.mark.parametrize("knot", ["T2_3", "T2_5", "T2_7", "4_1"])
def test_21_structural_identities(knot):
    """rCS_[2,1] のミラー対称性、q = ±1・t = ±1 の EGP、a = t² での U(2) への簡約（仮説）が厳密に成立。"""
    r = {c: ex(R[(knot, c)]) for c in ["1", "2", "11", "21"]}
    assert sp.expand(r["21"].subs({Q: 1 / T, T: 1 / Q}, simultaneous=True) - r["21"]) == 0
    for s in (1, -1):
        assert sp.expand((r["21"] - r["1"] * r["11"]).subs(Q, s)) == 0
        assert sp.expand((r["21"] - r["1"] * r["2"]).subs(T, s)) == 0
        assert sp.expand((r["21"] - r["1"] * r["11"]).subs(A, s * T ** 2)) == 0


# トーラス結び目は判定のコードの機能テスト（1 秒以下）。ケーブル（各 23〜35 万項）は結果の再検証
U4 = ["T3_4", "T3_5", "T2_11", "T3_7", "T4_5"] + [pytest.param(k, marks=pytest.mark.results) for k in ("K9", "K8", "K7")]


@pytest.mark.parametrize("knot", U4)
def test_u4_three_box_positivity(knot):
    """U4：|ρ| ≤ 3 の \\hat f が、G 展開の前後とも（同じ全体の符号で）非負整数。"""
    paths = [ROOT / f"data/lmov/{knot}_{n}.json" for n in (1, 2, 3)]
    if not all(p.exists() for p in paths):
        pytest.skip("run scripts/u4_run.py first")
    for n, p in zip((1, 2, 3), paths):
        for key, poly in json.loads(p.read_text())["f"].items():
            if not key.startswith("hat"):
                continue
            hat = {tuple(map(int, e.split(","))): c for e, c in poly.items()}
            bold = to_bold(hat, A_SHIFT["nontorus"][n])
            N = g_expand(bold)
            vals = list(bold.values()) + list(N.values())
            assert all(v >= 0 for v in vals) or all(v <= 0 for v in vals), (knot, key)
