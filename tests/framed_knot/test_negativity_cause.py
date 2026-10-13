"""非代数的なケーブルの負の係数の原因（scripts/framed_knot/negativity_cause.py）の検査。"""
from pathlib import Path

import pytest
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]


def braid_homfly(name):
    for line in (ROOT / "data/framed_knot/reference/braid_homfly.txt").read_text().splitlines():
        cols = line.split("\t")
        if cols[0] == name:
            return sp.sympify(cols[3])
    raise KeyError(name)


@pytest.mark.parametrize("name,obstructed", [("T2_3", False), ("K1", False), ("K2", False), ("K5c", False),
                                             ("K4", True), ("K5a", True)])
def test_uncolored_obstruction_from_braid_homfly(name, obstructed):
    """Independent of the refined pipeline: the z^0 part of the braid HOMFLY-PT (Sage, M = 0), with a = −L^{−2}
    (KN's a; P(a=1, z=0) = ±1), times (1 − a) gives the unrefined N_{□,0,β}.  Refined positivity with uniform parity
    requires these to alternate in sign; it fails for K4 and K5a."""
    L, M, a = sp.symbols("L M a")
    p0 = sp.expand(braid_homfly(name).subs(M, 0).subs(L, (-a) ** sp.Rational(-1, 2)))
    p0 = sp.expand(p0 * a ** 20)
    assert abs(p0.subs(a, 1)) == 1
    coeffs = [c for c in sp.Poly(sp.expand((1 - a) * p0), a).all_coeffs() if c != 0]
    deg = sp.Poly(sp.expand((1 - a) * p0), a).all_coeffs()
    signs = {(1 if c > 0 else -1) * (-1) ** i for i, c in enumerate(deg) if c != 0}
    assert (len(signs) > 1) == obstructed, coeffs


@pytest.mark.parametrize("name", ["T2_3", "K1", "K5b", "K5a", "K4", "K2", "K5c", "K6b"])
def test_constructed_superpolynomial_constraints(name):
    """data/framed_knot/superpoly の候補（scripts/framed_knot/construct_superpolynomial.py）：非負、𝐭 = −1 で HOMFLY-PT、DGR の対称性。"""
    import json
    import sys
    sys.path.insert(0, str(ROOT / "scripts/framed_knot"))
    from construct_superpolynomial import daha_bold, homfly
    path = ROOT / f"data/framed_knot/superpoly/{name}.json"
    if not path.exists():
        pytest.skip("not constructed")
    P = {tuple(map(int, k.split(","))): v for k, v in json.loads(path.read_text())["solution"].items()}
    assert all(v > 0 for v in P.values())
    assert homfly(P) == homfly(daha_bold(name))
    assert all(P.get((a, -q, t - q), 0) == v for (a, q, t), v in P.items())
    if name == "T2_3":
        assert P == {(2, -2, 0): 1, (2, 2, 2): 1, (4, 0, 3): 1}
