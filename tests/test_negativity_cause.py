"""非代数的なケーブルの負の係数の原因（scripts/negativity_cause.py）の検査。"""
from pathlib import Path

import pytest
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent


def braid_homfly(name):
    for line in (ROOT / "data/reference/braid_homfly.txt").read_text().splitlines():
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
