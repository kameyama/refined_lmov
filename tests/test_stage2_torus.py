"""Stage 2: torus knots from the elliptic Hall algebra (GN) vs. the 1703.05408 notebook."""
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest
import sympy as sp

from rcsknots.field import V, P31, make_conv
from rcsknots.gn import GN
from rcsknots.homfly import torus_reduced, a as ha, Q as hQ
from rcsknots.laurent import load_table
from rcsknots.macdonald import principal_specialization
from rcsknots.modified import solve_in_basis, multiply_p, phi_factor
from rcsknots.partitions import partitions

ROOT = Path(__file__).resolve().parent.parent
REF = load_table(ROOT / "data/reference/notebook_torus.txt")
COLORS = {"1": (1,), "2": (2,), "11": (1, 1), "3": (3,), "111": (1, 1, 1), "21": (2, 1)}
KNOTS = {"T2_3": (2, 3), "T2_5": (2, 5), "T2_7": (2, 7), "T2_9": (2, 9), "T3_4": (3, 4)}


@pytest.fixture(scope="module")
def G():
    rng = np.random.default_rng(42)
    q, t, a = (V(rng.integers(2, P31, 4), P31) for _ in range(3))
    return GN(q, t, a, make_conv(P31, 4))


def add(u, v):
    out = dict(u)
    for k, x in v.items():
        out[k] = out[k] + x if k in out else x
    return out


@pytest.mark.parametrize("k", [1, 2, 3])
def test_p_k0_is_multiplication(G, k):
    """P̃_{k,0} = multiplication by φ(p_k) = p_k/(1 − t_w^{−k}) in the H̃ basis (fixes the SYT orientation and c_k)."""
    for n in range(0, 6 - k):
        for lam in partitions(n):
            got = G.op(k, 1, 0, {lam: G.one})
            target = multiply_p({(k,): phi_factor((k,), G.tw) * G.one}, G.H(n)[lam])
            exp = solve_in_basis(target, G.H(n + k), G.zero)
            for mu in partitions(n + k):
                assert got.get(mu, G.zero) == exp.get(mu, G.zero)


@pytest.mark.parametrize("x,y,s", [((1, 2), (1, 1), (2, 3)), ((1, 1), (1, 0), (2, 1)), ((1, 3), (1, 2), (2, 5))])
def test_commutation_relation(G, x, y, s):
    """[P_x, P_y] = κ P_{x+y} for quasi-empty triangles (GN Fig. 1, Burban–Schiffmann).

    With our normalization (P̃_{k,0} = multiplication by φ(p_k)) the constant is κ = (1 − q)/t for every
    triangle; GN normalize the generators so that κ = 1 (they omit the central element, GN §3.3).
    A slope-independent κ is what SL(2,Z)-covariance requires.
    """
    kappa = (1 - G.q) * G.tw
    for n in range(0, 3):
        for lam in partitions(n):
            v = {lam: G.one}
            lhs = add(G.op(1, *x, G.op(1, *y, v)), {k: -c for k, c in G.op(1, *y, G.op(1, *x, v)).items()})
            rhs = G.op(1, *s, v)
            for mu in set(lhs) | set(rhs):
                assert lhs.get(mu, G.zero) == kappa * rhs.get(mu, G.zero)


@pytest.mark.parametrize("lam", [(1,), (2,), (1, 1), (2, 1), (3,)])
def test_unknot(G, lam):
    """Slope (1,0): ε(P_λ · 1) is the principal specialization."""
    assert G.torus(lam, 1, 0) == principal_specialization(lam, G.a, G.q, G.t)


@pytest.mark.parametrize("key", sorted(REF))
def test_notebook_t_equals_q_is_rosso_jones(key):
    """t = q limit of the notebook values = framing-0 colored HOMFLY-PT from the Rosso–Jones formula (DGR orientation)."""
    knot, cname = key
    nb = sum(c * sp.sqrt(ha) ** ea * hQ ** (eq + et) for (ea, eq, et), c in REF[key].items())
    assert sp.expand(nb - torus_reduced(COLORS[cname], *KNOTS[knot], sign=-1)) == 0


def stored(knot, cname):
    path = ROOT / f"data/torus/{knot}_{cname}.json"
    if not path.exists():
        pytest.skip("run scripts/stage2_torus.py first")
    return {tuple(map(int, k.split(","))): v for k, v in json.loads(path.read_text())["poly"].items()}


@pytest.mark.parametrize("key", sorted(REF))
def test_gn_reproduces_notebook(key):
    """GN construction + normalization (N1)(N2) reproduces every coefficient of the notebook values."""
    assert stored(*key) == REF[key]


@pytest.mark.parametrize("knot", sorted(KNOTS))
@pytest.mark.parametrize("pair", [("2", "11"), ("3", "111"), ("1", "1"), ("21", "21")])
def test_mirror_transposition(knot, pair):
    """rCS_{λ^T}(a,q,t) = rCS_λ(a,t^{-1},q^{-1}) (KN eq. (mirror-sym)): (eA,eQ,eT) ↦ (eA,−eT,−eQ)."""
    if (knot, pair[0]) not in REF or (knot, pair[1]) not in REF:
        pytest.skip("color not in notebook")
    lam, lamT = stored(knot, pair[0]), stored(knot, pair[1])
    assert {(ea, -et, -eq): c for (ea, eq, et), c in lam.items()} == lamT
