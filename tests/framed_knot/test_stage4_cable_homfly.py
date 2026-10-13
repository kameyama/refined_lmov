"""Stage 4: unrefined cabling formula vs. HOMFLY-PT from cable braids (Sage), and colored consistency checks."""
from functools import lru_cache
from pathlib import Path

import pytest
import sympy as sp

from rcsknots.exact import reconstruct_points
from rcsknots.homfly import a, Q, cable_unreduced, unknot
from rcsknots.homfly_num import homfly_evaluator

ROOT = Path(__file__).resolve().parents[2]
L, M, Abold = sp.symbols("L M A")
# research-plan knots: (p1, q1, p2, q2) for the (p2, q2)-cable of T(p1, q1)
KNOTS = {"K1": (2, 3, 2, 13), "K2": (2, 3, 2, 1), "K3": (2, 3, 2, 15), "K4": (2, 3, 2, 3),
         "K5a": (2, 3, 2, 5), "K5b": (2, 3, 2, 7), "K5c": (2, 3, 2, -1), "K6a": (2, 5, 2, 21), "K6b": (2, 5, 2, 1)}


def load_braids():
    out = {}
    for line in open(ROOT / "data/framed_knot/reference/braid_homfly.txt"):
        name, n, word, hom = line.rstrip("\n").split("\t")
        out[name] = sp.sympify(hom, locals={"L": L, "M": M})
    return out


BRAID = load_braids()


def sage_to_kn(expr):
    """Sage (L P+ + L^{-1} P- + M P0 = 0) → DGR (a P+ − a^{-1} P- = (q − q^{-1}) P0): L = i𝐚, M = −i(𝐪 − 𝐪^{-1}).

    Sage's positive braids give positive knots; DGR's (and 1703.05408's) torus knots are negative, so we take
    the mirror (𝐚, 𝐪) → (𝐚^{-1}, 𝐪^{-1}).  Finally KN's a = 𝐚², Q = 𝐪.
    """
    e = sp.expand(expr.subs({L: sp.I * Abold, M: -sp.I * (Q - 1 / Q)}, simultaneous=True))
    e = e.subs({Abold: 1 / Abold, Q: 1 / Q}, simultaneous=True)
    return sp.expand(e.subs(Abold, sp.sqrt(a)))


@pytest.mark.parametrize("q2", [1, 3, 5, -3])
def test_framing_correction_on_kinked_unknot(q2):
    """The unknot drawn as σ1 (writhe 1), 2-cabled with the correction (σ1)^{q2 − 2w}, is T(2, q2)."""
    assert sp.expand(BRAID[f"check_kinkunknot_2_{q2}"] - BRAID[f"check_T2_{q2}"]) == 0


def test_p_equal_one_gives_companion():
    assert sp.expand(BRAID["check_p1_trefoil"] - BRAID["T2_3"]) == 0


def test_sage_trefoil_is_dgr():
    assert sp.expand(sage_to_kn(BRAID["T2_3"]) - (a / Q**2 + a * Q**2 - a**2)) == 0


@lru_cache(maxsize=None)
def companion(p1, q1):
    @lru_cache(maxsize=None)
    def W(mu):
        return cable_unreduced(mu, p1, q1, unknot, sign=-1)
    return W


@lru_cache(maxsize=None)
def cable_reduced(lam, p1, q1, p2, q2):
    val = cable_unreduced(lam, p2, q2, companion(p1, q1), sign=-1) / unknot(lam)
    return sp.factor(sp.simplify(val))


@pytest.mark.parametrize("name", sorted(KNOTS))
def test_cable_formula_matches_braid(name):
    """Uncolored: the cabling formula (θ_λ = q^{κ}a^{|λ|/2}, framing correction θ^{−p q}) = braid HOMFLY-PT."""
    assert sp.simplify(cable_reduced((1,), *KNOTS[name]) - sage_to_kn(BRAID[name])) == 0


@lru_cache(maxsize=None)
def cable_exact(lam, name):
    """Exact reduced colored HOMFLY-PT {(e_A, e_Q): c} of the cable, from the numerical cabling formula + interpolation."""
    p1, q1, p2, q2 = KNOTS[name]
    raw = reconstruct_points(homfly_evaluator(lam, [(p1, q1), (p2, q2)]), log=lambda s: None)
    return {(ea, eq): c for (_, ea, eq), c in raw.items()}


@pytest.mark.parametrize("name", ["K1", "K2", "K3", "K4", "K6a", "K6b"])
def test_numeric_formula_matches_braid(name):
    """The numerical cabling formula reproduces the braid HOMFLY-PT (uncolored)."""
    expr = sum(c * sp.sqrt(a) ** ea * Q ** eq for (ea, eq), c in cable_exact((1,), name).items())
    assert sp.expand(expr - sage_to_kn(BRAID[name])) == 0


@pytest.mark.parametrize("name", ["K1", "K2", "K3", "K4", "K5a", "K5b", "K5c", "K6a", "K6b"])
def test_colored_transposition_symmetry(name):
    """Reduced colored HOMFLY-PT of a knot in framing 0: H_{λ^T}(a, q) = H_λ(a, q^{-1})."""
    for lam, lamT in [((2,), (1, 1)), ((3,), (1, 1, 1)), ((2, 1), (2, 1))]:
        h, hT = cable_exact(lam, name), cable_exact(lamT, name)
        assert {(ea, -eq): c for (ea, eq), c in h.items()} == hT


@pytest.mark.parametrize("name", ["K1", "K2", "K3", "K4", "K5a", "K5b", "K5c", "K6a", "K6b"])
def test_special_polynomial(name):
    """H_λ(a, q=1) = H_□(a, q=1)^{|λ|}."""
    def at_q1(h):
        out = {}
        for (ea, eq), c in h.items():
            out[ea] = out.get(ea, 0) + c
        return sp.expand(sum(c * Abold ** ea for ea, c in out.items()))
    h1 = at_q1(cable_exact((1,), name))
    for lam in [(2,), (1, 1), (3,), (2, 1), (1, 1, 1)]:
        assert sp.expand(at_q1(cable_exact(lam, name)) - h1 ** sum(lam)) == 0
