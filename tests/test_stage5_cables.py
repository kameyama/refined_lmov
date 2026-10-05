"""Stage 5: refined invariants of cables (GN iteration) vs. DAHA superpolynomials (CD14), and consistency checks."""
import json
from pathlib import Path

import pytest
import sympy as sp

from rcsknots.daha_compare import daha_to_kn, monomial_ratio, q, t, a
from rcsknots.exact import reconstruct_points
from rcsknots.homfly_num import homfly_evaluator

ROOT = Path(__file__).resolve().parent.parent
EXAMPLES = [l.rstrip("\n").split("\t") for l in open(ROOT / "data/reference/daha_examples.txt")]
KNOTS = {"K1": (2, 3, 2, 13), "K2": (2, 3, 2, 1), "K3": (2, 3, 2, 15), "K4": (2, 3, 2, 3),
         "K5a": (2, 3, 2, 5), "K5b": (2, 3, 2, 7), "K5c": (2, 3, 2, -1), "K6a": (2, 5, 2, 21), "K6b": (2, 5, 2, 1)}


def stored(knot, color):
    path = ROOT / f"data/cables/{knot}_{color}.json"
    if not path.exists():
        pytest.skip("run scripts/stage5_cables.py first")
    return {tuple(map(int, k.split(","))): v for k, v in json.loads(path.read_text())["poly"].items()}


@pytest.mark.parametrize("name,knot,color,expr", EXAMPLES, ids=[e[0] for e in EXAMPLES])
def test_matches_daha_example(name, knot, color, expr):
    """Every coefficient agrees with the CD14 superpolynomial up to an overall monomial (a_DAHA = −a_KN)."""
    theirs = daha_to_kn(sp.sympify(expr, locals={"q": q, "t": t, "a": a}))
    assert monomial_ratio(stored(knot, color), theirs) is not None


@pytest.mark.parametrize("knot", ["K1", "K2", "K3", "K4", "K5a", "K5b", "K5c", "K6a", "K6b"])
@pytest.mark.parametrize("pair", [("2", "11"), ("3", "111"), ("1", "1"), ("21", "21")])
def test_mirror_transposition(knot, pair):
    """rCS_{λ^T}(a,q,t) = rCS_λ(a,t^{-1},q^{-1}) exactly, with (N1)(N2) imposed separately on λ and λ^T."""
    lam, lamT = stored(knot, pair[0]), stored(knot, pair[1])
    assert {(ea, -et, -eq): c for (ea, eq, et), c in lam.items()} == lamT


@pytest.mark.parametrize("knot", ["K1", "K2", "K4", "K5c"])
@pytest.mark.parametrize("color,lam", [("1", (1,)), ("2", (2,)), ("21", (2, 1))])
def test_t_equals_q_is_homfly(knot, color, lam):
    """t = q limit = framing-0 colored HOMFLY-PT from the cabling formula (independent of the normalization step)."""
    poly = stored(knot, color)
    p1, q1, p2, q2 = KNOTS[knot]
    hom = {(ea, eq): c for (_, ea, eq), c in
           reconstruct_points(homfly_evaluator(lam, [(p1, q1), (p2, q2)]), log=lambda s: None).items()}
    lim = {}
    for (ea, eq, et), c in poly.items():
        lim[(ea, eq + et)] = lim.get((ea, eq + et), 0) + c
    assert {k: v for k, v in lim.items() if v} == hom


@pytest.mark.parametrize("knot", ["K1", "K2", "K3", "K4", "K5a", "K5b", "K5c", "K6a", "K6b"])
def test_min_bold_t_degree_zero(knot):
    """(N2) holds: minimal 𝐭-degree 0 in the bold variables (e_𝐭 = e_a + 2 e_q, here in half units)."""
    poly = stored(knot, "1")
    assert min((ea + 2 * eq) // 2 for (ea, eq, et) in poly) == 0


@pytest.mark.parametrize("knot", ["K1", "K2", "K3", "K4", "K5a", "K5b", "K5c", "K6a", "K6b"])
def test_normalization_N2prime(knot):
    """(N2'): M_λ = M_□^{|λ|} T_λ^{p₂q₂} for every color (docs/conventions.md §5)."""
    from fractions import Fraction
    from rcsknots.partitions import n_stat, conjugate
    p1, q1, p2, q2 = KNOTS[knot]
    path1 = ROOT / f"data/cables/{knot}_1.json"
    if not path1.exists():
        pytest.skip("not computed")
    m1 = json.loads(path1.read_text())["monomial_applied"]
    for cname, lam in {"2": (2,), "11": (1, 1), "3": (3,), "21": (2, 1), "111": (1, 1, 1)}.items():
        path = ROOT / f"data/cables/{knot}_{cname}.json"
        if not path.exists():
            continue
        m = json.loads(path.read_text())["monomial_applied"]
        k = sum(lam)
        assert Fraction(m["q"]) == Fraction(m1["q"]) * k + p2 * q2 * n_stat(conjugate(lam))
        assert Fraction(m["t"]) == Fraction(m1["t"]) * k - p2 * q2 * n_stat(lam)
        assert Fraction(m["a"]) == Fraction(m1["a"]) * k and m["sign"] == m1["sign"] ** k


@pytest.mark.parametrize("knot", ["K1", "K2"])
def test_bps_recompute_matches_saved(knot):
    """Stage 6 pipeline: \\hat N_{[1]} recomputed from the stored refined CS invariant equals the saved table."""
    from rcsknots.extract import hat_f_exact, bps_numbers
    from rcsknots.tables import normalize_sign
    path = ROOT / f"data/bps/{knot}_1.json"
    if not path.exists():
        pytest.skip("run scripts/stage6_run.py first")
    hf = hat_f_exact({(1,): stored(knot, "1")}, 1, log=lambda s: None)[(1,)]
    _, N = normalize_sign(bps_numbers(hf, (1,), "nontorus"))
    assert {f"{g},{b},{F}": v for (g, b, F), v in N.items()} == json.loads(path.read_text())["N"]


@pytest.mark.parametrize("knot", ["T2_3", "T2_5", "4_1", "K1", "K3", "K6a"])
def test_positive_before_and_after_G(knot):
    """Torus knots, 4_1 (KN) and the algebraic cables: \\hat f (before G) and \\hat N (after G) are non-negative
    with one common overall sign, for every |ρ| ≤ 3."""
    from rcsknots.positivity import judge
    for r in ["1", "2", "11", "3", "21", "111"]:
        path = ROOT / f"data/bps/{knot}_{r}.json"
        if not path.exists():
            pytest.skip("not computed")
        j = judge(path)
        assert j["both_positive_same_sign"] and j["integral"]


def test_K5c_before_and_after_G_differ():
    """K5c [1]: \\hat N (after G) is non-negative, while the coefficients before G are not."""
    from rcsknots.positivity import judge
    path = ROOT / "data/bps/K5c_1.json"
    if not path.exists():
        pytest.skip("not computed")
    j = judge(path)
    assert j["after_G"]["min_negative"] == 0 and j["before_G"]["min_negative"] == 5
