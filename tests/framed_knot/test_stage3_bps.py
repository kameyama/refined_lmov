"""Stage 3: refined BPS extraction (M matrix, f^q, f^t, \\hat f, g-expansion) vs. 1703.05408."""
import json
import random
from pathlib import Path

import pytest
import sympy as sp

from rcsknots.bps import f_formulas, C, C1, Q, T, M_matrix
from rcsknots.extract import hat_f_exact, bps_numbers
from rcsknots.laurent import load_table
from rcsknots.tables import normalize_sign

ROOT = Path(__file__).resolve().parents[2]
q, t = Q**2, T**2
F = f_formulas(3)

# ---- 1703.05408 App. B (verbatim transcription) --------------------------------------------------
pre=T-1/T
paper_q={
 (1,): pre*C((1,)),
 (2,): pre*sp.sqrt(q/t)*((q*t-1)/(q**2-1)*C((2,)) - (t-1)/(2*(q-1))*C((1,))**2 - (t+1)/(2*(q+1))*C1(2)),
 (1,1): pre*sp.sqrt(q/t)*((t-q)/(q**2-1)*C((2,)) + (t**2-1)/(q*t-1)*C((1,1)) - (t-1)/(2*(q-1))*C((1,))**2 + (t+1)/(2*(q+1))*C1(2)),
 (3,): pre*(q/t)*((q*t-1)*(q**2*t-1)/((q**2-1)*(q**3-1))*C((3,)) - (t-1)*(q*t-1)/((q-1)**2*(q+1))*C((2,))*C((1,)) + (t-1)**2/(3*(q-1)**2)*C((1,))**3 - (t**2+t+1)/(3*(q**2+q+1))*C1(3)),
 (2,1): pre*(q/t)*(-(q-t)*(q*t-1)/((q-1)*(q**3-1))*C((3,)) + (t-1)*(q*t**2-1)/((q-1)*(q**2*t-1))*C((2,1)) - ((t-1)**2/(q-1)**2*C((2,))+(t-1)**2*(t+1)/((q-1)*(q*t-1))*C((1,1)))*C((1,)) + 2*(t-1)**2/(3*(q-1)**2)*C((1,))**3 + (t**2+t+1)/(3*(q**2+q+1))*C1(3)),
 (1,1,1): pre*(q/t)*((q-t)*(q**2-t)/((q**2-1)*(q**3-1))*C((3,)) - (t**2-1)*(q-t)/((q-1)*(q**2*t-1))*C((2,1)) + (t**2-1)*(t**3-1)/((q*t-1)*(q*t**2-1))*C((1,1,1)) + ((t-1)*(q-t)/((q-1)**2*(q+1))*C((2,)) - (t-1)**2*(t+1)/((q-1)*(q*t-1))*C((1,1)))*C((1,)) + (t-1)**2/(3*(q-1)**2)*C((1,))**3 - (t**2+t+1)/(3*(q**2+q+1))*C1(3)),
}
paper_t={
 (1,): pre*C((1,)),
 (2,): -pre*(C((1,1)) + C1(2)/2 - C((1,))**2/2),
 (1,1): -pre*(C((2,)) + (q-t)/(q*t-1)*C((1,1)) - C((1,))**2/2 - C1(2)/2),
 (3,): pre*(C((1,1,1)) - C((1,))*C((1,1)) + C((1,))**3/3 - C1(3)/3),
 (2,1): pre*(C((2,1)) + (t+1)*(q-t)/(q*t**2-1)*C((1,1,1)) - (C((2,)) + (q-1)*(t+1)/(q*t-1)*C((1,1)))*C((1,)) + 2*C((1,))**3/3 + C1(3)/3),
 (1,1,1): pre*(C((3,)) + (q+1)*(q-t)/(q**2*t-1)*C((2,1)) + (q-t)*(q-t**2)/((q*t-1)*(q*t**2-1))*C((1,1,1)) + (-C((2,)) + (q-t)/(q*t-1)*C((1,1)))*C((1,)) + C((1,))**3/3 - C1(3)/3),
}

# notebook `ft111` (refined-LMOV-non-torus.nb); the printed App. B formula for f^t_[1,1,1] differs (typo)
C1_,C2,C11,C3,C21,C111,C13=C((1,)),C((2,)),C((1,1)),C((3,)),C((2,1)),C((1,1,1)),C1(3)
ft111_nb=((-1+t)*((-1+q**4*t**4+q*t*(1+t)-q**3*t**2*(1+t)+q**2*(t-t**3))*C1_**3+(1-q**4*t**4-q*t*(1+t)+q**3*t**2*(1+t)+q**2*t*(-1+t**2))*C13-3*(1-q**2*t-q*t**2+q**3*t**3)*C1_*((-1+q*t)*C2+(q-t)*C11)+3*((-1+q**4*t**4+q*t*(1+t)-q**3*t**2*(1+t)+q**2*(t-t**3))*C3+(q-t)*((1+q**3*t**3+q**2*t*(-1-t+t**2)-q*(-1+t+t**2))*C21+(-q+q**3*t+t**2-q**2*t**3)*C111))))/(3*T*(-1+q*t)*(-1+q**2*t)*(-1+q*t**2))


VALS = {Q: sp.Rational(3, 2), T: sp.Rational(5, 7)}
random.seed(1)
for _s in [C((1,)), C((2,)), C((1, 1)), C((3,)), C((2, 1)), C((1, 1, 1)), C1(2), C1(3)]:
    VALS[_s] = sp.Rational(random.randint(1, 50), random.randint(1, 50))


@pytest.mark.parametrize("mu", list(paper_q))
def test_fq_matches_appendix_b(mu):
    assert F[mu][0].subs(VALS) == paper_q[mu].subs(VALS)


@pytest.mark.parametrize("mu", [(1,), (2,), (1, 1), (3,), (2, 1)])
def test_ft_matches_appendix_b(mu):
    assert F[mu][1].subs(VALS) == paper_t[mu].subs(VALS)


def test_ft111_matches_notebook():
    assert F[(1, 1, 1)][1].subs(VALS) == ft111_nb.subs(VALS)
    assert F[(1, 1, 1)][1].subs(VALS) != paper_t[(1, 1, 1)].subs(VALS)   # misprint in App. B


def test_M_matrices_match_notebook():
    assert sp.simplify(M_matrix(2) - sp.Matrix([[1 / T, -T], [-T, 1 / T]])) == sp.zeros(2)
    assert sp.simplify(M_matrix(3) - sp.Matrix([[1 / t, -1, t], [-1, -1 + 1 / t + t, -1], [t, -1, 1 / t]])) == sp.zeros(3)


# ---- tables -------------------------------------------------------------------------------------
PAPER = json.load(open(ROOT / "data/framed_knot/reference/kn_tables.json"))
CASES = [("T2_3", r) for r in ["1", "2", "11", "3", "21", "111"]] + \
        [("4_1", r) for r in ["1", "2", "11", "3", "21", "111"]] + [("T2_5", "1")]


def paper_table(knot, rho):
    key = [",".join(rho)]
    return next(x for x in PAPER if x["knot"] == knot and x["rho"] == key)["data"]


@pytest.mark.parametrize("knot,rho", CASES)
def test_saved_tables_match_paper(knot, rho):
    path = ROOT / f"data/framed_knot/bps/{knot}_{rho}.json"
    if not path.exists():
        pytest.skip("run scripts/framed_knot/stage3_bps.py first")
    assert json.loads(path.read_text())["N"] == paper_table(knot, rho)


def test_pipeline_recomputes_trefoil():
    """Fresh run of the whole pipeline for T(2,3), |ρ| ≤ 2 (torus-notebook a-shift convention)."""
    ref = load_table(ROOT / "data/framed_knot/reference/notebook_torus.txt")
    colors = {"1": (1,), "2": (2,), "11": (1, 1)}
    reduced = {colors[c]: p for (k, c), p in ref.items() if k == "T2_3" and c in colors}
    for n in (1, 2):
        for rho, hf in hat_f_exact(reduced, n, log=lambda s: None).items():
            _, N = normalize_sign(bps_numbers(hf, rho, "torus"))
            assert {f"{g},{b},{F_}": v for (g, b, F_), v in N.items()} == paper_table("T2_3", "".join(map(str, rho)))
