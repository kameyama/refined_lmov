"""見直し：refined CS 不変量の小文字 (a,q,t) での一貫性と、太字への変換の検証。

太字への変換は refined LMOV を読み取る直前（src/rcsknots/lmov.py）でだけ行う。ここでは
(1) 新しく計算した rCS が、全体の因子を除いてトーラス結び目（KN）と同じ流儀になっていることを小文字で確かめ、
(2) 変換の規約を KN 本文に印刷された例で確かめる。
"""
import glob
import json
from pathlib import Path

import pytest
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
COLORS = ["1", "2", "11", "3", "21", "111"]
CABLES = ["K1", "K2", "K3", "K4", "K5a", "K5b", "K5c", "K6a", "K6b"]
TORUS = ["T2_3", "T2_5", "T2_7", "T2_9"]
A, Q, T = sp.symbols("A Q T")


def ld(path):
    return {tuple(map(int, k.split(","))): v for k, v in json.loads(Path(path).read_text())["poly"].items()}


def files():
    return sorted(glob.glob(str(ROOT / "data/torus/*.json"))) + sorted(glob.glob(str(ROOT / "data/cables/*.json")))


@pytest.mark.parametrize("path", files(), ids=lambda p: Path(p).stem)
def test_uniform_parity(path):
    """rCS = (monomial) × (polynomial in a, q, t with integer exponents): parities of (e_A, e_Q, e_T) are uniform.
    Hence the sign convention of q^{1/2}, t^{1/2}, a^{1/2} only affects the overall sign."""
    P = ld(path)
    assert len({(ea % 2, eq % 2, et % 2) for (ea, eq, et) in P}) == 1


def value(P, q1=False, t1=False):
    return sp.expand(sum(c * A ** ea * (1 if q1 else Q ** eq) * (1 if t1 else T ** et) for (ea, eq, et), c in P.items()))


@pytest.mark.parametrize("prefix", [f"torus/{k}" for k in TORUS] + [f"cables/{k}" for k in CABLES])
def test_exponential_growth(prefix):
    """KN eq. (EGP) and its transpose, in lowercase:
    rCS_{Σ b_i ω_i}(q=1) = Π rCS_{ω_i}(q=1)^{b_i}  and  rCS_λ(t=1) = Π over columns (transposed statement)."""
    R = {c: ld(ROOT / f"data/{prefix}_{c}.json") for c in COLORS}
    q1 = {c: value(R[c], q1=True) for c in COLORS}
    t1 = {c: value(R[c], t1=True) for c in COLORS}
    assert sp.expand(q1["2"] - q1["1"] ** 2) == 0
    assert sp.expand(q1["21"] - q1["1"] * q1["11"]) == 0
    assert sp.expand(q1["3"] - q1["1"] ** 3) == 0
    assert sp.expand(t1["11"] - t1["1"] ** 2) == 0
    assert sp.expand(t1["21"] - t1["1"] * t1["2"]) == 0
    assert sp.expand(t1["111"] - t1["1"] ** 3) == 0


@pytest.mark.parametrize("q2,name", [(3, "T2_3")])
def test_degenerate_cable_reproduces_torus(q2, name):
    """The cable pipeline (GN iteration + (N1)(N2')) with the unknot as companion reproduces KN's T(2,q2)."""
    from rcsknots.exact import reconstruct, reconstruct_points
    from rcsknots.homfly_num import homfly_evaluator
    from rcsknots.laurent import load_table
    from rcsknots.normalize import normalize_exact
    from rcsknots.partitions import n_stat, conjugate
    ref = load_table(ROOT / "data/reference/notebook_torus.txt")
    lams = {"1": (1,), "2": (2,), "11": (1, 1), "3": (3,), "21": (2, 1), "111": (1, 1, 1)}
    F, m1 = 2 * q2, None
    for c, lam in lams.items():
        raw = reconstruct(lambda G: G.cable_reduced_apoly(lam, (2, -q2), (1, 0)), log=lambda s: None)
        hom = {(ea, eq): v for (_, ea, eq), v in
               reconstruct_points(homfly_evaluator(lam, [(1, 0), (2, q2)]), log=lambda s: None).items()}
        norm, mono = normalize_exact(raw, hom)
        m1 = mono if lam == (1,) else m1
        k = sum(lam)
        assert norm == ref[(name, c)]
        assert tuple(mono) == (m1[0] ** k, m1[1] * k, m1[2] * k + F * n_stat(conjugate(lam)), m1[3] * k - F * n_stat(lam))


def kn_printed():
    a, q, t = sp.symbols("a q t")
    return (a, q, t), {
        ("hat", (1,)): a**2 / q**2 * (a**2 * t + 1) * (a**2 * q**2 * t**3 + q**4 * t**2 + 1),
        ("hat", (2,)): a**2 / q**4 * (a**2 * t + 1) * (a**2 * t**3 + 1) * (q**4 * t**2 + 1) * (a**2 * t + q**2) * (a**2 * q**2 * t**3 + 1),
        ("hat", (1, 1)): a**2 / q**6 * (a**2 * t + 1) * (a**2 * t + q**2) * (a**2 * q**2 * t**3 + 1)
        * (q**8 * (a**2 * t**7 + t**4) + q**4 * t**2 * (a**2 * (t**3 + t) + 2) + a**2 * t**3 + 1),
        ("fq", (2,)): -a**2 * t**2 / q * (a**2 * t + 1) * (q**8 * (a**4 * t**8 + a**2 * t**5) + q**6 * (a**6 * t**9 + a**4 * t**6 + a**2 * t**5 + t**2)
                                                     + a**2 * q**4 * t**3 * (a**2 * (t**3 + t) + 2) + q**2 * (a**6 * t**5 + a**4 * t**4 + a**2 * t + 1) + a**2 * t * (a**2 * t + 1)),
        ("ft", (2,)): -a**2 / (q**7 * t) * (a**2 * t + 1) * (a**2 * q**8 * t**5 * (a**2 * t + 1) + q**6 * t**2 * (a**6 * t**5 + a**4 * t**4 + a**2 * t + 1)
                                                        + a**2 * q**4 * t**3 * (a**2 * (t**3 + t) + 2) + q**2 * (a**6 * t**7 + a**4 * t**4 + a**2 * t**3 + 1) + a**4 * t**4 + a**2 * t)}


@pytest.mark.parametrize("key", list(kn_printed()[1]))
def test_change_of_variables_matches_kn_examples(key):
    """KN §4 examples for T(2,3): our lowercase f's, changed to bold variables just before reading them off,
    reproduce the printed polynomials.  KN's printed examples use q^{1/2} = +𝐪𝐭 (torus notebook `knotparam`)
    and the T(2,3)-table a-shift; with (cov) q^{1/2} = −𝐪𝐭 they agree up to an overall sign."""
    from rcsknots.lmov import load_lmov
    from rcsknots.extract import A_SHIFT
    (a, q, t), printed = kn_printed()
    path = ROOT / "data/lmov/T2_3_1.json"
    if not path.exists():
        pytest.skip("run scripts/stage8_lmov.py T2_3")
    F = load_lmov("T2_3")
    shift = A_SHIFT["torus"][sum(key[1])]
    for sq, expected_ratio in ((1, 1), (-1, None)):
        expr = 0
        for (eA, eQ, eT), c in F[key].items():
            eA += shift
            expr += c * (-a**2 * t) ** (eA // 2) * (sq * q * t) ** eQ * q ** eT
        ratio = sp.simplify(sp.expand(expr) / sp.expand(printed[key]))
        if expected_ratio is not None:
            assert ratio == expected_ratio
        else:
            assert ratio in (1, -1)


@pytest.mark.parametrize("name", ["T2_3", "T2_5", "4_1"] + CABLES)
def test_lmov_hat_f_matches_stage6(name):
    """\\hat f recomputed together with f^q, f^t (stage 8) equals the Stage-6 value."""
    for n in (1, 2, 3):
        path = ROOT / f"data/lmov/{name}_{n}.json"
        if not path.exists():
            pytest.skip("not computed yet")
        assert json.loads(path.read_text())["hat_f_matches_stage6"] is True


@pytest.mark.parametrize("name", ["T2_3", "T2_5", "K1", "K3"])
def test_N_q_t_symmetric_positive(name):
    """N^q, N^{\\bar t}: non-negative integers for symmetric [r] (KN §4 for torus knots; also the algebraic cables),
    negative coefficients for [1,1] (also for torus knots).  \\hat N: non-negative before and after G."""
    from rcsknots.lmov import judge
    if not (ROOT / f"data/lmov/{name}_3.json").exists():
        pytest.skip("not computed")
    import os
    cwd = os.getcwd()
    os.chdir(ROOT)
    try:
        J = judge(name)
    finally:
        os.chdir(cwd)
    for rho in [(1,), (2,), (3,)]:
        for kind in ("fq", "ft"):
            assert J[rho][kind]["min_negative"] == 0 and J[rho][kind]["integral"]
    assert J[(1, 1)]["fq"]["min_negative"] > 0
    for rho, r in J.items():
        assert r["hat_and_hatN_same_sign"]
