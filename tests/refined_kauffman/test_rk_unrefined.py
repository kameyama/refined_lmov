"""refined KauffmanK0：unrefined の Kauffman large N duality の再現と、refined unknot の検算。

- trefoil の \\hat g_[1,1] が Mariño 0904.1088 (4.8) と一致（\\hat g_[2] は z^9 の項が原論文の誤植で、他の項は一致）。
- 4_1 の \\hat g_[2], \\hat g_[1,1] が NRZ 1310.2240 (C.10) と一致し、N^{c=0} の表が NRZ の表と一致。
- AS 1202.4456 の D_N refined unknot は [1] で 2018 年の式と一致し、[2] では (q−t) に比例する差がある。
"""
import importlib.util
import pathlib

import pytest
import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parents[2]


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / "refined_kauffman" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def rk():
    return _load("rk_unrefined_check")


def test_trefoil_marino_48(rk):
    nu, t, z = rk.nu, rk.t, rk.z
    from rcsknots import homfly
    sub = {homfly.a: nu**2, homfly.Q: t}
    H = {lam: rk.simp(homfly.torus_reduced(lam, 2, 3, -1).subs(sub) * rk.H_unk[lam]) for lam in ((1,), (2,), (1, 1))}
    Hadj = rk.simp(rk.H_adj_unk * (4*nu**4 - 4*nu**6 + nu**8 + z**2*(4*nu**4 - 7*nu**6 + 2*nu**8 + nu**10) + z**4*(nu**4 - 2*nu**6 + nu**8)))
    F31 = 2*nu**2 - nu**4 + z*(-nu**3 + nu**5) + z**2*(nu**2 - nu**4)
    assert sp.simplify(rk.nrz_trefoil(1, nu, t, sp.Integer(-1)) - F31) == 0
    G = {(1,): rk.simp(F31 * rk.G_unk[(1,)]),
         (2,): rk.simp(rk.nrz_trefoil(2, nu, t, sp.Integer(-1)) * rk.G_unk[(2,)]),
         (1, 1): rk.simp(rk.nrz_trefoil(2, -nu, 1/t, sp.Integer(-1)) * rk.G_unk[(1, 1)])}
    g1, g2, g11, h2, h11 = rk.reformulated(H, Hadj, G)
    gh2, gh11 = rk.hat2(g2, g11)
    M11 = (-15*nu**11 + 53*nu**9 - 69*nu**7 + 39*nu**5 - 8*nu**3 + 15*z*(nu**12 - 3*nu**10 + 3*nu**8 - nu**6)
           + z**2*(-35*nu**11 + 126*nu**9 - 146*nu**7 + 61*nu**5 - 6*nu**3) + 5*z**3*(7*nu**12 - 24*nu**10 + 24*nu**8 - 7*nu**6)
           + z**4*(-28*nu**11 + 120*nu**9 - 128*nu**7 + 37*nu**5 - nu**3) + 7*z**5*(4*nu**12 - 17*nu**10 + 17*nu**8 - 4*nu**6)
           + z**6*(-9*nu**11 + 55*nu**9 - 56*nu**7 + 10*nu**5) + z**7*(9*nu**12 - 55*nu**10 + 55*nu**8 - 9*nu**6)
           + z**8*(-nu**11 + 12*nu**9 - 12*nu**7 + nu**5) + z**9*(nu**12 - 12*nu**10 + 12*nu**8 - nu**6)
           + z**10*(nu**9 - nu**7) + z**11*(nu**8 - nu**10))
    assert sp.simplify(gh11 - M11) == 0
    # \hat g_[2] の z^9 の項（Mariño の誤植の訂正値）と z^8 の項
    e = rk.z_expansion(gh2)
    assert {k: v for k, v in e.items() if k[0] == 9} == {(9, 6): -11, (9, 8): 78, (9, 10): -78, (9, 12): 11}
    assert {k: v for k, v in e.items() if k[0] == 8} == {(8, 5): 12, (8, 7): -79, (8, 9): 78, (8, 11): -11}
    hh2, hh11 = rk.hat2(h2, h11)
    assert rk.in_ring(gh2, "g") and rk.in_ring(hh2, "h")


def test_fig8_nrz_tables(rk):
    nu, t, z = rk.nu, rk.t, rk.z
    H41 = {}
    for lam, r in (((1,), 1), ((2,), 2)):
        H41[lam] = rk.simp(rk.supertwist(r+1, -1, nu**2, t**2, sp.Integer(-1), -1) * rk.H_unk[lam])
    H41[(1, 1)] = rk.simp(rk.supertwist(3, -1, nu**2, t**-2, sp.Integer(-1), -1) * rk.H_unk[(1, 1)])
    adj41 = (3 - 2*z**2 - 6*z**4 - 2*z**6 + (-2 - z**2 + 2*z**4 + z**6)*(nu**-2 + nu**2) + (1 + 2*z**2 + z**4)*(nu**-4 + nu**4))
    Hadj41 = rk.simp(adj41 * rk.H_adj_unk)
    G41 = {(1,): rk.simp(rk.nrz_fig8(1, nu, t, sp.Integer(-1)) * rk.G_unk[(1,)]),
           (2,): rk.simp(rk.nrz_fig8(2, nu, t, sp.Integer(-1)) * rk.G_unk[(2,)]),
           (1, 1): rk.simp(rk.nrz_fig8(2, -nu, 1/t, sp.Integer(-1)) * rk.G_unk[(1, 1)])}
    g1, g2, g11, h2, h11 = rk.reformulated(H41, Hadj41, G41)
    gh2, gh11 = rk.hat2(g2, g11)
    hh2, hh11 = rk.hat2(h2, h11)
    # NRZ の N^{c=0}_{[2]} の表（g = 0..3、β = −6..6 偶数）
    tab = {}
    for (k, Q), c in rk.z_expansion(sp.expand(sp.cancel(hh2*z))).items():
        tab[(k//2, Q)] = c
    expected = {0: [4, -10, 6, -5, 16, -18, 6], 1: [1, -4, -5, 10, 9, -14, 3], 2: [0, 0, -6, 10, -2, -2, 0], 3: [0, 0, -1, 2, -1, 0, 0]}
    for g, row in expected.items():
        for Q, c in zip(range(-6, 7, 2), row):
            assert tab.get((g, Q), 0) == c, (g, Q)
    # \hat g の c=1（偶数 z）と c=2（奇数 z）の先頭の表（NRZ の表の g=0 の行）
    e2 = rk.z_expansion(gh2)
    assert [e2.get((0, Q), 0) for Q in range(-5, 6, 2)] == [5, -19, 34, -38, 25, -7]
    assert [e2.get((1, Q), 0) for Q in range(-6, 7, 2)] == [10, -32, 40, -45, 65, -55, 17]
    e11 = rk.z_expansion(gh11)
    assert [e11.get((0, Q), 0) for Q in range(-5, 6, 2)] == [7, -25, 38, -34, 19, -5]


def test_refined_so_unknot_vs_2018():
    m = _load("rk_so_unknot_check")
    q, t, u = m.q, m.t, m.u
    N = 5
    for lam, nb, expect in (((1,), m.rkauf1, 0), ((2,), m.rkauf2, (q - t)*(q - t**4)/((q*t - 1)*(q*t**4 - 1)))):
        d = sp.simplify(m.eval_generic(lam, N) - nb.subs(u, t**(N-1)))
        assert sp.simplify(d - expect) == 0
    # AS (5.5) の t = q^β の積と一般 t の評価式の一致
    qv = sp.Rational(3, 7)
    for lam in ((2,), (1, 1)):
        assert sp.simplify(m.dim_AS(lam, 3, 2).subs({q: qv, t: qv**2}) - m.eval_generic(lam, 3).subs({q: qv, t: qv**2})) == 0
