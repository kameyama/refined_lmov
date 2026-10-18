"""タスク C の高速版（src/rcsknots/daha_fast.py）：現行版（src/rcsknots/daha_dn.py）と全係数一致すること。"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def test_symmetrizer_gives_orbit_sum():
    """Π_+(X_{ε_1}) = Σ_i (X_{ε_i} + X_{−ε_i})（D_4）、Π_−(X_{−ε_1}) = Σ_i X_{−ε_i}（GL_4）。H_2 = Π_+(X_{2ε_1}) は対称。"""
    from rcsknots.daha_dn import DN, GLN
    from rcsknots.daha_fast import Fast, _with_rho2
    D = _with_rho2(DN(4))
    F = Fast(D, 1)
    e = lambda i, x: tuple(x if j == i else 0 for j in range(4))
    S1 = F.H(1)
    assert set(S1) == {e(i, x) for i in range(4) for x in (1, -1)} and all(v == 1 for v in S1.values())
    H2 = F.H(2)
    for j in range(4):
        assert all(H2.get(D.reflect(j, a)) == v for a, v in H2.items())
    G = Fast(_with_rho2(GLN(4)), 1, sign=-1)
    P = G.H(1)
    assert set(P) == {e(i, -1) for i in range(4)} and all(v == 1 for v in P.values())


def test_dn_matches_current():
    """D_N の Z, K：N = 2 の 2 箱（m = 1, 2）、N = 3 の 3 箱 (2,1)（m = 1）と 2 箱 (1,1)（m = 2）。λ の間で途中結果を共有。"""
    from rcsknots.daha_dn import DN, daha_jones
    from rcsknots.daha_fast import Fast, daha_jones_fast, _with_rho2
    for N, m, lams in [(2, 1, [(1,), (2,), (1, 1)]), (2, 2, [(2,), (1, 1)]), (3, 1, [(2, 1), (1,)]), (3, 2, [(1, 1)])]:
        F = Fast(_with_rho2(DN(N)), m)
        for lam in lams:
            Zo, Ko = daha_jones(N, lam, m)
            Zn, Kn = daha_jones_fast(N, lam, m, F=F)
            assert (Zo - Zn).is_zero() and (Ko - Kn).is_zero()


def test_gl_matches_current():
    """GL_N：γ̂(P_λ)·1 の評価（N = 3、λ = (2,1)、m = 1；N = 2、λ = (2)、m = 2）と随伴（N = 3、m = 1, 2）。"""
    from rcsknots.daha_dn import GLN, gl_gamma_P, gl_adj_gamma
    from rcsknots.daha_fast import gl_jones_fast, gl_adj_fast
    for N, lam, m in [(3, (2, 1), 1), (2, (2,), 2)]:
        D = GLN(N)
        assert (D.evaluate(gl_gamma_P(D, lam, m, 2)) - gl_jones_fast(N, lam, m)).is_zero()
    for m in (1, 2):
        D = GLN(3)
        assert (D.evaluate(gl_adj_gamma(D, m, 2)) - gl_adj_fast(3, m)).is_zero()
