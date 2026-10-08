"""第 4 期 V0：高速化と並列・再開の版が、従来の計算と全係数で一致すること。"""
import shutil
from pathlib import Path

import numpy as np

from rcsknots.extract import reformulated_exact
from rcsknots.interp import newton_to_monomial, newton_to_monomial_direct
from rcsknots.laurent import load_table

ROOT = Path(__file__).resolve().parent.parent
C = {"1": (1,), "2": (2,), "11": (1, 1), "3": (3,), "111": (1, 1, 1), "21": (2, 1)}


def test_interpolation_matmul_equals_direct():
    p = 2**31 - 1
    rng = np.random.default_rng(0)
    for N, M in [(1, 1), (7, 5), (60, 300)]:
        xs = rng.choice(np.arange(2, 10**6), N, replace=False)
        vals = rng.integers(0, p, (N, M))
        assert np.array_equal(newton_to_monomial(xs, vals, p), newton_to_monomial_direct(xs, vals, p))


def test_parallel_and_resume_equal_serial(tmp_path):
    R = load_table(ROOT / "data/reference/notebook_torus.txt")
    red = {lam: R[("T2_3", c)] for c, lam in C.items()}
    serial = reformulated_exact(red, 3, log=lambda s: None)
    cache = tmp_path / "cache"
    par = reformulated_exact(red, 3, log=lambda s: None, workers=3, cache=cache)
    resumed = reformulated_exact(red, 3, log=lambda s: None, workers=3, cache=cache)
    assert serial == par == resumed
    assert any(cache.glob("*.npy"))
    shutil.rmtree(cache)


def test_single_prime_then_crt_fallback():
    """係数が小さければ素数 1 つで復元し、P31/2 を超える係数があれば 2 つ目の素数（CRT）に切り替える。"""
    from rcsknots.exact import reconstruct_points
    for poly in ({(0, -2, 1): 7, (1, 3, -1): -5, (0, 0, 0): 1},
                 {(0, -2, 1): 3_000_000_000, (1, 3, -1): -5_000_000_000, (0, 0, 0): 1}):
        def ev(xs, ys, p, poly=poly):
            out = np.zeros((len(xs), 2), dtype=np.int64)
            for (d, i, j), c in poly.items():
                out[:, d] = (out[:, d] + c % p * np.array([pow(int(x), i % (p - 1), p) * pow(int(y), j % (p - 1), p) % p
                                                             for x, y in zip(xs, ys)], dtype=np.int64)) % p
            return out
        assert reconstruct_points(ev, log=lambda s: None) == poly


def test_blocked_parallel_interpolation(monkeypatch, tmp_path):
    """補間を a の係数の列のまとまりに分けて並列にしても、素数 1 つ・2 つ（CRT）のどちらでも同じ多項式を復元する。"""
    from rcsknots import exact
    monkeypatch.setattr(exact, "INTERP_BLOCK", 2)
    for poly in ({(0, -2, 1): 7, (3, 3, -1): -5, (4, 0, 0): 1, (2, 1, 2): 11},
                 {(0, -2, 1): 3_000_000_000, (3, 3, -1): -5_000_000_000, (4, 0, 0): 1, (2, 1, 2): 11}):
        def ev(xs, ys, p, poly=poly):
            out = np.zeros((len(xs), 5), dtype=np.int64)
            for (d, i, j), c in poly.items():
                out[:, d] = (out[:, d] + c % p * np.array([pow(int(x), i % (p - 1), p) * pow(int(y), j % (p - 1), p) % p
                                                             for x, y in zip(xs, ys)], dtype=np.int64)) % p
            return out
        cache = tmp_path / str(len(str(max(poly.values()))))
        assert exact.reconstruct_points(ev, log=lambda s: None, workers=3, cache=cache) == poly
        assert exact.reconstruct_points(ev, log=lambda s: None, workers=3, cache=cache) == poly    # 保存した評価から再開


def test_g_expand_fast_equals_direct():
    """G 展開（python-flint の多項式で計算）が、旧来の実装 g_expand_direct と全係数で一致し、ランダムな Σ c_g(𝐭) G^g を復元する。"""
    import random
    import pytest
    from rcsknots.bps import g_expand, g_expand_direct
    G = {(2, 1): 1, (0, 1): -1, (0, -1): -1, (-2, -1): 1}
    rng = random.Random(1)
    for _ in range(20):
        coeffs = {(rng.randint(0, 6), rng.randint(0, 2), rng.randint(-4, 4)): rng.randint(-5, 5) for _ in range(12)}
        coeffs = {k: v for k, v in coeffs.items() if v}
        bold = {}
        for (g, beta, F), c in coeffs.items():
            P = {(0, 0): 1}
            for _ in range(g):
                Q = {}
                for (a1, b1), x in P.items():
                    for (a2, b2), y in G.items():
                        Q[(a1 + a2, b1 + b2)] = Q.get((a1 + a2, b1 + b2), 0) + x * y
                P = Q
            for (eq, et), x in P.items():
                k = (2 * beta, eq, et + F)
                bold[k] = bold.get(k, 0) + c * x
        bold = {k: v for k, v in bold.items() if v}
        assert g_expand(bold) == coeffs == g_expand_direct(bold)
    for bad in ({(0, 1, 0): 1}, {(0, 2, 1): 1}):        # 𝐪 の奇数冪、w ↔ 1/w で非対称
        with pytest.raises(ValueError):
            g_expand(bad)
        with pytest.raises(ValueError):
            g_expand_direct(bad)


def test_evaluate_closed_form_equals_via_p():
    """ε の評価：H̃_μ[1 − a] = Π(1 − a q^c t_w^r) による直接の評価が、p 基底を経由する旧来の評価と一致する（mod p、乱数の点）。"""
    from rcsknots.field import V, P31, make_conv
    from rcsknots.gn import GN
    n = 4
    rng = np.random.default_rng(3)
    G = GN(V(rng.integers(2, P31, n), P31), V(rng.integers(2, P31, n), P31), V(np.full(n, 3), P31), make_conv(P31, n))
    for lam, slope in [((1,), (2, -3)), ((2,), (2, -1)), ((1, 1), (3, -2)), ((2, 1), (2, -1)), ((3,), (2, -3))]:
        vec = G.apply_symmetric(G.P(sum(lam))[lam], *slope, {(): G.one})
        A, B = G.evaluate_apoly(vec), G.evaluate_apoly_via_p(vec)
        m = max(len(A), len(B))
        A, B = A + [G.zero] * (m - len(A)), B + [G.zero] * (m - len(B))
        assert all(np.array_equal(x.a % P31, y.a % P31) for x, y in zip(A, B))


def test_eha_apply_merged_equals_direct():
    """EHA の演算子：全ての λ をまとめた動的計画法（apply）が、λ ごとに回す旧来の実装（apply_direct）と一致する（mod p、乱数のベクトル）。"""
    from rcsknots.field import V, P31, make_conv
    from rcsknots.gn import GN
    from rcsknots.partitions import partitions
    n = 4
    rng = np.random.default_rng(5)
    E = GN(V(rng.integers(2, P31, n), P31), V(rng.integers(2, P31, n), P31), V(np.full(n, 3), P31), make_conv(P31, n)).eha
    for d in (3, 5):
        vec = {mu: V(rng.integers(1, P31, n), P31) for mu in partitions(d)}
        for k, nn, m in [(1, 2, -3), (2, 2, -3), (1, 3, 2), (3, 1, 0), (2, 1, -1)]:
            A, B = E.apply(k, nn, m, vec), E.apply_direct(k, nn, m, vec)
            assert set(A) == set(B) and all(np.array_equal(A[x].a % P31, B[x].a % P31) for x in A)
