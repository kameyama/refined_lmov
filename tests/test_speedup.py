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
