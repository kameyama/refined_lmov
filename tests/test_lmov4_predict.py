"""lmov4_predict の法 p の消去（ρ ごとに流し込む RREF の状態）の機能テスト。"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from lmov4_predict import PRIMES, eliminate, new_state, nullspace_of, solution_of  # noqa: E402


def test_eliminate_streaming_matches_rank_and_solves():
    rng = np.random.default_rng(1)
    ncols, rank = 40, 25
    G = rng.integers(-5, 6, size=(rank, ncols))
    C = rng.integers(-3, 4, size=(90, rank))
    A = C @ G                                                  # 階数 25 の 90 × 40
    ytrue = rng.integers(-9, 10, size=ncols)
    rhs = A @ ytrue
    rows = [{j: int(v) for j, v in enumerate(r) if v} for r in A]
    p = PRIMES[0]
    st = new_state(ncols)
    for start in (0, 30, 60):                                  # 3 つのまとまりに分けて流し込む
        eliminate(st, rows[start:start + 30], [int(b) for b in rhs[start:start + 30]], p, batch=7)
    y, nullity, bad = solution_of(st)
    assert not bad and nullity == ncols - rank
    yv = np.array([int(v) for v in y], dtype=object)
    assert all((A @ yv - rhs) % p == 0)
    basis = nullspace_of(st, p)
    assert len(basis) == nullity
    for vec in basis:
        v = np.array([vec.get(j, 0) for j in range(ncols)], dtype=object)
        assert all((A @ v) % p == 0)
    # 矛盾する右辺
    st2 = new_state(ncols)
    bad_rhs = [int(b) for b in rhs]; bad_rhs[5] += 1
    eliminate(st2, rows, bad_rhs, p)
    assert solution_of(st2)[2]
