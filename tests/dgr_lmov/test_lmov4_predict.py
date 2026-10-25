"""lmov4_predict の法 p の消去（ρ ごとに流し込む RREF の状態）の機能テスト。"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts/framed_knot"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts/dgr_lmov"))
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


def test_pickle_rows_to_csr_and_eliminate_csr(tmp_path):
    import pickle
    from lmov4_predict import eliminate_csr, load_csr, pickle_rows_to_csr
    rng = np.random.default_rng(2)
    ncols = 1300
    rows = [{int(j): int(rng.integers(-2**29, 2**29)) for j in rng.choice(ncols, size=1105, replace=False)},   # 1000 項を超える辞書
            {7: 3}, {}, {int(j): int(v) for j, v in zip(range(0, 1300, 13), rng.integers(1, 300, size=100))}]
    rows += [{int(j): int(rng.integers(-9, 10)) for j in rng.choice(ncols, size=50, replace=False)} for _ in range(1200)]  # APPENDS が分かれる
    rhs = [int(v) for v in rng.integers(-2**30, 2**30, size=len(rows))]
    src = tmp_path / "rows.pkl"
    src.write_bytes(pickle.dumps((rows, rhs)))
    n, nnz = pickle_rows_to_csr(src, tmp_path / "csr")
    assert (n, nnz) == (len(rows), sum(len(r) for r in rows))
    indptr, indices, vals, rhs2 = load_csr(tmp_path / "csr")
    back = [{int(c): int(v) for c, v in zip(indices[indptr[i]:indptr[i + 1]], vals[indptr[i]:indptr[i + 1]])} for i in range(n)]
    assert back == rows and list(rhs2) == rhs
    p = PRIMES[1]
    k = 40                                                     # 消去の比較は先頭の 40 行だけ（速さのため）
    st_a = eliminate(new_state(ncols), rows[:k], rhs[:k], p, batch=16)
    st_b = eliminate_csr(new_state(ncols), indptr[:k + 1], indices, vals, rhs2, p, batch=16)
    assert st_a["piv"] == st_b["piv"] and np.array_equal(st_a["B"], st_b["B"]) and st_a["bad"] == st_b["bad"]
