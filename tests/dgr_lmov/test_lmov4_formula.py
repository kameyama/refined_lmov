"""N 箱の \\hat f の公式（scripts/dgr_lmov/lmov4_formula.py）と 4 箱の Laurent 性の検査（scripts/dgr_lmov/lmov4_joint.py）の機能テスト。"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/framed_knot"))
sys.path.insert(0, str(ROOT / "scripts/dgr_lmov"))

from framing_check import load_reduced  # noqa: E402
from lmov4_formula import formula_n, hat_n, parse_symbol  # noqa: E402


def test_parse_symbols():
    assert parse_symbol("C_21") == ((2, 1), 1)
    assert parse_symbol("C1_3") == ((1,), 3)
    assert parse_symbol("C_11_2") == ((1, 1), 2)


def test_formula_n_3_reproduces_hat3():
    """3 箱：一般化した公式の \\hat f_ρ が dgr_box3.hat3（既存の公式）と全係数で一致する（T(2,3)）。"""
    from dgr_box3 import hat3
    small = load_reduced("T2_3", 3)
    h = hat_n(small, 3)
    h3 = hat3(small)
    assert set(h) == set(h3) and all(h[r] == h3[r] for r in h)


def test_formula_n_4_symbols_and_denominators():
    """4 箱の記号は 14 個（Adams 像 C_2(x²), C_11(x²), C_1(x²) を含む。C_1(x⁴) は Möbius の μ(4) = 0 で現れず、C_1(x³) は 3 ∤ 4 で現れない）。"""
    fo, keys = formula_n(4)
    assert len(keys) == 14 and ((2,), 2) in keys and ((1, 1), 2) in keys and ((1,), 2) in keys
    assert ((1,), 4) not in keys and ((1,), 3) not in keys
    assert set(fo) == {(4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1)}


@pytest.mark.results
def test_formula_n_4_matches_stored_hat_f():
    """4 箱：T(2,3) の \\hat f（data/framed_knot/highbox、本体のチェックアウト）と全係数で一致（GN からの補間を含めて約 30 秒）。"""
    import json
    import lzma
    import numpy as np
    from rcsknots.exact import reconstruct_points
    from rcsknots.extract import A_OFFSET, _arr
    from rcsknots.field import V
    from rcsknots.bpsgen import torus_rcs_factory
    from rcsknots.partitions import partitions
    small = load_reduced("T2_3", 3)
    fac = torus_rcs_factory(2, 3)
    for lam in partitions(4):
        def ev(xs, ys, p, lam=lam):
            k = len(xs)
            Qv, Tv = V(xs, p), V(ys, p)
            lo, cs = fac(lam, Qv, Tv)
            arr = np.zeros((k, 2 * A_OFFSET + 1), dtype=np.int64)
            for i, c in enumerate(cs):
                arr[:, lo + i + A_OFFSET] = _arr(c, k, p) % p
            return arr
        raw = reconstruct_points(ev, seed=7, log=lambda s: None)
        small[lam] = {(d - A_OFFSET, eq, et): c for (d, eq, et), c in raw.items()}
    ref = json.loads(lzma.open("/Users/masaya/projects/ai_research/refined_cs/data/framed_knot/highbox/T2_3_N4.json.xz").read())
    for rho, h in hat_n(small, 4).items():
        stored = {tuple(map(int, k.split(","))): v for k, v in ref["rho"]["".join(map(str, rho))]["hat_f"].items()}
        assert h == stored
