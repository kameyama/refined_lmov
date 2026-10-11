"""[2,1] の構成の線形解析（scripts/lmov21_linear.py）の機能テスト。"""
import random
import sys
from pathlib import Path

import pytest
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from lmov21_linear import _images, rank_nullity, solve_unique  # noqa: E402


def test_laurent_images_match_direct_substitution():
    """因子ごとの像（_images）が、乱数の Laurent 多項式を各因子の零点に代入したものと一致する（分母の各既約因子が割り切ることと同値）。"""
    A, Q, T = sp.symbols("A Q T")
    rng = random.Random(3)
    P = {(rng.randint(-2, 2), rng.randint(-5, 5), rng.randint(-5, 5)): rng.randint(-4, 4) for _ in range(25)}
    P = {k: v for k, v in P.items() if v}
    f = sum(c * A ** a * Q ** q * T ** t for (a, q, t), c in P.items())
    img = {}
    for e, c in P.items():
        for key, v in _images(e).items():
            img[key] = img.get(key, 0) + c * v
    img = {k: v for k, v in img.items() if v}

    def collect(expr, gens):
        e = sp.expand(expr * A ** 20 * Q ** 40 * T ** 40)
        return {tuple(x - s for x, s in zip(m, (20, 40, 40))): int(c) for m, c in sp.Poly(e, A, Q, T).terms()}
    # Q = 1 の値と微分、QT = 1（T = 1/Q）、Q²T = 1（T = Q^{-2}）
    assert {k[1]: v for k, v in img.items() if k[0] == "Q=1"} == {(a, t): c for (a, q, t), c in collect(f.subs(Q, 1), None).items()}
    assert {k[1]: v for k, v in img.items() if k[0] == "dQ=1"} == {(a, t): c for (a, q, t), c in collect(sp.diff(f, Q).subs(Q, 1), None).items()}
    assert {k[1]: v for k, v in img.items() if k[0] == "QT-1"} == {(a, q): c for (a, q, t), c in collect(f.subs(T, 1 / Q), None).items() if t == 0}
    assert {k[1]: v for k, v in img.items() if k[0] == "Q2T-1"} == {(a, q): c for (a, q, t), c in collect(f.subs(T, Q ** -2), None).items() if t == 0}
    # Q² + Q + 1 の剰余：Q = ω で 0 ⇔ 成分 (1, Q) が 0
    w = sp.Rational(-1, 2) + sp.sqrt(3) * sp.I / 2                       # 1 の原始 3 乗根（厳密）
    val = sp.expand(f.subs({A: 2, T: 3}).subs(Q, w))
    comp = {k[1]: v for k, v in img.items() if k[0] == "Q^2+Q+1"}
    direct = sum(v * sp.Integer(2) ** a * sp.Integer(3) ** t * (1 if c == 0 else w) for (a, t, c), v in comp.items())
    assert abs(complex(val - direct)) < 1e-9


def test_rank_nullity_and_solve_on_small_system():
    """x1 + x2 = 3, x1 − x2 = 1 → 階数 2、零空間 0、解 (2, 1)。従属な式を足しても変わらない。"""
    rows = [({0: 1, 1: 1}, 3), ({0: 1, 1: -1}, 1), ({0: 2, 1: 0}, 4)]
    assert rank_nullity(rows, 2) == (2, True)
    assert solve_unique(rows, 2) == {0: 2, 1: 1}
    assert rank_nullity(rows + [({0: 1, 1: 1}, 4)], 2) == (2, False)
    assert rank_nullity([({0: 1, 1: 1}, 3)], 3) == (1, True)


@pytest.mark.results
def test_trefoil_21_is_determined_linearly():
    """三葉結び目：U+T+E+U2+L（正値性なし）で rCS_{[2,1]} の零空間が 0 で、唯一解が EHA の値に一致する（約 30 秒）。"""
    import lmov21_linear as LL
    from dgr_box3 import build21
    LL.MARGIN = 0
    small, H21u, window, true21 = LL.known_inputs("T2_3")
    lin, mons = build21(small, H21u, window, (1, 1, 1), conds=("U", "T", "E", "U2"))
    rows = [(dict(r), b) for r, b in zip(lin.rows, lin.rhs)] + LL.laurent_rows(small, mons, lin.col)
    rank, feas = rank_nullity(rows, len(mons))
    assert feas and rank == len(mons)
    x = solve_unique(rows, len(mons))
    inv = {v: k for k, v in lin.col.items()}
    assert {inv[c][1]: v for c, v in x.items() if v} == true21
