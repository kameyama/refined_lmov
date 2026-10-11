"""refined LMOV の既約表示・種数表示の道具（src/rcsknots/lmovstruct.py）の機能テスト。小さな入力（数秒）。"""
import random
from pathlib import Path

import pytest

from rcsknots.bps import g_expand
from rcsknots.knot_table import hat_f1
from rcsknots.lmovstruct import (W, W_lower, delta_values, from_rows, from_spin, is_W_invariant, load_bps, rows,
                                 specialize_bold, specialize_lower, spin, spin_count, unimodal_rows)

ROOT = Path(__file__).resolve().parent.parent
TREFOIL = {(2, -2, 0): 1, (2, 2, 2): 1, (4, 0, 3): 1}        # DGR の T(2,3)（負の交差）


def test_rows_and_spin_of_trefoil():
    """手計算：P の行は (β, f) = (1, 1): {−1: 1, 1: 1}、(2, 3): {0: 1}。既約表示は (1, 1, 1), (2, 0, 3)。"""
    assert rows(TREFOIL) == {(1, 1): {-1: 1, 1: 1}, (2, 3): {0: 1}}
    assert spin(TREFOIL) == {(1, 1, 1): 1, (2, 0, 3): 1}
    f = hat_f1(TREFOIL)                                          # (1 + 𝐚^{−2}𝐭^{−1}) P
    assert spin(f) == {(0, 1, 0): 1, (1, 0, 2): 1, (1, 1, 1): 1, (2, 0, 3): 1}
    assert from_spin(spin(f)) == f and from_rows(rows(f)) == f


def test_W_is_dgr_symmetry_and_equals_g_closing():
    """W 不変 ⇔ G 展開が閉じる（乱数の多項式 200 例）。"""
    rng = random.Random(1)
    for _ in range(200):
        b = {}
        for _ in range(rng.randint(1, 6)):
            beta, k, f, c = rng.randint(-2, 2), rng.randint(-3, 3), rng.randint(-3, 3), rng.randint(-3, 3)
            if c:
                b[(2 * beta, 2 * k, f + k)] = c
        if rng.random() < 0.5:                                   # 半分は対称化する
            b = {k: v + W(b).get(k, 0) for k, v in {**b, **W(b)}.items()}
            b = {k: v for k, v in b.items() if v}
        closes = True
        try:
            g_expand(b)
        except ValueError:
            closes = False
        assert closes == is_W_invariant(b)
        if closes:
            assert from_spin(spin(b)) == b


def test_W_lower_matches_bold():
    """小文字の (A,Q,T) ↦ (A, 1/T, 1/Q) は太字の W と同じ（4_1 の [2]）。"""
    meta, lower, bold = load_bps(ROOT / "data/bps/4_1_2.json")
    assert W_lower(lower) == lower and is_W_invariant(bold)
    assert spin_count(bold)[1] == 0 and unimodal_rows(bold)


def test_spin_nonneg_is_not_implied_by_nonneg():
    """非負だが単峰でない行：𝐚⁰(w² + w^{−2}) 𝐭⁰（k = 0 が欠ける）。既約表示は χ_1 − χ_0。"""
    b = {(0, 4, 2): 1, (0, -4, -2): 1}
    assert is_W_invariant(b) and all(v > 0 for v in b.values())
    assert spin(b) == {(0, 2, 0): 1, (0, 0, 0): -1} and not unimodal_rows(b)


def test_specializations_of_trefoil():
    """rCS_[1](T(2,3)) = a(1 + qt − aq)/√(qt) は a = t で (t/q)^{1/2}（1 項）、a = 1 で \\hat f_[1] = 0。"""
    r = {(2, -1, -1): 1, (2, 1, 1): 1, (4, 1, -1): -1}
    assert specialize_lower(r, 1) == {(-1, 1): 1}
    meta, lower, bold = load_bps(ROOT / "data/bps/T2_3_1.json")
    assert specialize_lower(lower, 0) == {}
    # 太字：a = t^N ⇔ 𝐚² = −𝐪^{2N}𝐭^{−1}（d_N の対が打ち消す点）。[1] の \\hat f は a = t で 2 項
    assert len(specialize_bold(bold, 2, -1, -1)) == 2 == len(specialize_lower(lower, 1))
    assert specialize_bold(bold, 0, -1, -1) == {}
    # d_1 の対：三葉結び目の (4,0,3) と (2,2,2) は 𝐚² = −𝐪²/𝐭 で打ち消し、(2,−2,0) が残る
    assert specialize_bold(TREFOIL, 2, -1, -1) == {(0, -1): -1}


def test_delta_values():
    assert delta_values(TREFOIL) == [-1] and delta_values(hat_f1(TREFOIL)) == [-1, 0]


def test_colored_d1_vanishing_and_symmetric_monomials():
    """彩色 d_1：a = t で |ρ| ≥ 2 の \\hat f_ρ は 0。根拠は rCS_[r](a=t) = c^r（c = rCS_[1](a=t) = (t/q)^{δ/2}）。T(2,3)（δ = 1）で確認。"""
    import json
    for rho in ("2", "11", "3", "21", "111"):
        meta, lower, bold = load_bps(ROOT / f"data/bps/T2_3_{rho}.json")
        assert specialize_lower(lower, 1) == {}
    for color, r in (("1", 1), ("2", 2), ("3", 3)):
        poly = {tuple(map(int, k.split(","))): v for k, v in json.loads((ROOT / f"data/torus/T2_3_{color}.json").read_text())["poly"].items()}
        assert specialize_lower(poly, 1) == {(-r, r): 1}            # (T/Q)^r = (t/q)^{r/2}
    poly = {tuple(map(int, k.split(","))): v for k, v in json.loads((ROOT / "data/torus/T2_3_11.json").read_text())["poly"].items()}
    assert len(specialize_lower(poly, 1)) > 1                        # 列 [1,1] は a = t で単項式にならない（unknot_[1,1](a=t) = 0 が効く）
    assert len(specialize_lower(poly, 2)) == 1                       # [1,1] は a = t²（SU(2) の det）で単項式
