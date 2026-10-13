"""別研究 [2,1] の X0：ダブルツイスト結び目 K_{m,n} の t = q の彩色 HOMFLY-PT（文献の式とケーブルの逆解き）。

機能テストは小さな結び目（3_1, 4_1）で数秒〜十数秒。全ての結び目の再照合は results の印（scripts/framed_knot/x0_check.py の保存値）。
"""
import sys
from pathlib import Path

import pytest
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/framed_knot"))
from colored_from_cables import L, M, colored_fast, colored_p, load, sage_to_kn  # noqa: E402
from rcsknots.double_twist import (A, q, B21, Bkntz, Dprinc, F21_1612, F21_1903, F_kntz, H21_1612, H21_1903,  # noqa: E402
                                   H_rect, Z21_1612, Z21_1903, to_kn)
from rcsknots.homfly import Q, a  # noqa: E402

REF = ROOT / "data/framed_knot/reference"


def zero(x):
    return sp.expand(sp.cancel(sp.together(x))) == 0


# ---------------------------------------------------------------- 文献の式の内部の性質
@pytest.mark.parametrize("lam", [(1,), (2,), (1, 1), (2, 1), (2, 2), (3,)])
def test_kntz_B_empty_is_D(lam):
    """KNTZ :627「μ = ∅ では B_{λ,∅} = D_λ」（ガウスの二項係数の規約の確認）。"""
    for s in (len(lam), len(lam) + 1):
        assert zero(Bkntz(lam, (), s) - Dprinc(lam, s))


@pytest.mark.parametrize("lam", [(1,), (2,), (1, 1), (2, 1), (2, 2)])
def test_kntz_F_independent_of_s(lam):
    """KNTZ :641「F^{(m)}_λ は s（≥ length(λ)）によらない Laurent 多項式」。"""
    for m in (2, 3):
        f = F_kntz(lam, m, len(lam))
        assert zero(f - F_kntz(lam, m, len(lam) + 1))
        assert len(sp.Poly(sp.fraction(sp.cancel(f))[1], A, q).terms()) == 1   # 分母は単項式（Laurent 多項式）


def test_morozov_B21_properties():
    """1903.00259 :250-254：F^{(−1)} = 1（4_1）、F^{(0)} = 0（自明な結び目、∅ 以外）、F^{(1)} は単項式（3_1）。"""
    assert all(zero(f - 1) for f in F21_1903(-1))
    assert [sp.simplify(f) for f in F21_1903(0)] == [1, 0, 0, 0, 0, 0]
    assert [sp.simplify(f) for f in F21_1903(1)] == [1, -A**2, A**4 / q**2, q**2 * A**4, -A**6, -A**6]
    assert B21().is_lower


def test_morozov_two_papers_share_F_and_Z1():
    """1612.00422 と 1903.00259 の F^{(±1)} と Z^{[1]} の 2 つの書き方は一致する。"""
    for m in (-1, 1):
        assert all(zero(x - y) for x, y in zip(F21_1903(m), F21_1612(m)))
    assert zero(Z21_1903()[1] - Z21_1612()[1])


def test_morozov_F_from_B_equals_1612_F():
    """1903 の F^{(m)} = Σ_Y (𝓑^{m+1})_{XY} は 1612 の (Ftwist21) と一致（m = −3..3。固有値は 6 個なので 6 点で全ての m に従う）。"""
    for m in range(-3, 4):
        assert all(zero(x - y) for x, y in zip(F21_1903(m), F21_1612(m))), m


@pytest.mark.parametrize("mn", [(2, 2), (-2, 2), (1, 1)])
def test_morozov_1903_normalized_equals_1612(mn):
    """1903 の (difexpan) を Z/(F^{(1)}F^{(−1)}) で読むと 1612 の (H21de) に等しい。書かれたとおりの Σ Z F F は K_{m,2} で異なる。"""
    from rcsknots.double_twist import Z21_1903 as Zf
    assert zero(H21_1903(*mn) - H21_1612(*mn))
    Fm, Fn = F21_1903(mn[0]), F21_1903(mn[1])
    literal = sum(z * x * y for z, x, y in zip(Zf(), Fm, Fn))
    assert not zero(literal - H21_1612(*mn)) or mn[1] == -1


# ---------------------------------------------------------------- 結び目の同定（[1]）
def test_braids_are_kntz_knots():
    """data/framed_knot/reference/x0_double_twist.txt の各ブレイドの [1] が KNTZ の K_{m,n} の式と（鏡像でなく）一致。"""
    n_ok = 0
    for line in open(REF / "x0_double_twist.txt"):
        name, _, _, h = line.rstrip("\n").split("\t")
        m, n = (int(x.replace("m", "-")) for x in name[1:].split("_"))
        assert sp.expand(sage_to_kn(sp.sympify(h, locals={"L": L, "M": M})) - to_kn(H_rect(1, 1, m, n))) == 0, name
        n_ok += 1
    assert n_ok >= 24


# ---------------------------------------------------------------- 長方形での較正（KNTZ = ケーブルの逆解き）
@pytest.mark.parametrize("name,mn", [("D1_1", (1, 1)), ("Dm1_1", (-1, 1)), ("D2_1", (2, 1))])
def test_rect_2_and_11(name, mn):
    data = load(str(REF / "x0_cable_p2.txt"), [name])
    col, checked = colored_fast(name, data, 2)
    assert checked == 1
    assert sp.expand(col[(2,)] - to_kn(H_rect(2, 1, *mn))) == 0
    assert sp.expand(col[(1, 1)] - to_kn(H_rect(1, 2, *mn))) == 0


def test_rect_22_trefoil():
    """[2,2]：2 段の 2 本ケーブル（6 つの傾きの 5 つで解き 1 つで検証）。[4], [1⁴] も KNTZ と一致。"""
    data = load(str(REF / "x0_cable_p22.txt"), ["D1_1"])
    col, checked = colored_fast("D1_1", data, 4)
    assert checked == 1
    for lam, (r, s) in {(2, 2): (2, 2), (4,): (4, 1), (1, 1, 1, 1): (1, 4)}.items():
        assert sp.expand(col[lam] - to_kn(H_rect(r, s, 1, 1))) == 0


# ---------------------------------------------------------------- [2,1]
@pytest.mark.parametrize("name,mn", [("D1_1", (1, 1)), ("Dm1_1", (-1, 1))])
def test_21_trefoil_figure_eight(name, mn):
    data = load(str(REF / "x0_cable_p3.txt"), [name])
    col, checked = colored_fast(name, data, 3)
    assert checked == 1
    assert sp.expand(col[(2, 1)] - to_kn(H21_1612(*mn))) == 0
    assert sp.expand(col[(2, 1)] - to_kn(H21_1903(*mn))) == 0
    assert sp.expand(col[(3,)] - to_kn(H_rect(3, 1, *mn))) == 0


def test_fast_equals_old_trefoil():
    """高速版（逆行列を先に作る）と旧来の colored_p（sympy の solve）が 3_1 の 3 本ケーブルで全係数一致。"""
    old = colored_p("3_1", 3, str(REF / "cable_homfly_p3.txt"))
    new, _ = colored_fast("3_1", load(str(REF / "cable_homfly_p3.txt")), 3)
    assert all(sp.expand(sp.expand(old[mu]) - new[mu]) == 0 for mu in new)


def test_literature_gj_and_mmm_figure_eight():
    """Gu–Jockers 1407.5643 の 4_1, 5_2 と MMM 1408.3076 の k = −1, 2 が 1612.00422 の式に一致（conventions.md §14）。

    GJ は λ = a, q = Q² で、5_2 は KNTZ の K_{2,1} の鏡像（λ → 1/λ）。[2,1] は q ↔ q^{-1} で不変なので、変数の対応と鏡像は区別できない。"""
    from x0_literature import gj_21, mmm_21, GJ_MAPS
    gj = gj_21()
    assert sp.expand(gj["4_1"].subs(GJ_MAPS["λ=a, q=Q²"], simultaneous=True) - to_kn(H21_1612(-1, 1))) == 0
    assert sp.expand(gj["5_2"].subs(GJ_MAPS["λ=1/a, q=Q²"], simultaneous=True) - to_kn(H21_1612(2, 1))) == 0
    for k in (-1, 2):
        assert zero(mmm_21(k) - H21_1612(k, 1))


# ---------------------------------------------------------------- 保存値の再照合（重い）
@pytest.mark.results
def test_fast_equals_old_K3(tmp_path):
    """E の調査の確認（K_3 = 7_2、旧来の方法で約 2 分）。"""
    lines = [ln for ln in open(REF / "x0_cable_p3.txt") if ln.startswith("D3_1\t") and ln.split("\t")[1] in ("1", "2", "4")]
    p = tmp_path / "p3.txt"
    p.write_text("".join(lines))
    old = colored_p("D3_1", 3, str(p))
    new, _ = colored_fast("D3_1", load(str(p)), 3)
    assert all(sp.expand(sp.expand(old[mu]) - new[mu]) == 0 for mu in new)


@pytest.mark.results
def test_saved_values_match_formulas():
    """data/framed_knot/x0/colored_t_eq_q.txt の全ての値を文献の式と照合（scripts/framed_knot/x0_check.py の結果の再確認）。"""
    from x0_check import reference
    n = 0
    for line in open(ROOT / "data/framed_knot/x0/colored_t_eq_q.txt"):
        name, lam, val = line.rstrip("\n").split("\t")
        lam = tuple(int(x) for x in lam.strip("()").split(",") if x.strip())
        m, nn = (int(x.replace("m", "-")) for x in name[1:].split("_"))
        ref = reference(lam, m, nn)
        if ref is None:
            continue
        v = sp.sympify(val, locals={"a": a, "Q": Q})
        assert all(sp.expand(r - v) == 0 for r in ref), (name, lam)
        n += 1
    assert n > 0
