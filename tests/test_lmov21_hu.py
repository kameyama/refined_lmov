"""D12 段階 2（scripts/lmov21_hu.py）：仮説 H_u″（P 正規化の hook の Schur 係数の冪級数の非負性）を課した構成の機能テスト。"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import lmov21_hu as HU  # noqa: E402


def quiet(_):
    pass


def test_series_rows_divide_by_macdonald_denominators():
    """冪級数の行：N = 1 − 𝐪² を (1 − 𝐪²)²(1 − 𝐪⁶𝐭²) で割ると 1/((1 − 𝐪²)(1 − 𝐪⁶𝐭²)) の展開になる（最初の係数を手計算と比較）。"""
    rows = HU.series_rows({None: {(0, 0, 0): 1, (0, 2, 0): -1}}, (2, 1), order=12)
    s = {k: v[None] for k, v in rows.items()}
    assert s[(0, 0, 0)] == 1 and s[(0, 2, 0)] == 1 and s[(0, 4, 0)] == 1       # 1/(1 − 𝐪²)
    assert s[(0, 6, 2)] == 1 and s[(0, 8, 2)] == 1                             # × 𝐪⁶𝐭²
    assert (0, 2, 2) not in s


def test_t25_series_positivity_bounds_and_with_R_is_unique():
    """T(2,5)：S だけで多面体が有界（整数解 2 個、正解を含む）、S ＋ R で正解 1 個（約 10 秒）。"""
    r = HU.run("T2_5", use_R=False, log=quiet)
    assert r["true_feasible"] and r["bounded"] and r["count"] == 2 and r["contains_true"]
    r = HU.run("T2_5", use_R=True, log=quiet)
    assert r["bounded"] and r["count"] == 1 and r["contains_true"]
