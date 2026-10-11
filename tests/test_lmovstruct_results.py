"""保存済みの refined LMOV（data/bps/、|ρ| ≤ 3）の再判定：既約表示の非負性と彩色 d_1（results の印。既定では実行しない）。

主張（docs/REPORT_dgr_lmov_spin.md）：
- \\hat f_ρ は全件で W 不変（G 展開が閉じることと同値）。
- 既約表示（SU(2)_L の指標 χ_j(w)、w = 𝐪²𝐭）の係数は、KN の G 展開の前の正値性が成り立つ (結び目, ρ) で非負。
- a = t で |ρ| ≥ 2 の \\hat f_ρ は 0、a = 1 で全ての \\hat f_ρ は 0。
"""
import glob
from pathlib import Path

import pytest

from rcsknots.lmovstruct import is_W_invariant, load_bps, specialize_lower, spin_count

ROOT = Path(__file__).resolve().parent.parent
FILES = sorted(glob.glob(str(ROOT / "data/bps/*.json")))      # .xz（K7〜K10 の 3 箱）は除く


@pytest.mark.results
@pytest.mark.parametrize("path", FILES, ids=[Path(p).stem for p in FILES])
def test_spin_positivity_matches_kn_positivity_and_colored_d1(path):
    meta, lower, bold = load_bps(path)
    assert is_W_invariant(bold)
    n_terms, n_neg = spin_count(bold)
    kn_nonneg = all(v >= 0 for v in bold.values())
    if kn_nonneg:
        assert n_neg == 0
    assert specialize_lower(lower, 0) == {}
    if sum(meta["rho"]) >= 2:
        assert specialize_lower(lower, 1) == {}
    else:
        assert len(specialize_lower(lower, 1)) == 2
