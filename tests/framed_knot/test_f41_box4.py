"""4_1 の 4 箱の弱い判定（scripts/framed_knot/f41_box4.py）の定式化を、正解の分かる T(2,3) で確かめる機能テスト。"""
import json
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts/framed_knot"))
import f41_box4 as F  # noqa: E402
import highbox_torus as H  # noqa: E402
from rcsknots.exact import reconstruct  # noqa: E402
from rcsknots.partitions import conjugate, n_stat, partitions  # noqa: E402


def trefoil(lam):
    """T(2,3) の rCS_λ（GN の値 × conventions.md §13 の正規化）。"""
    raw = reconstruct(lambda G: G.torus_reduced_apoly(lam, 2, -3), log=lambda s: None)
    k = sum(lam)
    eq = Fraction(k) + Fraction(k, 2) + 6 * n_stat(conjugate(lam))
    et = -Fraction(3 * k, 2) - 6 * n_stat(lam)
    return {(int(2 * (x + k)), int(2 * (y + eq)), int(2 * (z + et))): (-1) ** k * c for (x, y, z), c in raw.items()}


def load(c):
    d = json.loads(Path(f"data/framed_knot/torus/T2_3_{c}.json").read_text())
    return {tuple(map(int, k.split(","))): v for k, v in d["poly"].items()}


def test_linear_decomposition_and_conditions_on_trefoil():
    red = {lam: load(c) for c, lam in {"1": (1,), "2": (2,), "11": (1, 1), "21": (2, 1), "3": (3,), "111": (1, 1, 1)}.items()}
    for lam in [(4,), (2, 2), (1, 1, 1, 1), (3, 1), (2, 1, 1)]:
        red[lam] = trefoil(lam)
    r1, r2, r11, r3, r31 = red[(1,)], red[(2,)], red[(1, 1)], red[(3,)], red[(3, 1)]
    # 条件：転置、指数増大、U(2)、窓
    assert F.tr(r31) == red[(2, 1, 1)]
    prods = [F.mul(r11, r1, r1), F.mul(r3, r1), F.mul(r2, r11)]
    for (which, prod) in zip("QTU", prods):
        assert all(F.spec(r31, which, s) == F.spec(prod, which, s) for s in (1, -1))
    assert set(r31) <= set(F.window(prods, F.t_eq_q(r31), 0))
    # 一次の分解：D·\hat f_ρ = Σ_λ (D α unk)_λρ r_λ + D·K_ρ が新しい経路の \hat f と一致
    al, Dp, _ = F.alpha()
    dk = F.dK({k: v for k, v in red.items() if k not in F.UNK}, Dp, path=None, log=lambda s: None, workers=1)
    hats, _ = H.run(2, 3, 4, workers=1, log=lambda s: None)
    for rho in partitions(4):
        num = dict(dk[rho])
        for P in (F.mul(al[(F.UNK[0], rho)], r31), F.mul(al[(F.UNK[1], rho)], F.tr(r31))):
            for k, v in P.items():
                num[k] = num.get(k, 0) + v
        assert {k: v for k, v in num.items() if v} == F.mul(hats[rho], Dp)
