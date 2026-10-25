"""段階 D8・D9（docs/dgr_lmov/research_plan.md）：\\hat f_ρ から GT1・GT2 を逆に解いた rCS_λ が、既存の rCS_λ と一致するか（往復の照合）。

点での厳密な比較（有理数の 3 点）。データは data/framed_knot/lmov（\\hat f）と、framing_check.load_reduced が読む既存の rCS。
記録：logs/dgr_lmov/dgr_inverse_check.log。使い方：PYTHONPATH=src python scripts/dgr_lmov/dgr_inverse_check.py [結び目...]
"""
import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "framed_knot"))
from framing_check import load_reduced  # noqa: E402
from rcsknots.dgr_lmov.large_n_inverse import ev_poly, rcs_from_hat, reduced  # noqa: E402
from rcsknots.lmov import RHOS, load_lmov  # noqa: E402

KNOTS = ["T2_3", "T2_5", "T2_11", "T3_4", "T3_5", "T3_7", "T4_5", "4_1", "K1", "K3", "K6a", "K7", "K8", "K9",
         "K2", "K4", "K5a", "K5b", "K5c", "K6b"]
POINTS = [(Fr(2, 3), Fr(5, 7), Fr(3, 11)), (Fr(7, 5), Fr(2, 13), Fr(9, 4)), (Fr(11, 3), Fr(4, 9), Fr(6, 5))]


def main():
    out = open(Path(__file__).resolve().parents[2] / "logs/dgr_lmov/dgr_inverse_check.log", "a")
    for K in sys.argv[1:] or KNOTS:
        t0 = time.time()
        F = load_lmov(K)
        hat = {rho: F[("hat", rho)] for rho in RHOS}
        R = load_reduced(K, 3)
        ok = {"GT1": True, "GT2": True}
        for A, Qv, Tv in POINTS:
            for which, name in (("q", "GT1"), ("t", "GT2")):
                red = reduced(rcs_from_hat(hat, A, Qv, Tv, which), A, Qv, Tv)
                ok[name] &= all(red[lam] == ev_poly(R[lam], A, Qv, Tv) for lam in red)
        line = f"{K}: GT1 から解いた rCS が既存と一致 {ok['GT1']}、GT2 から {ok['GT2']}（|λ| ≤ 3、3 点、{time.time() - t0:.0f} 秒）"
        print(line, flush=True)
        out.write(line + "\n")


if __name__ == "__main__":
    main()
