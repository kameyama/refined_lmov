"""Stage 3：reduced refined 不変量から \\hat f と refined BPS 不変量を求め、data/bps/ に保存し、1703.05408 の表と比べる。

使い方：python scripts/stage3_bps.py KNOT CONVENTION
  KNOT は data/reference/notebook_torus.txt の結び目名（T2_3, 4_1 など）
  CONVENTION は a 次数シフトの規約（torus / nontorus、src/rcsknots/extract.py）
"""
import json
import sys
import time

from rcsknots.extract import hat_f_exact, bps_numbers, A_SHIFT
from rcsknots.laurent import load_table
from rcsknots.tables import normalize_sign, save

COLORS = {"1": (1,), "2": (2,), "11": (1, 1), "3": (3,), "111": (1, 1, 1), "21": (2, 1)}


def run(knot, reduced, convention, out_dir="data/bps", paper=None):
    results = {}
    for n in (1, 2, 3):
        if not all(c in reduced for c in [(1,)] + [l for l in COLORS.values() if sum(l) == n]):
            continue
        t0 = time.time()
        H = hat_f_exact(reduced, n, log=lambda s: None)
        for rho, hf in H.items():
            raw = bps_numbers(hf, rho, convention)
            sigma, N = normalize_sign(raw)
            meta = {"knot": knot, "a_shift_convention": convention, "a_shift_A_units": A_SHIFT[convention][n],
                    "overall_sign_applied": sigma, "ft_consistency_checked": True,
                    "variables_hat_f": "exponents of (A=a^(1/2), Q=q^(1/2), T=t^(1/2)) of 1703.05408, before the a-shift",
                    "N_keys": "g,beta,2J_r", "seconds_for_size_n": round(time.time() - t0, 1)}
            stem = f"{out_dir}/{knot}_{''.join(map(str, rho))}"
            save(stem, rho, hf, N, meta)
            results[rho] = N
            msg = f"{knot} {rho}: {len(N)} terms, nonneg={all(v >= 0 for v in N.values())}"
            if paper is not None:
                tb = next((x for x in paper if x["knot"] == knot and x["rho"] == [",".join(map(str, rho))]), None)
                if tb:
                    msg += f", matches KN table: {({f'{g},{b},{F}': v for (g, b, F), v in N.items()} == tb['data'])}"
            print(msg, flush=True)
    return results


if __name__ == "__main__":
    knot, convention = sys.argv[1], sys.argv[2]
    ref = load_table("data/reference/notebook_torus.txt")
    reduced = {COLORS[c]: p for (k, c), p in ref.items() if k == knot}
    paper = json.load(open("data/reference/kn_tables.json"))
    run(knot, reduced, convention, paper=paper)
