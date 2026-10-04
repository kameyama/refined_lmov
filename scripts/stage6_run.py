"""Stage 6：ケーブル結び目の refined CS 不変量（全 6 彩色）と refined BPS 不変量を計算して保存する。

使い方：python scripts/stage6_run.py K1 K2 ...
  refined CS：data/cables/{K}_{color}.json（既にあれば読む）
  BPS       ：data/bps/{K}_{rho}.json と .txt（a 次数シフトは nontorus 規約）
"""
import resource
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage5_cables import compute, COLORS          # noqa: E402
from stage3_bps import run as run_bps              # noqa: E402
from renormalize_cables import main as renormalize  # noqa: E402
import json                                         # noqa: E402


def main(knots):
    for name in knots:
        t0 = time.time()
        for cname in COLORS:
            d = compute(name, cname, log=lambda s: None)
            print(f"{name} {cname}: rCS {d['n_terms']} terms ({d['seconds_total']} s)", flush=True)
        renormalize([name])                       # (N2'): M_λ = M_□^{|λ|} T_λ^{p₂q₂}
        reduced = {}
        for cname, lam in COLORS.items():
            d = json.loads(Path(f"data/cables/{name}_{cname}.json").read_text())
            reduced[lam] = {tuple(map(int, k.split(","))): v for k, v in d["poly"].items()}
        t1 = time.time()
        run_bps(name, reduced, "nontorus")
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20
        print(f"{name}: rCS {t1 - t0:.0f} s, BPS {time.time() - t1:.0f} s, max RSS {rss:.0f} MB", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
