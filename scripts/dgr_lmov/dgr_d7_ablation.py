"""D7（rCS_{[2]} ≥ 0）で、どの条件が一意性に効いているかの切り分け（docs/dgr_lmov/research_plan.md）。

基準：V1（scripts/dgr_lmov/dgr_d7_pos.py）。ここから外す：hook（変種 a）、EGP、genus 表示（σ 不変性・N ≥ 0）、\\hat f の正値性（h, N の符号の制約）、\\hat f の Laurent 性（2 箱の公式そのもの）。
各組で、全体の符号 σ の 4 通り × 乱数の目的 6 回で異なる rCS_{[2]} の解を数える。記録：logs/dgr_lmov/dgr_d7_ablation.log。
使い方：PYTHONPATH=src python scripts/dgr_lmov/dgr_d7_ablation.py 結び目...
"""
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "framed_knot"))
from dgr_d7_construct import build, parity_class, solution, solve, unrefined  # noqa: E402
from dgr_d7_pos import SIGMAS, load  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
CASES = {"V1（全て）": dict(variant="b", egp=True, genus=True, lmov=True),
         "hook なし": dict(variant="a", egp=True, genus=True, lmov=True),
         "EGP なし": dict(variant="b", egp=False, genus=True, lmov=True),
         "genus 表示なし": dict(variant="b", egp=True, genus=False, lmov=True),
         "LMOV の正値性なし": dict(variant="b", egp=True, genus=True, lmov=False),
         "EGP と LMOV の正値性なし": dict(variant="b", egp=False, genus=False, lmov=False),
         "Laurent なし（\\hat f を使わない）": dict(variant="b", egp=True, genus=False, lmov=False, laurent=False),
         "Laurent と EGP なし": dict(variant="b", egp=False, genus=False, lmov=False, laurent=False)}


def main():
    logf = open(ROOT / "logs/dgr_lmov/dgr_d7_ablation.log", "a")
    for name in sys.argv[1:]:
        P1, P2 = load(name)
        ts = [k[2] for k in P2]
        for case, kw in CASES.items():
            t0 = time.time(); sols = []
            for sig in SIGMAS:
                lin, mons = build(P1, unrefined(P2), (min(ts) - 2, max(ts) + 3), kw["variant"], sig, genus=kw["genus"], egp=kw["egp"],
                                  parity=parity_class(P1) if kw["genus"] else None, laurent=kw.get("laurent", True))
                if solve(lin, rcs_pos=True, lmov_pos=kw["lmov"]).x is None:
                    continue
                rng = np.random.default_rng(0)
                for _ in range(6):
                    obj = np.zeros(len(lin.col))
                    for m in mons:
                        obj[lin.col[("x", m)]] = rng.normal()
                    r = solve(lin, objective=obj, rcs_pos=True, lmov_pos=kw["lmov"])
                    if r.x is not None:
                        s = solution(lin, r, mons)
                        if s not in sols:
                            sols.append(s)
                if not kw["lmov"]:
                    break                                   # 正値性なしでは σ は意味をもたない
            line = f"{name} {case}: 異なる解 {len(sols)} 個、正解を含む {P2 in sols}（{time.time() - t0:.0f} 秒）"
            print(line, flush=True); logf.write(line + "\n"); logf.flush()


if __name__ == "__main__":
    main()
