"""第 4 期 V0：K10（2 段の反復ケーブル）の |ρ| ≤ 2 の refined LMOV 不変量の判定（[3] 以降は見積もり 24 時間超のため保留）。

使い方：python scripts/v0_K10_lmov.py [--workers N]。\\hat f は data/lmov/K10_{n}.json、判定は標準出力。
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from u4_run import judge  # noqa: E402
from rcsknots.extract import reformulated_exact  # noqa: E402

COLORS = {"1": (1,), "2": (2,), "11": (1, 1)}
workers = int(sys.argv[2]) if sys.argv[1:2] == ["--workers"] else 1
reduced = {lam: {tuple(map(int, k.split(","))): v for k, v in json.loads(Path(f"data/cables/K10_{c}.json").read_text())["poly"].items()}
           for c, lam in COLORS.items()}
for n in (1, 2):
    path = Path(f"data/lmov/K10_{n}.json")
    cache = Path(f"data/cache/lmov/K10_{n}") if workers > 1 else None
    F = reformulated_exact(reduced, n, log=lambda s: print(f"    {s}", flush=True), workers=workers, cache=cache)
    path.write_text(json.dumps({"knot": "K10", "n": n, "variables": "lowercase (A, Q, T); no a-shift",
                                "f": {f"{k[0]}_{''.join(map(str, k[1]))}": {",".join(map(str, e)): c for e, c in sorted(v.items())}
                                      for k, v in F.items()}}))
    for (kind, rho), poly in sorted(F.items()):
        if kind == "hat":
            print(f"K10 ρ={list(rho)}: {judge(poly, n)}", flush=True)
