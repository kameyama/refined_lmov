"""第 4 期 V0：K10（2 段の反復ケーブル）の |ρ| ≤ 3 の refined LMOV 不変量の判定（3 箱は 2026-10-07 の高速化の後に追加）。

使い方：python scripts/v0_K10_lmov.py [--workers N]。\\hat f は data/lmov/K10_{n}.json、判定は標準出力。
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from u4_run import judge  # noqa: E402
from rcsknots.extract import reformulated_exact  # noqa: E402

COLORS = {"1": (1,), "2": (2,), "11": (1, 1), "3": (3,), "111": (1, 1, 1), "21": (2, 1)}
workers = int(sys.argv[2]) if sys.argv[1:2] == ["--workers"] else 1
reduced = {lam: {tuple(map(int, k.split(","))): v for k, v in json.loads(Path(f"data/cables/K10_{c}.json").read_text())["poly"].items()}
           for c, lam in COLORS.items() if Path(f"data/cables/K10_{c}.json").exists()}
for n in (1, 2, 3):
    path = Path(f"data/lmov/K10_{n}.json")
    if path.exists():
        continue                                # 計算済み（1・2 箱は 2026-10-06）
    cache = Path(f"data/cache/lmov/K10_{n}") if workers > 1 else None
    F = reformulated_exact(reduced, n, log=lambda s: print(f"    {s}", flush=True), workers=workers, cache=cache)
    path.write_text(json.dumps({"knot": "K10", "n": n, "variables": "lowercase (A, Q, T); no a-shift",
                                "f": {f"{k[0]}_{''.join(map(str, k[1]))}": {",".join(map(str, e)): c for e, c in sorted(v.items())}
                                      for k, v in F.items()}}))
    for (kind, rho), poly in sorted(F.items()):
        if kind == "hat":
            print(f"K10 ρ={list(rho)}: {judge(poly, n)}", flush=True)
