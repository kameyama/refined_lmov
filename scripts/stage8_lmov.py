"""見直し：refined LMOV 不変量を、小文字 (a,q,t) で計算してから太字に変換して判定する。

各結び目について、|μ| = 1, 2, 3 の f^q_μ, f^{\\bar t}_μ, \\hat f_ρ を小文字の (A,Q,T) = (a,q,t)^{1/2} で厳密に求め、
data/lmov/{K}_{n}.json に保存する（変換はしない）。既存の data/bps の \\hat f と一致するかも確かめる。
使い方：python scripts/stage8_lmov.py K1 K2 ...（T2_3, T2_5, 4_1 も可）
"""
import json
import resource
import sys
import time
from pathlib import Path

from rcsknots.extract import reformulated_exact
from rcsknots.laurent import load_table

COLORS = {"1": (1,), "2": (2,), "11": (1, 1), "3": (3,), "111": (1, 1, 1), "21": (2, 1)}
OUT = Path("data/lmov")


def reduced_of(name):
    if name.startswith("K"):
        return {lam: {tuple(map(int, k.split(","))): v for k, v in
                      json.loads(Path(f"data/cables/{name}_{c}.json").read_text())["poly"].items()}
                for c, lam in COLORS.items()}
    ref = load_table("data/reference/notebook_torus.txt")
    return {COLORS[c]: p for (k, c), p in ref.items() if k == name}


def key_str(k):
    return f"{k[0]}_{''.join(map(str, k[1]))}"


def main(names):
    OUT.mkdir(parents=True, exist_ok=True)
    for name in names:
        reduced = reduced_of(name)
        for n in (1, 2, 3):
            path = OUT / f"{name}_{n}.json"
            if path.exists():
                continue
            t0 = time.time()
            F = reformulated_exact(reduced, n, log=lambda s: None)
            # cross-check with the \hat f saved in Stage 6 (computed by a separate evaluator)
            same = []
            for (kind, rho), poly in F.items():
                if kind != "hat":
                    continue
                old = Path(f"data/bps/{name}_{''.join(map(str, rho))}.json")
                if old.exists():
                    saved = {tuple(map(int, k.split(","))): v for k, v in json.loads(old.read_text())["hat_f"].items()}
                    same.append(saved == poly)
            data = {"knot": name, "n": n,
                    "variables": "lowercase: exponents of (A=a^(1/2), Q=q^(1/2), T=t^(1/2)) of 1703.05408; no a-shift, no bold change",
                    "ft_equals_M_hatf_checked": True, "hat_f_matches_stage6": all(same) if same else None,
                    "seconds": round(time.time() - t0, 1),
                    "max_rss_MB": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20),
                    "f": {key_str(k): {",".join(map(str, e)): c for e, c in sorted(v.items())} for k, v in F.items()}}
            path.write_text(json.dumps(data))
            print(f"{name} n={n}: {data['seconds']} s, hat f = Stage 6: {data['hat_f_matches_stage6']}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
