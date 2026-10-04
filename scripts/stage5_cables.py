"""Stage 5/6：ケーブル結び目の refined CS 不変量（reduced）を GN の構成で計算し、(N1)(N2) で正規化して data/cables/ に保存する。

使い方：python scripts/stage5_cables.py K1 [colors...]   （colors 省略時は 1 2 11 3 21 111）
"""
import json
import resource
import sys
import time
from pathlib import Path

from rcsknots.exact import reconstruct, reconstruct_points
from rcsknots.homfly_num import homfly_evaluator
from rcsknots.normalize import normalize_exact

COLORS = {"1": (1,), "2": (2,), "11": (1, 1), "3": (3,), "111": (1, 1, 1), "21": (2, 1)}
# (p2, q2)-cable of T(p1, q1); q2 in the Seifert framing of the companion
KNOTS = {"K1": (2, 3, 2, 13), "K2": (2, 3, 2, 1), "K3": (2, 3, 2, 15), "K4": (2, 3, 2, 3),
         "K5a": (2, 3, 2, 5), "K5b": (2, 3, 2, 7), "K5c": (2, 3, 2, -1), "K6a": (2, 5, 2, 21), "K6b": (2, 5, 2, 1)}
OUT = Path("data/cables")


def slopes(p1, q1, p2, q2):
    """GN slopes (DGR orientation): companion (p1, −q1); pattern (p2, −s2) with the DAHA Newton pair s2 = q2 − p1 q1 p2."""
    return (p2, -(q2 - p1 * q1 * p2)), (p1, -q1)


def compute(name, cname, log=print):
    path = OUT / f"{name}_{cname}.json"
    if path.exists():
        return json.loads(path.read_text())
    p1, q1, p2, q2 = KNOTS[name]
    lam = COLORS[cname]
    inner, outer = slopes(p1, q1, p2, q2)
    t0 = time.time()
    raw = reconstruct(lambda G: G.cable_reduced_apoly(lam, inner, outer), log=log)
    t1 = time.time()
    hraw = reconstruct_points(homfly_evaluator(lam, [(p1, q1), (p2, q2)]), log=lambda s: None)
    hom = {(ea, eq): c for (_, ea, eq), c in hraw.items()}
    norm, mono = normalize_exact(raw, hom)
    data = {"knot": name, "cable": {"companion": [p1, q1], "pattern": [p2, q2]},
            "daha_newton_pairs": {"r": [p1, p2], "s": [q1, q2 - p1 * q1 * p2]},
            "algebraic": q2 > p1 * q1 * p2, "gn_slopes": {"pattern": list(inner), "companion": list(outer)},
            "color": list(lam), "reduced": True,
            "variables": "exponents of (A=a^(1/2), Q=q^(1/2), T=t^(1/2)) of 1703.05408",
            "framing": "0; normalization (N1)+(N2) of docs/conventions.md §5 (hypothesis for cables)",
            "monomial_applied": {"sign": mono[0], "a": str(mono[1]), "q": str(mono[2]), "t": str(mono[3])},
            "seconds_refined": round(t1 - t0, 1), "seconds_total": round(time.time() - t0, 1),
            "max_rss_MB": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20, 1),
            "n_terms": len(norm),
            "poly": {",".join(map(str, k)): v for k, v in sorted(norm.items())}}
    OUT.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=1))
    return data


if __name__ == "__main__":
    name = sys.argv[1]
    for cname in (sys.argv[2:] or list(COLORS)):
        d = compute(name, cname)
        print(name, cname, d["n_terms"], "terms", d["seconds_total"], "s", d["max_rss_MB"], "MB", flush=True)
