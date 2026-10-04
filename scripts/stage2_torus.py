"""Stage 2：GN の構成でトーラス結び目の refined 不変量を厳密に求め、(N1)(N2) で正規化し、ノートブックと全係数で比較する。

結果を data/torus/{knot}_{color}.json に保存する（再実行時は読み込む）。
"""
import json
import sys
import time
from pathlib import Path

from rcsknots.exact import reconstruct
from rcsknots.homfly import torus_reduced
from rcsknots.laurent import load_table
from rcsknots.normalize import normalize

COLORS = {"1": (1,), "2": (2,), "11": (1, 1), "3": (3,), "111": (1, 1, 1), "21": (2, 1)}
# notebook T(m,n) (DGR orientation) ↔ GN slope (n_GN, m_GN)
KNOTS = {"T2_3": ((2, 3), (2, -3)), "T2_5": ((2, 5), (2, -5)), "T2_7": ((2, 7), (2, -7)),
         "T2_9": ((2, 9), (2, -9)), "T3_4": ((3, 4), (3, -4))}
OUT = Path("data/torus")


def compute(knot, cname):
    path = OUT / f"{knot}_{cname}.json"
    if path.exists():
        d = json.loads(path.read_text())
        return {tuple(map(int, k.split(","))): v for k, v in d["poly"].items()}, d
    (m, n), slope = KNOTS[knot]
    lam = COLORS[cname]
    t0 = time.time()
    raw = reconstruct(lambda G: G.torus_reduced_apoly(lam, *slope), log=lambda s: print(s, flush=True))
    norm, mono = normalize(raw, torus_reduced(lam, m, n, sign=-1))
    meta = {"knot": knot, "torus": [m, n], "gn_slope": list(slope), "color": list(lam),
            "variables": "A=a^(1/2), Q=q^(1/2), T=t^(1/2) of 1703.05408; exponents of (A,Q,T)",
            "reduced": True, "framing": "0 (normalization (N1)+(N2), docs/conventions.md §5)",
            "monomial_applied": {"sign": mono[0], "a": str(mono[1]), "q": str(mono[2]), "t": str(mono[3])},
            "seconds": round(time.time() - t0, 1),
            "poly": {",".join(map(str, k)): v for k, v in sorted(norm.items())}}
    OUT.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(meta, indent=1))
    return norm, meta


def main(knots):
    ref = load_table("data/reference/notebook_torus.txt")
    for knot in knots:
        for cname in COLORS:
            if (knot, cname) not in ref:
                continue
            norm, meta = compute(knot, cname)
            print(knot, cname, "match notebook:", norm == ref[(knot, cname)], meta["monomial_applied"], meta["seconds"], "s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or list(KNOTS))
