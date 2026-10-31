"""Stage 2：GN の構成でトーラス結び目の refined 不変量を厳密に求め、(N1)(N2) で正規化し、ノートブックと全係数で比較する。

結果を data/framed_knot/torus/{knot}_{color}.json に保存する（再実行時は読み込む）。
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
         "T2_9": ((2, 9), (2, -9)), "T3_4": ((3, 4), (3, -4)),
         # 第 3 期 U4（KN の表に無いもの）
         "T3_5": ((3, 5), (3, -5)), "T4_5": ((4, 5), (4, -5)), "T3_7": ((3, 7), (3, -7)), "T2_11": ((2, 11), (2, -11)),
         # refined Kauffman の研究（n = 1 の原理の検証、docs/refined_kauffman/progress.md 続き 99）で T(2,13) の 1, 2, 11 を追加（参照表には無いので compute を直接呼ぶ）
         "T2_13": ((2, 13), (2, -13)), "T2_15": ((2, 15), (2, -15))}
OUT = Path("data/framed_knot/torus")
WORKERS = 1        # 格子点の評価の並列数（scripts/framed_knot/u4_run.py が設定する）


def compute(knot, cname):
    path = OUT / f"{knot}_{cname}.json"
    if path.exists():
        d = json.loads(path.read_text())
        return {tuple(map(int, k.split(","))): v for k, v in d["poly"].items()}, d
    (m, n), slope = KNOTS[knot]
    lam = COLORS[cname]
    t0 = time.time()
    cache = Path(f"data/cache/rcs/{knot}_{cname}") if WORKERS > 1 else None
    raw = reconstruct(lambda G: G.torus_reduced_apoly(lam, *slope), log=lambda s: print(s, flush=True), workers=WORKERS, cache=cache)
    norm, mono = normalize(raw, torus_reduced(lam, m, n, sign=-1))
    meta = {"knot": knot, "torus": [m, n], "gn_slope": list(slope), "color": list(lam),
            "variables": "A=a^(1/2), Q=q^(1/2), T=t^(1/2) of 1703.05408; exponents of (A,Q,T)",
            "reduced": True, "framing": "0 (normalization (N1)+(N2), docs/framed_knot/conventions.md §5)",
            "monomial_applied": {"sign": mono[0], "a": str(mono[1]), "q": str(mono[2]), "t": str(mono[3])},
            "seconds": round(time.time() - t0, 1),
            "poly": {",".join(map(str, k)): v for k, v in sorted(norm.items())}}
    OUT.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(meta, indent=1))
    return norm, meta


def main(knots):
    ref = load_table("data/framed_knot/reference/notebook_torus.txt")
    for knot in knots:
        for cname in COLORS:
            if (knot, cname) not in ref:
                continue
            norm, meta = compute(knot, cname)
            print(knot, cname, "match notebook:", norm == ref[(knot, cname)], meta["monomial_applied"], meta["seconds"], "s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or list(KNOTS))
