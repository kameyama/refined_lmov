"""S4：\hat f_ρ の「微分の点」の走査。太字で 𝐚² → ±𝐪^m 𝐭^n（|m|, |n| ≤ 8）に特殊化したときの項数を全ての点で数え、縮退の強い点を列挙する。

出力：data/dgr_lmov/lmovstruct/loci.json、logs/dgr_lmov/lmovstruct_loci.log。使い方：PYTHONPATH=src python scripts/dgr_lmov/lmovstruct_loci.py [knot ...]
"""
import glob
import json
import sys
import time
from pathlib import Path

from rcsknots.dgr_lmov.lmovstruct import load_bps, specialize_bold

OUT = Path("data/dgr_lmov/lmovstruct")
RANGE = range(-8, 9)


def scan(bold):
    res = {}
    for sign in (1, -1):
        for m in RANGE:
            for n in RANGE:
                res[(sign, m, n)] = len(specialize_bold(bold, m, n, sign))
    return res


def main():
    knots = sys.argv[1:] or ["T2_3", "T2_5", "T3_4", "4_1", "K1", "K3"]
    OUT.mkdir(parents=True, exist_ok=True)
    log = open("logs/dgr_lmov/lmovstruct_loci.log", "a")
    out = {}
    for knot in knots:
        for path in sorted(glob.glob(f"data/framed_knot/bps/{knot}_*.json")):
            meta, lower, bold = load_bps(path)
            t0 = time.time()
            res = scan(bold)
            generic = max(res.values())
            # 縮退の強い点：項数が最大の 1/3 以下、または 0
            special = sorted(((v, s, m, n) for (s, m, n), v in res.items() if v <= generic / 3), key=lambda x: (x[0], x[2], x[3]))
            key = f"{meta['knot']} {meta['rho']}"
            out[key] = {"terms": len(bold), "generic_terms": generic,
                        "special": [{"sign": s, "m": m, "n": n, "terms": v, "a_lower": f"a = {'-' if s < 0 else ''}q^{m/2} t^{(n+1)/2}" if True else ""} for v, s, m, n in special]}
            line = f"{key:16s} terms={len(bold):7d} generic={generic:7d} special={[(s, m, n, v) for v, s, m, n in special[:12]]} ({time.time() - t0:.1f}s)"
            print(line); log.write(line + "\n"); log.flush()
    (OUT / "loci.json").write_text(json.dumps(out, ensure_ascii=False, indent=0))


if __name__ == "__main__":
    main()
