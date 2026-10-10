"""段階 D1（docs/research_plan_dgr_lmov.md）：トーラス結び目の 4 箱以上で、因子 D の割り切り・商 R_ρ の正値性・d_{1|1} 後の最大 genus を判定する。

入力：data/highbox/<結び目>_N<N>.json.xz（scripts/highbox_torus.py の出力、git の外。既定は本体のチェックアウトのもの）。
出力：data/dgr_lmov/highbox/<結び目>_N<N>.json（ρ ごとの判定）。保存済みは読み飛ばす。ログは logs/dgr_highbox.log に追記。
使い方：PYTHONPATH=src python scripts/dgr_highbox.py [結び目_N ...]（省略時は入力の全て、小さい順）
"""
import json
import lzma
import os
import resource
import sys
import time
from pathlib import Path

from rcsknots.bps import g_expand, to_bold
from rcsknots.dgr import d11_top_genus, divides, negatives, quotient

ROOT = Path(__file__).resolve().parents[1]
SRC = Path(os.environ.get("HIGHBOX_DIR", "/Users/masaya/projects/ai_research/refined_cs/data/highbox"))
OUT = ROOT / "data/dgr_lmov/highbox"
PAIRS = [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2), (2, 1), (1, 2)]


def a_shift(N):
    """conventions.md §18：偶数 N は 0、N ≡ 1 (mod 4) は −1、N ≡ 3 (mod 4) は +1（A の単位）。"""
    return 0 if N % 2 == 0 else (-1 if N % 4 == 1 else 1)


def run(stem, log):
    out = OUT / f"{stem}.json"
    if out.exists():
        return
    t0 = time.time()
    d = json.loads(lzma.open(SRC / f"{stem}.json.xz").read())
    N = d["N"]
    res = {}
    for rho, r in d["rho"].items():
        hat = {tuple(map(int, k.split(","))): int(v) for k, v in r["hat_f"].items()}
        b = to_bold(hat, a_shift(N))
        divs = {f"{n}|{m}": divides(b, n, m) for n, m in PAIRS}
        R = quotient(b, N)
        NR = g_expand(R)
        sR = sum(1 for v in R.values() if v < 0) <= len(R) / 2
        sN = sum(1 for v in NR.values() if v < 0) <= len(NR) / 2
        res[rho] = {"divides": divs, "R_spin_neg": negatives(R.values()), "R_spin_terms": len(R),
                    "R_genus_neg": negatives(NR.values()), "R_genus_terms": len(NR), "same_sign": sR == sN,
                    "R_genus_max": max(k[0] for k in NR), "R_d11_top_genus": d11_top_genus(R)}
        log(f"{stem} {rho:>5}: 因子 {''.join('1' if v else '0' for v in divs.values())}  "
            f"R spin 負 {res[rho]['R_spin_neg']}/{len(R)}  genus 負 {res[rho]['R_genus_neg']}/{len(NR)}  "
            f"符号一致 {sR == sN}  genus 最大 {res[rho]['R_genus_max']}  d11 後 {res[rho]['R_d11_top_genus']}")
    sec = time.time() - t0
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20
    OUT.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"knot": d["knot"], "N": N, "seconds": round(sec, 1), "max_rss_MB": round(rss),
                               "rho": res}, ensure_ascii=False, indent=1))
    log(f"{stem}: {sec:.0f} 秒、最大 RSS {rss:.0f} MB")


def main():
    stems = sys.argv[1:] or sorted((p.name[:-8] for p in SRC.glob("*.json.xz")), key=lambda s: (int(s.split("_N")[1]), s))
    (ROOT / "logs").mkdir(exist_ok=True)
    logf = open(ROOT / "logs/dgr_highbox.log", "a")

    def log(s):
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    for s in stems:
        run(s, log)


if __name__ == "__main__":
    main()
