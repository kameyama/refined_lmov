"""S2：保存済みの \\hat f_ρ（data/bps/ の |ρ| ≤ 3 と data/highbox/ の 4〜9 箱）で、既約表示の非負性・δ の台・特殊化の縮退を調べる。

出力：data/lmovstruct/survey.jsonl（1 行 1 件）、data/lmovstruct/survey.md（表）、logs/lmovstruct_survey.log。
使い方：PYTHONPATH=src python scripts/lmovstruct_survey.py [--highbox DIR] [--skip-xz]
"""
import argparse
import glob
import json
import sys
import time
from pathlib import Path

from rcsknots.lmovstruct import (delta_values, is_W_invariant, load_bps, load_highbox, specialize_lower, spin_count)

OUT = Path("data/lmovstruct")
LOCI = [(N, M) for N in range(-3, 4) for M in (0,)] + [(0, M) for M in (-2, -1, 1, 2)] + [(1, -1), (1, 1), (2, -1), (2, 1), (-1, 1), (-1, -1)]


def analyze(name, rho, lower, bold, extra=None):
    t0 = time.time()
    rec = {"knot": name, "rho": rho, "terms": len(bold), "W_invariant": is_W_invariant(bold),
           "nonneg": all(v >= 0 for v in bold.values())}
    if rec["W_invariant"]:
        rec["spin_terms"], rec["spin_neg"] = spin_count(bold)
    ds = delta_values(bold)
    rec["delta_min"], rec["delta_max"] = ds[0], ds[-1]
    rec["spec"] = {f"t^{N}q^{M}": len(specialize_lower(lower, N, M)) for N, M in LOCI}
    if extra:
        rec.update(extra)
    rec["seconds"] = round(time.time() - t0, 2)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--highbox", default="/Users/masaya/projects/ai_research/refined_cs/data/highbox")
    ap.add_argument("--skip-xz", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    log = open("logs/lmovstruct_survey.log", "a")

    def say(s):
        print(s); log.write(s + "\n"); log.flush()

    say(f"=== lmovstruct survey {time.strftime('%Y-%m-%d %H:%M:%S')} ===")
    recs = []
    files = sorted(glob.glob("data/bps/*.json")) + ([] if args.skip_xz else sorted(glob.glob("data/bps/*.json.xz")))
    for path in files:
        meta, lower, bold = load_bps(path)
        rec = analyze(meta["knot"], meta["rho"], lower, bold, {"source": path})
        recs.append(rec)
        say(f"{rec['knot']:6s} {str(rec['rho']):12s} terms={rec['terms']:8d} W={rec['W_invariant']} nonneg={rec['nonneg']} "
            f"spin_neg={rec.get('spin_neg')} delta=[{rec['delta_min']},{rec['delta_max']}] a=t:{rec['spec']['t^1q^0']} a=t^2:{rec['spec']['t^2q^0']} ({rec['seconds']}s)")
    hb = sorted(glob.glob(f"{args.highbox}/*.json.xz")) if args.highbox else []
    for path in hb:
        d, polys = load_highbox(path)
        for rho, (lower, bold) in polys.items():
            rec = analyze(d["knot"], [int(c) for c in rho.strip("[]").split(",")] if "[" in rho else rho, lower, bold,
                          {"source": path, "N": d["N"]})
            recs.append(rec)
            say(f"{rec['knot']:6s} {str(rec['rho']):12s} terms={rec['terms']:8d} W={rec['W_invariant']} nonneg={rec['nonneg']} "
                f"spin_neg={rec.get('spin_neg')} delta=[{rec['delta_min']},{rec['delta_max']}] a=t:{rec['spec']['t^1q^0']} a=t^2:{rec['spec']['t^2q^0']} ({rec['seconds']}s)")
    with open(OUT / "survey.jsonl", "w") as fp:
        for r in recs:
            fp.write(json.dumps(r, ensure_ascii=False) + "\n")
    # 表
    lines = ["# 既約表示・特殊化の調査（scripts/lmovstruct_survey.py）", "",
             "| 結び目 | ρ | 項数 | W 不変 | 非負 | 既約項数 | 既約負 | δ の範囲 | a=1 | a=t | a=t² | a=t³ | a=t⁻¹ | a=q | a=q² |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in recs:
        s = r["spec"]
        lines.append(f"| {r['knot']} | {r['rho']} | {r['terms']} | {r['W_invariant']} | {r['nonneg']} | {r.get('spin_terms', '')} | {r.get('spin_neg', '')} | "
                     f"[{r['delta_min']}, {r['delta_max']}] | {s['t^0q^0']} | {s['t^1q^0']} | {s['t^2q^0']} | {s['t^3q^0']} | {s['t^-1q^0']} | {s['t^0q^1']} | {s['t^0q^2']} |")
    (OUT / "survey.md").write_text("\n".join(lines) + "\n")
    n_pos = sum(1 for r in recs if r["nonneg"])
    n_spin = sum(1 for r in recs if r.get("spin_neg") == 0)
    both = sum(1 for r in recs if r["nonneg"] and r.get("spin_neg") == 0)
    say(f"件数 {len(recs)}、W 不変 {sum(r['W_invariant'] for r in recs)}、G 展開前に非負 {n_pos}、既約表示で非負 {n_spin}、両方 {both}")
    say(f"a=t で 0（|ρ|≥2）: {sum(1 for r in recs if sum(r['rho']) >= 2 and r['spec']['t^1q^0'] == 0)} / {sum(1 for r in recs if sum(r['rho']) >= 2)}")


if __name__ == "__main__":
    sys.exit(main())
