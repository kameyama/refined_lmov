"""段階 D1 の集計：data/dgr_lmov/highbox/*.json → data/dgr_lmov/highbox_summary.md。

d_{1|1} 後の最大 genus を、予想 |ρ|(g−1) − [ρ = [n]]（トーラス結び目は Σ = g なので行だけが 1 落ちる、docs/dgr_lmov/report.md D2）と比べる。
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GENUS = {"T2_3": 1, "T2_5": 2, "T2_7": 3, "T2_9": 4, "T3_4": 3, "T3_5": 4, "T3_7": 6, "T4_5": 6}


def main():
    rows = ["# トーラス結び目の 4 箱以上での DGR 型の構造（scripts/dgr_highbox.py、【計算】）", "",
            "因子は d_{0|0}, d_{1|0}, d_{0|1}, d_{1|1}, d_{2|0}, d_{0|2}, d_{2|1}, d_{1|2} の順に 1 = 割り切れる。"
            "予想は d_{1|1} 後の最大 genus = |ρ|(g−1) − [ρ が 1 行]（0 未満は特殊化が 0）。", "",
            "| 結び目 | N | ρ の数 | 因子の型 | R の負（spin / genus） | 予想との不一致 | 秒 |", "|---|---|---|---|---|---|---|"]
    files = sorted((ROOT / "data/dgr_lmov/highbox").glob("*.json"), key=lambda p: (json.loads(p.read_text())["N"], p.name))
    for p in files:
        d = json.loads(p.read_text())
        K, N, g = d["knot"], d["N"], GENUS[d["knot"]]
        types = sorted({"".join("1" if v else "0" for v in r["divides"].values()) for r in d["rho"].values()})
        neg = (sum(r["R_spin_neg"] for r in d["rho"].values()), sum(r["R_genus_neg"] for r in d["rho"].values()))
        bad = []
        for rho, r in d["rho"].items():
            pred = N * (g - 1) - (1 if len(rho) == 1 else 0)
            got = r["R_d11_top_genus"]
            if (got is None and pred >= 0) or (got is not None and got != pred):
                bad.append(f"{rho}:{got}≠{pred}")
        rows.append(f"| {K} | {N} | {len(d['rho'])} | {', '.join(types)} | {neg[0]} / {neg[1]} | {' '.join(bad) or 'なし'} | {d['seconds']} |")
    (ROOT / "data/dgr_lmov/highbox_summary.md").write_text("\n".join(rows) + "\n")
    print("\n".join(rows))


if __name__ == "__main__":
    main()
