"""ステージ T1：結び目表（yui-kr、11 交点以下の全素結び目）の [1] の refined LMOV 不変量を判定する。

出力：data/framed_knot/knot_table/judge1.json（結び目ごとの判定と \\hat N）、data/framed_knot/knot_table_summary.md（集計）。
使い方：python scripts/framed_knot/knot_table_t1.py
"""
import csv
import json
import time
from pathlib import Path

from rcsknots.knot_table import DATA, delta_gradings, judge1, load_homology, superpoly

# Rolfsen の表の非交代結び目（10 交点以下）。11 交点は名前の a/n で分かる。
NONALT = {"8_19", "8_20", "8_21"} | {f"9_{i}" for i in range(42, 50)} | {f"10_{i}" for i in range(124, 166)}


def crossings(name):
    return int(name.split("_")[0].rstrip("an"))


def alternating(name):
    if name.startswith("11"):
        return name.startswith("11a")
    return name not in NONALT


def main():
    t0 = time.time()
    braids = {r["name"]: json.loads(r["braid"]) for r in csv.DictReader((DATA / "braid.csv").open())}
    out = {}
    for name, hom in load_homology().items():
        P = superpoly(hom)
        r = judge1(P)
        w = braids[name]
        rec = {
            "crossings": crossings(name),
            "alternating": alternating(name),
            "positive_braid": all(g > 0 for g in w) or all(g < 0 for g in w),
            "rank": sum(P.values()),
            "delta": delta_gradings(P),
            "before_terms": len(r["before"]),
            "before_neg": r["before_neg"],
            "closes": r["closes"],
        }
        if r["closes"]:
            rec.update(after_terms=len(r["after"]), after_neg=r["after_neg"], parity_aligned=r["parity_aligned"],
                       N={f"{g},{b},{F}": v for (g, b, F), v in sorted(r["after"].items())})
        else:
            rec["error"] = r["error"]
        out[name] = rec
    Path("data/framed_knot/knot_table").mkdir(exist_ok=True)
    Path("data/framed_knot/knot_table/judge1.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    Path("data/framed_knot/knot_table_summary.md").write_text(summary(out))
    print(f"{len(out)} knots, {time.time() - t0:.1f} s")


def summary(out):
    rows = []
    for cls, pred in [("thin・交代", lambda r: len(r["delta"]) == 1 and r["alternating"]),
                      ("thin・非交代", lambda r: len(r["delta"]) == 1 and not r["alternating"]),
                      ("thick", lambda r: len(r["delta"]) > 1)]:
        rs = [r for r in out.values() if pred(r)]
        rows.append(f"| {cls} | {len(rs)} | {sum(r['closes'] for r in rs)} | {sum(r['before_neg'] == 0 for r in rs)} | "
                    f"{sum(r.get('after_neg') == 0 for r in rs)} | {sum(bool(r.get('parity_aligned')) for r in rs)} |")
    aligned = ", ".join(n for n, r in out.items() if r.get("parity_aligned"))
    thick = ", ".join(n for n, r in out.items() if len(r["delta"]) > 1)
    return "\n".join([
        "# 結び目表の [1] の refined LMOV 不変量の判定（ステージ T1）",
        "",
        "入力：yui-kr（Nakagane–Sano arXiv:2111.00388）の reduced HOMFLY-PT ホモロジー、11 交点以下の全素結び目。生成：`scripts/framed_knot/knot_table_t1.py`。",
        "判定：\\hat f_{[1]} = (1 + 𝐚²𝐭)/(𝐚²𝐭)·P の係数（G 展開の前）と、G 展開の後の \\hat N。全体の符号は反転していない。",
        "",
        "| 分類 | 結び目数 | G 展開が閉じる | 前：負なし | 後：負なし | 2J_r の偶奇がそろう |",
        "|---|---|---|---|---|---|",
        *rows,
        "",
        f"2J_r の偶奇がそろう結び目：{aligned}",
        "",
        f"thick な結び目：{thick}",
        "",
    ])


if __name__ == "__main__":
    main()
