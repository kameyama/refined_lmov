"""data/framed/ の判定（scripts/framing_big.py）を (結び目, n) ごとにまとめた表を data/framed_summary.md に書く。"""
import json
import re
from collections import defaultdict
from pathlib import Path

ORDER = ["T3_5", "T2_11", "T3_7", "T4_5", "K6a", "K7", "K8", "K9", "K10"]
rows = defaultdict(list)
for f in Path("data/framed").glob("*.json"):
    d = json.loads(f.read_text())
    rows[(d["knot"], d["n"])].append(d)

lines = ["# framed refined LMOV の判定のまとめ（scripts/framing_big.py、scripts/framing_big_summary.py が生成）", "",
         "負の項の数は、全体の符号を良い方に選んだ後の数。前 = G 展開の前（太字の係数）、後 = G 展開の後（\\hat N）。", "",
         "| 結び目 | \\|ρ\\| | framing x | 件数（x × ρ） | 前の負 | 後の負 | 前の最大項数 | 1 framing の最大秒数 |",
         "|---|---|---|---|---|---|---|---|"]
for knot in ORDER + sorted({k for k, _ in rows} - set(ORDER)):
    for n in (1, 2, 3):
        ds = sorted(rows.get((knot, n), []), key=lambda d: d["framing"])
        if not ds:
            continue
        neg_b = neg_a = terms = 0
        bad = []
        for d in ds:
            for rho, msg in d["judge"].items():
                m = re.match(r"前 (\d+)/(\d+) 負、後 (\d+)/\d+ 負", msg)
                if not m:
                    bad.append(f"x={d['framing']} ρ={rho}：{msg}")
                    continue
                neg_b += int(m[1])
                neg_a += int(m[3])
                terms = max(terms, int(m[2]))
        xs = [d["framing"] for d in ds]
        xr = f"{min(xs)}..{max(xs)}" if xs == list(range(min(xs), max(xs) + 1)) else ",".join(map(str, xs))
        cnt = sum(len(d["judge"]) for d in ds)
        lines.append(f"| {knot} | {n} | {xr} | {cnt} | {neg_b} | {neg_a} | {terms} | {max(d['seconds'] for d in ds)} |")
        lines += [f"  - 判定できない：{b}" for b in bad]
Path("data/framed_summary.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
