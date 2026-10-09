"""論文の framing の表の素材：framing を振った判定を集計する。計算はしない。
入力：logs/framing_batch.log（scripts/framing_batch.sh、framing_check.py）と data/framed/（scripts/framing_big.py）。
出力：docs/paper/generated/framing_table.tex。結び目ごとに framing の範囲・最大の箱数・判定件数・genus 展開の前後の負の項の合計・\\hat f = 0 の件数。
同じ (結び目, x, ρ) が両方にあるときは 1 件と数える（data/framed/ の K1 の x = −2 は framing_batch.log と重なる）。
使い方：python scripts/paper_framing_table.py
"""
import json
import re
from collections import defaultdict
from pathlib import Path

LINE = re.compile(r"^(\S+) framing ([+-]?\d+) ρ=\[([\d, ]+)\]：(.*)$")
JUDGE = re.compile(r"前 (\d+)/\d+ 負、後 (\d+)/\d+ 負")
ORDER = ["unknot", "T2_3", "T2_5", "T2_7", "T3_4", "T3_5", "T2_11", "T3_7", "T4_5", "4_1",
         "K1", "K3", "K6a", "K7", "K8", "K9", "K10"] + [f"twist:{p}" for p in (1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6)]


def label(name):
    if name == "unknot":
        return "自明な結び目"
    if name == "4_1":
        return "$\\mathbf{4_1}$"
    if name.startswith("twist:"):
        return f"$K_{{{name.split(':')[1]}}}$"
    m = re.fullmatch(r"T(\d+)_(\d+)", name)
    return f"$T_{{{m[1]},{m[2]}}}$" if m else f"${name}$"


def ranges(xs_by_box):
    """{箱の数: x の集合} を「n 箱まで x の範囲」の文字列に。"""
    parts = []
    for n in sorted(xs_by_box):
        xs = sorted(xs_by_box[n])
        r = f"{xs[0]}..{xs[-1]}" if xs == list(range(xs[0], xs[-1] + 1)) else ",".join(map(str, xs))
        parts.append(f"{n}: ${r}$")
    return "; ".join(parts)


def main():
    seen = {}                    # (結び目, x, ρ) → (前の負, 後の負, \hat f = 0)
    for line in Path("logs/framing_batch.log").read_text().splitlines():
        m = LINE.match(line)
        if not m:
            continue
        name, x, rho, msg = m[1], int(m[2]), tuple(int(c) for c in m[3].split(",")), m[4]
        if "\\hat f = 0" in msg:
            seen[(name, x, rho)] = (0, 0, 1)
            continue
        j = JUDGE.search(msg)
        seen[(name, x, rho)] = (int(j[1]), int(j[2]), 0) if j else (10**9, 10**9, 0)
    for f in Path("data/framed").glob("*.json"):
        d = json.loads(f.read_text())
        for rho, msg in d["judge"].items():
            j = JUDGE.search(msg)
            seen[(d["knot"], d["framing"], tuple(int(c) for c in rho))] = (int(j[1]), int(j[2]), 0) if j else (10**9, 10**9, 0)
    stats = defaultdict(lambda: {"xs": defaultdict(set), "n": 0, "nb": 0, "na": 0, "zero": 0})
    for (name, x, rho), (nb, na, z) in seen.items():
        s = stats[name]
        s["xs"][sum(rho)].add(x)
        s["n"] += 1
        s["nb"] += nb
        s["na"] += na
        s["zero"] += z
    rows, total = [], 0
    for name in ORDER + sorted(set(stats) - set(ORDER)):
        if name not in stats:
            continue
        s = stats[name]
        total += s["n"]
        rows.append(f"{label(name)} & {ranges(s['xs'])} & {s['n']} & {s['nb']} & {s['na']} & {s['zero']} \\\\")
    tex = ("% scripts/paper_framing_table.py が logs/framing_batch.log と data/framed/ から生成。\n"
           f"% 判定件数の合計 {total}（重なりは 1 件と数える）。\n"
           "\\begin{tabular}{llcccc}\n\\hline\n結び目 & 箱の数: $x$ の範囲 & 件数 & 前の負 & 後の負 & $\\wh f=0$ \\\\\n\\hline\n"
           + "\n".join(rows) + "\n\\hline\n\\end{tabular}\n")
    Path("docs/paper/generated/framing_table.tex").write_text(tex)
    print(tex)


if __name__ == "__main__":
    main()
