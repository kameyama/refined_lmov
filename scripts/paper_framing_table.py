"""論文の framing の表の素材：logs/framing_batch.log（scripts/framing_batch.sh の判定ログ）を集計する。計算はしない。
出力：docs/paper/generated/framing_table.tex。結び目ごとに framing の範囲・最大の箱数・判定件数・G 展開の前後の負の項の合計・\\hat f = 0 の件数。
使い方：python scripts/paper_framing_table.py
"""
import re
from collections import defaultdict
from pathlib import Path

LINE = re.compile(r"^(\S+) framing ([+-]?\d+) ρ=\[([\d, ]+)\]：(.*)$")


def main():
    stats = defaultdict(lambda: {"x": set(), "n": 0, "neg_before": 0, "neg_after": 0, "zero": 0, "boxes": 0})
    for line in Path("logs/framing_batch.log").read_text().splitlines():
        m = LINE.match(line)
        if not m:
            continue
        name, x, rho, msg = m.group(1), int(m.group(2)), [int(c) for c in m.group(3).split(",")], m.group(4)
        s = stats[name]
        s["x"].add(x)
        s["boxes"] = max(s["boxes"], sum(rho))
        s["n"] += 1
        if "\\hat f = 0" in msg:
            s["zero"] += 1
            continue
        b, a = re.search(r"前 (\d+)/\d+ 負", msg), re.search(r"後 (\d+)/\d+ 負", msg)
        s["neg_before"] += int(b.group(1))
        s["neg_after"] += int(a.group(1)) if a else 10**9          # G 展開が閉じない場合は目立たせる
    rows, total = [], 0
    for name, s in stats.items():
        total += s["n"]
        rows.append(f"{name.replace('_', '\\_')} & ${min(s['x'])}..{max(s['x'])}$ & {s['boxes']} & {s['n']} & "
                    f"{s['neg_before']} & {s['neg_after']} & {s['zero']} \\\\")
    tex = ("% scripts/paper_framing_table.py が logs/framing_batch.log から生成。\n"
           f"% 判定件数の合計 {total}。\n"
           "\\begin{tabular}{lcccccc}\n\\hline\n結び目 & $x$ & 箱数 & 件数 & 前の負 & 後の負 & $\\hat f=0$ \\\\\n\\hline\n"
           + "\n".join(rows) + "\n\\hline\n\\end{tabular}\n")
    Path("docs/paper/generated/framing_table.tex").write_text(tex)
    print(tex)


if __name__ == "__main__":
    main()
