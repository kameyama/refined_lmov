"""data/framed_knot/highbox/ の判定をまとめた表 data/framed_knot/highbox_summary.md を作る（保存済みの判定を読むだけ）。"""
import json
import lzma
import re
from pathlib import Path


def main():
    rows = []
    for f in Path("data/framed_knot/highbox").glob("T*_N*.json.xz"):
        with lzma.open(f, "rt") as fh:
            d = json.load(fh)
        J = [v["judge"] for v in d["rho"].values()]
        m, n = d["torus"]
        rows.append((m, n, d["N"], len(J), sum(j["before_terms"] for j in J), sum(j["before_neg"] for j in J),
                     sum(j["after_terms"] or 0 for j in J), sum(j["after_neg"] if j["after_neg"] is not None else 10**9 for j in J),
                     sorted({j["sign"] for j in J}), d["seconds"]))
    rows.sort()
    out = ["# トーラス結び目の 4 箱以上の refined LMOV 不変量（scripts/framed_knot/highbox_torus.py）", "",
           "全て【計算】。framing 0、a 次数のずらしは nontorus 規約の延長（偶数 N は 0）。負の数は全体の符号を良い方に選んだ後の数。",
           "", "| 結び目 | 箱数 N | ρ の数 | G 展開の前：項数 / 負 | 後：項数 / 負 | 全体の符号 | 計算時間（秒） |", "|---|---|---|---|---|---|---|"]
    for m, n, N, k, bt, bn, at, an, sg, sec in rows:
        out.append(f"| T({m},{n}) | {N} | {k} | {bt} / {bn} | {at} / {an} | {','.join(f'{s:+d}' for s in sg)} | {sec:.0f} |")
    Path("data/framed_knot/highbox_summary.md").write_text("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
