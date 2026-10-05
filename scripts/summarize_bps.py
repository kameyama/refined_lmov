"""refined LMOV 不変量の正値性・整数性の判定表（Markdown）を作る。出力：data/bps_summary.md

判定は src/rcsknots/lmov.py（小文字で求めた f を、読み取る直前に太字へ変換）による。
"""
from pathlib import Path

from rcsknots.lmov import judge

KNOTS = {"T2_3": "T(2,3)", "T2_5": "T(2,5)", "4_1": "4_1",
         "K1": "K1 (2,13)-cable 代数的", "K3": "K3 (2,15)-cable 代数的", "K6a": "K6a T(2,5) の (2,21)-cable 代数的",
         "K2": "K2 (2,1)-cable 非代数的", "K4": "K4 (2,3)-cable 非代数的", "K5a": "K5a (2,5)-cable 非代数的",
         "K5b": "K5b (2,7)-cable 非代数的", "K5c": "K5c (2,−1)-cable 非代数的", "K6b": "K6b T(2,5) の (2,1)-cable 非代数的"}
RHOS = [(1,), (2,), (1, 1), (3,), (2, 1), (1, 1, 1)]
TITLES = {"fq": "N^q（f^q_μ の係数）", "ft": "N^{\\bar t}（f^{\\bar t}_μ の係数）",
          "hat": "\\hat N（G 展開の前：\\hat f_ρ の係数）", "hatN": "\\hat N_{ρ,g,β,J_r}（G 展開の後）"}


def main():
    J = {k: judge(k) for k in KNOTS}
    lines = ["小文字の (a,q,t) で求めた f^q, f^{\\bar t}, \\hat f を、読み取る直前に太字へ変換して判定した（a 次数シフトは nontorus 規約、√q = −𝐪𝐭）。",
             "✓：全体の符号を選べば全項が非負。m/N：N 項中、全体の符号を選んでも m 項が負。†：(g,β) ごとの 2J_r の偶奇がそろわない。",
             "全ての係数は整数（確認済み）。", ""]
    for which, title in TITLES.items():
        lines += [f"## {title}", "", "| 結び目 | " + " | ".join(f"[{','.join(map(str, r))}]" for r in RHOS) + " |",
                  "|---|" + "---|" * len(RHOS)]
        for k, label in KNOTS.items():
            cells = []
            for rho in RHOS:
                r = J[k][rho][which]
                assert r["integral"]
                c = "✓" if r["min_negative"] == 0 else f"{r['min_negative']}/{r['terms']}"
                if which == "hatN" and not r["parity_aligned"]:
                    c += " †"
                cells.append(c)
            lines.append(f"| {label} | " + " | ".join(cells) + " |")
        lines.append("")
    both = [k for k in KNOTS if all(J[k][rho]["hat_and_hatN_same_sign"] for rho in RHOS)]
    lines.append("G 展開の前と後の両方が同じ全体の符号で全彩色非負：" + ", ".join(both))
    Path("data/bps_summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[-3:]))


if __name__ == "__main__":
    main()
