"""data/bps/*.json から、正値性・整数性・2J_r の偶奇の判定表を作る（Markdown）。

出力：data/bps_summary.md
"""
import json
from pathlib import Path

KNOTS = {"T2_3": "T(2,3)", "T2_5": "T(2,5)", "4_1": "4_1",
         "K1": "K1 (2,13)-cable 代数的", "K3": "K3 (2,15)-cable 代数的", "K6a": "K6a T(2,5) の (2,21)-cable 代数的",
         "K2": "K2 (2,1)-cable 非代数的", "K4": "K4 (2,3)-cable 非代数的", "K5a": "K5a (2,5)-cable 非代数的",
         "K5b": "K5b (2,7)-cable 非代数的", "K5c": "K5c (2,−1)-cable 非代数的", "K6b": "K6b T(2,5) の (2,1)-cable 非代数的"}
RHOS = ["1", "2", "11", "3", "21", "111"]


def cell(d):
    p = d["positivity"]
    nmin = min(p["n_negative"], p["n_terms"] - p["n_negative"])     # best over the overall sign
    s = "✓" if nmin == 0 else f"✗({p['n_negative']}/{p['n_terms']}, 最小 {nmin})"
    if p["mixed_parity_g_beta"]:
        s += f" 偶奇混在{len(p['mixed_parity_g_beta'])}"
    return s


def main():
    lines = ["| 結び目 | " + " | ".join(f"[{','.join(r)}]" for r in RHOS) + " |",
             "|---|" + "---|" * len(RHOS)]
    negatives = []
    for k, label in KNOTS.items():
        row = []
        for r in RHOS:
            f = Path(f"data/bps/{k}_{r}.json")
            if not f.exists():
                row.append("—")
                continue
            d = json.loads(f.read_text())
            if d.get("status") == "hat_f_not_polynomial":
                row.append("f̂ 非多項式")
                continue
            assert all(isinstance(v, int) for v in d["N"].values())       # integrality
            row.append(cell(d))
            if "positivity" not in d:
                continue
        lines.append(f"| {label} | " + " | ".join(row) + " |")
    lines += ["", "（上の表は G 展開の後の \\hat N_{ρ,g,β,J_r}。）✓：全体の符号を適当に選べば全項が非負。✗(n/N, 最小 m)：規約の符号で N 項中 n 項が負、全体の符号を選び直しても m 項が負。偶奇混在 m：(g,β) を固定して 2J_r の偶奇がそろわない組が m 個。",
              "全体の符号は最低次の項（g=0、最小 β、最小 2J_r）を正にする規約（docs/conventions.md §8）。", "",
              "## 負の係数（結び目, ρ, g, β, 2J_r, 値）", ""]
    # list only the minority-sign entries of small tables; the full lists are in data/bps/*.json
    lines += ["（各表の負の項は data/bps/{結び目}_{ρ}.json の positivity.negative に全件ある）"]
    from rcsknots.positivity import judge
    lines += ["", "## G 展開の前（\\hat f_ρ の太字変数での係数、KN の素朴な変数変換）", "",
              "| 結び目 | " + " | ".join(f"[{','.join(r)}]" for r in RHOS) + " |", "|---|" + "---|" * len(RHOS)]
    for k, label in KNOTS.items():
        cells = []
        for r in RHOS:
            j = judge(f"data/bps/{k}_{r}.json")
            m, n = j["before_G"]["min_negative"], j["before_G"]["terms"]
            cells.append("✓" if m == 0 else f"✗(最小 {m}/{n})")
        lines.append(f"| {label} | " + " | ".join(cells) + " |")
    lines += ["", "両方の表で ✓ の場合（トーラス結び目、4_1、K1、K3、K6a）は、同じ全体の符号の下で両方が非負。"]
    Path("data/bps_summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:20]))


if __name__ == "__main__":
    main()
