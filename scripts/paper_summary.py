"""論文の表の素材：data/bps/ の保存値（framing 0）から、結び目・ρ ごとの genus 展開の前後の項数と負の項の数を集める。

計算はしない（保存済みの \\hat f と \\hat N を読むだけ）。genus 展開の前は太字の変数（a 次数のずらしは nontorus 規約）、
全体の符号は保存時に適用済みのもの（overall_sign_applied）。
出力：data/paper_summary.json、docs/paper/generated/summary_table{,_main,_nonalg}.tex（全体、本文用、付録用）。
使い方：PYTHONPATH=src python scripts/paper_summary.py
"""
import json
import lzma
import re
from pathlib import Path

from rcsknots.bps import to_bold
from rcsknots.extract import A_SHIFT

RHOS = ["1", "2", "11", "3", "21", "111"]
RHO_TEX = {"1": "\\yng(1)", "2": "\\yng(2)", "11": "\\yng(1,1)", "3": "\\yng(3)", "21": "\\yng(2,1)", "111": "\\yng(1,1,1)"}   # ヤング図は KN と同じ youngtab の \yng
CLASSES = [  # (名前, 表示, 分類)。分類は論文の表の行のまとまり
    ("T2_3", "$T_{2,3}$", "torus"), ("T2_5", "$T_{2,5}$", "torus"), ("T3_4", "$T_{3,4}$", "torus"),
    ("T3_5", "$T_{3,5}$", "torus"), ("T2_11", "$T_{2,11}$", "torus"), ("T3_7", "$T_{3,7}$", "torus"),
    ("T4_5", "$T_{4,5}$", "torus"),
    ("K1", "K1", "alg"), ("K3", "K3", "alg"), ("K6a", "K6a", "alg"), ("K9", "K9", "alg"),
    ("K8", "K8", "alg"), ("K7", "K7", "alg"), ("K10", "K10", "alg"),
    ("4_1", "$4_1$", "kn"),
    ("K2", "K2", "non"), ("K4", "K4", "non"), ("K5a", "K5a", "non"), ("K5b", "K5b", "non"),
    ("K5c", "K5c", "non"), ("K6b", "K6b", "non"),
]


def load(name, rho):
    p = Path(f"data/bps/{name}_{rho}.json")
    if p.exists():
        return json.loads(p.read_text())
    px = Path(f"{p}.xz")
    return json.loads(lzma.open(px, "rt").read()) if px.exists() else None


def summarize(d, rho):
    hat = {tuple(map(int, k.split(","))): v for k, v in d["hat_f"].items()}
    bold = to_bold(hat, A_SHIFT["nontorus"][sum(map(int, rho))])
    pos = d["positivity"]
    neg = sum(v < 0 for v in bold.values())
    return {"before_terms": len(bold), "before_neg": min(neg, len(bold) - neg),     # 全体の符号は良い方（判定の規則）
            "after_terms": pos["n_terms"], "after_neg": min(pos["n_negative"], pos["n_terms"] - pos["n_negative"])}


def main():
    out_json = Path("data/paper_summary.json")
    summary = json.loads(out_json.read_text()) if out_json.exists() else {}
    for name, _, _ in CLASSES:
        for rho in RHOS:
            key = f"{name}_{rho}"
            if key in summary:
                continue
            d = load(name, rho)
            if d is None:
                continue
            summary[key] = summarize(d, rho)
            print(key, summary[key], flush=True)
            out_json.write_text(json.dumps(summary, indent=1, sort_keys=True))   # 1 件ごとに保存（再開できる）
    for suffix, kinds in (("", None), ("_main", {"torus", "alg", "kn"}), ("_nonalg", {"non"})):
        write_table(summary, suffix, kinds)
    print("done")


def write_table(summary, suffix, kinds):
    rows = []
    for name, label, kind in CLASSES:
        if kinds is not None and kind not in kinds:
            continue
        cells = []
        for rho in RHOS:
            s = summary.get(f"{name}_{rho}")
            if s is None:
                cells.append("--")
            elif s["before_neg"] == 0 and s["after_neg"] == 0:
                cells.append(f"\\checkmark\\,({s['before_terms']})")
            else:
                cells.append(f"{s['before_neg']}/{s['before_terms']}; {s['after_neg']}/{s['after_terms']}")
        rows.append(f"{label} & {kind} & " + " & ".join(cells) + " \\\\")
    head = " & ".join(f"${RHO_TEX[r]}$" for r in RHOS)
    tex = ("% scripts/paper_summary.py が生成。\\checkmark(n) は genus 展開の前後とも負の項 0（n は前の項数）。"
           "それ以外は 前の負/前の項数; 後の負/後の項数（前後それぞれで全体の符号は良い方）。framing 0。\n"
           "\\begin{tabular}{ll" + "c" * len(RHOS) + "}\n\\hline\n結び目 & 分類 & " + head + " \\\\\n\\hline\n"
           + "\n".join(rows) + "\n\\hline\n\\end{tabular}\n")
    Path(f"docs/paper/generated/summary_table{suffix}.tex").write_text(tex)


if __name__ == "__main__":
    main()
