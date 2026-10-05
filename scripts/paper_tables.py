"""論文（docs/paper/）に載せる表を data/ から LaTeX で生成する。出力：docs/paper/generated/*.tex"""
import json
from pathlib import Path

from rcsknots.bps import to_bold

OUT = Path("docs/paper/generated")
RHOS = ["1", "2", "11", "3", "21", "111"]
RHO_TEX = {"1": "[1]", "2": "[2]", "11": "[1,1]", "3": "[3]", "21": "[2,1]", "111": "[1,1,1]"}
KNOTS = [("T2_3", "$T(2,3)$", "torus"), ("T2_5", "$T(2,5)$", "torus"), ("4_1", "$\\mathbf{4_1}$", "--"),
         ("K1", "K1", "alg"), ("K3", "K3", "alg"), ("K6a", "K6a", "alg"),
         ("K2", "K2", "non"), ("K4", "K4", "non"), ("K5a", "K5a", "non"),
         ("K5b", "K5b", "non"), ("K5c", "K5c", "non"), ("K6b", "K6b", "non")]


def load_bps(k, r):
    d = json.loads(Path(f"data/bps/{k}_{r}.json").read_text())
    return {tuple(map(int, x.split(","))): v for x, v in d["N"].items()}, d


def minority(N):
    neg = sum(1 for v in N.values() if v < 0)
    s = 1 if neg <= len(N) - neg else -1
    return s, sorted((k, s * v) for k, v in N.items() if s * v < 0)


def bps_table():
    rows = []
    for k, label, kind in KNOTS:
        cells = []
        for r in RHOS:
            N, d = load_bps(k, r)
            _, mino = minority(N)
            c = "\\checkmark" if not mino else f"{len(mino)}/{len(N)}"
            if d["positivity"]["mixed_parity_g_beta"]:
                c += "$^\\dagger$"
            cells.append(c)
        rows.append(f"{label} & " + " & ".join(cells) + " \\\\")
    return "\n".join(rows)


def superpoly_table():
    rows = []
    for k, label, kind in KNOTS:
        if k == "4_1":
            continue
        src = Path(f"data/cables/{k}_1.json") if k.startswith("K") else Path(f"data/torus/{k}_1.json")
        cells = []
        for r in RHOS:
            src = Path(f"data/cables/{k}_{r}.json") if k.startswith("K") else Path(f"data/torus/{k}_{r}.json")
            poly = {tuple(map(int, x.split(","))): v for x, v in json.loads(src.read_text())["poly"].items()}
            b = to_bold(poly, 0)
            neg = sum(1 for v in b.values() if v < 0)
            m = min(neg, len(b) - neg)
            cells.append("$+$" if m == 0 else f"{m}/{len(b)}")
        rows.append(f"{label} & " + " & ".join(cells) + " \\\\")
    return "\n".join(rows)


def tpoly(terms):
    """Σ c 𝐭^F as LaTeX."""
    out = []
    for F, c in sorted(terms.items()):
        mono = "" if F == 0 else ("\\mathbf{t}" if F == 1 else f"\\mathbf{{t}}^{{{F}}}")
        coef = ("" if abs(c) == 1 and mono else str(abs(c)))
        out.append(("-" if c < 0 else "+") + coef + mono)
    s = "".join(out)
    return s[1:] if s.startswith("+") else s


def polytable_by_g(k, r):
    """Rows g, columns β."""
    N, _ = load_bps(k, r)
    s, _ = minority(N)
    N = {key: s * v for key, v in N.items()}
    betas = sorted({b for _, b, _ in N})
    gs = sorted({g for g, _, _ in N})
    lines = []
    for g in gs:
        cells = []
        for b in betas:
            terms = {F: v for (gg, bb, F), v in N.items() if gg == g and bb == b}
            cells.append(f"${tpoly(terms)}$" if terms else "$0$")
        lines.append(f"{g} & " + " & ".join(cells) + " \\\\")
    header = "$g \\backslash \\beta$ & " + " & ".join(map(str, betas)) + " \\\\"
    return header, "\n".join(lines), len(betas)


def polytable(k, r):
    N, _ = load_bps(k, r)
    s, _ = minority(N)
    N = {key: s * v for key, v in N.items()}
    betas = sorted({b for _, b, _ in N})
    gs = sorted({g for g, _, _ in N})
    lines = []
    for b in betas:
        cells = []
        for g in gs:
            terms = {F: v for (gg, bb, F), v in N.items() if gg == g and bb == b}
            cells.append(f"${tpoly(terms)}$" if terms else "$0$")
        lines.append(f"{b} & " + " & ".join(cells) + " \\\\")
    header = "$\\beta \\backslash g$ & " + " & ".join(map(str, gs)) + " \\\\"
    return header, "\n".join(lines), len(gs)


def minority_list():
    out = []
    for k, r in [("K2", "1"), ("K5a", "1"), ("K6b", "1"), ("K5c", "3"), ("K5c", "21")]:
        N, _ = load_bps(k, r)
        _, mino = minority(N)
        cells = ", ".join(f"$({g},{b},{F})\\!:{v}$" for (g, b, F), v in mino)
        diag = sorted({F - b - g for (g, b, F), _ in mino})
        out.append(f"{k} & ${RHO_TEX[r]}$ & {cells} & ${', '.join(map(str, diag))}$ \\\\")
    return "\n".join(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    head = "\\begin{tabular}{lcccccc}\n\\toprule\n & $[1]$ & $[2]$ & $[1,1]$ & $[3]$ & $[2,1]$ & $[1,1,1]$ \\\\\n\\midrule\n"
    tail = "\\bottomrule\n\\end{tabular}\n"
    (OUT / "bps_table.tex").write_text(head + bps_table() + "\n" + tail)
    (OUT / "superpoly_table.tex").write_text(head + superpoly_table() + "\n" + tail)
    for k, fn, col in [("K1", polytable_by_g, ">{\\raggedright\\arraybackslash}p{4.0cm}"), ("K2", polytable, "l")]:
        h, body, n = fn(k, "1")
        (OUT / f"{k}_1_poly.tex").write_text(
            f"\\begin{{tabular}}{{c|{col * n}}}\n\\hline\n{h}\n\\hline\n{body}\n\\hline\n\\end{{tabular}}\n")
    (OUT / "minority.tex").write_text("\\begin{tabular}{llp{9.5cm}c}\n\\toprule\n結び目 & $\\rho$ & 項 & $2J_r-\\beta-g$ \\\\\n\\midrule\n"
                                      + minority_list() + "\n" + tail)
    print("written to", OUT)


if __name__ == "__main__":
    main()
