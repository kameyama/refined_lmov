r"""docs/refined_kauffman/tex/gen/*.tex（rk_report_tex.py が生成した日本語の見出し付きの表）を英語の見出しに置き換えて
docs/refined_kauffman/tex/gen_en/ に書く（論文の英語版 paper/en が \input する。続き 142）。数式は変えない。
置き換えの後に日本語が残っていれば止まる。"""
import re
import sys
from pathlib import Path

SRC = Path("docs/refined_kauffman/tex/gen")
OUT = Path("docs/refined_kauffman/tex/gen_en")

RULES = [
    # 見出し（全体の形）
    (r"^\$m = (\d+)\$（\$T\((\d+),(\d+)\)\$、(\d+) 項）：\$g\^\{\\mathrm\{true\}\}_\{\\yone\}\$ の \$\\ba\$ の冪ごとの係数",
     r"$m = \1$ ($T(\2,\3)$, \4 terms): coefficients of $g^{\mathrm{true}}_{\yone}$ by powers of $\ba$"),
    (r"^\$m = (\d+)\$：\$c = 1\$ 部分の \$G\$ 展開（\$\\lambda\^\\beta\$ ごと、\$\\qr = \(q/t\)\^\{1/2\}\$）",
     r"$m = \1$: $G$-expansion of the $c = 1$ part (per $\lambda^\beta$, $\qr = (q/t)^{1/2}$)"),
    (r"^\$m = (\d+)\$：同じものを \$\\lambda\^\\beta\$ ごとに因数分解した形（",
     r"$m = \1$: the same, factorized per $\lambda^\beta$ ("),
    (r"、同じものを \$\\lambda\^\\beta\$ ごとに因数分解した形（", r": the same, factorized per $\lambda^\beta$ ("),
    (r"、境界因子をすべて外した \$N = ", r": $N$ with all boundary factors removed, $N = "),
    (r"\s*の \$G\$ 展開（\$\\lambda\^\\beta\$ ごと、(\d+) 項）", r": $G$-expansion (per $\lambda^\beta$, \1 terms)"),
    (r"^\$\\Rred_\{(\\y\w+)\}\(t\^N\)\$（\$\(t\^N\)\^j\$ の係数）：", r"$\Rred_{\1}(t^N)$ (coefficients of $(t^N)^j$):"),
    (r"^被覆の随伴 \$\(1;1\)\$ の \$\\Rred_\{\\mathrm\{adj\}\}\$（\$\\lambda = t\^\{N_A/2\}\$ の冪の係数）：",
     r"$\Rred_{\mathrm{adj}}$ of the adjoint $(1;1)$ of the cover (coefficients of powers of $\lambda = t^{N_A/2}$):"),
    (r"（整数 Laurent、(\d+) 項）：", r" (integer Laurent, \1 terms):"),
    (r"^\$\\hat g_\{\\yone\}\$ の商：", r"quotient of $\hat g_{\yone}$:"),
    (r"^和の商 (\$[^$]*\$)（(\d+) 項）：", r"quotient of the sum, \1 (\2 terms):"),
    (r"（(\d+) 項、\$\\tlam\$ の冪ごと）：", r" (\1 terms, by powers of $\tlam$):"),
    (r"（(\d+) 項）：", r" (\1 terms):"),
    (r"（\$z = t_M - t_M\^\{-1\}\$ の冪ごと、\$\\lambda\$ は Mariño の \$\\nu\$）：",
     r" (by powers of $z = t_M - t_M^{-1}$; $\lambda$ is the $\nu$ of \cite{Marino09}):"),
    (r"：\$t = q\$ の部分 \$=", r": $t = q$ part $="),
    (r"：refined だけの部分 \$=", r": refined part $="),
    (r"（refined だけの部分は 0）", r" (the refined part is 0)"),
    (r"\\text\{既約（上の形のまま）\}", r"\text{irreducible (as above)}"),
    # 表
    (r"^成分 & \$t = q\$ の部分（\$\\tlam\$ の冪ごとの符号） & refined だけの部分（商の符号） \\\\",
     r"component & $t = q$ part (sign by power of $\tlam$) & refined part (sign of the quotient) \\"),
    (r"^結び目 & セクター & \$\\lambda\^\\beta\$ ごとの符号 \\\\", r"knot & sector & sign per $\lambda^\beta$ \\"),
    (r"境界因子を外した brane--image \$N\$", r"brane--image $N$ with boundary factors removed"),
    (r"^(\$T\(\d+,\d+\)\$)、", r"\1, "),
    (r"\\textbf\{混在\}", r"\textbf{mixed}"),
    (r"（混在）", r" (mixed)"),
    (r"（\+）", r" (+)"),
    (r"（−）", r" (−)"),
    (r"そろう", r"aligned"),
    (r"、", r", "),
    (r"（", r" ("),
    (r"）", r")"),
    (r"\\begin\{tabular\}\{llp\{9\.2cm\}\}", r"\begin{tabular}{llp{7cm}}"),      # 英語の見出しは長いので列を狭める
    # 英語版の記号に合わせる（paper/en：fixed sector の母関数は \Nfix、brane--image は \Bim、refined part は X_1、随伴の冪は λ）
    (r"\$N\$ with all boundary factors removed, \$N = \(1-t\)x/", r"$\Nfix$ with all boundary factors removed, $\Nfix = (1-t)\Bim/"),
    (r"brane--image \$N\$ with boundary factors removed", r"brane--image $\Nfix$ with boundary factors removed"),
    (r"\\,x/z_q\$", r"\,\Bim/z_q$"),
    (r"quotient of the sum, \$X/", r"quotient of the sum, $X_1/"),
    (r"\[A\^\{", r"[\lambda^{"),
]


def _expand(m, rep):
    """置き換え文字列の \1, \2 … だけを群に置き換え、ほかのバックスラッシュは TeX の命令としてそのまま残す。"""
    return re.sub(r"\\(\d)", lambda g: m.group(int(g.group(1))), rep)


def convert(text):
    out = []
    for line in text.split("\n"):
        for pat, rep in RULES:
            line = re.sub(pat, lambda m, rep=rep: _expand(m, rep), line)
        out.append(line)
    return "\n".join(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    bad = []
    for f in sorted(SRC.glob("*.tex")):
        t = convert(f.read_text(encoding="utf-8"))
        (OUT / f.name).write_text(t, encoding="utf-8")
        for i, line in enumerate(t.split("\n"), 1):
            if re.search(r"[ぁ-んァ-ヶ一-龥（）：、]", line):
                bad.append(f"{f.name}:{i}: {line[:100]}")
    if bad:
        print("\n".join(bad))
        sys.exit(1)
    print(f"wrote {len(list(OUT.glob('*.tex')))} files to {OUT}")


if __name__ == "__main__":
    main()
