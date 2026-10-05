"""CD14（1408.4348）の明示例の DAHA superpolynomial を TeX ソースから読み取り、data/reference/daha_examples.txt に保存する。

各行：名前<TAB>結び目<TAB>彩色<TAB>多項式（DAHA 変数 q, t, a、sympy の式）
"""
import re
from pathlib import Path

import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application

SRC = Path("references/1408.4348/src/cable.tex")
if not SRC.exists():                              # worktree では references/ が空なので、メインチェックアウトのものを読む
    SRC = Path(__file__).resolve().parents[4] / "references/1408.4348/src/cable.tex"
OUT = Path("data/reference/daha_examples.txt")
# label in the TeX source -> (our name, knot, color)
EXAMPLES = {r"(\mathbf 1)": ("CD14_1", "K1", "1"),
            r"(\mathbf 5)": ("CD14_5", "K1", "2"),
            r"({\bf c})": ("CD14_c", "K4", "1"),
            r"({\bf d})": ("CD14_d", "K2", "1"),
            r"({\bf b})": ("CD14_b", "K5a", "1"),
            r"({\bf a})": ("CD14_a", "K5b", "1"),
            r"({\bf e})": ("CD14_e", "K5c", "1"),
            r"(\mathbf 2)": ("CD14_2", "K3", "1"),
            r"(\mathbf 3)": ("CD14_3", "K7", "1"),      # 第 3 期 U4
            r"(\mathbf 4)": ("CD14_4", "K8", "1")}
q, t, a = sp.symbols("q t a")


def latex_to_sympy(s):
    s = re.sub(r"\\[Bb]igl|\\[Bb]igr|\\left|\\right|\\!|\\,|\\;|\\\\|\\smallskip|%[^\n]*", " ", s)
    s = re.sub(r"\\frac\{([^{}]*)\}\{([^{}]*)\}", r"((\1)/(\2))", s)
    s = re.sub(r"\^\{([^}]*)\}", r"^(\1)", s)
    s = s.replace("{", "(").replace("}", ")").replace("^", "**").replace("\n", " ")
    s = s.strip().rstrip(".").rstrip(",")
    return parse_expr(s, local_dict={"q": q, "t": t, "a": a},
                      transformations=standard_transformations + (implicit_multiplication_application,))


def main():
    text = SRC.read_text(encoding="latin-1").replace("\r", "")
    rows = []
    for label, (name, knot, color) in EXAMPLES.items():
        i = text.index(label)
        seg = text[i:text.index(r"\renewcommand{\baselinestretch}{1.2}", i)]
        blocks = re.split(r"\\\)|\\\(", seg)[1::2]      # the polynomial may be split over several \( \) blocks
        poly = sp.expand(latex_to_sympy(" ".join(blocks)))
        rows.append(f"{name}\t{knot}\t{color}\t{poly}")
        print(name, knot, color, len(poly.as_ordered_terms()), "terms; at q=t=1:", poly.subs({q: 1, t: 1}))
    OUT.write_text("\n".join(rows) + "\n")


if __name__ == "__main__":
    main()
