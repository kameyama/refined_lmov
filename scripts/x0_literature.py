"""別研究 [2,1] の X0：[2,1] 彩色 HOMFLY-PT（t = q）の文献値を TeX から読み、本研究の逆解きの値と照合する。

- Gu–Jockers 1407.5643（GJ.v3.arXiv.tex :1601-2000）：4_1, 5_2, 6_1, 7_2, 7_4, 8_1, 8_3 など。変数 q = e^{2πi/(k+N)}, λ = q^N（:143）、framing 0、自明な結び目で 1（:167, :1598）。
- Mironov–Morozov–Morozov 1408.3076（main.lf.tex）：ツイスト結び目の evolution の式 (anza21)（:553-555）と係数 (u's21)（:603-629）。u_2 は (anza21) の u_{20}。
変数の対応（本研究の (a, Q)）は chiral な結び目で 1 通りに決まる（conventions.md §14）。
使い方：python scripts/x0_literature.py（data/x0/colored_t_eq_q.txt の [2,1] と照合）
"""
import os
import re
import sys

import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application

from rcsknots.homfly import a, Q

lam, qg = sp.symbols("lambda q")
GJ_TEX = "references/1407.5643/src/GJ.v3.arXiv.tex"
MMM_TEX = "references/1408.3076/src/main.tex"
# KNTZ の K_{m,n}（Cyclotomic-Rectangular.tex の図 1）
GJ_KNOTS = {"4_1": (-1, 1), "5_2": (2, 1), "6_1": (-2, 1), "7_2": (3, 1), "7_4": (2, 2), "8_1": (-3, 1), "8_3": (-2, 2)}


def _frac(s):
    """\\frac{x}{y} → ((x)/(y))（入れ子に対応）。"""
    while "\\frac" in s:
        i = s.index("\\frac")
        args, j = [], i + 5
        for _ in range(2):
            assert s[j] == "{"
            depth, k = 0, j
            while True:
                depth += {"{": 1, "}": -1}.get(s[k], 0)
                if depth == 0:
                    break
                k += 1
            args.append(s[j + 1:k])
            j = k + 1
        s = s[:i] + f"(({args[0]})/({args[1]}))" + s[j:]
    return s


def tex_to_sympy(s, symbols):
    s = _frac(s)
    for x in ("\\nonumber", "\\left.", "\\right.", "\\left", "\\right", "\\\\", "&"):
        s = s.replace(x, " ")
    s = " ".join(s.split()).replace("\\lambda", "lam").replace("{", "(").replace("}", ")").replace("^", "**")
    tr = standard_transformations + (implicit_multiplication_application,)
    return parse_expr(s, local_dict=symbols, transformations=tr)


def gj_21():
    text = open(GJ_TEX, encoding="latin-1").read()
    out = {}
    for m in re.finditer(r"\\bar\{H\}&_\{\\ydiagram\{2,1\}\}\(\\mathbf\{(\d+)_\{?(\d+)\}?\}\)\s*=(.*?)\\end\{align\}", text, re.S):
        name = f"{m.group(1)}_{m.group(2)}"
        out[name] = sp.expand(tex_to_sympy(m.group(3), {"lam": lam, "q": qg}))
    return out


def mmm_21(k):
    """1408.3076 の H^{(k)}_{[21]}(A, q)（(anza21), (u's21)）。k は KNTZ の K_{k,1} の k（k = 1 が 3_1、−1 が 4_1）。"""
    from rcsknots.double_twist import A, q, br, qnum
    q2, q3, q5, q7, q10, q14 = (qnum(n) for n in (2, 3, 5, 7, 10, 14))
    u3 = -A**3 * br(A * q**3) * br(A / q**3) * br(A * q) * br(A / q) / br(A)
    u2p = q3 / q2**2 * A**3 * br(A * q**3) * br(A) / br(A * q) * (q2 * A**2 / q**3 - (q**4 + 1 - q**-2 + q**-4))
    u20 = -2 * br(q)**2 * A**3 * (q3 / q2)**2 * br(A * q**2) * br(A / q**2) / br(A)
    u2m = q3 / q2**2 * A**3 * br(A / q**3) * br(A) / br(A / q) * (q2 * A**2 * q**3 - (q**4 + 1 - q**2 + q**-4))
    u1 = -q3 * A**3 / br(A) * (A**4 - (q**6 + q**-6) * A**2 + (2 * q**6 - 4 * q**4 + 4 * q**2 - 3 + 4 * q**-2 - 4 * q**-4 + 2 * q**-6))
    u0 = A**3 / (br(A * q) * br(A) * br(A / q)) * (A**6 - q3 * q10 / (q2 * q5) * A**4 + q3 * q10**2 / (q2**2 * q5**2) * A**2
                                                    - q10 * q14 / (q2**2 * q5 * q7))
    H = u3 * A**(6 * k) + (u2p * q**(4 * k) + u20 + u2m * q**(-4 * k)) * A**(4 * k) + u1 * A**(2 * k) + u0
    return sp.expand(sp.cancel(sp.together(H)))


GJ_MAPS = {"λ=a, q=Q²": {lam: a, qg: Q**2}, "λ=1/a, q=Q^-2": {lam: 1 / a, qg: Q**-2},
           "λ=a, q=Q^-2": {lam: a, qg: Q**-2}, "λ=1/a, q=Q²": {lam: 1 / a, qg: Q**2}}


def load_ours(path="data/x0/colored_t_eq_q.txt"):
    out = {}
    for line in open(path):
        name, mu, val = line.rstrip("\n").split("\t")
        if mu == "(2, 1)":
            m, n = (int(x.replace("m", "-")) for x in name[1:].split("_"))
            out[(m, n)] = sp.sympify(val, locals={"a": a, "Q": Q})
    return out


if __name__ == "__main__":
    from rcsknots.double_twist import to_kn
    ours = load_ours()
    gj = gj_21()
    for name, (m, n) in GJ_KNOTS.items():
        res = [k for k, sub in GJ_MAPS.items() if sp.expand(gj[name].subs(sub, simultaneous=True) - ours[(m, n)]) == 0]
        print(f"GJ {name} = K_{{{m},{n}}}：一致する変数の対応 {res}")
    for k in sorted(m for (m, n) in ours if n == 1):
        v = to_kn(mmm_21(k))
        mir = sp.expand(v.subs({a: 1 / a, Q: 1 / Q}, simultaneous=True))
        print(f"MMM k={k}：A²=a, q=Q で {'一致' if sp.expand(v - ours[(k, 1)]) == 0 else '不一致'}、鏡像で {'一致' if sp.expand(mir - ours[(k, 1)]) == 0 else '不一致'}")
    sys.stdout.flush()
