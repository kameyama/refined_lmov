"""B4：8_19 = T(3,4) の S² 四重次数ホモロジーの t_c 次数が、KN の入力（EHA の正解）の t_c 次数と一致するかを確かめる。

入力：
- GGS（arXiv:1304.3481）の \\tilde P^{[2]}(T_{3,4})(a, Q, t_r, t_c)：references/1304.3481/src/higher_v17.tex 2058-2072 行。
  TeX から直接読む。a^{16} の行（2069 行）の小文字 q^6, q^8 は Q の誤植とみなす（次元 121、refined 指数増大性、self-symmetry で確認）。
- 正解：data/torus/T3_4_2.json, T3_4_11.json（KN の (A, Q, T) = (a, q, t)^{1/2} の指数、正規化 (N1)(N2')、conventions.md §5）。
変換：
- 四重次数 → t_c 次数の三重次数（DGR の太字）：P(𝐚, 𝐪, 𝐭) = \\tilde P(a=𝐚, Q=𝐪, t_r=𝐪⁻¹, t_c=𝐭𝐪)
  （GGS 2076 行の (a,q,t_r,t_c) への変換と、t_r = 1 の置き換え。KRSS 1707.04017 (quad-super) 560 行は t_r 次数側）。
  すなわち生成元 a^A Q^q t_r^r t_c^c → (e_𝐚, e_𝐪, e_𝐭) = (A, q − r + c, c)。
- 太字 → KN の小文字：scripts/t4_colored.py の bold_to_small（(cov)、conventions.md §1.1）。
- [1,1]：ミラー (e_𝐚, e_𝐪, e_𝐭) ↦ (e_𝐚, −e_𝐪, e_𝐭 − e_𝐪) を [2] に施す。
照合は単項式倍（±A^α Q^β T^γ）を除いて全項で行う。ずれる場合は生成元ごとの T 次数のずれを出す。
使い方：python scripts/b4_819_quad.py
"""
import json
import re
import sys
from pathlib import Path

import sympy as sp
from sympy.parsing.sympy_parser import implicit_multiplication_application, parse_expr, standard_transformations

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from rcsknots.knot_table import load_homology, superpoly  # noqa: E402
from t4_colored import bold_to_small  # noqa: E402

a, Q, R, C = sp.symbols("a Q R C")


def ggs_s2_T34():
    L = (ROOT / "references/1304.3481/src/higher_v17.tex").read_text(encoding="latin-1").split("\n")
    s = " ".join(L[2057:2072]).split("=", 1)[1].rstrip(". ")
    s = s.replace("\\a", "a").replace("\\!", "").replace("\\\\", "").replace("t_r", "R").replace("t_c", "C")
    s = s.replace("q^", "Q^")  # 2069 行の誤植
    s = re.sub(r"\^\{(-?\d+)\}", r"**(\1)", s)
    s = re.sub(r"\^(-?\d)", r"**\1", s).replace("{", "(").replace("}", ")")
    e = parse_expr(s, local_dict=dict(a=a, Q=Q, R=R, C=C),
                   transformations=standard_transformations + (implicit_multiplication_application,))
    return sp.expand(e)


def quad_terms(e):
    out = {}
    for (ea, eq, er, ec), c in sp.Poly(sp.expand(e * (a * Q * R * C) ** 60), a, Q, R, C).terms():
        out[(ea - 60, eq - 60, er - 60, ec - 60)] = int(c)
    return out


def quad_to_tc(quad):
    out = {}
    for (A, q, r, c), m in quad.items():
        k = (A, q - r + c, c)
        out[k] = out.get(k, 0) + m
    return out


def mirror11(P):
    return {(ea, -eq, et - eq): c for (ea, eq, et), c in P.items()}


def load_eha(name):
    d = json.loads((ROOT / f"data/torus/{name}.json").read_text())
    return {tuple(int(x) for x in k.split(",")): v for k, v in d["poly"].items()}


def compare(mine, ref):
    """単項式倍を除いて比べる。最低次の項どうしでずらし量を決める。"""
    ka, kr = min(mine), min(ref)
    sh = tuple(x - y for x, y in zip(kr, ka))
    for sign in (1, -1):
        moved = {tuple(x + y for x, y in zip(k, sh)): sign * v for k, v in mine.items()}
        if moved == ref:
            return f"一致（ずらし {sh}、符号 {sign}）"
    # (A, Q) ごとの T 次数の比較
    def byAQ(P):
        D = {}
        for (A, q, t), v in P.items():
            D.setdefault((A, q), []).append((t, v))
        return {k: sorted(v) for k, v in D.items()}
    mm, rr = byAQ({tuple(x + y for x, y in zip(k, sh)): v for k, v in mine.items()}), byAQ(ref)
    diffs = [(k, mm.get(k), rr.get(k)) for k in sorted(set(mm) | set(rr)) if mm.get(k) != rr.get(k)]
    return f"不一致（最低次でずらし {sh}）。(A,Q) ごとの (T, 係数) が違う箇所 {len(diffs)} 個：" + \
        "".join(f"\n      {k}: 本式 {x} / 正解 {y}" for k, x, y in diffs[:12])


def main():
    e = ggs_s2_T34()
    P1 = (a**6*Q**-6 + a**8*Q**-4*R**3 + a**6*Q**-2*R**2 + a**8*R**5 + a**6*Q**2*R**4 + a**8*Q**4*R**7 + a**6*Q**6*R**6
          + a**10*R**8 + a**8*Q**-2*R**5 + a**8*Q**2*R**7 + a**6*R**4)  # GGS 2047-2048 行（DGR の [1]）
    print("GGS S²(T34)：次元", e.subs({a: 1, Q: 1, R: 1, C: 1}),
          "、refined 指数増大性", sp.expand(e.subs(C, 1) - P1**2) == 0,
          "、self-symmetry", sp.expand(e.subs(Q, 1 / (Q * R * C**2)) - e) == 0)
    # [1] の照合（yui-kr → 小文字 → EHA）
    print("[1]：yui-kr 8_19 と T3_4_1：", compare(bold_to_small(superpoly(load_homology()["8_19"])), load_eha("T3_4_1")))
    P2 = quad_to_tc(quad_terms(e))
    for label, P in (("S² → [2]", P2), ("S² → [1,1]（ミラー）", mirror11(P2))):
        for ref in ("T3_4_2", "T3_4_11"):
            print(f"{label} vs {ref}：", compare(bold_to_small(P), load_eha(ref)), flush=True)


if __name__ == "__main__":
    main()
