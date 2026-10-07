"""B4 の前段：Stošić arXiv:2402.03066 の quiver 型の式から S² 彩色 HOMFLY-PT 多項式を作り、
thick な結び目で次元が二乗になり得るかを、ホモロジーを使わずに判定する。

入力：references/2402.03066/src/SuperExponential-arxiv.tex の §3（式 (specred2)、323-327 行）の C, n, a, q, l, r。
TeX から直接読み、書き写しを data/reference/stosic_2402.03066/quivers.json に行番号付きで保存する。
Stošić の (a, q) は GS の old 規約（DGR の太字の半分）：a = 𝐚², q = 𝐪²（三葉結び目で確認）。

判定（命題 2 からの帰結）：KN の形に乗るには非負の入力で rCS_{[2]}(q=1) = rCS_{[1]}(q=1)²（打ち消し無し）が必要なので、
H^{[2]} の a 次数の集合は H^{[1]} の a 次数の対の和の集合に一致しなければならない。
彩色 HOMFLY-PT P_2 の a 次数は H^{[2]} の a 次数に含まれるので、P_2 に対の和の外の a 次数があれば、
ホモロジーと t_c 次数の選び方によらず乗らない。また P_2 の単項式の数 > (dim H^{[1]})² なら次元は二乗でない。
使い方：python scripts/b4_stosic_quiver.py [--refine]（--refine で [1] の refine を yui-kr と照合）
"""
import json
import re
import sys
from itertools import product
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from rcsknots.knot_table import homfly, load_homology, mirror, superpoly  # noqa: E402

TEX = ROOT / "references/2402.03066/src/SuperExponential-arxiv.tex"
OUT = ROOT / "data/reference/stosic_2402.03066/quivers.json"
# 結び目名 → 節の見出しの行（TeX）
SECTIONS = {"9_42": 330, "9_46": 410, "8_20": 481, "10_132": 554, "10_145": 645, "10_139": 749,
            "10_152": 838, "8_19": 2236, "10_124": 2300}
a, q = sp.symbols("a q")


def parse():
    lines = TEX.read_text().split("\n")
    data = {}
    for name, start in SECTIONS.items():
        i = start - 1
        rec = {"section_line": start}
        while "C=\\left(" not in lines[i]:
            i += 1
        j = i
        while "\\end{array}" not in lines[j]:
            j += 1
        body = " ".join(lines[i:j]).split("\\begin{array}", 1)[1].split("}", 1)[1]
        rows = [r for r in body.replace("\\hline", "").split("\\\\") if r.strip()]
        rec["C"] = [[int(x) for x in re.findall(r"-?\d+", r)] for r in rows]
        rec["C_lines"] = [i + 1, j + 1]
        for key in "naqlr":
            while not re.search(r"\\mathbf\{%s\}=" % key, lines[i]):
                i += 1
            j = i
            while "\\end{array}" not in lines[j]:
                j += 1
            body = " ".join(lines[i + 1:j]).replace("\\hline", "")
            rec[key] = [int(x) for x in re.findall(r"-?\d+", body.split("}", 1)[1] if "begin{array}" in body else body)]
            rec[key + "_lines"] = [i + 1, j + 1]
            i = j
        k = len(rec["C"])
        assert all(len(r) == k for r in rec["C"]) and all(len(rec[x]) == k for x in "naqlr"), name
        data[name] = rec
    return data


def qp(x, k, base):
    out = sp.Integer(1)
    for i in range(k):
        out *= 1 - x * base ** i
    return out


def P_r(rec, r):
    C, n, av, qv, lv, rv = (rec[x] for x in ("C", "n", "a", "q", "l", "r"))
    k = len(n)
    tot = sp.Integer(0)
    for d in product(*[range(r // ni + 1) for ni in n]):
        if sum(ni * di for ni, di in zip(n, d)) != r:
            continue
        quad = sum(C[i][j] * d[i] * d[j] for i in range(k) for j in range(k))
        diag = sum(C[i][i] * d[i] for i in range(k))
        term = qp(q, r, q) / sp.prod([qp(q, di, q) for di in d])
        term *= (-1) ** diag * a ** sum(x * y for x, y in zip(av, d)) * q ** sum(x * y for x, y in zip(qv, d))
        term *= q ** sp.Rational(quad - diag, 2)
        term *= qp(q / a, sum(x * y for x, y in zip(rv, d)), 1 / q) * qp(q ** (-r) / a, sum(x * y for x, y in zip(lv, d)), 1 / q)
        tot += sp.cancel(term)
    return sp.expand(sp.cancel(tot))


def monos(expr):
    """{(e_a, e_q): 係数}。半整数の q 次数も扱う。"""
    expr = sp.expand(expr)
    out = {}
    for t in sp.Add.make_args(expr):
        c, m = t.as_coeff_Mul()
        pw = m.as_powers_dict()
        ea, eq = pw.get(a, 0), pw.get(q, 0)
        out[(ea, eq)] = out.get((ea, eq), 0) + c
    return {k: v for k, v in out.items() if v}


def main():
    data = parse()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"source": "arXiv:2402.03066v1 TeX SuperExponential-arxiv.tex（行番号はこのファイル）",
                               "formula": "(specred2) 323-327 行", "knots": data}, indent=1))
    hom = load_homology()
    for name, rec in data.items():
        P1 = monos(P_r(rec, 1))
        P2s = P_r(rec, 2)
        P2 = monos(P2s)
        sp_ok = sp.expand(P2s.subs(q, 1) - P_r(rec, 1).subs(q, 1) ** 2) == 0
        # yui-kr と照合：Stošić (i, j) → 太字 (2i, 2j)。鏡像も試す。
        H = superpoly(hom[name])
        match = None
        for lab, HH in (("same", H), ("mirror", mirror(H))):
            if {(2 * i, 2 * j): c for (i, j), c in P1.items()} == homfly(HH):
                match, Hm = lab, HH
        dim1 = sum(H.values())
        line = f"{name}: n={rec['n']} 節点の水準 {sorted(set(rec['n']))}、dim H^[1]={dim1}、(dim)^2={dim1 ** 2}、" \
               f"P_1 と yui-kr の HOMFLY-PT：{match}、special polynomial (q=1) P_2=P_1^2：{sp_ok}、P_2 の単項式数={len(P2)}"
        if match:
            A1 = sorted({e for e, _, _ in Hm})
            pairs = sorted({x + y for x in A1 for y in A1})
            A2 = sorted({2 * i for i, _ in P2})
            out_ = [x for x in A2 if x not in pairs]
            line += f"\n    H^[1] の 𝐚 次数 {A1}、対の和 {pairs}、P_2 の 𝐚 次数 {A2}、対の和の外：{out_ or 'なし'}"
        print(line, flush=True)


if __name__ == "__main__":
    main()


def refined_P1_check():
    """[1] の refine の試み：節点 i（水準 1）→ a^{a_i} q^{q_i} t^{C_ii}（Stošić 2551-2553 行：対角成分は t 次数）×
    (1 + d)^{r_i} (1 + d')^{l_i}。d, d' は (a^{-1}q;q^{-1}), (a^{-1}q^{-r};q^{-1}) の refine の候補（DGR の太字の次数）。
    8_19 の生成元から読み取った (−2, 2, −1), (−2, −2, −3) と、入れ替えた割り当ての両方を試し、yui-kr の H^{[1]} と比べる。【計算、補足】"""
    data = json.loads(OUT.read_text())["knots"]
    hom = load_homology()
    D1, Dm1 = (-2, 2, -1), (-2, -2, -3)
    for name in ("8_19", "10_124", "10_139", "10_152"):
        rec = data[name]
        H = superpoly(hom[name])
        for lab, (dr, dl) in (("r→(−2,2,−1), l→(−2,−2,−3)", (D1, Dm1)), ("入れ替え", (Dm1, D1))):
            P = {}
            for i, ni in enumerate(rec["n"]):
                if ni != 1:
                    continue
                gens = [(2 * rec["a"][i], 2 * rec["q"][i], rec["C"][i][i])]
                for d, m in ((dr, rec["r"][i]), (dl, rec["l"][i])):
                    for _ in range(m):
                        gens = gens + [tuple(x + y for x, y in zip(g, d)) for g in gens]
                for g in gens:
                    P[g] = P.get(g, 0) + 1
            for lab2, HH in (("same", H), ("mirror", mirror(H))):
                if P == HH:
                    print(f"{name}: refine した [1] が yui-kr の H^[1]（{lab2}）と一致（{lab}）")


if __name__ == "__main__" and "--refine" in sys.argv:
    refined_P1_check()
