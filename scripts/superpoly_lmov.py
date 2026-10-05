"""構成した superpolynomial の候補（data/superpoly）から [1] の refined LMOV 不変量を作り、G 展開の前後の正整数係数性を判定する。

\\hat f_{[1]} = (a^{1/2} − a^{−1/2}) rCS_{[1]} なので、太字の変数では \\hat f_{[1]} = F · P と書ける。
因子 F（単項式 × (1 + 𝐚²𝐭)、符号込み）は、DAHA の rCS_□ と \\hat f_{[1]}（data/lmov）の比から決める。
使い方：python scripts/superpoly_lmov.py T2_3 K1 K2 ...
"""
import json
import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
from construct_superpolynomial import daha_bold  # noqa: E402
from rcsknots.bps import g_expand, to_bold  # noqa: E402
from rcsknots.extract import A_SHIFT  # noqa: E402
from rcsknots.lmov import load_lmov  # noqa: E402

a, q, t = sp.symbols("a q t")


def expr(P):
    return sum(c * a ** i * q ** j * t ** k for (i, j, k), c in P.items())


def poly_dict(e):
    e = sp.expand(e)
    out = {}
    for term in sp.Add.make_args(e):
        c, m = term.as_coeff_Mul()
        d = m.as_powers_dict()
        out[(int(d.get(a, 0)), int(d.get(q, 0)), int(d.get(t, 0)))] = int(c)
    return out


def judge(name, P):
    """[1] の \\hat f = F·P の G 展開の前後の判定を文字列で返す。"""
    D = daha_bold(name)
    hat = to_bold(load_lmov(name)[("hat", (1,))], A_SHIFT["nontorus"][1])
    F = sp.factor(sp.cancel(expr(hat) / expr(D)))
    f = poly_dict(sp.expand(F * expr(P)))
    neg = sum(v < 0 for v in f.values())
    try:
        N = g_expand(f)
        gneg = sum(v < 0 for v in N.values())
        par = {}
        for (g, b, Fj) in N:
            par.setdefault((g, b), set()).add(Fj % 2)
        g_res = f"G 展開の後：{len(N)} 項中 負 {gneg}" + ("" if all(len(s) == 1 for s in par.values()) else "（2J_r の偶奇がそろわない）")
    except ValueError as e:
        g_res = f"G 展開が閉じない（{str(e)[:60]}）"
    return f"因子 {F}  G 展開の前：{len(f)} 項中 負 {neg}  {g_res}"


def main(names):
    for name in names:
        sol = json.loads(Path(f"data/superpoly/{name}.json").read_text())["solution"]
        if sol is None:
            print(name, "候補なし", flush=True)
            continue
        P = {tuple(map(int, k.split(","))): v for k, v in sol.items()}
        print(name, judge(name, P), flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
