"""refined framing のずれ δ に対する refined LMOV 不変量の正値性。

正規化 (N2') の rCS_λ に T_λ^δ（T_λ = q^{n(λ')} t^{−n(λ)}）を掛けてから \\hat f_ρ（|ρ| = n）を求め、
G 展開の前（太字の係数）と後（\\hat N）の負の項の数（全体の符号を選んだ後の少数派）を出す。
t = q では T_λ^2 = q^{κ_λ} なので、δ = 2 が framing 1 単位に当たる。δ が奇数では \\hat f が閉じない。
使い方：python scripts/framing_scan.py K5b 2 -16,-12,0,12   （出力は標準出力）
"""
import json
import sys
from pathlib import Path

from rcsknots.bps import g_expand, to_bold
from rcsknots.extract import A_SHIFT, reformulated_exact
from rcsknots.laurent import load_table
from rcsknots.partitions import conjugate, n_stat

COLORS = {"1": (1,), "2": (2,), "11": (1, 1), "3": (3,), "21": (2, 1), "111": (1, 1, 1)}


def ld(path):
    return {tuple(map(int, k.split(","))): v for k, v in json.loads(Path(path).read_text())["poly"].items()}


def minority(values):
    vals = list(values)
    neg = sum(v < 0 for v in vals)
    return min(neg, len(vals) - neg)


def reduced(knot):
    if knot.startswith("K"):
        return {lam: ld(f"data/cables/{knot}_{c}.json") for c, lam in COLORS.items()}
    ref = load_table("data/reference/notebook_torus.txt")
    return {lam: ref[(knot, c)] for c, lam in COLORS.items()}


def framed(base, d):
    return {lam: {(ea, eq + 2 * d * n_stat(conjugate(lam)), et - 2 * d * n_stat(lam)): v for (ea, eq, et), v in P.items()}
            for lam, P in base.items()}


def main(knot, n, deltas):
    base = reduced(knot)
    for d in deltas:
        F = reformulated_exact(framed(base, d), n, log=lambda s: None)
        row = []
        for (kind, rho), P in F.items():
            if kind != "hat":
                continue
            b = to_bold(P, A_SHIFT["nontorus"][n])
            row.append(f"{''.join(map(str, rho))}: 前 {minority(b.values())}/{len(b)} 後 {minority(g_expand(b).values())}")
        print(knot, n, f"δ={d}", " | ".join(row), flush=True)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), [int(x) for x in sys.argv[3].split(",")])
