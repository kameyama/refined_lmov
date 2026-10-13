"""非代数的なケーブルの refined LMOV 不変量が負になる原因の切り分け。出力：data/framed_knot/negativity_cause.md

1. rCS_λ 自体（太字の変数）の負の項の数。
2. 非 refined の障害：\\hat f の (A,Q,T) 指数の偶奇がそろっているので、太字の 𝐭 次数の偶奇は (𝐚,𝐪) ごと（展開の前）、
   2J_r の偶奇は (g,β) ごと（展開の後）に決まる。したがって refined の係数が全て非負なら、
   𝐭 = −1（t = q）に落とした非 refined の係数 S は、その偶奇で決まる符号 (−1)^k を一様にもつ。
   これが破れていれば、t = q で同じ彩色 HOMFLY-PT を与えるどんな refinement（偶奇がそろうもの）でも正値性は不可能。
使い方：python scripts/framed_knot/negativity_cause.py
"""
import json
from pathlib import Path

from rcsknots.bps import to_bold
from rcsknots.extract import A_SHIFT
from rcsknots.laurent import load_table
from rcsknots.lmov import load_lmov

KNOTS = ["T2_3", "T2_5", "K1", "K3", "K6a", "K2", "K4", "K5a", "K5b", "K5c", "K6b"]
COLORS = {"1": (1,), "2": (2,), "11": (1, 1), "3": (3,), "21": (2, 1), "111": (1, 1, 1)}


def ld(path):
    return {tuple(map(int, k.split(","))): v for k, v in json.loads(Path(path).read_text())["poly"].items()}


def rcs(knot, c, ref):
    return ld(f"data/framed_knot/cables/{knot}_{c}.json") if knot.startswith("K") else ref[(knot, c)]


def minority(values):
    vals = list(values)
    neg = sum(v < 0 for v in vals)
    m = min(neg, len(vals) - neg)
    return "✓" if m == 0 else f"{m}/{len(vals)}"


def compatible(coeffs):
    """coeffs: {(column, k): N}.  Positivity is possible only if sign(Σ_k (−1)^k N)·(−1)^k is the same for all columns."""
    S, parity = {}, {}
    for (col, k), v in coeffs.items():
        S[col] = S.get(col, 0) + (-1) ** k * v
        parity.setdefault(col, set()).add(k % 2)
    assert all(len(p) == 1 for p in parity.values()), "parity not uniform"
    return len({(1 if v > 0 else -1) * (-1) ** next(iter(parity[c])) for c, v in S.items() if v}) == 1


def bps_N(knot, c):
    return {tuple(map(int, k.split(","))): v for k, v in json.loads(Path(f"data/framed_knot/bps/{knot}_{c}.json").read_text())["N"].items()}


def main():
    ref = load_table("data/framed_knot/reference/notebook_torus.txt")
    head = "| 結び目 | " + " | ".join(f"[{','.join(c)}]" for c in COLORS) + " |\n|---|" + "---|" * len(COLORS)
    lines = ["# 非代数的なケーブルの負の係数の原因（scripts/framed_knot/negativity_cause.py）", "",
             "## 1. rCS_λ 自体の負の項（太字の変数、全体の符号を選んだ後の少数派／全項数）", "", head]
    for k in KNOTS:
        cells = []
        for c in COLORS:
            P = rcs(k, c, ref)
            cells.append(minority(to_bold(P, next(iter(P))[0] % 2).values()))
        lines.append(f"| {k} | " + " | ".join(cells) + " |")
    lines += ["", "## 2. 非 refined（t = q）の障害", "",
              "○：t = q の係数の符号が正値性と両立。×：t = q の彩色 HOMFLY-PT だけで正値性が不可能。「展開の前 / 展開の後」の順。", "", head]
    for k in KNOTS:
        F = load_lmov(k)
        cells = []
        for c, rho in COLORS.items():
            b = to_bold(F[("hat", rho)], A_SHIFT["nontorus"][sum(rho)])
            before = compatible({((i, j), kk): v for (i, j, kk), v in b.items()})
            after = compatible({((g, beta), Fj): v for (g, beta, Fj), v in bps_N(k, c).items()})
            cells.append(("○" if before else "×") + " / " + ("○" if after else "×"))
        lines.append(f"| {k} | " + " | ".join(cells) + " |")
    Path("data/framed_knot/negativity_cause.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
