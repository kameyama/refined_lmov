"""4_1 の [3,1]（4 箱）の線形解析：scripts/f41_box4.py の一次の系（t = q、転置、両手型、EGP、U(2)、\\hat f_ρ の Laurent 性）に正値性を課さず、
零空間の次元を法 p で数える。さらに hook 性質（a = t/q で rCS_{[3,1]} = c³·A(𝐪⁴𝐭², 𝐭)、docs/dgr_lmov/report.md §5 の [3,1] への適用）を加え、
一意になれば厳密解を求めて検算する。
使い方：PYTHONPATH=src:scripts python scripts/lmov31_linear.py。出力：data/lmovstruct/lmov31_linear.json、logs/lmov31_linear.log。
"""
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import f41_box4 as F  # noqa: E402
from lmov21_linear import rank_nullity, solve_unique  # noqa: E402
from rcsknots.bps import to_bold  # noqa: E402
from rcsknots.dgr import specialize  # noqa: E402


def rows_from(Amat, b):
    A = Amat.tocsr()
    rows = []
    for i in range(A.shape[0]):
        s, e = A.indptr[i], A.indptr[i + 1]
        row = {int(c): int(round(v)) for c, v in zip(A.indices[s:e], A.data[s:e]) if round(v)}
        rows.append((row, int(round(b[i]))))
    return rows


def hook_rows(xs, r1_small):
    """a = t/q（A → T/Q）での rCS_{[3,1]} = c³·A(𝐪⁴𝐭², 𝐭)。太字で作って小文字 (Q, T) に直す：𝐪 = T、𝐭 = −Q/T。"""
    r1 = to_bold(r1_small, 1 if any(k[0] % 2 for k in r1_small) else 0)
    (cq, ct), cs = next(iter(specialize(r1, 1, 0).items()))
    assert cq == 0
    Ab = specialize(r1, 1, 1)
    target_bold = defaultdict(int)
    for (eq, et), v in Ab.items():
        target_bold[(4 * eq, et + 2 * eq + 3 * ct)] += cs ** 3 * v
    target = defaultdict(int)
    for (eq, et), v in target_bold.items():                # 𝐪^{eq}𝐭^{et} → (−1)^{et} Q^{et} T^{eq − et}
        target[(et, eq - et)] += v * (-1) ** (et % 2)
    img = defaultdict(dict)
    for n, (ea, eq, et) in enumerate(xs):
        img[(eq - ea, et + ea)][n] = img[(eq - ea, et + ea)].get(n, 0) + 1
    rows = []
    for key in set(img) | set(target):
        rows.append((dict(img.get(key, {})), target.get(key, 0)))
    return rows


def main():
    Path("data/lmovstruct").mkdir(parents=True, exist_ok=True)
    lg = open("logs/lmov31_linear.log", "a")

    def log(s):
        print(s, flush=True)
        lg.write(s + "\n")
        lg.flush()

    log(f"=== lmov31_linear {time.strftime('%Y-%m-%d %H:%M:%S')} ===")
    t0 = time.time()
    red, h31 = F.inputs(log)
    al, Dp, dstr = F.alpha()
    dk = F.dK(red, Dp, log)
    res = {}
    for margin in (0, 2):
        xs, Amat, b, lb, ub = F.build(red, h31, al, Dp, dk, sigma=1, positive=False, margin=margin, log=log, with_G=False)
        rows = rows_from(Amat, b)
        ncols = Amat.shape[1]
        rank, feas = rank_nullity(rows, ncols)
        log(f"  margin {margin}：未知数 x {len(xs)}、列 {ncols}、式 {len(rows)}、階数 {rank}、零空間 {ncols - rank}、実行可能 {feas}（{time.time() - t0:.0f}s）")
        res[f"margin{margin}"] = {"x": len(xs), "cols": ncols, "rows": len(rows), "rank": rank, "nullity": ncols - rank, "feasible": feas}
        hr = hook_rows(xs, red[(1,)])
        rank2, feas2 = rank_nullity(rows + hr, ncols)
        log(f"  margin {margin} ＋ hook：式 {len(rows) + len(hr)}、階数 {rank2}、零空間 {ncols - rank2}、実行可能 {feas2}（{time.time() - t0:.0f}s）")
        res[f"margin{margin}_hook"] = {"rows": len(rows) + len(hr), "rank": rank2, "nullity": ncols - rank2, "feasible": feas2}
        if ncols - rank2 == 0 and feas2:
            x = solve_unique(rows + hr, ncols)
            if x is not None:
                sol = {",".join(map(str, xs[c])): v for c, v in x.items() if c < len(xs) and v}
                res[f"margin{margin}_solution"] = sol
                log(f"  一意解 rCS_[3,1]：{len(sol)} 項、負 {sum(1 for v in sol.values() if v < 0)}")
                Path("data/lmovstruct/lmov31_41_solution.json").write_text(json.dumps(sol))
            else:
                log("  一意解の整数での検算に失敗")
        Path("data/lmovstruct/lmov31_linear.json").write_text(json.dumps(res, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
