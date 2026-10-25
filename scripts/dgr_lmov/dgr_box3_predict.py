"""3 箱の予言：ツイスト結び目 K_p の rCS_{[2,1]}（open_problems.md §F の対象）を、見つけた構造から構成する（docs/dgr_lmov/research_plan.md）。

入力：
  refined の [1], [2], [3]：GNSSS の予想式（scripts/framed_knot/t3_pretzel.py の tc_graded(1, p, r)）。[1,1], [1,1,1] は転置（太字の σ）。
  非 refined の [2,1]：2・3 本ケーブルからの値（data/framed_knot/x0/colored_t_eq_q.txt、別研究 [2,1] の X0）。
条件（scripts/dgr_lmov/dgr_box3.py の build21）：t = q、転置で不変、hook（a = t/q で c²·A(𝐪³𝐭, 𝐭)）、EGP（q = 1 で r₁r₁₁）、U(2)（a = t² で r₁r₁₁）、
  \\hat f_{[3]}, \\hat f_{[2,1]}, \\hat f_{[1,1,1]} の Laurent 性と spin 表示の非負性。rCS_{[2,1]} には正値性を課さない。
加えて予想 R（rCS_{[2,1]} の 𝐭 の両端の線 = [2] と [1] の線の min-plus・max-plus の積。𝐭 の窓もこれで決まる）。
全体の符号の 8 通りのうち実行可能なものごとに、一意性（乱数の方向 c で c·x > c·x₀ と < c·x₀ の実行可能性）を判定し、
一意でなければ【仮説】L1 最小（rCS_{[2,1]} の係数の絶対値の和が最小）の解を選ぶ。2026-10-11 に旧版（乱数の目的 10 回、窓を広げる）から変更。
既知の値との照合：K_1（三葉結び目）は EHA の rCS_{[2,1]}、K_{−1}（4_1）は KN の候補。
出力：data/dgr_lmov/box3/K_<p>.json（旧版の K_1, K_-1 は K_<p>_old.json に移した）、記録 logs/dgr_lmov/dgr_box3_predict.log。使い方：PYTHONPATH=src python scripts/dgr_lmov/dgr_box3_predict.py 1 -1 2 ...
"""
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "framed_knot"))
from dgr_box3 import build21, hat3, l1_min, other_solution, rcs_window21  # noqa: E402
from dgr_d10_predict import to_bold_unref  # noqa: E402
from dgr_d7_construct import solution, solve  # noqa: E402
from t3_pretzel import tc_graded  # noqa: E402
from t4_colored import bold_to_small  # noqa: E402
from rcsknots.bps import g_expand, to_bold  # noqa: E402
from rcsknots.homfly import Q, a  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
CONDS = ("U", "T", "H", "E", "U2", "L", "R")       # 仮説 B は 4_1 = K_{−1} で成り立たないので使わない
mir = lambda P: {(x, -y, z - y): c for (x, y, z), c in P.items()}  # noqa: E731   太字の転置


def x0(p, lam):
    nm = f"D{p}_1" if p > 0 else f"Dm{-p}_1"
    for line in open(ROOT / "data/framed_knot/x0/colored_t_eq_q.txt"):
        n, l_, e = line.rstrip("\n").split("\t")
        if n == nm and l_ == str(lam):
            return to_bold_unref(sp.sympify(e, locals={"a": a, "Q": Q}))
    raise KeyError((p, lam))


def inputs3(p):
    P1, P2, P3 = ({k: v for k, v in tc_graded(1, p, r).items() if v} for r in (1, 2, 3))
    bold = {(1,): P1, (2,): P2, (1, 1): mir(P2), (3,): P3, (1, 1, 1): mir(P3)}
    return bold, {lam: bold_to_small(P) for lam, P in bold.items()}


def known21(p):
    from framing_check import load_reduced
    if p == 1:
        R = load_reduced("T2_3", 3)
    elif p == -1:
        R = load_reduced("4_1", 3)
    else:
        return None
    return {k: v for k, v in to_bold(R[(2, 1)], 0 if all(k[0] % 2 == 0 for k in R[(2, 1)]) else 1).items() if v}


def run(p, log):
    t0 = time.time()
    bold, small = inputs3(p)
    H21u = x0(p, (2, 1))
    win = rcs_window21(bold[(1,)], bold[(2,)])                 # 予想 R から 𝐭 の窓が決まる
    window = (min(lo for lo, _ in win.values()), max(hi for _, hi in win.values()) + 1)
    rng = np.random.default_rng(0)
    rec = {"p": p, "conds": list(CONDS), "window_t": window, "by_sigma": {}}
    for sig in itertools.product((1, -1), repeat=3):
        lin, mons = build21(small, H21u, window, sig, CONDS)
        r = solve(lin, xbound=400)
        if r.x is None:
            continue
        sol = solution(lin, r, mons)
        c = np.zeros(len(lin.col))
        for m in mons:
            c[lin.col[("x", m)]] = rng.integers(1, 10 ** 6)
        others = [other_solution(lin, mons, sol, c, bd) for bd in (lambda v: (v + 1, np.inf), lambda v: (-np.inf, v - 1))]
        unique = all(o.status == 2 for o in others)
        x1 = l1_min(lin, mons)
        rec["by_sigma"][",".join(map(str, sig))] = {"unique": unique, "l1_solution": {",".join(map(str, k)): v for k, v in sorted(x1.items())}}
        log(f"  K_{p} σ={sig}: 一意 {unique}、L1 最小の解の項数 {len(x1)}（負 {sum(v < 0 for v in x1.values())}、Σ|係数| {sum(abs(v) for v in x1.values())}）"
            f"（{time.time() - t0:.0f} 秒）")
    sols = [{tuple(map(int, k.split(","))): v for k, v in d["l1_solution"].items()} for d in rec["by_sigma"].values()]
    kn = known21(p)
    if kn is not None:
        rec["known_match"] = [s_ == kn for s_ in sols]
    if len(sols) == 1:
        s_ = sols[0]
        sm = dict(small); sm[(2, 1)] = bold_to_small(s_)
        H = hat3(sm)
        gpos = {}
        for rho, b in H.items():
            sh = 0 if all(k[0] % 2 == 0 for k in b) else 1
            N = g_expand(to_bold(b, sh))
            gpos["".join(map(str, rho))] = min(N.values()) >= 0 or max(N.values()) <= 0
        rec.update({"n_terms": len(s_), "n_negative_coeffs": sum(1 for v in s_.values() if v < 0), "genus_uniform_sign": gpos})
    rec["seconds"] = round(time.time() - t0)
    log(f"K_{p}: 実行可能な σ {list(rec['by_sigma'])}、窓 {window}"
        + (f"、既知と一致 {rec['known_match']}" if kn is not None else "")
        + (f"、項数 {rec['n_terms']}（負 {rec['n_negative_coeffs']}）、genus 表示の一様な符号 {rec['genus_uniform_sign']}" if len(sols) == 1 else "")
        + f"（{rec['seconds']} 秒）")
    return rec


def main():
    logf = open(ROOT / "logs/dgr_lmov/dgr_box3_predict.log", "a")

    def log(s):
        print(s, flush=True); logf.write(s + "\n"); logf.flush()
    (ROOT / "data/dgr_lmov/box3").mkdir(parents=True, exist_ok=True)
    for p in [int(x) for x in sys.argv[1:]]:
        out = ROOT / f"data/dgr_lmov/box3/K_{p}.json"
        if out.exists():
            continue
        out.write_text(json.dumps(run(p, log), ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
