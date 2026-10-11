"""3 箱の予言：ツイスト結び目 K_p の rCS_{[2,1]}（open_problems.md §F の対象）を、見つけた構造から構成する（docs/dgr_lmov/research_plan.md）。

入力：
  refined の [1], [2], [3]：GNSSS の予想式（scripts/t3_pretzel.py の tc_graded(1, p, r)）。[1,1], [1,1,1] は転置（太字の σ）。
  非 refined の [2,1]：2・3 本ケーブルからの値（data/x0/colored_t_eq_q.txt、別研究 [2,1] の X0）。
条件（scripts/dgr_box3.py の build21）：t = q、転置で不変、hook（a = t/q で c²·A(𝐪³𝐭, 𝐭)）、EGP（q = 1 で r₁r₁₁）、U(2)（a = t² で r₁r₁₁）、
  \\hat f_{[3]}, \\hat f_{[2,1]}, \\hat f_{[1,1,1]} の Laurent 性と spin 表示の非負性。rCS_{[2,1]} には正値性を課さない。
全体の符号の 8 通り × 乱数の目的 10 回で異なる解を数える。窓は [1] の 𝐭 の範囲の 6 倍で、解が端に接すれば広げる。
既知の値との照合：K_1（三葉結び目）は EHA の rCS_{[2,1]}、K_{−1}（4_1）は KN の候補。
出力：data/dgr_lmov/box3/K_<p>.json、記録 logs/dgr_box3_predict.log。使い方：PYTHONPATH=src python scripts/dgr_box3_predict.py 1 -1 2 ...
"""
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dgr_box3 import build21, hat3  # noqa: E402
from dgr_d10_predict import to_bold_unref  # noqa: E402
from dgr_d7_construct import solution, solve  # noqa: E402
from t3_pretzel import tc_graded  # noqa: E402
from t4_colored import bold_to_small  # noqa: E402
from rcsknots.bps import g_expand, to_bold  # noqa: E402
from rcsknots.homfly import Q, a  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
mir = lambda P: {(x, -y, z - y): c for (x, y, z), c in P.items()}  # noqa: E731   太字の転置


def x0(p, lam):
    nm = f"D{p}_1" if p > 0 else f"Dm{-p}_1"
    for line in open(ROOT / "data/x0/colored_t_eq_q.txt"):
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


def run(p, log, samples=10):
    t0 = time.time()
    bold, small = inputs3(p)
    H21u = x0(p, (2, 1))
    ts = [e[2] for e in bold[(1,)]]
    window = (min(6 * min(ts), 0) - 4, max(6 * max(ts), 0) + 6)
    for _ in range(4):
        sols, feas = [], []
        for sig in itertools.product((1, -1), repeat=3):
            lin, mons = build21(small, H21u, window, sig)
            if solve(lin, xbound=400).x is None:
                continue
            feas.append(sig)
            rng = np.random.default_rng(0)
            for _ in range(samples):
                obj = np.zeros(len(lin.col))
                for m in mons:
                    obj[lin.col[("x", m)]] = rng.normal()
                r = solve(lin, objective=obj, xbound=400)
                if r.x is not None:
                    s = solution(lin, r, mons)
                    if s not in sols:
                        sols.append(s)
        if not sols:
            break
        lo = min(e[2] for s in sols for e in s); hi = max(e[2] for s in sols for e in s)
        if lo > window[0] and hi < window[1] - 1:
            break
        window = (window[0] - 4, window[1] + 4)
    rec = {"p": p, "window_t": window, "feasible_sigmas": feas, "n_solutions": len(sols),
           "solutions": [{",".join(map(str, k)): v for k, v in sorted(s.items())} for s in sols]}
    kn = known21(p)
    if kn is not None:
        # 既知の値は正規化が違う可能性があるので、単項式倍を除いて比べる
        def same_up_to_mono(A_, B_):
            if len(A_) != len(B_):
                return False
            ka, kb = min(A_), min(B_)
            s_ = tuple(ka[i] - kb[i] for i in range(3))
            return all(A_.get(tuple(k[i] + s_[i] for i in range(3))) == v for k, v in B_.items())
        rec["known_match"] = [s == kn for s in sols]
        rec["known_match_up_to_monomial"] = [same_up_to_mono(s, kn) for s in sols]
    if len(sols) == 1:
        s = sols[0]
        sm = dict(small); sm[(2, 1)] = bold_to_small(s)
        H = hat3(sm)
        gpos = {}
        for rho, b in H.items():
            sh = 0 if all(k[0] % 2 == 0 for k in b) else 1
            N = g_expand(to_bold(b, sh))
            gpos["".join(map(str, rho))] = min(N.values()) >= 0 or max(N.values()) <= 0
        rec.update({"n_terms": len(s), "n_negative_coeffs": sum(1 for v in s.values() if v < 0), "genus_uniform_sign": gpos})
    rec["seconds"] = round(time.time() - t0)
    log(f"K_{p}: 実行可能な σ {feas}、解 {len(sols)} 個、窓 {window}"
        + (f"、既知と一致 {rec['known_match']}（単項式倍を除いて {rec['known_match_up_to_monomial']}）" if kn is not None else "")
        + (f"、項数 {rec['n_terms']}（負 {rec['n_negative_coeffs']}）、genus 表示の一様な符号 {rec['genus_uniform_sign']}" if len(sols) == 1 else "")
        + f"（{rec['seconds']} 秒）")
    return rec


def main():
    logf = open(ROOT / "logs/dgr_box3_predict.log", "a")

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
