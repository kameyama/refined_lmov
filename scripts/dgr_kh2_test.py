"""D10 の補足（docs/dgr_lmov/research_plan.md）：6_2 で、色付き Khovanov（d₂ のホモロジーを A = Q² で写したもの）が解を区別するか、一意にするか。

色付き Khovanov の実際の値は無いので、文献の正解 P₂（Nawata–Oblomkov 1510.01795 §2.3）を d₂ の残り ＋ 次数 (−2, 4, −1) の対に分けて、
残りの総数が最小になるものを代理とする（DGR の「d_N が単独でホモロジーを与える」仮定のもとでの、最小の sl(2) ホモロジー）。
  1. 各解（文献の解と別解）について、最小の残りの像（A = Q²）と、それが一意か。互いの像を作れるか。
  2. 文献の解の像を D10 の整数計画（V2 ＋ rCS_{[2]} ≥ 0）に (K) として加えたときの解の数。
使い方：PYTHONPATH=src python scripts/dgr_kh2_test.py [結び目]
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import lil_matrix

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dgr_d10_predict as D  # noqa: E402
from dgr_d7_construct import Lin, add_kh2, build, solution, solve  # noqa: E402
from dgr_no_6263 import tc_bold  # noqa: E402
from rcsknots.knot_table import load_homology, superpoly  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def kh_images(P2, samples=12, kh_fixed=None):
    """P₂ を残り ＋ 対に分ける。kh_fixed が None なら残りの総数の最小値と、最小の分け方の像（異なるもの）を返す。
    kh_fixed を与えると、その像を作れるか（実行可能か）を返す。"""
    lin = Lin()
    mons = sorted(P2)
    ts = [m[2] for m in mons]
    box = sorted({(a, q, t) for (a, q, _) in mons for t in range(min(ts) - 2, max(ts) + 3)})
    for m in box:
        lin.eq({("x", m): 1}, P2.get(m, 0))
    add_kh2(lin, box, kh_fixed)                  # kh_fixed が None なら分解だけ
    n = len(lin.col)
    M = lil_matrix((len(lin.rows), n))
    for i, row in enumerate(lin.rows):
        for j, v in row.items():
            M[i, j] = v
    rhs = np.array(lin.rhs, float)
    cons = LinearConstraint(M.tocsr(), rhs, rhs)
    cost = np.array([1.0 if nm[0] == "d2r" else 0.0 for nm in lin.col])
    res = milp(cost, constraints=[cons], integrality=np.ones(n), bounds=Bounds(0, np.inf))
    if kh_fixed is not None:
        return res.x is not None
    mn = round(res.fun)
    row = cost.reshape(1, -1)
    imgs = []
    rng = np.random.default_rng(0)
    for _ in range(samples):
        c2 = np.array([rng.normal() if nm[0] == "d2r" else 0.0 for nm in lin.col])
        r = milp(c2, constraints=[cons, LinearConstraint(row, mn, mn)], integrality=np.ones(n), bounds=Bounds(0, np.inf))
        img = {}
        for nm, v in zip(lin.col, r.x):
            if nm[0] == "d2r" and round(v):
                k = (2 * nm[1][0] + nm[1][1], nm[1][2])
                img[k] = img.get(k, 0) + int(round(v))
        if img not in imgs:
            imgs.append(img)
    return mn, imgs


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "6_2"
    out = open(ROOT / "logs/dgr_kh2_test.log", "a")

    def log(s):
        print(s, flush=True); out.write(s + "\n"); out.flush()
    rec = json.loads((ROOT / f"data/dgr_lmov/d10/{name}.json").read_text())
    sols = [{tuple(map(int, k.split(","))): v for k, v in s.items()} for s in rec["all_solutions"]]
    lit = tc_bold(name, 2)
    labels = ["文献" if s == lit else f"別解 {i}" for i, s in enumerate(sols)]
    log(f"--- {name}（解 {len(sols)} 個） ---")
    images = {}
    for lab, s in zip(labels, sols):
        mn, imgs = kh_images(s)
        images[lab] = imgs
        log(f"1. {lab}：残りの最小の総数 {mn}、最小の分け方の像 {len(imgs)} 通り（乱数 12 回）")
    lit_img = images["文献"][0]
    for lab, s in zip(labels, sols):
        if lab != "文献":
            log(f"1. 文献の像を {lab} から作れるか {kh_images(s, kh_fixed=lit_img)}、{lab} の像を文献の解から作れるか "
                f"{kh_images(lit, kh_fixed=images[lab][0])}、像が等しい {images[lab][0] == lit_img}")
    # 2. D10 の整数計画に (K) を加える
    P1 = {k: int(v) for k, v in superpoly(load_homology()[name]).items()}
    H2u = D.to_bold_unref(D.colored2(D.load_cables(name))[(2,)])
    w = tuple(rec["window_t"])
    found = []
    for sig in D.SIGMAS:
        lin, mons = build(P1, H2u, w, "b", sig, egp=True, extra_special=D.targets(P1), kh2=lit_img)
        if solve(lin, rcs_pos=True).x is None:
            continue
        rng = np.random.default_rng(0)
        for _ in range(10):
            obj = np.zeros(len(lin.col))
            for m in mons:
                obj[lin.col[("x", m)]] = rng.normal()
            r = solve(lin, objective=obj, rcs_pos=True)
            if r.x is not None:
                s = solution(lin, r, mons)
                if s not in found:
                    found.append(s)
    log(f"2. 文献の像を (K) として加えた解の数 {len(found)}、文献の解に一致 {[s == lit for s in found]}")


if __name__ == "__main__":
    main()
