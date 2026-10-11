"""D7 の再実行（docs/dgr_lmov/research_plan.md）：rCS_{[2]} ≥ 0 を課して、非 refined のデータと構造から 2 箱の refined LMOV が一意に決まるか。

条件の組（全て scripts/dgr_d7_construct.py の build）：
  V1：(U) t = q、(S) a = t、(H) hook、(E) EGP、(L) Laurent、(P) \\hat f の spin 表示の非負性 ＋ rCS_{[2]} ≥ 0
      （genus 表示の非負性・σ 不変性・偶奇は、変数が膨らみソルバーが異常終了するため条件にせず、得られた解で検査する。2026-10-11）
  V2：V1 ＋ d_{0|2}（単項式）と d_{0|3}（P₁ の d_{0|3} の像の単項式倍）の特殊化。係数は正解から（S² の追跡と同じ扱い）
  V3：V2 ＋ sl(2) の像（d_{2|0} の特殊化、正解から）
全体の符号 σ の 4 通りを試し、実行可能なもの全てで乱数の目的関数を使って解を集める（rCS_{[2]} として異なるもの）。
記録：logs/dgr_d7_pos.log、data/dgr_lmov/d7_pos.json。使い方：PYTHONPATH=src python scripts/dgr_d7_pos.py [結び目...]
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dgr_d7_construct import build, parity_class, solution, solve, unrefined  # noqa: E402
from dgr_s2_follow import ratio  # noqa: E402
from v4_closed import hat_from_inputs, inputs  # noqa: E402
from rcsknots.bps import g_expand  # noqa: E402
from rcsknots.dgr import specialize  # noqa: E402
from rcsknots.knot_table import load_homology, superpoly  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SIGMAS = ((1, 1), (-1, -1), (1, -1), (-1, 1))
GENUS = False      # genus 表示の条件は変数が多く HiGHS が異常終了する（6_1 で SIGBUS）。切り分けで一意性への寄与が小さいので外し、解で検査する


def load(name):
    if name.startswith("K_"):
        P1, P2 = inputs(int(name[2:]))
    else:
        d = json.loads((ROOT / f"data/b4/{name}_s2_tc.json").read_text())
        P1 = superpoly(load_homology()[name])
        P2 = {tuple(map(int, k.split(","))): v for k, v in d["[2]"].items()}
    return {k: int(v) for k, v in P1.items() if v}, {k: int(v) for k, v in P2.items() if v}


def extras(P1, P2, level):
    out = []
    if level >= 2:
        out.append((0, 2, specialize(P2, 0, 2)))
        S1 = specialize(P1, 0, 3); S2 = specialize(P2, 0, 3)
        (s, c) = ratio(S2, S1)
        tgt = {(q + s[0], t + s[1]): c * v for (q, t), v in S1.items()}
        assert tgt == S2
        out.append((0, 3, tgt))
    if level >= 3:
        out.append((2, 0, specialize(P2, 2, 0)))
    return out


def run(name, log, samples=10):
    P1, P2 = load(name)
    ts = [k[2] for k in P2]
    th = hat_from_inputs(P1, P2)
    true_sig = tuple(1 if sum(1 for v in th[r].values() if v < 0) <= len(th[r]) / 2 else -1 for r in th)
    res = {}
    for level in (1, 2, 3):
        t0 = time.time()
        sols, feas = [], []
        for sig in SIGMAS:
            lin, mons = build(P1, unrefined(P2), (min(ts) - 2, max(ts) + 3), "b", sig, genus=GENUS, egp=True,
                              parity=parity_class(P1) if GENUS else None, extra_special=extras(P1, P2, level))
            first = solve(lin, rcs_pos=True)
            if first.x is None:
                continue
            feas.append(sig)
            rng = np.random.default_rng(0)
            for k in range(samples):
                obj = np.zeros(len(lin.col))
                for m in mons:
                    obj[lin.col[("x", m)]] = rng.normal()
                r = solve(lin, objective=obj, rcs_pos=True)
                if r.x is None:
                    continue
                s = solution(lin, r, mons)
                if s not in sols:
                    sols.append(s)
        ok = P2 in sols
        # 解の \\hat f の genus 表示の正値性を後から検査（genus を条件にしない場合）
        gpos = []
        for s_ in sols:
            hh = hat_from_inputs(P1, s_)
            try:
                gpos.append(all(min(g_expand(b).values()) >= 0 or max(g_expand(b).values()) <= 0 for b in hh.values()))
            except ValueError:
                gpos.append(False)
        res[f"V{level}"] = {"feasible_sigmas": feas, "n_solutions_found": len(sols), "true_found": ok, "true_sigma": true_sig, "genus_positive": gpos,
                            "seconds": round(time.time() - t0)}
        log(f"{name} V{level}: 実行可能な σ {feas}（正解の σ {true_sig}）、異なる rCS_[2] の解 {len(sols)} 個（各 σ で乱数 {samples} 回）、"
            f"正解を含む {ok}、各解の genus 表示が一様な符号 {gpos}（{time.time() - t0:.0f} 秒）")
        if len(sols) == 1 and ok:
            break
    return res


def main():
    names = sys.argv[1:] or ["K_1", "K_-1", "K_2", "K_-2", "K_3", "K_-3", "10_124", "10_139", "10_152"]
    logf = open(ROOT / "logs/dgr_d7_pos.log", "a")

    def log(s):
        print(s, flush=True); logf.write(s + "\n"); logf.flush()
    out_p = ROOT / "data/dgr_lmov/d7_pos.json"
    out = json.loads(out_p.read_text()) if out_p.exists() else {}
    for n in names:
        if n in out:
            continue
        out[n] = run(n, log)
        out_p.write_text(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
