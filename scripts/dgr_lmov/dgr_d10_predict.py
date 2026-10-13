"""D10（docs/dgr_lmov/research_plan.md）：thin な結び目の refined の 2 箱（rCS_{[2]} と \\hat f_{[2]}, \\hat f_{[1,1]}）の予言。

入力（正解は使わない）：
  refined の [1]：yui-kr の HOMFLY-PT ホモロジーの Poincaré 多項式 P₁（太字、KnotInfo のブレイドの向き）。
  非 refined の [2]：同じブレイドの (2, s)-ケーブル（s = 1, 3）の HOMFLY-PT からケーブル公式を逆に解いた値（scripts/framed_knot/colored_from_cables.py の colored と同じ式。
    data/dgr_lmov/d10/homfly.txt、scripts/dgr_lmov/dgr_d10_cables.sage）。
構成：scripts/dgr_lmov/dgr_d7_construct.py の build（V2）＋ rCS_{[2]} ≥ 0。追加の特殊化の目標は c = P₁|_{d_{1|0}} から作る（D7 の再実行で 9 例全てで確認した規則）：
  d_{0|2} の像 = c²、d_{0|3} の像 = c·P₁|_{d_{0|3}}（d_{1|0} の像 = c² は条件 (S)）。
𝐭 の窓は P₁ の 𝐭 の範囲の 4 倍に余白 4 を足したもの。解が窓の端に接していれば窓を 4 ずつ広げて解き直す（最大 3 回）。全体の符号 σ の 4 通り × 乱数の目的 10 回で、異なる解を数える。
確認（予言に使わない独立な検査）：
  (a) 結び目自身の HOMFLY-PT（Sage）が P₁ の 𝐭 = −1 と一致（向き）
  (b) 答えが既知の結び目では、予言が既知の値と一致（または解の族に含まれる）。ツイスト結び目・ダブルツイスト結び目 7_4, 8_3：GNSSS（yui-kr の 6_1, 8_1 は鏡像）、
      T(2,5), T(2,7)：EHA、6_2, 6_3：Nawata–Oblomkov 1510.01795 §2.3（scripts/dgr_lmov/dgr_no_6263.py）
  (c) \\hat f の genus 表示の正値性（条件には使っていない）
  （t_c 次数の S² は thin な結び目でも δ が一つにならない（4_1 で確認）ので、thin 性の検査はしない）
出力：data/dgr_lmov/d10/<結び目>.json、記録 logs/dgr_lmov/dgr_d10_predict.log。
使い方：PYTHONPATH=src python scripts/dgr_lmov/dgr_d10_predict.py 5_2 6_2 ...
"""
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "framed_knot"))
from colored_from_cables import L as SL, M as SM, sage_to_kn  # noqa: E402
from dgr_d7_construct import build, solution, solve  # noqa: E402
from v4_closed import hat_from_inputs  # noqa: E402
from rcsknots.bps import g_expand  # noqa: E402
from rcsknots.dgr_lmov.dgr import specialize  # noqa: E402
from rcsknots.homfly import Q, a, adams_coeffs, theta, unknot  # noqa: E402
from rcsknots.knot_table import homfly, load_homology, superpoly  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
SIGMAS = ((1, 1), (-1, -1), (1, -1), (-1, 1))


def load_cables(name):
    out = {}
    for line in open(ROOT / "data/dgr_lmov/d10/homfly.txt"):
        nm, s, h = line.rstrip("\n").split("\t")
        if nm == name:
            if max(int(x) for x in re.findall(r"\d+", h)) >= 2 ** 62:
                raise OverflowError(f"{name} {s}: libhomfly の係数あふれの疑い")
            out[s] = sage_to_kn(sp.sympify(h, locals={"L": SL, "M": SM}))
    return out


def to_bold_unref(expr):
    """KN の (a, Q) の式 → 太字の 𝐭 = −1：a^i Q^j ↦ 𝐚^{2i} 𝐪^j。"""
    e = sp.expand(expr * a ** 40 * Q ** 80)
    P = sp.Poly(e, a, Q)
    return {(2 * (i - 40), j - 80): int(c) for (i, j), c in P.terms()}


def colored2(data):
    c = adams_coeffs((1,), 2)
    mus = [(2,), (1, 1)]
    W = sp.symbols("W2 W11")
    eqs = [sp.Eq(data[str(s)] * unknot((1,)), theta((1,), 2 * s) * sum(c[mu] * theta(mu, sp.Rational(-s, 2)) * w for mu, w in zip(mus, W)))
           for s in (1, 3)]
    sol = sp.solve(eqs, W, dict=True)[0]
    return {mu: sp.factor(sp.simplify(sol[w] / unknot(mu))) for mu, w in zip(mus, W)}


def targets(P1):
    (cq, ct), cs = next(iter(specialize(P1, 1, 0).items()))
    assert cq == 0
    c2 = {(0, 2 * ct): cs * cs}
    d03 = {(q, t + ct): cs * v for (q, t), v in specialize(P1, 0, 3).items()}
    return [(0, 2, c2), (0, 3, d03)]


def delta_values(P):
    return sorted({e[0] + e[1] / 2 - e[2] for e in P})


def predict(name, log, samples=10, known=None):
    t0 = time.time()
    P1 = {k: int(v) for k, v in superpoly(load_homology()[name]).items()}
    data = load_cables(name)
    H1 = to_bold_unref(data["knot"])
    orient = H1 == homfly(P1)
    if not orient:
        log(f"{name}: 結び目の HOMFLY-PT が P₁ の 𝐭 = −1 と一致しない（向きの確認に失敗）。中止")
        return None
    H2u = to_bold_unref(colored2(data)[(2,)])
    ts = [e[2] for e in P1]
    window = (min(2 * min(ts), 4 * min(ts)) - 4, max(2 * max(ts), 4 * max(ts)) + 5)
    for _ in range(4):                      # 解が窓の端に接していれば窓を広げて解き直す
        sols, feas = solve_all(P1, H2u, window, samples)
        if not sols:
            break
        lo_ = min(e[2] for s_ in sols for e in s_); hi_ = max(e[2] for s_ in sols for e in s_)
        if lo_ > window[0] and hi_ < window[1] - 1:
            break
        window = (window[0] - 4, window[1] + 4)
    rec = {"knot": name, "orientation_ok": orient, "window_t": window, "feasible_sigmas": feas, "n_solutions": len(sols),
           "delta_P1": delta_values(P1)}
    rec_tail(rec, P1, sols, window, known)
    rec["seconds"] = round(time.time() - t0)
    msg = (f"{name}: 向き {orient}、実行可能な σ {feas}、解 {len(sols)} 個、全解で genus 表示が一様な符号 {all(rec.get('genus_uniform_sign_all', [False]))}、"
           f"共通の係数 {rec.get('n_coeffs_common')}/{rec.get('n_coeffs_total')}"
           + (f"、既知の値を含む {rec.get('known_among_solutions')}" if known is not None else "")
           + (f"、総階数 {rec['rank_P2']}（dim P₁² = {rec['rank_P1_squared']}）、"
              f"genus 表示の一様な符号 {rec['genus_uniform_sign']}、窓 {window}、P₂ の 𝐭 の範囲 {rec['t_range_P2']}"
              + (f"、既知の値と一致 {rec['matches_known']}" if known is not None else "") if len(sols) == 1 else "")
           + f"（{rec['seconds']} 秒）")
    log(msg)
    return rec


CANCEL_HOM = False      # --hom：打ち消す微分 d_{1|0}, d_{0|2} をホモロジーの形で課す


def solve_all(P1, H2u, window, samples):
    sols, feas = [], []
    for sig in SIGMAS:
        tg = targets(P1)
        ch = [(1, 0, tg[0][2]), (0, 2, tg[0][2])] if CANCEL_HOM else ()
        lin, mons = build(P1, H2u, window, "b", sig, egp=True, extra_special=tg, cancel_hom=ch)
        if solve(lin, rcs_pos=True).x is None:
            continue
        feas.append(sig)
        rng = np.random.default_rng(0)
        for _ in range(samples):
            obj = np.zeros(len(lin.col))
            for m in mons:
                obj[lin.col[("x", m)]] = rng.normal()
            r = solve(lin, objective=obj, rcs_pos=True)
            if r.x is not None:
                s = solution(lin, r, mons)
                if s not in sols:
                    sols.append(s)
    return sols, feas


def rec_tail(rec, P1, sols, window, known):
    # 全ての解：共通の係数（部分的な予言）と、各解の genus 表示の符号
    keys = set().union(*sols) if sols else set()
    common = {k: sols[0].get(k, 0) for k in keys if all(s_.get(k, 0) == sols[0].get(k, 0) for s_ in sols)}
    rec["all_solutions"] = [{",".join(map(str, k)): v for k, v in sorted(s_.items())} for s_ in sols]
    rec["n_coeffs_common"] = sum(1 for v in common.values() if v)
    rec["n_coeffs_total"] = len(keys)
    gall = []
    for s_ in sols:
        hh = hat_from_inputs(P1, s_)
        ok = True
        for b in hh.values():
            N = g_expand(b)
            ok &= min(N.values()) >= 0 or max(N.values()) <= 0
        gall.append(ok)
    rec["genus_uniform_sign_all"] = gall
    if known is not None:
        rec["known_among_solutions"] = known in sols
    if len(sols) == 1:
        P2 = sols[0]
        hats = hat_from_inputs(P1, P2)
        gpos = {}
        for rho, b in hats.items():
            N = g_expand(b)
            gpos["".join(map(str, rho))] = min(N.values()) >= 0 or max(N.values()) <= 0
        edges = [min(e[2] for e in P2), max(e[2] for e in P2)]
        rec.update({"P2": {",".join(map(str, k)): v for k, v in sorted(P2.items())}, "rank_P2": sum(P2.values()),
                    "rank_P1_squared": sum(P1.values()) ** 2, "delta_P2": delta_values(P2), "thin_P2": len(delta_values(P2)) == 1,
                    "genus_uniform_sign": gpos, "t_range_P2": edges,
                    "hat_f": {"".join(map(str, rho)): {",".join(map(str, k)): v for k, v in sorted(b.items())} for rho, b in hats.items()}})
        if known is not None:
            rec["matches_known"] = P2 == known


KNOWN_TORUS = {"5_1": "T2_5", "7_1": "T2_7"}
KNOWN_DOUBLE = {"7_4": (2, 2), "8_3": (2, -2)}          # GNSSS の pretzel K_{s,t}（KNTZ の表で 7_4 = K_{2,2}, 8_3 = K_{−2,2}）


def torus_known(name):
    """EHA の rCS_{[2]}（data/framed_knot/torus、reduced、framing 0、KN の小文字）を太字に。yui-kr の向きと合わなければ鏡像。"""
    from framing_check import load_reduced
    from rcsknots.bps import to_bold
    R = load_reduced(KNOWN_TORUS[name], 2)
    P1t = to_bold(R[(1,)], 0 if all(k[0] % 2 == 0 for k in R[(1,)]) else 1)
    P2t = to_bold(R[(2,)], 0)
    P1y = superpoly(load_homology()[name])
    mir = lambda P: {(-x, -y, -z): c for (x, y, z), c in P.items()}  # noqa: E731
    sh = lambda A, B: (lambda ka, kb: (ka[0] - kb[0], ka[1] - kb[1], ka[2] - kb[2]))(min(A), min(B))  # noqa: E731
    for P1c, P2c in ((P1t, P2t), (mir(P1t), mir(P2t))):
        d = sh(P1y, P1c)
        if {(k[0] + d[0], k[1] + d[1], k[2] + d[2]): v for k, v in P1c.items()} == P1y:
            return {(k[0] + 2 * d[0], k[1] + 2 * d[1], k[2] + 2 * d[2]): v for k, v in P2c.items()}   # [1] の単項式のずれの 2 乗（仮の比較）
    return None


def main():
    from v4_closed import inputs
    KNOWN_TWIST = {"5_2": 2, "6_1": -2, "7_2": 3, "8_1": -3, "4_1": -1, "3_1": 1}
    logf = open(ROOT / "logs/dgr_lmov/dgr_d10_predict.log", "a")

    def log(s):
        print(s, flush=True); logf.write(s + "\n"); logf.flush()
    global CANCEL_HOM
    if "--hom" in sys.argv:
        CANCEL_HOM = True
    for name in [x for x in sys.argv[1:] if not x.startswith("--")]:
        known = None
        if name in KNOWN_TWIST:
            P1k, P2k = (({k: int(v) for k, v in P.items() if v}) for P in inputs(KNOWN_TWIST[name]))
            P1y = superpoly(load_homology()[name])
            mir = lambda P: {(-x, -y, -z): c for (x, y, z), c in P.items()}  # noqa: E731
            known = P2k if P1y == P1k else (mir(P2k) if P1y == mir(P1k) else None)   # yui-kr の 6_1, 8_1 は GNSSS の鏡像
        elif name in KNOWN_TORUS:
            known = torus_known(name)
        elif name in KNOWN_DOUBLE:
            from t3_pretzel import tc_graded
            known = {k: v for k, v in tc_graded(*KNOWN_DOUBLE[name], 2).items() if v}
            P1k = {k: v for k, v in tc_graded(*KNOWN_DOUBLE[name], 1).items() if v}
            mir = lambda P: {(-x, -y, -z): c for (x, y, z), c in P.items()}  # noqa: E731
            P1y = superpoly(load_homology()[name])
            known = known if P1y == P1k else (mir(known) if P1y == mir(P1k) else None)
        elif name in ("6_2", "6_3"):
            from dgr_no_6263 import tc_bold
            known = tc_bold(name, 2)                     # Nawata–Oblomkov 1510.01795 §2.3（予想式）
        rec = predict(name, log, known=known)
        if rec is not None:
            (ROOT / f"data/dgr_lmov/d10/{name}.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
