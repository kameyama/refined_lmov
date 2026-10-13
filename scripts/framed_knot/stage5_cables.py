"""Stage 5/6：ケーブル結び目の refined CS 不変量（reduced）を GN の構成で計算し、(N1)(N2) で正規化して data/framed_knot/cables/ に保存する。

使い方：python scripts/framed_knot/stage5_cables.py K1 [colors...]   （colors 省略時は 1 2 11 3 21 111）
"""
import json
import resource
import sys
import time
from fractions import Fraction
from pathlib import Path

from rcsknots.exact import gn_evaluator, reconstruct, reconstruct_points, verify
from rcsknots.homfly_num import homfly_evaluator
from rcsknots.normalize import normalize_exact

COLORS = {"1": (1,), "2": (2,), "11": (1, 1), "3": (3,), "111": (1, 1, 1), "21": (2, 1)}
# (p2, q2)-cable of T(p1, q1); q2 in the Seifert framing of the companion
KNOTS = {"K1": (2, 3, 2, 13), "K2": (2, 3, 2, 1), "K3": (2, 3, 2, 15), "K4": (2, 3, 2, 3),
         "K5a": (2, 3, 2, 5), "K5b": (2, 3, 2, 7), "K5c": (2, 3, 2, -1), "K6a": (2, 5, 2, 21), "K6b": (2, 5, 2, 1),
         # 第 3 期 U4：代数的なケーブル（CD14・CP16 で代数的と確認済み、notes_DAHA.md の辞書表）
         "K7": (2, 3, 3, 19), "K8": (3, 4, 2, 25), "K9": (2, 7, 2, 31)}
OUT = Path("data/framed_knot/cables")
WORKERS = 1        # 格子点の評価の並列数（scripts/framed_knot/u4_run.py が設定する）


def slopes(p1, q1, p2, q2):
    """GN slopes (DGR orientation): companion (p1, −q1); pattern (p2, −s2) with the DAHA Newton pair s2 = q2 − p1 q1 p2."""
    return (p2, -(q2 - p1 * q1 * p2)), (p1, -q1)


TRANSPOSE = {"111": "3"}     # 転置の彩色から作る（転置対称性、補間に使っていない点で GN と照合）


def from_transpose(name, cname, func):
    """rCS_{λ^T}(a, q, t) = rCS_λ(a, 1/t, 1/q)（正規化済み、docs/framed_knot/conventions.md §7）から λ^T の reduced 値（正規化前）を作る。

    正規化の単項式は (N2') の予測 M_λ = M_□^{|λ|} T_λ^{p₂q₂}。得た式は独立な素数の 16 点で GN の値と照合する（全係数の検証）。
    """
    from renormalize_cables import predicted
    src = json.loads((OUT / f"{name}_{TRANSPOSE[cname]}.json").read_text())
    m1 = json.loads((OUT / f"{name}_1.json").read_text())["monomial_applied"]
    p1, q1, p2, q2 = KNOTS[name]
    m = predicted(m1, COLORS[cname], p2 * q2)
    raw = {}
    for k, v in src["poly"].items():
        ea, eq, et = map(int, k.split(","))
        key = (Fraction(ea, 2) - m["a"], Fraction(-et, 2) - m["q"], Fraction(-eq, 2) - m["t"])   # (a, q, t) ↦ (a, 1/t, 1/q)、単項式を外す
        assert all(x.denominator == 1 for x in key)
        raw[tuple(int(x) for x in key)] = m["sign"] * v
    verify(gn_evaluator(func), raw, seed=2000, npts=16)
    return raw


def compute(name, cname, log=print):
    path = OUT / f"{name}_{cname}.json"
    if path.exists():
        return json.loads(path.read_text())
    p1, q1, p2, q2 = KNOTS[name]
    lam = COLORS[cname]
    inner, outer = slopes(p1, q1, p2, q2)
    t0 = time.time()
    func = lambda G: G.cable_reduced_apoly(lam, inner, outer)      # noqa: E731
    method = "GN の値の補間"
    if cname in TRANSPOSE and (OUT / f"{name}_{TRANSPOSE[cname]}.json").exists():
        raw = from_transpose(name, cname, func)
        method = f"[{TRANSPOSE[cname]}] から転置対称性で作り、独立な素数の 16 点で GN の値と照合"
    else:
        cache = Path(f"data/cache/rcs/{name}_{cname}") if WORKERS > 1 else None     # まとまりごとの保存（再開用）
        raw = reconstruct(func, log=log, workers=WORKERS, cache=cache)
    t1 = time.time()
    hraw = reconstruct_points(homfly_evaluator(lam, [(p1, q1), (p2, q2)]), log=lambda s: None)
    hom = {(ea, eq): c for (_, ea, eq), c in hraw.items()}
    norm, mono = normalize_exact(raw, hom)
    data = {"knot": name, "cable": {"companion": [p1, q1], "pattern": [p2, q2]},
            "daha_newton_pairs": {"r": [p1, p2], "s": [q1, q2 - p1 * q1 * p2]},
            "algebraic": q2 > p1 * q1 * p2, "gn_slopes": {"pattern": list(inner), "companion": list(outer)},
            "color": list(lam), "reduced": True, "method": method,
            "variables": "exponents of (A=a^(1/2), Q=q^(1/2), T=t^(1/2)) of 1703.05408",
            "framing": "0; normalization (N1)+(N2) of docs/framed_knot/conventions.md §5 (hypothesis for cables)",
            "monomial_applied": {"sign": mono[0], "a": str(mono[1]), "q": str(mono[2]), "t": str(mono[3])},
            "seconds_refined": round(t1 - t0, 1), "seconds_total": round(time.time() - t0, 1),
            "max_rss_MB": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20, 1),
            "n_terms": len(norm),
            "poly": {",".join(map(str, k)): v for k, v in sorted(norm.items())}}
    OUT.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=1))
    return data


if __name__ == "__main__":
    name = sys.argv[1]
    for cname in (sys.argv[2:] or list(COLORS)):
        d = compute(name, cname)
        print(name, cname, d["n_terms"], "terms", d["seconds_total"], "s", d["max_rss_MB"], "MB", flush=True)
