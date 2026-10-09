"""大きい結び目の framed refined LMOV の判定：framing x を掛けた rCS から \\hat f_ρ（|ρ| = n）だけを補間し、G 展開の前後の正整数性を判定する。

scripts/framing_check.py と同じ framing の因子・判定（framing_check.framed, framing_check.judge）。違いは次の 2 点。
- 補間するのは \\hat f_ρ だけ（rcsknots.extract.reformulated_exact は f^q, f^t も補間するので、列の数が約 3 倍）。全ての ρ を 1 回の補間でまとめて求める。
- 格子点の評価を並列にし（--workers）、評価のまとまりと補間の係数を data/cache/framing/{K}_x{x}_n{n}/ に保存して再開できるようにする（判定を保存したら消す）。
f^t = M(q^{-1}) \\hat f の照合は全ての評価点で行う（reformulated_exact と同じ）。\\hat f が恒等的に 0 になる ρ は扱わない（自明な結び目は framing_check.py で）。
判定の結果は (結び目, x, n) ごとに data/framed/{K}_x{x}_n{n}.json に保存し、保存済みの単位は飛ばす。
使い方：python scripts/framing_big.py [--workers N] 結び目 n x1,x2,...
"""
import json
import resource
import shutil
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from framing_check import framed, judge, load_reduced  # noqa: E402
from rcsknots.bpsnum import Reformulated  # noqa: E402
from rcsknots.exact import reconstruct_points  # noqa: E402
from rcsknots.extract import _arr  # noqa: E402
from rcsknots.field import V, make_conv  # noqa: E402
from rcsknots.partitions import partitions  # noqa: E402


def hat_f_all(reduced, n, log, workers, cache):
    """{ρ: \\hat f_ρ}（KN の (A, Q, T)、a 次数シフトなし）を 1 回の補間で。"""
    rhos = partitions(n)

    def values(xs, ys, p):
        k = len(xs)
        hat, resid = Reformulated(reduced, V(xs, p), V(ys, p), make_conv(p, k)).hat_f(n)
        for r in resid.values():
            if any(np.any(_arr(c, k, p) != 0) for c in r.terms.values()):
                raise ArithmeticError("f^t ≠ M(q^{-1}) \\hat f at some point")
        return hat, k

    rng = np.random.default_rng(7)
    hat, _ = values(rng.integers(2, 10**6, 4), rng.integers(2, 10**6, 4), 2**31 - 1)
    lo = min(min(hat[r].terms) for r in rhos) - 2
    width = max(max(hat[r].terms) for r in rhos) + 2 - lo + 1

    def ev(xs, ys, p):
        hat, k = values(xs, ys, p)
        arr = np.zeros((k, width * len(rhos)), dtype=np.int64)
        for i, r in enumerate(rhos):
            for e, c in hat[r].terms.items():
                arr[:, i * width + e - lo] = (arr[:, i * width + e - lo] + _arr(c, k, p)) % p
        return arr

    raw = reconstruct_points(ev, seed=11 + n, log=log, workers=workers, cache=cache)
    out = {r: {} for r in rhos}
    for (d, eq, et), c in raw.items():
        i, e = divmod(d, width)
        out[rhos[i]][(e + lo, eq, et)] = c
    return out


def run_unit(name, n, x, workers=1):
    path = Path(f"data/framed/{name}_x{x}_n{n}.json")
    if path.exists():
        print(f"{name} framing {x:+d} n={n}：保存済み（{path}）", flush=True)
        return
    t0 = time.time()
    R = framed(load_reduced(name, n), x)
    cache = Path(f"data/cache/framing/{name}_x{x}_n{n}") if workers > 1 else None
    hats = hat_f_all(R, n, lambda s: print(f"    {s}", flush=True), workers, cache)
    res = {}
    for rho, hat in hats.items():
        msg, N = judge(hat, n)
        res["".join(map(str, rho))] = msg
        print(f"{name} framing {x:+d} ρ={list(rho)}：{msg}", flush=True)
    sec = time.time() - t0
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 2**20      # macOS はバイト単位
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"knot": name, "framing": x, "n": n, "judge": res,
                                "seconds": round(sec), "max_rss_mb": rss}, ensure_ascii=False, indent=1))
    if cache is not None:
        shutil.rmtree(cache, ignore_errors=True)      # 判定を保存した単位の途中結果は再開に要らない
    print(f"--- {name} framing {x:+d} n={n}：{sec:.0f} 秒、最大 RSS {rss} MB", flush=True)


if __name__ == "__main__":
    args = sys.argv[1:]
    workers = 1
    if args[0] == "--workers":
        workers, args = int(args[1]), args[2:]
    name, n, xs = args[0], int(args[1]), [int(v) for v in args[2].split(",")]
    for x in xs:
        run_unit(name, n, x, workers)
