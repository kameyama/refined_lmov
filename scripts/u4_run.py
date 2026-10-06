"""第 3 期 U4：トーラス結び目・代数的なケーブルの全 6 彩色の rCS と、|ρ| ≤ 3 の refined LMOV 不変量の判定。

使い方：python scripts/u4_run.py [--workers N] T3_5 K7 ...
  --workers N：格子点の評価の並列数（既定 1）。N > 1 のとき、まとまりごとの結果を data/cache/ に保存し、再実行時は再利用する。
  rCS：data/torus/ または data/cables/（stage2_torus.compute、stage5_cables.compute + 正規化 (N2')）
  \\hat f：data/lmov/{K}_{n}.json（小文字 (A, Q, T)、a 次数シフトなし）
  判定は太字（a 次数シフトは nontorus 規約）で、G 展開の前と後。全体の符号は良い方。
"""
import json
import resource
import sys
import time
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stage2_torus  # noqa: E402
import stage5_cables  # noqa: E402
from renormalize_cables import main as renormalize, predicted  # noqa: E402
from rcsknots.bps import g_expand, to_bold  # noqa: E402
from rcsknots.extract import A_SHIFT, reformulated_exact  # noqa: E402

COLORS = {"1": (1,), "2": (2,), "11": (1, 1), "3": (3,), "111": (1, 1, 1), "21": (2, 1)}


def reduced_of(name):
    if name.startswith("K"):
        for cname in COLORS:
            d = stage5_cables.compute(name, cname, log=lambda s: None)
            print(f"  {name} [{cname}] rCS {d['n_terms']} 項 {d['seconds_total']} s", flush=True)
        renormalize([name])
        return {lam: {tuple(map(int, k.split(","))): v for k, v in
                      json.loads(Path(f"data/cables/{name}_{c}.json").read_text())["poly"].items()}
                for c, lam in COLORS.items()}
    out, metas = {}, {}
    for cname, lam in COLORS.items():
        out[lam], metas[cname] = stage2_torus.compute(name, cname)
        print(f"  {name} [{cname}] rCS {len(out[lam])} 項 {metas[cname]['seconds']} s", flush=True)
    # トーラス結び目では (N2) と (N2')（F = m n）が一致するか
    m, n = stage2_torus.KNOTS[name][0]
    m1 = metas["1"]["monomial_applied"]
    for cname, lam in COLORS.items():
        pr, got = predicted(m1, lam, m * n), metas[cname]["monomial_applied"]
        dq = pr["q"] - Fraction(got["q"])
        print(f"  {name} [{cname}] (N2') − (N2) の q 指数の差 {dq}", flush=True)
    return out


def judge(hat, n):
    bold = to_bold(hat, A_SHIFT["nontorus"][n])
    neg_b, pos_b = sum(v < 0 for v in bold.values()), sum(v > 0 for v in bold.values())
    try:
        N = g_expand(bold)
    except ValueError:
        return f"前 {min(neg_b, pos_b)}/{len(bold)} 負、G 展開が閉じない"
    neg_a, pos_a = sum(v < 0 for v in N.values()), sum(v > 0 for v in N.values())
    par = {}
    for g, b, F in N:
        par.setdefault((g, b), set()).add(F % 2)
    common = (neg_b == 0 and neg_a == 0) or (pos_b == 0 and pos_a == 0)
    return (f"前 {min(neg_b, pos_b)}/{len(bold)} 負、後 {min(neg_a, pos_a)}/{len(N)} 負"
            f"{'' if common else '（共通の符号なし）'}{'' if all(len(s) == 1 for s in par.values()) else '、2J_r の偶奇がそろわない'}")


def main(names):
    Path("data/lmov").mkdir(exist_ok=True)
    for name in names:
        t0 = time.time()
        reduced = reduced_of(name)
        for n in (1, 2, 3):
            path = Path(f"data/lmov/{name}_{n}.json")
            if path.exists():
                F = {tuple([k.split("_")[0], tuple(int(c) for c in k.split("_")[1])]): {tuple(map(int, e.split(","))): v for e, v in d.items()}
                     for k, d in json.loads(path.read_text())["f"].items()}
            else:
                t1 = time.time()
                cache = Path(f"data/cache/lmov/{name}_{n}") if WORKERS > 1 else None
                F = reformulated_exact(reduced, n, log=lambda s: print(f"    {s}", flush=True), workers=WORKERS, cache=cache)
                path.write_text(json.dumps({"knot": name, "n": n, "seconds": round(time.time() - t1, 1),
                                            "variables": "lowercase (A, Q, T); no a-shift",
                                            "f": {f"{k[0]}_{''.join(map(str, k[1]))}": {",".join(map(str, e)): c for e, c in sorted(v.items())}
                                                  for k, v in F.items()}}))
            for (kind, rho), poly in sorted(F.items()):
                if kind == "hat":
                    print(f"{name} ρ={list(rho)}: {judge(poly, n)}", flush=True)
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20
        print(f"{name}: {time.time() - t0:.0f} s, 最大 RSS {rss:.0f} MB", flush=True)


WORKERS = 1

if __name__ == "__main__":
    args = sys.argv[1:]
    if args[:1] == ["--workers"]:
        WORKERS = int(args[1])
        args = args[2:]
    stage2_torus.WORKERS = stage5_cables.WORKERS = WORKERS
    main(args)
