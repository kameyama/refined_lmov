"""第 4 期 V0：2 段の反復ケーブル K10 = T(2,3) の (2,13)-ケーブルの (2,53)-ケーブル（代数的、CD14 の DAHA [r, s] = [(3,2,2), (2,1,1)]）。

rCS を楕円 Hall 代数の反復で計算し、(N1)（t = q で framing 0 の HOMFLY-PT）と (N2')（M_λ = M_□^{|λ|}·T_λ^{p₃q₃}）で正規化して data/cables/K10_{彩色}.json に保存する。
使い方：python scripts/v0_iterated.py [--workers N] 彩色 ...（例：1 2）
"""
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from renormalize_cables import predicted  # noqa: E402
from rcsknots.exact import reconstruct, reconstruct_points  # noqa: E402
from rcsknots.homfly_num import homfly_evaluator  # noqa: E402
from rcsknots.normalize import normalize_exact  # noqa: E402

COLORS = {"1": (1,), "2": (2,), "11": (1, 1), "3": (3,), "111": (1, 1, 1), "21": (2, 1)}
STEPS = [(2, 3), (2, 13), (2, 53)]          # (p, q)：companion T(2,3) から外側へ、q は直前の結び目の Seifert framing


def gn_slopes(steps):
    """[外側のパターン, ..., companion] の GN の傾き。s_i = q_i − q_{i−1} p_{i−1} p_i（CD14 の a_i = a_{i−1} r_{i−1} r_i + s_i）。"""
    (p1, q1), out = steps[0], [(steps[0][0], -steps[0][1])]
    for (pp, qp), (p, q) in zip(steps, steps[1:]):
        out.append((p, -(q - qp * pp * p)))
    return out[::-1]


def main(args):
    workers = 1
    if args[:1] == ["--workers"]:
        workers, args = int(args[1]), args[2:]
    slopes = gn_slopes(STEPS)
    for cname in args:
        lam = COLORS[cname]
        path = Path(f"data/cables/K10_{cname}.json")
        if path.exists():
            continue
        t0 = time.time()
        cache = Path(f"data/cache/rcs/K10_{cname}") if workers > 1 else None
        raw = reconstruct(lambda G: G.iterated_reduced_apoly(lam, slopes), log=lambda s: print(f"  {s}", flush=True),
                          workers=workers, cache=cache)
        hraw = reconstruct_points(homfly_evaluator(lam, STEPS), log=lambda s: None)
        hom = {(ea, eq): c for (_, ea, eq), c in hraw.items()}
        norm, mono = normalize_exact(raw, hom)
        data = {"knot": "K10", "steps": STEPS, "gn_slopes": slopes, "color": list(lam), "reduced": True,
                "variables": "exponents of (A=a^(1/2), Q=q^(1/2), T=t^(1/2)) of 1703.05408",
                "framing": "0; normalization (N1) + (N2) (then checked against (N2'))",
                "monomial_applied": {"sign": mono[0], "a": str(mono[1]), "q": str(mono[2]), "t": str(mono[3])},
                "seconds": round(time.time() - t0, 1), "n_terms": len(norm),
                "poly": {",".join(map(str, k)): v for k, v in sorted(norm.items())}}
        path.write_text(json.dumps(data, indent=1))
        print(f"K10 [{cname}] {len(norm)} 項 {data['seconds']} s", flush=True)
    # (N2') との比較（F = p₃q₃）
    m1p = Path("data/cables/K10_1.json")
    if m1p.exists():
        m1 = json.loads(m1p.read_text())["monomial_applied"]
        p3, q3 = STEPS[-1]
        for cname, lam in COLORS.items():
            f = Path(f"data/cables/K10_{cname}.json")
            if f.exists():
                got = json.loads(f.read_text())["monomial_applied"]
                print(f"K10 [{cname}] (N2') − (N2) の q 指数の差 {predicted(m1, lam, p3 * q3)['q'] - Fraction(got['q'])}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
