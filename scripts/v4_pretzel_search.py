"""V4 pretzel：2 変数のニュートン形の c_kl が一様な符号になる (根の並び_s, 根の並び_t, 基点) を探す。
使い方：python scripts/v4_pretzel_search.py REGION N0 [N0T]   （基点 |s| = N0, |t| = N0T（省略時 N0）、その先 5×5 の 25 点を使う）
"""
import sys
from itertools import permutations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from v4_closed import M_ROOTS, RHOS, newton_coeffs  # noqa: E402
from v4_pretzel import REGIONS, hat_direct  # noqa: E402


def search(region, n0, n0t=None):
    n0t = n0 if n0t is None else n0t
    ss, st = REGIONS[region]
    R = lambda sg: {k: v if sg == 1 else tuple(-e for e in v) for k, v in M_ROOTS.items()}  # noqa: E731
    Rs, Rt = R(ss), R(st)
    grid = {(i, j): hat_direct(ss * (n0 + i), st * (n0t + j)) for i in range(5) for j in range(5)}
    perms = list(permutations(M_ROOTS))
    for rho in RHOS:
        best, good = None, []
        for ps in perms:
            rows = [newton_coeffs([grid[(i, j)][rho] for i in range(5)], [Rs[k] for k in ps]) for j in range(5)]
            for pt in perms:
                cs = [newton_coeffs([rows[j][k] for j in range(5)], [Rt[x] for x in pt]) for k in range(5)]
                vals = [v for row in cs for c in row for v in c.values()]
                bad = min(sum(v > 0 for v in vals), sum(v < 0 for v in vals))
                if best is None or bad < best[0]:
                    best = (bad, ps, pt)
                if bad == 0:
                    good.append((ps, pt))
        print(f"{region} N0=({n0},{n0t}) ρ={rho}: 一様な符号の並び {len(good)} 個 {good[:4]}、最小の逆符号の数 {best}", flush=True)


if __name__ == "__main__":
    search(sys.argv[1], *map(int, sys.argv[2:]))
