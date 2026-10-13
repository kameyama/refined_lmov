"""V4 pretzel：外側の変数の並びを固定し、内側の変数の並びを外側の添字ごとに選べる形のニュートン形
  F = Σ_l [Σ_k c^{(l)}_{kl} h_{n−k}(π^{(l)}_{0..k})] h_{m−l}(π'_{0..l})
で、全ての c が一様な符号（≤ 0）になる並びを探す（d_l(s) = Σ_k ... は s の指数多項式なので、l ごとに別の並びで書ける）。
使い方：python scripts/framed_knot/v4_pretzel_search2.py REGION OUTER   （OUTER = t：t を外側（並び共通）、s：s を外側）
"""
import sys
from itertools import permutations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from v4_closed import M_ROOTS, RHOS, newton_coeffs  # noqa: E402
from v4_pretzel import REGIONS, hat_direct  # noqa: E402


def search(region, outer):
    ss, st = REGIONS[region]
    R = lambda sg: {k: v if sg == 1 else tuple(-e for e in v) for k, v in M_ROOTS.items()}  # noqa: E731
    Rs, Rt = R(ss), R(st)
    grid = {(i, j): hat_direct(ss * (1 + i), st * (1 + j)) for i in range(5) for j in range(5)}
    if outer == "s":                                   # 役割を入れ替える
        grid = {(j, i): v for (i, j), v in grid.items()}
        Rs, Rt = Rt, Rs
    perms = list(permutations(M_ROOTS))
    for rho in RHOS:
        found = []
        for pt in perms:
            d = [newton_coeffs([grid[(i, j)][rho] for j in range(5)], [Rt[x] for x in pt]) for i in range(5)]   # d[i][l]
            choice = []
            for l in range(5):
                ok = [ps for ps in perms
                      if all(v < 0 for c in newton_coeffs([d[i][l] for i in range(5)], [Rs[x] for x in ps]) for v in c.values())]
                if not ok:
                    break
                choice.append(ok[0])
            if len(choice) == 5:
                found.append((pt, choice))
        print(f"{region} 外側={outer} ρ={rho}: 見つかった外側の並び {len(found)} 個、例 {found[:2]}", flush=True)


if __name__ == "__main__":
    search(sys.argv[1], sys.argv[2])
