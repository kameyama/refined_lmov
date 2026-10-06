"""V4 pretzel：符号が異なる領域 np（s < 0 < t）で、直線ごとの 1 変数のニュートン形による証明。

np では 2 変数のニュートン形の c_kl が一様な符号にならない（logs/v4_pretzel_search_np*.log）。原因は t についての 4 番目の係数が |s| ≳ t の側で正の項を持つこと
（Y^{t} と Y^{−|s|} の型の打ち消し）。そこで、直線ごとに示す：
  (i)  s = s0 を固定：K_{s0,t} は t の指数多項式（根 {1, X, X², Y, Ȳ}）。基点 t_b から 5 点で c_l を作り、全て ≤ 0 となる (t_b, 並び) を探す。t < t_b は直接の判定。
  (ii) t = t0 を固定：K_{s,t0} は |s| の指数多項式（根は逆数）。同様。
値は 2 変数の閉じた形（v4_pretzel.closed_form("np")。25 点から作り、使っていない点で直接計算と一致済み）から取る。
使い方：python scripts/v4_pretzel_np.py NMAX BMAX   （|s0|, t0 = 1..NMAX、基点 ≤ BMAX）
"""
import sys
from itertools import permutations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from v4_closed import M_ROOTS, RHOS, newton_coeffs  # noqa: E402
from v4_pretzel import closed_eval, closed_form, fmt  # noqa: E402

PERMS = list(permutations(M_ROOTS))


def line_proof(F, sign, bmax):
    """F(n)（n ≥ 1）について、c_k が全て < 0 となる最小の基点 b と並び。見つからなければ None。"""
    R = {k: v if sign == 1 else tuple(-e for e in v) for k, v in M_ROOTS.items()}
    vals = {}
    for perms in (PERMS[:1], PERMS[1:]):               # まず既定の並び (1, X, X², Y, Ȳ) で全ての基点を試し、だめなら他の並び
        for b in range(1, bmax + 1):
            for n in range(b, b + 5):
                if n not in vals:
                    vals[n] = F(n)
            for pm in perms:
                cs = newton_coeffs([vals[n] for n in range(b, b + 5)], [R[x] for x in pm])
                if all(v < 0 for c in cs for v in c.values()):
                    return b, pm, vals
    return None, None, vals


def main(nmax, bmax):
    CS = closed_form("np")
    for kind in ("s0", "t0"):
        for x in range(1, nmax + 1):
            res = []
            for rho in RHOS:
                if kind == "s0":
                    F = lambda n: closed_eval(CS[rho], "np", -x, n)   # noqa: E731
                    b, pm, vals = line_proof(F, 1, bmax)
                else:
                    F = lambda n: closed_eval(CS[rho], "np", -n, x)   # noqa: E731
                    b, pm, vals = line_proof(F, -1, bmax)
                if b is None:
                    res.append(f"{rho}: 基点 ≤ {bmax} で見つからない")
                    continue
                pre = all("前 0/" in fmt(vals[n]) and "後 0/" in fmt(vals[n]) and "全体 −" in fmt(vals[n]) for n in range(1, b))
                res.append(f"{rho}: 基点 {b}、並び {pm}、基点より前（{b - 1} 点）の直接の判定 {pre}")
            name = f"K_{{{-x},t}}（t ≥ 1）" if kind == "s0" else f"K_{{s,{x}}}（s ≤ −1）"
            print(f"{name}: " + "；".join(res), flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))
