"""V4 pretzel：符号が異なる領域 np（s = −u < 0 < t）を 2 つのくさびに分けたニュートン形による証明。

np の指数多項式は F(−u, t) = Σ_{(j,k)} A_jk m_j^{−u} m_k^t（19 対、logs/v4_pretzel_wedge.log の冒頭）。
  くさび A（t = u + dA + n、n ≥ 0）：項は (m_k/m_j)^u m_k^{dA+n}。n の根は m_k（5 個）、u の根は比 m_k/m_j（13 個）。
  くさび B（u = t + eB + n、n ≥ 0）：項は (m_k/m_j)^t (m_j^{−1})^{eB+n}。n の根は m_j^{−1}、t の根は同じ 13 個の比。
各くさびで、まず n についてニュートン形の係数 d_l(·)（並び (1, X, X², Y, Ȳ)、B では逆数）を作り、
各 d_l を外側の変数（u または t）の 13 根のニュートン形に書き、全係数が ≤ 0 となる根の並びを深さ優先で探す。
c_k は最初の k 個の根の集合だけで決まる（(E − r_0)⋯(E − r_{k−1}) は対称）ので、集合で枝刈りできる。
外側の変数は基点から 15 点を使い、13 個の根を当てた後の残り 2 個が 0 になることで、根の集合が足りていることを確かめる。
使い方：python scripts/v4_pretzel_wedge.py A DA    /    python scripts/v4_pretzel_wedge.py B EB
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from v4_closed import RHOS, L, M_ROOTS, newton_coeffs  # noqa: E402
from v4_pretzel import closed_eval, closed_form  # noqa: E402

NAMES = ["1", "X", "X2", "Y", "Yb"]
PAIRS = [(j, k) for j in NAMES for k in NAMES
         if not ({j, k} <= {"X2", "Y", "Yb"} and j != k)]       # 19 対（(X², Y) などの 6 対が無い）
NPTS = 15


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


RATIOS = sorted({sub(M_ROOTS[k], M_ROOTS[j]) for j, k in PAIRS})   # m_k / m_j（13 個）


def nonpos(d):
    return all(v <= 0 for v in d.values())


def dfs(seq, roots):
    """seq（L の列）を roots の並びでニュートン形にし、全係数 ≤ 0 となる並びを返す（無ければ None）。"""
    failed = set()

    def go(cur, chosen):
        if not nonpos(cur[0].to_dict()):
            return None
        if len(chosen) == len(roots):
            return list(chosen) if all(c.is_zero() for c in cur) else "残りが 0 でない"
        key = frozenset(chosen)
        if key in failed:
            return None
        for r in roots:
            if r in chosen:
                continue
            m = L.from_dict({r: 1})
            nxt = [cur[i + 1] - m * cur[i] for i in range(len(cur) - 1)]
            res = go(nxt, chosen + [r])
            if res is not None:
                return res
        failed.add(key)
        return None

    return go([x.normalized() for x in seq], [])


def wedge(kind, shift):
    CS = closed_form("np")
    inner = [M_ROOTS[x] if kind == "A" else tuple(-e for e in M_ROOTS[x]) for x in NAMES]
    for rho in RHOS:
        if kind == "A":                                         # 外側 u ≥ 1、t = u + shift + n
            G = lambda o, n: closed_eval(CS[rho], "np", -o, o + shift + n)   # noqa: E731
        else:                                                   # 外側 t ≥ 1、u = t + shift + n
            G = lambda o, n: closed_eval(CS[rho], "np", -(o + shift + n), o)  # noqa: E731
        d = {o: newton_coeffs([G(o, n) for n in range(5)], inner) for o in range(1, 1 + NPTS)}
        for l in range(5):
            res = dfs([L.from_dict(d[o][l]) for o in range(1, 1 + NPTS)], RATIOS)
            print(f"くさび {kind}（ずれ {shift}）ρ={rho} l={l}: " + ("並びなし" if res is None else
                  ("根の集合が不足（" + res + "）" if isinstance(res, str) else f"全係数 ≤ 0 の並び {res}")), flush=True)


if __name__ == "__main__":
    print(f"19 対 {PAIRS}", flush=True)
    print(f"外側の根（比 m_k/m_j、13 個）{RATIOS}", flush=True)
    wedge(sys.argv[1], int(sys.argv[2]))
