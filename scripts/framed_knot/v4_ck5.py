"""第 4 期 V4：ツイスト結び目の \\hat f(K_p) の、根の部分集合ごとの完全斉次対称式による非負表示（LP 証明書）。

s(n) = \\hat f_ρ(K_{±n})（n ≥ 1）は根 R = {1, X, X², Y, Ȳ}^{±1} の指数多項式。母関数 Σ_{n≥1} s(n) z^{n−1} = N(z)/Π_R(1 − r z) について
    −N(z) = Σ_{S⊂R} P_S(z) Π_{R∖S}(1 − r z)、P_S の係数は全て ≥ 0（deg_z P_S < |S|）
を満たす P を LP（全ての部分集合を列に入れる。HiGHS）で探し、scripts/framed_knot/v4_pretzel_np2.py の exact_vertex で有理数に解き直し、check_identity で厳密に検証する。
P_S の台は N の台から始め、不能なら根の ±1 個分ずつ 2 回まで広げる（不能の判定はこの台の範囲でのもの）。
これは −s(n) = Σ_{S,k} p_{S,k} h_{n−1−k}(S)（p ≥ 0）と同値で、全ての n ≥ 1 での正値性の証明になる。
単一の並びのニュートン形は S = {μ_0..μ_k}（入れ子の鎖）だけを使う特別な場合。
ニュートン形（並び 1, X, X², Y, Ȳ）では c_2 は S = {1, X, X²} の係数。mode "nochain" はこの部分集合を禁じた版で、
c_2 の正値性を使わずに済む表示があるかを調べる。mode "maxK" は |S| ≤ K に制限した版。
使い方：python scripts/framed_knot/v4_ck5.py cert SIGN RHO [MODE]   （SIGN = 1 / -1、RHO = 2 / 11、MODE = all / nochain / maxK）
        python scripts/framed_knot/v4_ck5.py verify
"""
import pickle
import sys
import time
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from v4_closed import M_ROOTS, hat_direct, load_inputs  # noqa: E402
from v4_pretzel_np2 import cert_eval, charpoly, check_identity, exact_vertex, numerator  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
ORDER = ["1", "X", "X2", "Y", "Yb"]


def roots(sign):
    return [tuple(sign * e for e in M_ROOTS[k]) for k in ORDER]


def seq(sign, rho, nmax, cache=None):
    cache = cache or load_inputs()
    return [hat_direct(sign * n, cache)[rho] for n in range(1, nmax + 1)]


def cert_path(sign, rho, mode):
    return ROOT / f"data/framed_knot/v4/v4_ck5_cert_{'p' if sign == 1 else 'm'}_{''.join(map(str, rho))}_{mode}.pkl"


def allowed(S, R, mode):
    names = {ORDER[R.index(r)] for r in S}
    if mode == "nochain":
        return names != {"1", "X", "X2"}
    if mode.startswith("max"):
        return len(S) <= int(mode[3:])
    return True


def lp(N, R, mode):
    """全ての部分集合を列に入れた LP（部分集合は高々 31 個なので列生成は使わない）。P_S の台は N の台に制限、不能なら広げる。"""
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import csc_matrix
    subsets = [S for k in range(1, len(R) + 1) for S in combinations(R, k) if allowed(S, R, mode)]
    M = sorted({e for x in N for e in x})
    for grow in range(3):
        cols, ri, ci, data, rows = [], [], [], [], {}
        b = {(k,) + e: -v for k, x in enumerate(N) for e, v in x.items()}
        for key in b:
            rows.setdefault(key, len(rows))
        for S in subsets:
            cp = charpoly([r for r in R if r not in S])
            for k in range(len(S)):
                for mono in M:
                    c = len(cols)
                    cols.append((S, k, mono))
                    for i, d in enumerate(cp):
                        for e, v in d.items():
                            key = (k + i,) + tuple(x + y for x, y in zip(e, mono))
                            ri.append(rows.setdefault(key, len(rows)))
                            ci.append(c)
                            data.append(v)
        nr, nc = len(rows), len(cols)
        A = csc_matrix((data, (ri, ci)), shape=(nr, nc))
        bv = np.zeros(nr)
        for key, v in b.items():
            bv[rows[key]] = v
        t0 = time.time()
        res = linprog(np.zeros(nc), A_eq=A, b_eq=bv, bounds=(0, None), method="highs")
        print(f"  台 {len(M)} 点：列 {nc}、行 {nr}、状態 {res.status}（{time.time() - t0:.1f} 秒）", flush=True)
        if res.status == 0:
            return {cols[c]: res.x[c] for c in range(nc) if res.x[c] > 1e-12}
        M = sorted(set(M) | {tuple(a + s_ * b_ for a, b_ in zip(e, r)) for e in M for r in R for s_ in (1, -1)})
    return None


def cert(sign, rho, mode="all"):
    R = roots(sign)
    N = numerator(seq(sign, rho, len(R) + 2), R)
    print(f"符号 {sign} ρ={rho} mode={mode}: N の項数 {[len(x) for x in N]}", flush=True)
    sol = lp(N, R, mode)
    if sol is None:
        print("  証明書なし（LP が不能）", flush=True)
        return
    ex = exact_vertex(N, R, sol)
    ok = ex is not None and check_identity(N, R, ex)
    used = sorted({tuple(ORDER[R.index(r)] for r in S) for S, _, _ in (ex or {})}, key=len)
    print(f"  厳密な証明書 {ok}、項数 {len(ex or {})}、使った部分集合 {used}", flush=True)
    if ok:
        pickle.dump(ex, open(cert_path(sign, rho, mode), "wb"))


def verify(nmax=12):
    cache = load_inputs()
    for path in sorted((ROOT / "data/framed_knot/v4").glob("v4_ck5_cert_*.pkl")):
        _, _, _, sg, rr, mode = path.stem.split("_")
        sign, rho = (1 if sg == "p" else -1), ((2,) if rr == "2" else (1, 1))
        R = roots(sign)
        C = pickle.load(open(path, "rb"))
        s = seq(sign, rho, nmax, cache)
        ok = check_identity(numerator(s[:len(R) + 2], R), R, C)
        direct = all(cert_eval(C, n) == {k: -v for k, v in s[n - 1].items()} for n in range(1, nmax + 1))
        print(f"{path.name}: 厳密な検証 {ok}、n = 1..{nmax} で表示 = −\\hat f {direct}、非負 {all(v >= 0 for v in C.values())}", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "cert":
        cert(int(sys.argv[2]), (2,) if sys.argv[3] == "2" else (1, 1), *(sys.argv[4:5] or ["all"]))
    else:
        verify()
