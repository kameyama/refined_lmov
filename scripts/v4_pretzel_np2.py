"""V4 pretzel の異符号の領域 np（s = −u < 0 < t）：くさび ＋ 外側の LP 証明書による 2 箱の正値性の証明。

方針（2026-10-07）
  くさび A（t = u + d、d ≥ 0）：G(u, d) = F(−u, u + d) を d について根の並び (1, X, X², Y, Ȳ) のニュートン形にする：
      G(u, d) = Σ_l d_l(u) h_{d−l}(m_{π0..πl})、d_l(u) = [Π_{i<l}(E_d − m_{πi}) G(u, ·)](d = 0)。
  d_l(u) は u の指数多項式で、根は比 m_k/m_j（(j, k) は 19 対のうち k ∈ {π_l..π_4}）の集合 R_l（l = 0..4 で 13, 10, 7, 5, 3 個）。
  全ての d_l(u) は共通因子 K = (1 + 𝐚²𝐭)(𝐚²𝐭 + 𝐪²)(1 + 𝐚²𝐪²𝐭³)（l ≥ 3 はさらに (1 + 𝐚²𝐭³)）で割り切れる。
  e_l(u) = d_l(u)/K について、母関数の分解（証明書）
      −Σ_{u≥1} e_l(u) z^{u−1} = Σ_S P_S(z) / Π_{r∈S}(1 − r z)、  P_S の係数は全て ≥ 0（S ⊂ R_l、deg_z P_S < |S|）
  を LP で探し、有理数で厳密に検証する。これは −e_l(u) = Σ_{S,k} p_{S,k} h_{u−1−k}(S)（p ≥ 0）と同値なので、全ての u ≥ 1 で d_l(u) ≤ 0。
  したがって G(u, d) ≤ 0（全ての u ≥ 1, d ≥ 0）。
  くさび B（u ≥ t）：鏡像型の恒等式 F_ρ(−u, t) = τ F_ρ'(−t, u)（τ:(e_𝐚, e_𝐪, e_𝐭) ↦ (−e_𝐚, e_𝐪, e_𝐪 − e_𝐭)、ρ' は ρ の転置）を
  u, t = 1..5 の 25 点で確かめる（両辺とも (u, t) の指数多項式で、根は u に {m_j^{−1}}、t に {m_k}。25 点で一致すれば全体で一致）。
  τ は係数の符号を保つので、くさび B はくさび A（ρ' ）に帰着する。

使い方：
  python scripts/v4_pretzel_np2.py struct              （d_l の符号（u ≤ 20）、K での割り算、根の集合、鏡像型の恒等式）
  python scripts/v4_pretzel_np2.py cert L RHO [SMIN PER_ROUND]   （e_L(u)（ρ = 2 または 11）の証明書を列生成 LP で探し、data/cache に保存。
                                                   SMIN：使う部分集合の最小の大きさ、PER_ROUND：1 回に加える部分集合の数）
  python scripts/v4_pretzel_np2.py verify              （保存した証明書を厳密に検証し、u = 1..40 で表示と直接の値を照合）
  python scripts/v4_pretzel_np2.py check NMIN NMAX [UMIN UMAX]   （NMIN ≤ t ≤ NMAX、UMIN ≤ |s| ≤ UMAX（省略時は NMIN..NMAX）の np の点を閉じた形で判定）
"""
import pickle
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from v4_closed import RHOS, L, M_ROOTS, h_complete, newton_coeffs  # noqa: E402
from v4_pretzel import closed_eval, closed_form, fmt, hat_direct  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ORDER = ["1", "X", "X2", "Y", "Yb"]
PAIRS = [(j, k) for j in ORDER for k in ORDER if not ({j, k} <= {"X2", "Y", "Yb"} and j != k)]   # 19 対
K0 = L.from_dict({(0, 0, 0): 1, (2, 0, 1): 1}) * L.from_dict({(2, 0, 1): 1, (0, 2, 0): 1}) * L.from_dict({(0, 0, 0): 1, (2, 2, 3): 1})
K3 = L.from_dict({(0, 0, 0): 1, (2, 0, 3): 1})
NAMES = {(0, 0, 0): "1", (2, 0, 2): "X", (4, 0, 4): "X2", (4, 4, 8): "Y", (4, -4, 4): "Yb"}


def cert_path(l, rho):
    return ROOT / f"data/v4/v4_pretzel_np2_cert_{l}_{''.join(map(str, rho))}.pkl"


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def rname(r):
    for a, an in NAMES.items():
        for b, bn in NAMES.items():
            if sub(a, b) == r:
                return an if bn == "1" else (f"1/{bn}" if an == "1" else f"{an}/{bn}")
    return str(r)


def roots_l(l):
    """R_l = {m_k/m_j : (j, k) ∈ PAIRS, k ∈ ORDER[l:]}【導出】。"""
    return sorted({sub(M_ROOTS[k], M_ROOTS[j]) for j, k in PAIRS if k in ORDER[l:]})


def factor_l(l):
    return K0 * K3 if l >= 3 else K0


_CS = None


def d_seq(rho, umax, cs=None):
    """{u: [d_0(u)..d_4(u)]}（くさび A、L）。"""
    global _CS
    if cs is None:
        _CS = _CS or closed_form("np")
        cs = _CS
    inner = [M_ROOTS[x] for x in ORDER]
    return {u: [L.from_dict(c) for c in newton_coeffs([closed_eval(cs[rho], "np", -u, u + n) for n in range(5)], inner)]
            for u in range(1, umax + 1)}


def e_seq(rho, l, umax, ds=None):
    ds = ds or d_seq(rho, umax)
    return [(ds[u][l] / factor_l(l)).to_dict() for u in range(1, umax + 1)]   # 割り切れなければ例外


# ---- 母関数の算術 ----------------------------------------------------------------------------------------

def charpoly(R):
    """Π_{r∈R}(1 − r z) の z の係数（dict のリスト）。"""
    P = [{(0, 0, 0): 1}]
    for r in R:
        new = [dict(x) for x in P] + [{}]
        for i, c in enumerate(P):
            for e, v in c.items():
                k = tuple(x + y for x, y in zip(e, r))
                new[i + 1][k] = new[i + 1].get(k, 0) - v
        P = [{k: v for k, v in x.items() if v} for x in new]
    return P


def numerator(seq, R):
    """N(z) = Π_R(1 − r z) Σ_{u≥1} e(u) z^{u−1} の係数（次数 < |R|）。seq は e(1), e(2), …（長さ ≥ |R| + 2）。それより先が 0 になることも確かめる。"""
    P = charpoly(R)
    N = []
    for k in range(len(seq)):
        acc = {}
        for i in range(min(k, len(R)) + 1):
            for e1, c1 in P[i].items():
                for e2, c2 in seq[k - i].items():
                    kk = tuple(x + y for x, y in zip(e1, e2))
                    acc[kk] = acc.get(kk, 0) + c1 * c2
        N.append({k2: v for k2, v in acc.items() if v})
    assert all(not x for x in N[len(R):]), "根の集合が足りない"
    return N[:len(R)]


def check_identity(N, R, cert):
    """cert = {(S, k, mono): Fraction ≥ 0} について Σ_S P_S(z) Π_{R∖S}(1 − r z) = −N(z) を厳密に確かめる。"""
    if any(v < 0 for v in cert.values()):
        return False
    acc, cps = {}, {}
    for (S, k, mono), x in cert.items():
        if S not in cps:
            cps[S] = charpoly([r for r in R if r not in S])
        for i, c in enumerate(cps[S]):
            for e, v in c.items():
                key = (k + i,) + tuple(a + b for a, b in zip(e, mono))
                acc[key] = acc.get(key, 0) + v * x
    tgt = {(k,) + e: -v for k, x in enumerate(N) for e, v in x.items()}
    return all(acc.get(key, 0) == tgt.get(key, 0) for key in set(acc) | set(tgt))


def cert_eval(cert, u):
    """Σ_{S,k} p_{S,k} 𝐦 h_{u−1−k}(S)（= −e(u)）。"""
    tot = {}
    for (S, k, mono), x in cert.items():
        for e, c in h_complete(u - 1 - k, list(S)).items():
            key = tuple(a + b for a, b in zip(e, mono))
            tot[key] = tot.get(key, 0) + c * x
    return {k: v for k, v in tot.items() if v}


# ---- 列生成 LP（部分集合を単位に追加する。P_S の台は N の台の和集合に制限） --------------------------

def colgen(N, R, M=None, per_round=10, max_rounds=500, tol=1e-7, method="highs-ipm", smin=1, save=None, prune=True):
    from scipy.optimize import linprog
    from scipy.sparse import csc_matrix
    red = lambda e: (e[0] // 2, e[1] // 2, e[2])                      # noqa: E731（𝐚, 𝐪 の指数は偶数）
    m, Rr = len(R), [red(r) for r in R]
    M = sorted({red(e) for e in (M or {e for x in N for e in x})})
    lo = np.array([min(e[i] for e in M) + sum(min(r[i], 0) for r in Rr) for i in range(3)])
    hi = np.array([max(e[i] for e in M) + sum(max(r[i], 0) for r in Rr) for i in range(3)])
    shape = (m,) + tuple(int(x) for x in hi - lo + 1)
    Mi = np.array([[e[i] - lo[i] for i in range(3)] for e in M])
    b = {(k,) + tuple(int(red(e)[i] - lo[i]) for i in range(3)): -v for k, x in enumerate(N) for e, v in x.items()}
    rows = {key: i for i, key in enumerate(b)}
    chosen, cps = (pickle.load(open(save, "rb")) if save and Path(save).exists() else []), {}   # 再開：選んだ部分集合を読み込む

    def cp(S):
        if S not in cps:
            cps[S] = [(i, e, v) for i, d in enumerate(charpoly([Rr[i] for i in range(m) if i not in S])) for e, v in d.items()]
        return cps[S]

    def shift_sub(Z, r):                                              # 相関 Z[k, x] − Z[k+1, x + r]
        out = Z.copy()
        src, dst = [slice(1, None)], [slice(0, -1)]
        for j in range(3):
            n, s = Z.shape[j + 1], r[j]
            src.append(slice(s, n) if s >= 0 else slice(0, n + s))
            dst.append(slice(0, n - s) if s >= 0 else slice(-s, n))
        out[tuple(dst)] -= Z[tuple(src)]
        return out

    for rnd in range(max_rounds):
        cols = [(S, k, j) for S in chosen for k in range(len(S)) for j in range(len(M))]
        ri, ci, data = [], [], []
        for c, (S, k0, j) in enumerate(cols):
            for i, e, v in cp(S):
                key = (k0 + i,) + tuple(int(Mi[j][x] + e[x]) for x in range(3))
                ri.append(rows.setdefault(key, len(rows)))
                ci.append(c)
                data.append(v)
        nr, nc = len(rows), len(cols)
        tgt = [rows[key] for key in b]                                # 人工変数（Phase I）は N の台の行だけ（他の行は x = 0 で満たされる）
        for i2, i in enumerate(tgt):
            ri += [i, i]
            ci += [nc + 2 * i2, nc + 2 * i2 + 1]
            data += [1.0, -1.0]
        A = csc_matrix((data, (ri, ci)), shape=(nr, nc + 2 * len(tgt)))
        bv = np.zeros(nr)
        for key, v in b.items():
            bv[rows[key]] = v
        t0 = time.time()
        res = linprog(np.concatenate([np.zeros(nc), np.ones(2 * len(tgt))]), A_eq=A, b_eq=bv, bounds=(0, None), method=method)
        import resource
        mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30                 # macOS はバイト
        print(f"  round {rnd}: 部分集合 {len(chosen)}、列 {nc}、行 {nr}、残差 {res.fun:.5g}（{time.time() - t0:.1f} 秒、最大メモリ {mem:.2f} GB）", flush=True)
        if res.fun < tol:
            back = lambda e: (2 * (e[0] + int(lo[0])), 2 * (e[1] + int(lo[1])), e[2] + int(lo[2]))  # noqa: E731
            return {(tuple(R[i] for i in S), k, back(tuple(Mi[j]))): res.x[c] for c, (S, k, j) in enumerate(cols) if res.x[c] > 1e-12}
        y = np.zeros(shape)
        for key, i in rows.items():
            y[key] = res.eqlin.marginals[i]
        if prune:                                                     # 解で使われていない部分集合を外す（必要なら価格付けで戻る）
            used = {cols[c][0] for c in range(nc) if res.x[c] > 1e-12}
            chosen = [S for S in chosen if S in used]
        scores, cset = [], set(chosen)

        def rec(i, Z, T):
            if i == m:
                S = tuple(x for x in range(m) if x not in T)
                if len(S) >= smin and S not in cset:
                    scores.append((float(Z[:len(S), Mi[:, 0], Mi[:, 1], Mi[:, 2]].max()), S))
                return
            rec(i + 1, Z, T)
            rec(i + 1, shift_sub(Z, Rr[i]), T | {i})
        rec(0, y, frozenset())
        scores.sort(reverse=True)
        new = [S for v, S in scores[:per_round] if v > tol]
        if not new:
            return None
        chosen += new
        if save:
            pickle.dump(chosen, open(save, "wb"))
    return None


def exact_vertex(N, R, sol, p=2 ** 61 - 1):
    """頂点解の台（列は一次独立）の上で、有理数の厳密解を解き直す。"""
    import flint
    cols = [c for c, x in sol.items() if x > 1e-9]
    rows, ent, cps = {}, {}, {}
    for j, (S, k, mono) in enumerate(cols):
        if S not in cps:
            cps[S] = charpoly([r for r in R if r not in S])
        for i, c in enumerate(cps[S]):
            for e, v in c.items():
                key = (k + i,) + tuple(x + y for x, y in zip(e, mono))
                r = rows.setdefault(key, len(rows))
                ent[(r, j)] = ent.get((r, j), 0) + v
    b = [0] * len(rows)
    for k, x in enumerate(N):
        for e, v in x.items():
            b[rows[(k,) + e]] = -v
    nr, nc = len(rows), len(cols)
    At = flint.nmod_mat(nc, nr, p)
    for (i, j), v in ent.items():
        At[j, i] = v % p
    Rf, rank = At.rref()
    if rank < nc:
        return None
    piv, r = [], 0
    for c in range(nr):
        if r < rank and Rf[r, c] != 0:
            piv.append(c)
            r += 1
    pos = {c: i for i, c in enumerate(piv)}
    Mq, bq = flint.fmpq_mat(nc, nc), flint.fmpq_mat(nc, 1)
    for (i, j), v in ent.items():
        if i in pos:
            Mq[pos[i], j] = v
    for c, i in pos.items():
        bq[i, 0] = b[c]
    x = Mq.solve(bq)
    return {cols[j]: Fraction(int(x[j, 0].p), int(x[j, 0].q)) for j in range(nc) if x[j, 0] != 0}


# ---- 実行 --------------------------------------------------------------------------------------------------

def struct():
    for rho in RHOS:
        ds = d_seq(rho, 20)
        for l in range(5):
            npos = sum(v > 0 for u in ds for v in ds[u][l].to_dict().values())
            e = e_seq(rho, l, 20, ds)                                 # K で割り切れること
            R = roots_l(l)
            numerator(e, R)                                           # R_l で足りること（u = 1..20）
            print(f"ρ={rho} l={l}: d_l(u) の正の係数 {npos} 個（u ≤ 20）、K で割り切れる、根 {len(R)} 個 {[rname(r) for r in R]}", flush=True)
    tau = lambda P: {(-a, q, q - t): c for (a, q, t), c in P.items()}   # noqa: E731
    other = {(2,): (1, 1), (1, 1): (2,)}
    ok = all(hat_direct(-u, t)[rho] == tau(hat_direct(-t, u)[other[rho]]) for u in range(1, 6) for t in range(1, 6) for rho in RHOS)
    print(f"鏡像型の恒等式 F_ρ(−u, t) = τ F_ρ'(−t, u)（u, t = 1..5 の 25 点、両 ρ）: {ok}", flush=True)


def cert(l, rho, smin=1, per_round=10):
    R = roots_l(l)
    N = numerator(e_seq(rho, l, len(R) + 2), R)
    print(f"l={l} ρ={rho}: 根 {len(R)} 個、N の項数 {[len(x) for x in N]}", flush=True)
    t0 = time.time()
    M = {e for x in N for e in x}
    for grow in range(3):                                             # P_S の台：N の台、だめなら根の ±1 個分ずつ広げる
        sol = colgen(N, R, M, per_round=per_round, smin=smin, save=ROOT / f"data/cache/v4_pretzel_np2_cg_{l}_{''.join(map(str, rho))}_{grow}.pkl")
        if sol is not None:
            break
        print(f"  P_S の台 {len(M)} 点では LP が不能。広げる", flush=True)
        M = M | {tuple(a + s_ * b for a, b in zip(e, r)) for e in M for r in R for s_ in (1, -1)}
    if sol is None:
        print("証明書なし", flush=True)
        return
    ex = exact_vertex(N, R, sol)
    ok = ex is not None and check_identity(N, R, ex)
    print(f"l={l} ρ={rho}: 厳密な証明書 {ok}、項数 {len(ex or {})}、部分集合 {len({S for S, _, _ in ex or {}})} 個、"
          f"分母の最大 {max((v.denominator for v in (ex or {}).values()), default=0)}（{time.time() - t0:.0f} 秒）", flush=True)
    if ok:
        pickle.dump(ex, open(cert_path(l, rho), "wb"))


def verify(umax=25):
    for rho in RHOS:
        ds = d_seq(rho, umax)
        for l in range(5):
            R = roots_l(l)
            e = e_seq(rho, l, umax, ds)
            path = cert_path(l, rho)
            if not path.exists():
                print(f"ρ={rho} l={l}: 証明書なし。u ≤ {umax} で e_l(u) ≤ 0: {all(v <= 0 for x in e for v in x.values())}", flush=True)
                continue
            C = pickle.load(open(path, "rb"))
            ok = check_identity(numerator(e[:len(R) + 2], R), R, C)
            direct = all(cert_eval(C, u) == {k: -v for k, v in e[u - 1].items()} for u in range(1, umax + 1))
            print(f"ρ={rho} l={l}: 証明書の厳密な検証 {ok}、u = 1..{umax} で表示 = −e_l(u) {direct}、"
                  f"部分集合 {len({S for S, _, _ in C})} 個 {sorted({tuple(rname(r) for r in S) for S, _, _ in C}, key=len)[:6]}…", flush=True)


def check(nmin, nmax, umin=None, umax=None):
    CS = closed_form("np")
    bad = 0
    for u in range(umin or nmin, (umax or nmax) + 1):
        for t in range(nmin, nmax + 1):
            t0 = time.time()
            res = [fmt(closed_eval(CS[rho], "np", -u, t)) for rho in RHOS]
            bad += sum(not (r.startswith("前 0/") and "後 0/" in r and "全体 −" in r) for r in res)
            msg = "；".join(f"{rho}: {r}" for rho, r in zip(RHOS, res))
            print(f"K_{{{-u},{t}}}: {msg}（{time.time() - t0:.1f} 秒）", flush=True)
    print(f"負の項を含むもの {bad} 件", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "struct":
        struct()
    elif mode == "cert":
        cert(int(sys.argv[2]), (2,) if sys.argv[3] == "2" else (1, 1), *map(int, sys.argv[4:6]))
    elif mode == "verify":
        verify()
    elif mode == "check":
        check(*map(int, sys.argv[2:6]))
