r"""U = t^N の冪の係数 c_j(q, t) の s の次数の範囲を押さえて、少ない N から R_λ(U) = Σ_j c_j U^j を決める（progress.md 続き 71）。
  R(N) は U = s^{2N} を代入した (r, s) の Laurent 多項式（rk_c_interp.R_true と同じ規格化。約分は flint で行う）。
  r の次数 α ごとに、未知数 x_{j,β}（c_j の r^α s^β の係数）について
    Σ_j x_{j, d − 2Nj} = [R(N)]_{r^α s^d}（使う N すべて、すべての d）
  を解く。U の冪の範囲 [j_min, j_max] は s の次数の最小・最大の N に対する傾きから（rk_c_read_largeN と同じ）。
  x_{j,β} の範囲【仮定】（打ち消しで外にはみ出す項があれば、解が無いか、使わない N での照合で外れる）：
    mode "N0"：使う最小の N₀ で、各項 c_j U^j が R(N₀) の s の次数の範囲に収まる（β + 2N₀ j ∈ [m(N₀), M(N₀)]）。
    mode "allS"：使う N すべてで、c_j U^j が R(N) の s の次数の範囲に収まる（帯の共通部分）。
    mode "alpha"：N₀ の帯に加えて、r の次数 α ごとの U の冪の範囲を、使う N のうち大きい二つでの r^α の部分の s の次数の傾きから読む。
    mode "alphaL"：mode "alpha" に加えて、使う最大の N での r^α の部分の s の次数の範囲にも収まる。
    mode "all"：さらに r の次数 α ごとに、c_j U^j の r^α の部分が R(N) の r^α の部分の s の次数の範囲に収まる。
  解が一意で、使わない N すべてで全係数一致すれば採用する。
使い方：python rk_c_bounded_interp.py m λ Nmin Nmax_used [mode [Nstrip]]（例：1 111 4 9 alpha 4）"""
import sys, os, pickle
from fractions import Fraction
from collections import defaultdict
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import flint
from rcsknots.refined_kauffman.flint_rf import from_sympy, RF, CTX

OUT = "data/refined_kauffman/daha"
UNIT = 2       # U = s^{UNIT·N}（SO の 2 箱・4 箱は U = t^N = s^{2N}。随伴は A^{1/2} = s^N なので 1。rk_c_bounded_save が切り替える）
PAR = 2        # c_j の s の冪の偶奇（SO のデータは一定。随伴は 1 にして偶奇の制限を外す）


def _mirror(p):
    """多項式 p(r, s) → p(1/r, 1/s) × r^A s^B（A, B は p の最大次数）。(多項式, A, B) を返す。"""
    d = p.to_dict()
    A = max(e[0] for e in d); B = max(e[1] for e in d)
    return CTX.from_dict({(A - e[0], B - e[1]): c for e, c in d.items()}), A, B


CACHE = "data/cache/rk_Rterms"


def terms_true(N, m, lam):
    """R_λ(N) = T_λ^{−2n} × 鏡像(Z) / K を {(α, β): Fraction} で返す（保存値が無ければ None）。data/cache にキャッシュする。"""
    key = f"N{N}_m{m}_{''.join(map(str, lam))}"
    path = os.path.join(OUT, f"Z_{key}.pkl")
    if not os.path.exists(path):
        return None
    cpath = os.path.join(CACHE, f"R_{key}.pkl")
    if os.path.exists(cpath) and os.path.getmtime(cpath) >= os.path.getmtime(path):
        return pickle.load(open(cpath, "rb"))
    out = _terms_true(N, m, lam, path)
    os.makedirs(CACHE, exist_ok=True)
    pickle.dump(out, open(cpath, "wb"))
    return out


def _terms_true(N, m, lam, path):
    dat = pickle.load(open(path, "rb"))
    Z, K = from_sympy(dat["Z"]), from_sympy(dat["K"])
    if Z.is_zero():
        return None                             # λ の行の数が N を超える
    zn, A1, B1 = _mirror(Z.n)
    zd, A2, B2 = _mirror(Z.d)
    lamp = list(lam) + [0] * (N - len(lam))
    ll = sum(x * x for x in lamp)
    lr = sum(x * (N - 1 - i) for i, x in enumerate(lamp))
    n = 2 * m + 1
    sh_r = -2 * n * ll - A1 + A2               # 鏡像の分子・分母の補正と T_λ^{−2n} = r^{−2n ll} s^{−4n lr}
    sh_s = -4 * n * lr - B1 + B2
    W = RF(zn * K.d, zd * K.n)
    assert W.is_laurent(), (N, lam)
    (dm, dc), = W.d.to_dict().items()
    return {(e[0] - dm[0] + sh_r, e[1] - dm[1] + sh_s): Fraction(int(c), int(dc)) for e, c in W.n.to_dict().items()}


def s_range(T):
    return min(b for _, b in T), max(b for _, b in T)


def slopes(T):
    Ns = sorted(T)
    rng = {N: s_range(T[N]) for N in Ns}
    d = [((rng[b][0] - rng[a][0]) // (UNIT * (b - a)), (rng[b][1] - rng[a][1]) // (UNIT * (b - a))) for a, b in zip(Ns, Ns[1:])]
    return rng, d


def _alpha_jrange(T, used, al, jmin, jmax):
    """r^α の部分の s の次数の範囲の、使う N のうち大きい二つの間の傾き（なければ全体の範囲）。"""
    N1, N2 = sorted(used)[-2:]
    b1 = [b for a, b in T[N1] if a == al]; b2 = [b for a, b in T[N2] if a == al]
    if not b1 or not b2:
        return jmin, jmax
    d = UNIT * (N2 - N1)
    return max(jmin, (min(b2) - min(b1)) // d), min(jmax, (max(b2) - max(b1)) // d)


def _bounds(T, used, al, mode, Ns=None):
    """mode に応じた、各 N の帯 [lo_N, hi_N]（c_j U^j の s の次数 β + 2Nj が収まる範囲）のリスト。"""
    if mode in ("N0", "alpha", "alphaL"):
        N0 = Ns or min(used)
        out = [(N0, *s_range(T[N0]))]
        bs = [b for a, b in T[max(used)] if a == al]
        if mode == "alphaL" and bs:
            out.append((max(used), min(bs), max(bs)))
        return out
    if mode == "allS":
        return [(N, *s_range(T[N])) for N in used]
    out = []
    for N in used:
        bs = [b for a, b in T[N] if a == al]
        if bs:
            out.append((N, min(bs), max(bs)))
    return out


def solve_bounded(T, used, jmin, jmax, mode="all", info=None, Ns=None, partial=False):
    """T：{N: R(N) の項}。used の N で解く。戻り値 ({j: {(α, β): 係数}} または None, 状態, 未知数の総数)。
    Ns：帯を取る N（既定は使う最小の N。小さい N は打ち消しで s の次数の範囲が狭く、帯の仮定が破れることがある）。
    partial = True なら、一意でないときも一意に決まった α の部分を返す。
    info にリストを渡すと、α ごとに (α, U の冪の範囲, 未知数, 階数) を追記する。"""
    alphas = sorted({a for N in used for a, _ in T[N]})
    sol = defaultdict(dict)
    unique, nvar = True, 0
    for al in alphas:
        bd = _bounds(T, used, al, mode, Ns)
        N0, lo0, hi0 = bd[0]
        ja, jb = _alpha_jrange(T, used, al, jmin, jmax) if mode in ("alpha", "alphaL") else (jmin, jmax)
        var = []
        for j in range(ja, jb + 1):
            for be in range(lo0 - UNIT * N0 * j, hi0 - UNIT * N0 * j + 1):
                if (be - lo0) % PAR == 0 and all(lo <= be + UNIT * N * j <= hi for N, lo, hi in bd):
                    var.append((j, be))
        nvar += len(var)
        idx = {v: i for i, v in enumerate(var)}
        rows, rhs = [], []
        for N in used:
            ds = defaultdict(list)
            for (j, be), i in idx.items():
                ds[be + UNIT * N * j].append(i)
            target = {b: c for (a, b), c in T[N].items() if a == al}
            if any(d not in ds for d in target):
                return None, f"矛盾（r^{al}、N={N} で範囲の外に項）", nvar                # 範囲の外に項がある
            for d, ids in ds.items():
                row = [0] * len(var)
                for i in ids:
                    row[i] = 1
                rows.append(row); rhs.append(target.get(d, Fraction(0)))
        nv = len(var)
        if nv == 0:
            continue
        Ab = flint.fmpq_mat(len(rows), nv + 1, [x for row, r_ in zip(rows, rhs) for x in row + [flint.fmpq(r_.numerator, r_.denominator)]])
        R_, rk = Ab.rref()
        piv = []
        for i in range(rk):
            c = next(c for c in range(nv + 1) if R_[i, c] != 0)
            if c == nv:
                return None, f"矛盾（r^{al} の連立方程式に解が無い）", nvar                 # 範囲の仮定か冪の範囲が違う
            piv.append(c)
        if info is not None:
            info.append((al, (ja, jb), nv, rk))
        if rk < nv:
            unique = False
            continue                                    # 一意でない：この α は解かない
        for i, c in enumerate(piv):
            v = R_[i, nv]
            if v != 0:
                j, be = var[c]
                sol[j][(al, be)] = Fraction(int(v.p), int(v.q))
    return (sol, "一意", nvar) if unique else ((sol if partial else None), "一意でない", nvar)


def predict(sol, N):
    out = defaultdict(Fraction)
    for j, c in sol.items():
        for (a, b), v in c.items():
            out[(a, b + UNIT * N * j)] += v
    return {k: v for k, v in out.items() if v != 0}


if __name__ == "__main__":
    m, lam = int(sys.argv[1]), tuple(int(c) for c in sys.argv[2])
    Nmin, Nmax = int(sys.argv[3]), int(sys.argv[4])
    mode = sys.argv[5] if len(sys.argv) > 5 else "all"
    Ns = int(sys.argv[6]) if len(sys.argv) > 6 else None
    T = {}
    for N in range(2, 40):
        t = terms_true(N, m, lam)
        if t is not None:
            T[N] = t
    rng, d = slopes(T)
    jmin, jmax = d[-1]
    print(f"λ={lam}：N = {sorted(T)}、s の次数の範囲 {[rng[N] for N in sorted(T)]}、隣の N との傾き（U の冪の最小・最大）{d}", flush=True)
    used = [N for N in sorted(T) if Nmin <= N <= Nmax]
    sol, st, nvar = solve_bounded(T, used, jmin, jmax, mode, Ns=Ns)
    if sol is None:
        print(f"   [{mode}{'' if Ns is None else f', 帯 N={Ns}'}] 使う N = {used}、U^{jmin}〜U^{jmax}、未知数 {nvar}：{st}", flush=True)
    else:
        ok = {N: predict(sol, N) == T[N] for N in T if N not in used}
        print(f"   [{mode}{'' if Ns is None else f', 帯 N={Ns}'}] 使う N = {used}、U^{jmin}〜U^{jmax}、未知数 {nvar}：{st}、使わない N での照合 {ok}", flush=True)
