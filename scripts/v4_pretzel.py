"""第 4 期 V4 の続き：pretzel K_{s,t}（(2s−1, 1, 2t−1)-pretzel）の \\hat f_[2], \\hat f_[1,1] の 2 変数の閉じた形（ニュートン形）と正値性。

入力：GNSSS（arXiv:1512.07883）の予想式（scripts/t3_pretzel.py の tc_graded。2026-10-07 に浮動小数点の混入を修正した版）。
構造【導出】：varpi(p, i) は p について根 1, X, Y（X = 𝐚²𝐭², Y = 𝐚⁴𝐪⁴𝐭⁸）の指数多項式。P_r(s, t) は前因子 (−𝐭)^{−r·sgn} を除いて
  Σ_i varpi(s, i) varpi(t, i)(s, t によらない) なので、s, t それぞれについて根 {1, X, Y} の指数多項式。
  \\hat f は [2], [1,1]（ミラー：Y → Ȳ）, [1]², [1](x²)（Adams：X → X²）の一次結合なので、各変数について根 {1, X, X², Y, Ȳ}。
領域：pp（s, t > 0）、nn（s, t < 0）、np（s < 0 < t）。各領域で n = |s| − 1, m = |t| − 1 とし、根は p > 0 で m_j、p < 0 で m_j^{-1}。
2 変数のニュートン形：F(1+n, 1+m) = Σ_{k,l} c_{kl} h_{n−k}(m_{π0..πk}) h_{m−l}(m_{π0..πl})（c_{kl} は |s|, |t| = 1..5 の 25 点の差分）。
使い方：python scripts/v4_pretzel.py table SMAX      （|s|, |t| ≤ SMAX の直接計算：次元、special polynomial、正値性）
        python scripts/v4_pretzel.py newton          （各領域・ρ の c_{kl} の符号と 𝐪 の偶奇、σ 不変性）
        python scripts/v4_pretzel.py check NMAX      （|s|, |t| ≤ NMAX の閉じた形での正値性と、直接計算との照合）
        python scripts/v4_pretzel.py line NMAX       （K_{±1,t} と K_{s,s} 等の 1 パラメータ族を |t| ≤ NMAX で）
"""
import pickle
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from t3_pretzel import tc_graded  # noqa: E402
from v4_closed import RHOS, L, M_ROOTS, h_complete, hat_from_inputs, judge, newton_coeffs  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
INPUT_CACHE = ROOT / "data/cache/v4_pretzel_inputs.pkl"
ORDER = ["1", "X", "X2", "Y", "Yb"]                  # ツイスト結び目（v4_closed.ORDER）と同じ並び
REGIONS = {"pp": (1, 1), "nn": (-1, -1), "np": (-1, 1)}
# c_kl が全て ≤ 0 になる (s の並び, t の並び)（scripts/v4_pretzel_search.py で 120 × 120 通りから探索、pp・nn・両 ρ で 216 通り見つかったうちの一つ）。
# np では見つからない（logs/v4_pretzel_search_np*.log）。
PROOF_ORDERS = {(r, rho): (("1", "X", "X2", "Y", "Yb"), ("1", "X", "Yb", "Y", "X2")) for r in ("pp", "nn") for rho in [(2,), (1, 1)]}


def key(s, t):
    return (min(s, t), max(s, t))                     # K_{s,t} = K_{t,s}


_cache = None


def inputs(s, t):
    """GNSSS の [1], [2]（太字、t_c 次数）。data/cache/v4_pretzel_inputs.pkl に保存（再開用）。"""
    global _cache
    if _cache is None:
        _cache = pickle.load(open(INPUT_CACHE, "rb")) if INPUT_CACHE.exists() else {}
    k = key(s, t)
    if k not in _cache:
        _cache[k] = (tc_graded(*k, 1), tc_graded(*k, 2))
        pickle.dump(_cache, open(INPUT_CACHE, "wb"))
    return _cache[k]


def hat_direct(s, t):
    return hat_from_inputs(*inputs(s, t))


def euler_q1(P):
    """太字の Poincaré 多項式の 𝐭 = −1, 𝐪 = 1：{e_𝐚: c}。"""
    out = {}
    for (ea, eq, et), c in P.items():
        out[ea] = out.get(ea, 0) + c * (-1) ** (et % 2)
    return {k: v for k, v in out.items() if v}


def special_ok(P1, P2):
    e1, e2 = euler_q1(P1), euler_q1(P2)
    sq = {}
    for i, x in e1.items():
        for j, y in e1.items():
            sq[i + j] = sq.get(i + j, 0) + x * y
    return {k: v for k, v in sq.items() if v} == e2


def sigma(P):
    return {(a, -q, t - q): c for (a, q, t), c in P.items()}


# ---- 2 変数のニュートン形 ---------------------------------------------------------------------------

def roots(sign, names=ORDER):
    return [M_ROOTS[k] if sign == 1 else tuple(-e for e in M_ROOTS[k]) for k in names]


def orders_of(region, rho, orders):
    """orders が None なら (ORDER, ORDER)、dict なら ORDERS[(region, ρ)]。"""
    return (ORDER, ORDER) if orders is None else orders[(region, rho)]


def closed_form(region, orders=None, hat=None):
    """{ρ: c[k][l]}（太字 dict）。|s|, |t| = 1..5 の 25 点の直接計算から作る。"""
    hat = hat or hat_direct
    ss, st = REGIONS[region]
    grid = {(i, j): hat(ss * (1 + i), st * (1 + j)) for i in range(5) for j in range(5)}
    out = {}
    for rho in RHOS:
        os_, ot = orders_of(region, rho, orders)
        rs, rt = roots(ss, os_), roots(st, ot)
        rows = [newton_coeffs([grid[(i, j)][rho] for i in range(5)], rs) for j in range(5)]   # rows[j][k]
        out[rho] = [newton_coeffs([rows[j][k] for j in range(5)], rt) for k in range(5)]      # c[k][l]
    return out


_h = {}


def hL(n, ms):
    k = (n, tuple(ms))
    if k not in _h:
        _h[k] = L.from_dict(h_complete(n, list(ms)))
    return _h[k]


def closed_eval(cs, region, s, t, order_pair=(ORDER, ORDER)):
    ss, st = REGIONS[region]
    rs, rt = roots(ss, order_pair[0]), roots(st, order_pair[1])
    n, m = abs(s) - 1, abs(t) - 1
    tot = L.from_dict({})
    for k in range(5):
        if n - k < 0:
            continue
        for l in range(5):
            if m - l < 0 or not cs[k][l]:
                continue
            tot = tot + L.from_dict(cs[k][l]) * hL(n - k, rs[:k + 1]) * hL(m - l, rt[:l + 1])
    return tot.to_dict()


def region_of(s, t):
    if s > 0 and t > 0:
        return "pp", s, t
    if s < 0 and t < 0:
        return "nn", s, t
    return "np", min(s, t), max(s, t)


# ---- 実行 ----------------------------------------------------------------------------------------------

def fmt(P):
    b, a_, nb, na = judge(P)
    flip = b > nb // 2                                 # 全体の符号を良い方に
    if flip:
        b = nb - b
        a_ = None if a_ is None else na - a_
    return f"前 {b}/{nb} 負、" + ("G 展開が閉じない" if a_ is None else f"後 {a_}/{na} 負") + ("（全体 +）" if flip else "（全体 −）")


def table(smax):
    for s in range(-smax, smax + 1):
        for t in range(s, smax + 1):
            if not (s and t):
                continue
            t0 = time.time()
            P1, P2 = inputs(s, t)
            H = hat_from_inputs(P1, P2)
            r1, r2 = sum(P1.values()), sum(P2.values())
            print(f"K_{{{s},{t}}}: 次元 {r1}, {r2}（二乗 {r2 == r1 ** 2}）、special polynomial {special_ok(P1, P2)}、"
                  + "；".join(f"{rho}: {fmt(H[rho])}" for rho in RHOS) + f"（{time.time() - t0:.1f} 秒）", flush=True)


def newton():
    for region in REGIONS:
        orders = PROOF_ORDERS if region in ("pp", "nn") else None
        print(f"{region}: 並び " + ("PROOF_ORDERS " + str(PROOF_ORDERS[(region, (2,))]) if orders else f"({ORDER}, {ORDER})"), flush=True)
        CS = closed_form(region, orders)
        for rho in RHOS:
            c = CS[rho]
            pos = [[sum(v > 0 for v in c[k][l].values()) for l in range(5)] for k in range(5)]
            neg = [[sum(v < 0 for v in c[k][l].values()) for l in range(5)] for k in range(5)]
            even = all(e[1] % 2 == 0 for k in range(5) for l in range(5) for e in c[k][l])
            print(f"{region} ρ={rho}: c_kl の正の係数の数 {pos}", flush=True)
            print(f"{region} ρ={rho}: c_kl の負の係数の数 {neg}、𝐪 の指数が偶数 {even}", flush=True)
        ss, st = REGIONS[region]
        sym = all(hat_direct(ss * i, st * j)[rho] == sigma(hat_direct(ss * i, st * j)[rho])
                  for i in range(1, 6) for j in range(1, 6) for rho in RHOS)
        print(f"{region}: 25 点で σ 不変 {sym}", flush=True)
        pickle.dump(CS, open(ROOT / f"data/cache/v4_pretzel_newton_{region}.pkl", "wb"))


def check(nmax, direct_pts=()):
    CS = {r: closed_form(r) for r in REGIONS}
    for s in range(-nmax, nmax + 1):
        for t in range(s, nmax + 1):
            if not (s and t):
                continue
            t0 = time.time()
            region, s2, t2 = region_of(s, t)
            H = {rho: closed_eval(CS[region][rho], region, s2, t2) for rho in RHOS}
            msg = "；".join(f"{rho}: {fmt(H[rho])}" for rho in RHOS)
            if (s, t) in direct_pts or (max(abs(s), abs(t)) > 5 and max(abs(s), abs(t)) <= 7):
                D = hat_direct(s, t)
                msg += "；直接計算と一致 " + str(all(D[rho] == H[rho] for rho in RHOS))
            print(f"K_{{{s},{t}}}: {msg}（{time.time() - t0:.1f} 秒）", flush=True)


# ---- 領域の間の関係：pp の指数多項式の負の s, t への延長 ----------------------------------------------

def char_poly(ms):
    """Π (x − m_j) の係数（低次から、L）。"""
    char = [L.from_dict({(0, 0, 0): 1})]
    for m in ms:
        new = [L.from_dict({}) for _ in range(len(char) + 1)]
        for i, c in enumerate(char):
            new[i + 1] = new[i + 1] + c
            new[i] = new[i] - L.from_dict({m: 1}) * c
        char = new
    return char


def extend_back(seq, lo):
    """seq = {n: L}（n = 1..5 を含む）を、根 roots(1) の線形漸化式で n = lo..0 へ後ろ向きに延長する。"""
    char = char_poly(roots(1))
    for n in range(0, lo - 1, -1):
        acc = L.from_dict({})
        for i in range(1, 6):
            acc = acc + char[i] * seq[n + i]
        seq[n] = (L.from_dict({}) - acc) / char[0]
    return seq


def continuation(lo=-5):
    """pp の 25 点から、F(s, t)（s, t = lo..5）を漸化式で延長する。{(s, t): {ρ: dict}}。"""
    F = {}
    for rho in RHOS:
        cols = {}
        for t in range(1, 6):
            cols[t] = extend_back({s: L.from_dict(hat_direct(s, t)[rho]) for s in range(1, 6)}, lo)
        for s in range(lo, 6):
            row = extend_back({t: cols[t][s] for t in range(1, 6)}, lo)
            for t in range(lo, 6):
                F.setdefault((s, t), {})[rho] = row[t].to_dict()
    return F


def regions_check():
    """\\hat f(K_{s,t}) = 𝐭^{2(1 − sgn)} F(s, t)（sgn = (sign s + sign t)/2）と、F(0, t) = F(s, 0) = 0 を確かめる。"""
    F = continuation()
    zero = all(not F[(0, t)][rho] and not F[(t, 0)][rho] for t in range(-5, 6) for rho in RHOS)
    print(f"F(0, t) = F(s, 0) = 0（|s|, |t| ≤ 5）: {zero}", flush=True)
    for s in range(-5, 6):
        for t in range(-5, 6):
            if s and t and not (s > 0 and t > 0):
                sgn = ((s > 0) - (s < 0) + (t > 0) - (t < 0)) // 2
                eps = 2 * (1 - sgn)
                D = hat_direct(s, t)
                ok = all(D[rho] == {(a, q, x + eps): c for (a, q, x), c in F[(s, t)][rho].items()} for rho in RHOS)
                if not ok:
                    print(f"K_{{{s},{t}}}: 不一致", flush=True)
    print("負の s または t の 75 点の照合が終了（不一致の行が無ければ全て一致）", flush=True)


def lines(nmax):
    """1 パラメータ族 K_{−1,t}（|t| ≤ nmax）と対角 K_{t,t}, K_{−t,t}, K_{−t,−t}（t ≤ nmax）を閉じた形で判定。"""
    CS = {r: closed_form(r) for r in REGIONS}
    pts = [(-1, t) for t in range(-nmax, nmax + 1) if t] + [(x * t, t) for t in range(1, nmax + 1) for x in (1, -1)] \
        + [(-t, -t) for t in range(1, nmax + 1)]
    for s, t in pts:
        region, s2, t2 = region_of(s, t)
        H = {rho: closed_eval(CS[region][rho], region, s2, t2) for rho in RHOS}
        print(f"K_{{{s},{t}}}: " + "；".join(f"{rho}: {fmt(H[rho])}" for rho in RHOS), flush=True)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "table":
        table(int(sys.argv[2]))
    elif mode == "newton":
        newton()
    elif mode == "regions":
        regions_check()
    elif mode == "line":
        lines(int(sys.argv[2]))
    elif mode == "check":
        check(int(sys.argv[2]), direct_pts={(-10, 9), (8, 10), (-9, -8), (-12, 12)})
