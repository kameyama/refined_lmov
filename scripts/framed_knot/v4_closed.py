"""第 4 期 V4：ツイスト結び目 K_p の \\hat f_[2], \\hat f_[1,1] の厳密な直接計算と、p についての閉じた形（ニュートン形）。

直接計算：GNSSS の予想式（scripts/framed_knot/t3_pretzel.py の tc_graded、s = 1）の [1], [2] と、ミラーで作った [1,1] から、
  \\hat f_ρ = M(t)^{-1} f^q を 2 箱の明示式（src/rcsknots/bps.py の f_formulas, M_matrix から sympy で一度だけ導く）で厳密に求める。
  多項式演算は python-flint（分母は厳密な割り算で消す。割り切れなければ例外）。
閉じた形：p の符号ごとに \\hat f(K_p) = Σ_j A_j m_j^p（m_j は 5 つの単項式）。これをニュートン形
  s_p = Σ_{k=0}^{4} c_k h_{n−k}(m_{π0}, …, m_{πk})   （n = |p| − n0、h は完全斉次対称式）
  で表す。c_k = [(E − m_{π0})⋯(E − m_{π(k−1)}) s]_{n0} は s_{n0}, …, s_{n0+k} の整数多項式係数の一次結合なので、Laurent 多項式。
変数はすべて太字 (𝐚, 𝐪, 𝐭)。a 次数シフトは 2 箱なので 0（rcsknots.extract.A_SHIFT）。
"""
import pickle
import sys
from functools import lru_cache
from pathlib import Path

import flint
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from t3_pretzel import tc_graded  # noqa: E402
from t4_colored import bold_to_small  # noqa: E402
from rcsknots.bps import C, C1, M_matrix, Q, T, f_formulas, to_bold  # noqa: E402
from rcsknots.partitions import arm, boxes, leg  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
INPUT_CACHE = ROOT / "data/cache/v4_inputs.pkl"
RHOS = [(2,), (1, 1)]
CTX = flint.fmpz_mpoly_ctx.get(("a", "q", "t"), "lex")    # 小文字 (A, Q, T) と太字 (𝐚, 𝐪, 𝐭) の両方に使う（3 変数）


# ---- Laurent 多項式（flint の多項式と指数のずらし） -------------------------------------------------

class L:
    def __init__(self, poly, shift=(0, 0, 0)):
        self.p, self.s = poly, tuple(shift)

    @staticmethod
    def from_dict(d):
        if not d:
            return L(CTX.from_dict({}))
        lo = tuple(min(k[i] for k in d) for i in range(3))
        return L(CTX.from_dict({tuple(k[i] - lo[i] for i in range(3)): int(c) for k, c in d.items() if c}), lo)

    def to_dict(self):
        return {tuple(int(e[i]) + int(self.s[i]) for i in range(3)): int(c) for e, c in self.p.to_dict().items() if c}

    def _align(self, o):
        lo = tuple(min(self.s[i], o.s[i]) for i in range(3))
        up = lambda x: x.p * CTX.from_dict({tuple(x.s[i] - lo[i] for i in range(3)): 1})  # noqa: E731
        return up(self), up(o), lo

    def __add__(self, o):
        a, b, lo = self._align(o)
        return L(a + b, lo)

    def __sub__(self, o):
        a, b, lo = self._align(o)
        return L(a - b, lo)

    def __mul__(self, o):
        if isinstance(o, int):
            return L(self.p * o, self.s)
        return L(self.p * o.p, tuple(self.s[i] + o.s[i] for i in range(3)))

    def __truediv__(self, o):
        """Laurent 多項式としての厳密な割り算（割り切れなければ例外）。単項式の因子をずらしに移してから割る。"""
        a, b = self.normalized(), o.normalized()
        quo = a.p / b.p
        if quo * b.p != a.p:
            raise ArithmeticError("not divisible")
        return L(quo, tuple(a.s[i] - b.s[i] for i in range(3)))

    def normalized(self):
        return L.from_dict(self.to_dict())

    def is_zero(self):
        return self.p.is_zero()


def mono(e):
    return L(CTX.from_dict({(0, 0, 0): 1}), e)


# ---- 直接計算 ------------------------------------------------------------------------------------

@lru_cache(maxsize=None)
def formula():
    """\\hat f_ρ = (Σ_i n_{ρ,i} R_i) / D の (n_{ρ,i}, D)。R = (R2, R11, R1², R1(x²))。小文字 (A, Q, T) の多項式（L）。"""
    Asym = sp.Symbol("A")

    def unknot(lam, d=1):
        k = sum(lam)
        out = T ** (d * k) * Asym ** (-d * k)
        for i, j in boxes(lam):
            out *= (T ** (2 * d * i) - Asym ** (2 * d) * Q ** (2 * d * j)) / (1 - Q ** (2 * d * arm(lam, i, j)) * T ** (2 * d * (leg(lam, i, j) + 1)))
        return out

    F = f_formulas(2)
    fq = sp.Matrix([F[(2,)][0], F[(1, 1)][0]])
    hat = M_matrix(2).inv() * fq
    syms = [C((2,)), C((1, 1)), C((1,)), C1(2)]
    out = {}
    for rho, h in zip(RHOS, hat):
        h = sp.factor(sp.together(h))
        # h は C について 2 次（C_1² を含む）。各項の係数を個別に約分してから共通分母にそろえる。
        P = sp.Poly(sp.numer(h), *syms)
        den0 = sp.denom(h)
        parts = {}
        for mon, c in P.terms():
            key = {(1, 0, 0, 0): 0, (0, 1, 0, 0): 1, (0, 0, 2, 0): 2, (0, 0, 0, 1): 3}[mon]
            u = [unknot((2,)), unknot((1, 1)), unknot((1,)) ** 2, unknot((1,), 2)][key]
            parts[key] = sp.cancel(sp.together(c * u / den0))
        dens = [sp.denom(parts[k]) for k in range(4)]
        D = sp.lcm_list(dens) if hasattr(sp, "lcm_list") else sp.lcm(dens)
        coeffs = [sp.cancel(parts[k] * D) for k in range(4)]

        def to_L(e):
            P = sp.Poly(sp.expand(e * (Asym * Q * T) ** 40), Asym, Q, T)
            return L.from_dict({(i - 40, j - 40, k - 40): int(c) for (i, j, k), c in P.terms()})
        out[rho] = ([to_L(c) for c in coeffs], to_L(D))
    return out


def load_inputs():
    return pickle.load(open(INPUT_CACHE, "rb")) if INPUT_CACHE.exists() else {}


def inputs(p, cache=None):
    """GNSSS の [1], [2]（太字、t_c 次数）。data/cache/v4_inputs.pkl に p ごとに保存（再開用）。"""
    cache = load_inputs() if cache is None else cache
    if p not in cache:
        cache[p] = (tc_graded(1, p, 1), tc_graded(1, p, 2))
        pickle.dump(cache, open(INPUT_CACHE, "wb"))
    return cache[p]


def mir(P):
    return {(x, -y, z - y): c for (x, y, z), c in P.items()}


def hat_from_inputs(P1, P2):
    """太字の [1], [2] から \\hat f_[2], \\hat f_[1,1]（太字、{(e_𝐚, e_𝐪, e_𝐭): c}）。"""
    r1, r2, r11 = (L.from_dict(bold_to_small(P)) for P in (P1, P2, mir(P2)))
    r1d = L.from_dict({(2 * i, 2 * j, 2 * k): c for (i, j, k), c in bold_to_small(P1).items()})
    R = [r2, r11, r1 * r1, r1d]
    out = {}
    for rho, (ns, den) in formula().items():
        num = ns[0] * R[0] + ns[1] * R[1] + ns[2] * R[2] + ns[3] * R[3]
        out[rho] = to_bold((num / den).to_dict(), 0)
    return out


def hat_direct(p, cache=None):
    return hat_from_inputs(*inputs(p, cache))


# ---- 閉じた形（ニュートン形） ------------------------------------------------------------------------

M_ROOTS = {"1": (0, 0, 0), "X": (2, 0, 2), "X2": (4, 0, 4), "Y": (4, 4, 8), "Yb": (4, -4, 4)}   # 太字 (e_𝐚, e_𝐪, e_𝐭)


def h_complete(n, ms):
    """完全斉次対称式 h_n(m_0, …, m_k)（m は単項式の指数）を太字の dict で。h_n = Σ_i m_k^i h_{n−i}(m_0..m_{k−1}) で逐次に。"""
    if n < 0:
        return {}
    table = [{(0, 0, 0): 1}] + [{} for _ in range(n)]
    for m in ms:
        new = []
        for d in range(n + 1):
            acc = {}
            for i in range(d + 1):
                for e, c in table[d - i].items():
                    k = tuple(e[x] + i * m[x] for x in range(3))
                    acc[k] = acc.get(k, 0) + c
            new.append(acc)
        table = new
    return table[n]


def newton_coeffs(seq, order):
    """seq = [s_{n0}, …, s_{n0+4}]（太字 dict）、order = 根の並び → [c_0..c_4]（太字 dict）。"""
    cur = [L.from_dict(s) for s in seq]
    cs = []
    for k in range(5):
        cs.append(cur[0].to_dict())
        m = mono(order[k])
        cur = [cur[i + 1] - m * cur[i] for i in range(len(cur) - 1)]
    return cs


def newton_eval(cs, order, n):
    """s_{n0+n} = Σ_k c_k h_{n−k}(m_{π0..πk})。"""
    tot = L.from_dict({})
    for k, c in enumerate(cs):
        h = h_complete(n - k, order[:k + 1])
        if h and c:
            tot = tot + L.from_dict(c) * L.from_dict(h)
    return tot.to_dict()


def neg_count(P):
    """全体の符号を良い方に選んだときの負の係数の数。"""
    return min(sum(v < 0 for v in P.values()), sum(v > 0 for v in P.values()))


# ---- 実行：|p| ≤ 30 の閉じた形と、正値性の判定 -------------------------------------------------------

ORDER = ["1", "X", "X2", "Y", "Yb"]                 # ニュートン形の根の並び（scripts/framed_knot/v4_newton.py で c_k が一様な符号になる並びの一つ）
N0 = 1                                              # c_k は |p| = 1..5 の直接計算から作る


def roots_for(sign):
    """|p| = n の関数としての特性根：p > 0 は m_j、p < 0 は m_j^{-1}。"""
    return [M_ROOTS[k] if sign == 1 else tuple(-e for e in M_ROOTS[k]) for k in ORDER]


def closed_form(sign, cache=None):
    """{ρ: [c_0..c_4]}（p の符号ごと）。"""
    S = {n: hat_direct(sign * n, cache) for n in range(N0, N0 + 5)}
    return {rho: newton_coeffs([S[N0 + i][rho] for i in range(5)], roots_for(sign)) for rho in RHOS}


def closed_eval(cs, sign, p):
    return newton_eval(cs, roots_for(sign), abs(p) - N0)


def judge(P):
    """G 展開の前後の負の項の数（全体の符号は −1 に固定：c_k の符号）と、G 展開が閉じるか。"""
    from rcsknots.bps import g_expand
    before = sum(v > 0 for v in P.values())
    try:
        N = g_expand(P)
    except ValueError:
        return before, None, len(P), None
    return before, sum(v > 0 for v in N.values()), len(P), len(N)


def main(pmax=30, direct_extra=(15, 20, 30)):
    import time
    cache = load_inputs()
    for sign in (1, -1):
        t0 = time.time()
        CS = closed_form(sign, cache)
        for rho in RHOS:
            print(f"p の符号 {sign} ρ={rho}: c_k の正の係数の数 {[sum(v > 0 for v in c.values()) for c in CS[rho]]}、"
                  f"項数 {[len(c) for c in CS[rho]]}", flush=True)
        for n in range(1, pmax + 1):
            p = sign * n
            H = {rho: closed_eval(CS[rho], sign, p) for rho in RHOS}
            msg = []
            for rho in RHOS:
                b, a_, nb, na = judge(H[rho])
                msg.append(f"{rho}: 前 {b}/{nb} 負、" + ("G 展開が閉じない" if a_ is None else f"後 {a_}/{na} 負"))
            if n <= 12 or n in direct_extra:
                D = hat_direct(p, cache)
                msg.append("直接計算と一致 " + str(all(D[rho] == H[rho] for rho in RHOS)) + ("（c_k の作成に使用）" if n < N0 + 5 else ""))
            print(f"p={p}: " + "；".join(msg), flush=True)
        print(f"p の符号 {sign}: {time.time() - t0:.0f} 秒", flush=True)


if __name__ == "__main__":
    main()
