"""DGR / GGS の微分 d_{n|m} に当たる構造を refined LMOV 不変量（spin 表示・genus 表示）で調べる道具（docs/dgr_lmov/report.md）。

変数は太字 (𝐚, 𝐪, 𝐭)、多項式は {(e_𝐚, e_𝐪, e_𝐭): c}。GGS 1304.3481 の d_{n|m} は次数 (a,q,t_c) = (−2, 2(n−m), −2m−1) で、
対になる生成元の因子は 1 + u_{n|m}、u_{n|m} = 𝐚²𝐪^{−2(n−m)}𝐭^{2m+1}。KN の小文字では u_{n|m} = −1 ⇔ a = t^n q^{−m}（U(n|m)）。
"""
from flint import fmpz_mpoly_ctx

from .bps import g_expand

CTX = fmpz_mpoly_ctx.get(("a", "q", "t"), "lex")


def u(n, m):
    """d_{n|m} の対の単項式 u_{n|m} の指数。"""
    return (2, -2 * (n - m), 2 * m + 1)


def sigma(b):
    """DGR の対称性 (e_𝐚, e_𝐪, e_𝐭) ↦ (e_𝐚, −e_𝐪, e_𝐭 − e_𝐪)。d_{n|m} と d_{m|n} を入れ替える。"""
    return {(ea, -eq, et - eq): c for (ea, eq, et), c in b.items()}


def _poly(b):
    m = [min(k[i] for k in b) for i in range(3)]
    return CTX.from_dict({tuple(k[i] - m[i] for i in range(3)): int(v) for k, v in b.items()}), m


def _dict(P, m):
    return {tuple(e[i] + m[i] for i in range(3)): int(c) for e, c in P.to_dict().items()}


def factor_poly(pairs):
    """Π_{(n,m)} (1 + u_{n|m}) を多項式（指数は非負）で返す。u の 𝐪 の指数が負なら全体に 𝐪 の冪を掛けた形で、そのずれも返す。"""
    P, shift = CTX.from_dict({(0, 0, 0): 1}), 0
    for n, m in pairs:
        ea, eq, et = u(n, m)
        lo = min(0, eq)
        P *= CTX.from_dict({(0, -lo, 0): 1, (ea, eq - lo, et): 1})
        shift += lo
    return P, shift


def divides(b, n, m):
    """1 + u_{n|m} が b を割り切るか。"""
    P, _ = _poly(b)
    F, _ = factor_poly([(n, m)])
    return divmod(P, F)[1] == 0


D_PAIRS = ((0, 0), (1, 0), (0, 1))   # |ρ| ≥ 2 で \hat f_ρ を割り切る因子（framing 0）


def quotient(b, size):
    """\\hat f_ρ → R_ρ。|ρ| = 1 は (1 + u_{0|0}) で、|ρ| ≥ 2 は D = Π_{D_PAIRS} (1 + u) で割る。割り切れなければ ValueError。"""
    P, m = _poly(b)
    F, shift = factor_poly(D_PAIRS[:1] if size == 1 else D_PAIRS)
    R, rem = divmod(P, F)
    if rem != 0:
        raise ValueError("not divisible")
    m[1] -= shift
    return _dict(R, m)


def specialize(b, n, m):
    """1 + u_{n|m} = 0、すなわち 𝐚² = −𝐪^{2(n−m)}𝐭^{−2m−1} を代入した {(e_𝐪, e_𝐭): c}。"""
    out = {}
    for (ea, eq, et), c in b.items():
        be = int(ea) // 2
        k = (int(eq) + 2 * (n - m) * be, int(et) - (2 * m + 1) * be)
        out[k] = out.get(k, 0) + (-1) ** (be % 2) * int(c)
    return {k: v for k, v in out.items() if v}


def negatives(values):
    """全体の符号を良い方に選んだときの負の項の数。"""
    v = list(values)
    k = sum(1 for x in v if x < 0)
    return min(k, len(v) - k)


def judge_quotient(b, size):
    """(R の spin 表示の負の数, 項数, genus 表示の負の数, 項数, 同じ符号か)。"""
    R = quotient(b, size)
    N = g_expand(R)
    sR = sum(1 for x in R.values() if x < 0) <= len(R) / 2
    sN = sum(1 for x in N.values() if x < 0) <= len(N) / 2
    return negatives(R.values()), len(R), negatives(N.values()), len(N), sR == sN


def d11_top_genus(R):
    """genus 表示の R を d_{1|1}（𝐚² = −𝐭^{−3}、𝐪 を含まないので genus ごとに定義される）で特殊化し、0 でない最大の genus。"""
    acc = {}
    for (g, b, F), c in g_expand(R).items():
        d = acc.setdefault(int(g), {})
        k = int(F) - 3 * int(b)
        d[k] = d.get(k, 0) + (-1) ** (int(b) % 2) * int(c)
    gs = [g for g, d in acc.items() if any(d.values())]
    return max(gs) if gs else None


def hook_shift(lam):
    """hook λ の (|λ|, arm)。hook でなければ None。"""
    if len(lam) > 1 and lam[1] > 1:
        return None
    return sum(lam), lam[0] - 1


def hook_rule(A1, Alam, lam):
    """rCS_λ|_{d11} = c · A1(𝐪^{|λ|}𝐭^{arm}, 𝐭) となる単項式 c（{(e_𝐪, e_𝐭)}: 係数 の 1 項）を返す。成り立たなければ None。"""
    n, arm = hook_shift(lam)
    sub = {(n * eq, et + arm * eq): c for (eq, et), c in A1.items()}
    ks, ka = min(sub), min(Alam)
    if len(sub) != len(Alam) or Alam[ka] % sub[ks]:
        return None
    s, c = (ka[0] - ks[0], ka[1] - ks[1]), Alam[ka] // sub[ks]
    ok = all(Alam.get((k[0] + s[0], k[1] + s[1])) == c * v for k, v in sub.items())
    return (s, c) if ok else None
