"""open_problems D1 (I)：BGLX（arXiv:1404.4616）の演算子 Q_{m,n} と GN（1304.3328）の式 (fix3)/(adriano) の対応を数値で確かめる。

BGLX：D_k f[X] = f[X + M/z] Σ_i (−z)^i e_i[X] |_{z^k}（M = (1−q)(1−t)、式 (defD)）、Q_{1,n} = D_n、
Q_{m,n} = (1/M)[Q_{c,d}, Q_{a,b}]（split(m,n) = (a,b)+(c,d)、Def. defnQ、Prop. 3.3）。Q_{m,n} の次数は n。
BGLX の有理 Catalan の定義 (5.26)：C_{m,n} = ⟨Q_{m,n}·(−1)^n, e_n⟩_*、⟨p_λ,p_μ⟩_* = (−1)^{|μ|−ℓ(μ)} Π(1−t^{μ_i})(1−q^{μ_i}) z_μ δ。
GN：P̃_{n,m}·1 = Σ_λ c_{n,m}(λ) H̃_λ/g_λ（式 (adriano)、(fix3) を固定点基底 I_λ = H̃_λ の係数と読む）。

照合（有理数の点で厳密）：
  B1  Q_{m,n}·(−1)^n = GN の P̃_{n,m}·1（(adriano) をそのまま H̃_λ(q,t) で展開したもの）
  B2  ⟨Q_{m,n}·(−1)^n, e_n⟩_* と Hall の ⟨·, h_n⟩ = f[1] と ⟨·, e_n⟩ のうち、どれが Dyck 路の和 C_{m,n}(q,t) に等しいか
使い方：PYTHONPATH=src:scripts python scripts/d1_bglx.py（ログ logs/d1_bglx.log）。
"""
from fractions import Fraction as Fr
from math import factorial, gcd
from collections import Counter

from rcsknots.ehall import EHA
from rcsknots.modified import modified_H
from rcsknots.partitions import partitions
from d1_catalan import catalan, q as Q, t as T


def z(rho):
    out = 1
    for k, c in Counter(rho).items():
        out *= k ** c * factorial(c)
    return out


def add(f, g, s=1):
    out = dict(f)
    for k, v in g.items():
        out[k] = out.get(k, 0) + s * v
    return {k: v for k, v in out.items() if v != 0}


def mul(f, g):
    out = {}
    for a, x in f.items():
        for b, y in g.items():
            k = tuple(sorted(a + b, reverse=True))
            out[k] = out.get(k, 0) + x * y
    return out


def e(n):
    return {rho: Fr((-1) ** (n - len(rho)), z(rho)) for rho in partitions(n)} if n else {(): Fr(1)}


def D(k, f, q, t):
    """BGLX (defD)。"""
    out = {}
    for rho, c in f.items():
        # f[X + M/z] の p_ρ の項：Π_i (p_{ρ_i} + M_{ρ_i} z^{−ρ_i}) を展開、{(σ, j): 係数}（z^{−j}）
        terms = {((), 0): Fr(1)}
        for r in rho:
            new = {}
            Mr = (1 - q ** r) * (1 - t ** r)
            for (sig, j), x in terms.items():
                k1 = (tuple(sorted(sig + (r,), reverse=True)), j)
                new[k1] = new.get(k1, 0) + x
                k2 = (sig, j + r)
                new[k2] = new.get(k2, 0) + x * Mr
            terms = new
        for (sig, j), x in terms.items():
            deg = k + j                       # z^{-j} · (−z)^{deg} e_deg が z^k を与える
            if deg < 0:
                continue
            g = mul({sig: x * c}, e(deg))
            out = add(out, {kk: (-1) ** deg * v for kk, v in g.items()})
    return out


def split(m, n):
    if m == 1:
        return (1, n - 1), (0, 1)
    if n == 1:
        return (1, 0), (m - 1, 1)
    for a in range(1, m):
        for b in range(1, n):
            if m * b + 1 == n * a:
                return (a, b), (m - a, n - b)
    raise ValueError


def Qop(m, n, f, q, t):
    if m == 1:
        return D(n, f, q, t)
    if (m, n) == (1, 0):
        raise ValueError
    (a, b), (c, d) = split(m, n)
    M = (1 - q) * (1 - t)
    x = Qop(c, d, Qop_any(a, b, f, q, t), q, t)
    y = Qop_any(a, b, Qop(c, d, f, q, t), q, t)
    return {k: v / M for k, v in add(x, y, -1).items()}


def Qop_any(a, b, f, q, t):
    if (a, b) == (1, 0):
        # Q_{1,0} = D_0
        return D(0, f, q, t)
    if (a, b) == (0, 1):
        raise ValueError("Q_{0,1} は使わない")
    return Qop(a, b, f, q, t)


def gn_vector(n, m, q, t):
    """GN (adriano)：Σ_λ coef_λ H̃_λ(q,t) を p 基底で。"""
    coef = EHA(q, t, Fr(1)).apply(1, n, m, {(): Fr(1)})
    H = modified_H(n, q, t, lambda x: x)
    out = {}
    for lam, c in coef.items():
        out = add(out, {r: c * v for r, v in H[lam].items()})
    return out


def star_e(f, n, q, t):
    s = 0
    for rho, c in f.items():
        w = 1
        for r in rho:
            w *= (1 - q ** r) * (1 - t ** r)
        s += c * w            # ⟨p_ρ, e_n⟩_* = Π(1−q^{ρ_i})(1−t^{ρ_i})
    return s


def hall_h(f):
    return sum(f.values())                     # ⟨p_ρ, h_n⟩ = 1


def hall_e(f, n):
    return sum(c * (-1) ** (n - len(rho)) for rho, c in f.items())


def main():
    log = open("logs/d1_bglx.log", "w")
    pts = [(Fr(3, 7), Fr(5, 11)), (Fr(2, 9), Fr(7, 5))]
    for m, n in [(3, 2), (2, 3), (5, 2), (2, 5), (4, 3), (3, 4), (5, 3), (3, 5), (5, 4), (4, 5)]:
        res = {"B1": True, "star_e": True, "hall_h": True, "hall_e": True}
        for q, t in pts:
            Qv = {k: (-1) ** n * v for k, v in Qop(m, n, {(): Fr(1)}, q, t).items()}
            G = gn_vector(n, m, q, t)
            res["B1"] &= add(Qv, G, -1) == {}
            C = catalan(m, n).subs({Q: q, T: t})
            res["star_e"] &= star_e(Qv, n, q, t) == C
            res["hall_h"] &= hall_h(Qv) == C
            res["hall_e"] &= hall_e(Qv, n) == C
        line = f"BGLX Q_({m},{n})(−1)^{n} vs GN P̃_({n},{m})·1: " + "  ".join(f"{k} {v}" for k, v in res.items())
        print(line, flush=True)
        log.write(line + "\n")


if __name__ == "__main__":
    main()
