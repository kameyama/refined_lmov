"""open_problems D1 (I)：R_{m,n}(q,t) = (−t/q)^{m−1} C_{m,n}(1/q,t) を GN の規約から段ごとに確かめる。

記号：X_{n,s}(u; q, τ) = Σ_μ x_μ Π_{□∈μ}(1 − u χ(□))、x_μ は GN 式 (fix3) を λ = ∅ に使い固定点基底の係数と読んだもの
（GN 式 (sumsyt) そのもの。τ は GN の t）。GN の Conjecture catalan は X_{n,m}(0; q, τ) = C_{n,m}(q, τ)（m > 0）。
本研究の値 U(a) = G.torus((1,), m, −n)（gn.py、t_w = 1/t）、R = reduced 値の a⁰ の係数。

段（有理数の点で厳密に比較）：
  S1  reduced の a⁰：R = (1 − t)·U(0)                         （unknot = (1 − a)/(1 − t)）
  S2  gn.py の読み方：U(a) = X_{m,−n−m}(a; q, 1/t)/(1 − t)       （∇ による読み替えで傾き −n → −n−m、c_1 = 1 − t）
  S3  鏡像：X_{m,−n−m}(0; q, τ) = (−qτ)^{−(m−1)} X_{m,n}(0; 1/q, 1/τ)
  S4  GN の Conjecture catalan：X_{m,n}(0; q, τ) = C_{m,n}(q, τ)
  全体 R = (−t/q)^{m−1} C_{m,n}(1/q, t)
使い方：PYTHONPATH=src:scripts/framed_knot:scripts/dgr_lmov python scripts/framed_knot/d1_gn_chain.py（ログ logs/framed_knot/d1_gn_chain.log）。
"""
from fractions import Fraction as Fr
from math import gcd

from rcsknots.ehall import EHA, boxes_of
from rcsknots.gn import GN
from d1_catalan import catalan, q as Q, t as T

PTS = [(Fr(3, 7), Fr(5, 11), Fr(2, 13)), (Fr(2, 9), Fr(7, 5), Fr(-3, 4))]


def X(n, s, u, q, tau):
    w = EHA(q, tau, Fr(1)).apply(1, n, s, {(): Fr(1)})
    out = Fr(0)
    for mu, c in w.items():
        pr = Fr(1)
        for r, col in boxes_of(mu):
            pr *= 1 - u * q ** col * tau ** r
        out += c * pr
    return out


def main():
    log = open("logs/framed_knot/d1_gn_chain.log", "w")
    cases = [(m, n) for m in range(2, 6) for n in range(2, 6) if gcd(m, n) == 1]
    for m, n in cases:
        C = catalan(m, n)
        ok = {k: True for k in ("S1", "S2", "S3", "S4", "all")}
        for q, t, a in PTS:
            G = GN(q, t, a, lambda x: x)
            Ua = G.torus((1,), m, -n)
            U0 = GN(q, t, Fr(0), lambda x: x).torus((1,), m, -n)
            R = G.torus_reduced_apoly((1,), m, -n)[0]
            ok["S1"] &= R == (1 - t) * U0
            ok["S2"] &= Ua == X(m, -n - m, a, q, 1 / t) / (1 - t)
            for tau in (t, 1 / t):
                ok["S3"] &= X(m, -n - m, 0, q, tau) == (-q * tau) ** (1 - m) * X(m, n, 0, 1 / q, 1 / tau)
                ok["S4"] &= X(m, n, 0, q, tau) == C.subs({Q: q, T: tau})
            ok["all"] &= R == (-t / q) ** (m - 1) * C.subs({Q: 1 / q, T: t})
        line = f"T({m},{n}): " + "  ".join(f"{k} {v}" for k, v in ok.items())
        print(line, flush=True)
        log.write(line + "\n")


if __name__ == "__main__":
    main()
