"""open_problems D1：本研究の値 → BGLX の Q → Mellit の Dyck 路の式、の連鎖を段ごとに有理数の点で厳密に確かめる。

記号は scripts/framed_knot/d1_gn_chain.py と同じ（X_{n,s}(u; q, τ) は GN (sumsyt) を一般の u で）。V_{n,m}(q,t) = Q_{n,m}·1（BGLX、次数 m）。
  R1  u つきの鏡像：X_{m,s}(u; 1/q, 1/τ) = (−qτ)^{m−1} (−u)^m X_{m,−s−2m}(1/u; q, τ)
  R2  X_{m,n−m}(u; q, τ) = (−1)^m V_{n,m}(q,τ)[1−u]（(fix3) の ∇ の補題と Q_{m',n'} = (−1)^{n'} O_{n',m'}、scripts/framed_knot/d1_ops.py）
  R3  f[a−1] = a^d f[1−1/a]（次数 d の対称関数、フック展開から）
  R4  V_{n,m}(q,t) = V_{n,m}(t,q)（GN :1247 の対称性と H̃ の転置対称性）
  R5  まとめ：本研究の正規化前の reduced 値 R(a) = −(−q/t)^{1−m} V_{n,m}(1/q, t)[a−1]/(a−1)
  R6  Mellit の Cor (daha)（1704.07630 khro.tex:662）を P_{m,n}(1) = Q·1 と読んだときの照合：Thm (dyck paths) の式 Dyck_{n,m}(a,q,t) と
      −(a(qt)^{−1/2})^δ/(1−t) · ((−1)^m V_{n,m}(q,t))[a/q − 1]/(a/q − 1) の比
使い方：PYTHONPATH=src:scripts/framed_knot:scripts/dgr_lmov python scripts/framed_knot/d1_ops2.py（ログ logs/framed_knot/d1_ops2.log）。
"""
from fractions import Fraction as Fr

from rcsknots.gn import GN
from d1_bglx import Qop
from d1_gn_chain import X
from d1_dyck import terms


def pleth(f, x):
    """f[x − 1]（p_k ↦ x^k − 1）。"""
    s = Fr(0)
    for rho, c in f.items():
        w = c
        for k in rho:
            w *= x ** k - 1
        s += w
    return s


def pleth1(f, u):
    """f[1 − u]。"""
    s = Fr(0)
    for rho, c in f.items():
        w = c
        for k in rho:
            w *= 1 - u ** k
        s += w
    return s


def dyck(m, n, a, q, t, sqrt_qt):
    d = (m - 1) * (n - 1) // 2
    s = sum(t ** area * q ** (hp - sk) * (-a) ** S for area, hp, S, sk in terms(m, n))
    return (a / sqrt_qt) ** d * s / (1 - t)


def main():
    log = open("logs/framed_knot/d1_ops2.log", "w")

    def out(s):
        print(s, flush=True)
        log.write(s + "\n")

    q, t, a, u = Fr(4, 9), Fr(9, 25), Fr(2, 7), Fr(3, 11)
    sq = Fr(2, 5)                                    # √(qt)
    for m, n in [(2, 3), (3, 2), (2, 5), (3, 4), (4, 3), (3, 5), (5, 3)]:
        tau = 1 / t
        r1 = all(X(m, s, u, 1 / q, 1 / tau) == (-q * tau) ** (m - 1) * (-u) ** m * X(m, -s - 2 * m, 1 / u, q, tau)
                 for s in (n, n - m, -n))
        V = Qop(n, m, {(): Fr(1)}, q, tau)
        r2 = X(m, n - m, u, q, tau) == (-1) ** m * pleth1(V, u)
        r3 = pleth(V, a) == a ** m * pleth1(V, 1 / a)
        r4 = pleth(V, a) == pleth(Qop(n, m, {(): Fr(1)}, tau, q), a)
        G = GN(q, t, a, lambda x: x)
        R = sum(c * a ** i for i, c in enumerate(G.torus_reduced_apoly((1,), m, -n)))
        Vb = Qop(n, m, {(): Fr(1)}, 1 / q, t)
        r5 = R == -(-q / t) ** (1 - m) * pleth(Vb, a) / (a - 1)
        d = (m - 1) * (n - 1) // 2
        Vm = {k: (-1) ** m * v for k, v in Qop(n, m, {(): Fr(1)}, q, t).items()}
        cor = -(a / sq) ** d / (1 - t) * pleth(Vm, a / q) / (a / q - 1)
        r6 = dyck(n, m, a, q, t, sq) / cor
        out(f"T({m},{n})：R1 {r1}、R2 {r2}、R3 {r3}、R4 {r4}、R5 {r5}、R6 Dyck/Cor = {r6}")


if __name__ == "__main__":
    main()
