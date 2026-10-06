"""open_problems D1：トーラス結び目の正規化の単項式 M_□ を有理 q,t-Catalan 数から導くための照合。

GN の正規化前の reduced 値 ε(P_{(m,−n)}·1)/(unknot) の a⁰ の係数 R_{m,n}(q,t) が
(−t/q)^{m−1} C_{m,n}(1/q, t) に等しいかを調べる（C は GN 1304.3328 §7 の Dyck 路の和 Σ_D q^{δ−|D|} t^{h₊(D)}）。
あわせて (N1)+(N2) の単項式が b = m（傾きの第 1 成分）の式 M_□ = (−q/t)^{m−1}(a(q/t)^{1/2})^δ に一致するかを、m > n の向きも含めて調べる。
使い方：python scripts/d1_catalan.py。結果は logs/d1_catalan.log。
"""
import sys
from fractions import Fraction
from math import factorial

import sympy as sp

from rcsknots.exact import reconstruct
from rcsknots.homfly import torus_reduced
from rcsknots.normalize import normalize

q, t = sp.symbols("q t")


def dyck_paths(m, n):
    """m×n の長方形の対角線の下に入るヤング図（行の長さの列、下の行から）。"""
    def rows(j, prev):
        if j == n:
            yield ()
            return
        # 行 j（0 始まり）の箱 (i, j) は (i+1)/m + (j+1)/n ≤ 1 を満たす
        top = min(prev, max(0, (m * (n - j - 1)) // n))
        for r in range(top, -1, -1):
            for rest in rows(j + 1, r):
                yield (r,) + rest
    return [tuple(x for x in D if x) for D in rows(0, m)]


def catalan(m, n):
    delta = (m - 1) * (n - 1) // 2
    out = 0
    for D in dyck_paths(m, n):
        cols = [sum(1 for r in D if r > i) for i in range(D[0] if D else 0)]
        hplus = 0
        for j, r in enumerate(D):
            for i in range(r):
                arm, leg = r - i - 1, cols[i] - j - 1
                if Fraction(arm, leg + 1) < Fraction(m, n) < (Fraction(arm + 1, leg) if leg else Fraction(10**9)):
                    hplus += 1
        out += q ** (delta - sum(D)) * t ** hplus
    assert len(dyck_paths(m, n)) == factorial(m + n - 1) // (factorial(m) * factorial(n))
    return sp.expand(out)


def main():
    log = open("logs/d1_catalan.log", "w")

    def out(s):
        print(s, flush=True)
        log.write(s + "\n")

    for m, n in [(2, 3), (3, 2), (2, 5), (5, 2), (3, 4), (4, 3), (3, 5), (5, 3), (2, 7), (4, 5), (5, 4)]:
        delta = (m - 1) * (n - 1) // 2
        C = catalan(m, n)
        sym = sp.expand(C - C.subs({q: t, t: q}, simultaneous=True)) == 0
        raw = reconstruct(lambda G: G.torus_reduced_apoly((1,), m, -n), log=lambda s: None)
        R = sp.expand(sum(c * q ** eq * t ** et for (ea, eq, et), c in raw.items() if ea == 0))
        pred = sp.expand((-t / q) ** (m - 1) * C.subs(q, 1 / q))
        _, mono = normalize(raw, torus_reduced((1,), min(m, n), max(m, n), sign=-1))
        s, al, be, ga = mono
        mono_ok = (s, al, be, ga) == ((-1) ** (m - 1), Fraction(delta), Fraction(delta, 2) + m - 1, -Fraction(delta, 2) - m + 1)
        out(f"T({m},{n})：C の q,t 対称性 {sym}、a⁰ の係数 = (−t/q)^(m−1) C(1/q,t) {sp.expand(R - pred) == 0}、"
            f"単項式 {mono} が b = m の式に一致 {mono_ok}")


if __name__ == "__main__":
    sys.exit(main())
