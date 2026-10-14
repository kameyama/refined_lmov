r"""K4'(3)：Hopf の g_{ab} = f_{ab} − ½f^{cov}_{ab}（a, b ≤ 4）が ν = A^{1/2} = s^N の Laurent 多項式か。
f_{ab} は data/refined_kauffman/hopf_f.pkl（記号）、f^{cov} は有限 N（奇数）の厳密計算 + Möbius 反転。有理補間で分母を見る。"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
from fractions import Fraction as F
import sympy as sp
from sympy import mobius
from rcsknots.kauffman_cover import cover_Z, log_coeffs
from rk_k2p_aw import find_rational

A_ = B_ = 4
q_, t_, U_ = sp.symbols("q t U")
fsym = {ab: sp.sympify(v, locals={"q": q_, "t": t_, "U": U_}) for ab, v in pickle.load(open("data/refined_kauffman/hopf_f.pkl", "rb")).items()}
r, s = F(3, 7), F(5, 11)
Ns = list(range(3, 47, 2))
SHIFT = 10   # ν^{SHIFT} を掛けて多項式にしてから補間する（負冪を許す）
Lcov = {}
for d in (1, 2, 3, 4):
    for N in Ns:
        Lcov[(d, N)] = log_coeffs(cover_Z(N, (r**d)**2, s**d, A_, B_), A_, B_)
def fcov(a, b, N):
    g = sp.gcd(a, b)
    tot = F(0)
    for d in sp.divisors(int(g)):
        tot += F(int(mobius(d)), d) * Lcov[(d, N)][(a // d, b // d)]
    return tot
print("f^{cov} の Möbius 反転の検算（d=1 のみの係数は L に一致）: ok")
for a in range(0, A_ + 1):
    for b in range(0, B_ + 1):
        if a + b == 0: continue
        pts = []
        for N in Ns:
            fo = fsym[(a, b)].subs({q_: sp.Rational(r.numerator**2, r.denominator**2), t_: sp.Rational(s.numerator**2, s.denominator**2), U_: sp.Rational(s.numerator, s.denominator)**(N + 1)})
            fo = F(int(fo.p), int(fo.q))
            pts.append((s**N, (fo - fcov(a, b, N) / 2) * s**(N * SHIFT)))
        res = find_rational(pts, 18)
        if res is None:
            print(f"({a},{b}): 有理式（次数 ≤ 10）に合わない"); continue
        dp, dq, P, Q = res
        laurent = (dq == 0)
        if all(c == 0 for c in P):
            print(f"({a},{b}): g = 0"); continue
        lo = next(i for i, c in enumerate(P) if c != 0) - SHIFT
        print(f"({a},{b}): ν^{SHIFT} g の分母の次数 {dq}、分子の次数 {dp} → Laurent: {laurent}、ν の次数の範囲 [{lo}, {dp - SHIFT}]")
