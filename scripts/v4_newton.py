"""第 4 期 V4：p の符号ごとに、\\hat f(K_p) の漸化式の特性根を求め、ニュートン形の係数 c_k の正値性を全ての根の並びで調べる。

使い方：python scripts/v4_newton.py roots SIGN      （SIGN = 1 / -1、|p| = 1..12 を直接計算して特性根を求める）
        python scripts/v4_newton.py search SIGN     （ニュートン形の c_k が一様な符号で非負になる (n0, 並び) を探す）
"""
import sys
from fractions import Fraction
from itertools import permutations
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from v4_closed import M_ROOTS, RHOS, hat_direct, load_inputs, newton_coeffs  # noqa: E402


def ev(P, pt):
    return sum(Fraction(int(c)) * Fraction(pt[0]) ** i * Fraction(pt[1]) ** j * Fraction(pt[2]) ** k for (i, j, k), c in P.items())


def roots(sign):
    cache = load_inputs()
    S = {n: hat_direct(sign * n, cache) for n in range(1, 13)}
    x = sp.symbols("x")
    for rho in RHOS:
        for pt in [(2, 3, 5), (7, 11, 13)]:
            s = [sp.Rational(ev(S[n][rho], pt)) for n in range(1, 13)]
            for order in range(1, 7):
                H = sp.Matrix([[s[i + j] for j in range(order + 1)] for i in range(order + 1)])
                if H.det() == 0:
                    break
            n = order
            Mt = sp.Matrix([[s[i + j] for j in range(n)] for i in range(n)])
            c = Mt.solve(sp.Matrix([s[i + n] for i in range(n)]))
            ok = all(sum(c[j] * s[i + j] for j in range(n)) == s[i + n] for i in range(len(s) - n))
            char = x ** n - sum(c[j] * x ** j for j in range(n))
            out = []
            for r, mult in sp.roots(sp.Poly(char, x)).items():
                f = sp.factorint(sp.Rational(r).p) | {k: -v for k, v in sp.factorint(sp.Rational(r).q).items()}
                out.append(((sp.sign(r),) + tuple(f.get(k, 0) for k in pt), mult))
            print(f"符号 {sign} ρ={rho} 点 {pt}: 次数 {n}（|p| = 1..12 で成立 {ok}）、根（符号, 𝐚, 𝐪, 𝐭 の指数）{out}", flush=True)


def search(sign, roots_exp, n0s=(1, 2, 3)):
    cache = load_inputs()
    S = {n: hat_direct(sign * n, cache) for n in range(1, max(n0s) + 5)}
    names = list(roots_exp)
    found = {}
    for rho in RHOS:
        for n0 in n0s:
            good = []
            for perm in permutations(names):
                order = [roots_exp[k] for k in perm]
                cs = newton_coeffs([S[n0 + i][rho] for i in range(5)], order)
                for sg in (1, -1):
                    if all(sg * v >= 0 for c in cs for v in c.values()):
                        good.append((perm, sg))
            print(f"符号 {sign} ρ={rho} n0={n0}: c_k が一様な符号で非負になる並び {len(good)}/120 {good[:6]}", flush=True)
            found[(rho, n0)] = good
    return found


if __name__ == "__main__":
    mode, sign = sys.argv[1], int(sys.argv[2])
    if mode == "roots":
        roots(sign)
    else:
        R = M_ROOTS if sign == 1 else {k: tuple(-e for e in v) for k, v in M_ROOTS.items()}
        search(sign, R)
