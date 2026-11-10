"""第 4 期 V4：\\hat f_ρ(K_p) = Σ_j A_j m_j^p（p > 0）の分解を求め、各 A_j を調べる。

m_j = 1, X, X², Y, Ȳ（太字：X = 𝐚²𝐭², Y = 𝐚⁴𝐪⁴𝐭⁸, Ȳ = 𝐚⁴𝐪⁻⁴𝐭⁴）。p = 1..5 で A_j を解き、p = 6..12 で確かめる。
入力の \\hat f は scripts/framed_knot/v4_twist_structure.py の hats（data/cache に無ければ計算する）。
"""
import pickle

import sympy as sp

from v4_twist_structure import hats

a, q, t = sp.symbols("a q t")
M = [sp.Integer(1), a**2 * t**2, a**4 * t**4, a**4 * q**4 * t**8, a**4 * q**-4 * t**4]   # 特性根（scripts/framed_knot/v4_roots.py）
F = hats(1)
for rho in [(2,), (1, 1)]:
    Asym = sp.symbols("A0:5")
    eqs = [sp.Eq(sp.expand(F[p][rho]), sum(Asym[j] * M[j] ** p for j in range(5))) for p in range(1, 6)]
    sol = sp.solve(eqs, Asym, dict=True)[0]
    A = [sp.factor(sp.cancel(sol[x])) for x in Asym]
    pts = [{a: sp.Rational(2, 7), q: sp.Rational(3, 11), t: sp.Rational(5, 13)}, {a: 17, q: sp.Rational(19, 3), t: 23}]
    ok = all((F[p][rho] - sum(A[j] * M[j] ** p for j in range(5))).subs(v) == 0 for p in range(6, 13) if p in F for v in pts)
    print(f"ρ={rho}: p = 6..12 で成立 {ok}", flush=True)
    for j in range(5):
        num, den = sp.fraction(A[j])
        print(f"  m_{j} = {M[j]}: 分母 {sp.factor(den)}", flush=True)
    pickle.dump(A, open(f"data/cache/v4_A_{''.join(map(str, rho))}.pkl", "wb"))
