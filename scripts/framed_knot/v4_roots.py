"""第 4 期 V4：\\hat f_ρ(K_p)（p > 0）の p についての漸化式の特性根（＝ m_j）をデータから求める。

(a, q, t) を素数の点 (2, 3, 5)、(7, 11, 13) で評価し、次数 5 の漸化式の係数を厳密に解き、特性多項式の有理根を因数分解して単項式の指数を読む。
"""
import sympy as sp

from v4_twist_structure import hats

a, q, t = sp.symbols("a q t")
F = hats(1)
x = sp.symbols("x")
for rho in [(2,), (1, 1)]:
    for pt in [(2, 3, 5), (7, 11, 13)]:
        s = [sp.Rational(F[p][rho].subs({a: pt[0], q: pt[1], t: pt[2]})) for p in range(1, 13)]
        n = 5
        Mt = sp.Matrix([[s[i + j] for j in range(n)] for i in range(n)])
        rhs = sp.Matrix([s[i + n] for i in range(n)])
        c = Mt.solve(rhs)
        assert all(sum(c[j] * s[i + j] for j in range(n)) == s[i + n] for i in range(len(s) - n)), "recurrence fails"
        char = x**n - sum(c[j] * x**j for j in range(n))
        roots = sp.roots(sp.Poly(char, x))
        out = []
        for r, mult in roots.items():
            f = sp.factorint(sp.Rational(r).p) | {k: -v for k, v in sp.factorint(sp.Rational(r).q).items()}
            out.append((r, {k: f.get(k, 0) for k in pt}, mult))
        print(f"ρ={rho} 点 {pt}: 根 {[(o[1], o[2]) for o in out]}（全ての根が有理数：{sum(o[2] for o in out) == n}）", flush=True)
