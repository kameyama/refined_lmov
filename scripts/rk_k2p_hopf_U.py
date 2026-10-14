r"""K2'(a)：双対 Cauchy 型の Hopf 絡み目（probe 各 1 個）の log Z~ を U = t^N の形式変数で計算し、
Möbius 反転した f_{ab} の分母・分子の整数性・構造を出す。A = B = 4 で約 1 分。
出力：data/refined_kauffman/hopf_f.pkl（{(a,b): f_{ab} の sympy 式}）。
"""
import sys, pickle, time
from pathlib import Path
sys.path.insert(0, "src")
import sympy as sp
from rcsknots.kauffman_hopf import hopf_log_coeffs, mobius_f, qs as q, ts as t, Us as U

A = B = int(sys.argv[1]) if len(sys.argv) > 1 else 4
t0 = time.time()
L = hopf_log_coeffs(A, B)
print(f"L_ab 計算 {time.time()-t0:.0f}s")
f = {ab: mobius_f(L, *ab) for ab in L if ab != (0, 0)}
out = Path("data/refined_kauffman"); out.mkdir(parents=True, exist_ok=True)
pickle.dump({ab: str(v) for ab, v in f.items()}, open(out / "hopf_f.pkl", "wb"))

f10 = f[(1, 0)]
zq = (q - 1) / sp.sqrt(q)
nut = U * sp.sqrt(q) / t
print("f_10 =", sp.factor(f10), "  f_a0 (a≥2) =", [f[(a, 0)] for a in range(2, A + 1)])
print("f_11 + f_10 z_q(ν~ − 1/ν~) =", sp.simplify(f[(1, 1)] + f10 * zq * (nut - 1 / nut)), "  (ν~ = U q^{1/2}/t)")
chi = lambda m, x: sum(x**(m - 2*i) for i in range(m + 1))
for b in range(2, B + 1):
    print(f"f_1{b}/f_11 − χ_{b-1}(U/t) =", sp.cancel(f[(1, b)] / f[(1, 1)] - chi(b - 1, U / t)))
print("== f_ab の分母と分子の整数性")
for ab in sorted(f):
    if ab[0] >= 1 and ab[1] >= 1:
        num, den = sp.fraction(sp.cancel(f[ab]))
        P = sp.Poly(sp.expand(num), q, t, U)
        print(ab, " 分母 =", sp.factor(den), " 分子の係数が整数:", all(c.is_integer for c in P.coeffs()),
              " (q−1) で割り切れる:", sp.rem(P, sp.Poly(q - 1, q, t, U)).is_zero)
print("== t = q での f_ab（Laurent 多項式か）")
for ab in sorted(f):
    if ab[0] >= 1 and ab[1] >= 1 and ab[0] + ab[1] <= 5:
        num, den = sp.fraction(sp.cancel(f[ab].subs(t, q)))
        print(ab, " 分母 =", sp.factor(den))
