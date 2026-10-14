r"""K4'(3)：交差帽部分を Σ_{d 奇数} の構造で定義する。
  F = ½F^{cov} + Σ_{d 奇数} (1/d) Σ_{R,S} g_{R,S}(x^d) s_R(y^d) s_S(z^d)
  ⇒ g_{ab} = L_{ab} − ½L^{cov}_{ab} − Σ_{d 奇数 > 1, d|(a,b)} (1/d) g_{a/d,b/d}(x^d)。
判定：(i) q, t, ν の Laurent 多項式か、(ii) ν~ = ν (q/t)^{1/2} を固定すると t に依らないか、(iii) z_q = q^{1/2} − q^{−1/2} の冪で展開したときの係数（ν~, q^{1/2} の Laurent 多項式）の整数性と符号。
出力：data/refined_kauffman/hopf_g_odd.pkl。"""
import sys, pickle, time
sys.path.insert(0, "src")
import sympy as sp
from sympy import mobius
from rcsknots.kauffman_cover import cover_Z_sym, log_coeffs_sym, qs_, ss_, nus_

A_ = int(sys.argv[1]) if len(sys.argv) > 1 else 4
r, s, nu = sp.symbols("r s nu", positive=True)      # q = r², t = s², ν
nut, zq = sp.symbols("nut zq")                       # ν~、z_q
t0 = time.time()
Lc = {ab: v.as_expr().subs({qs_: r**2, ss_: s, nus_: nu}) for ab, v in log_coeffs_sym(cover_Z_sym(A_, A_), A_, A_).items()}
print(f"L^cov {time.time()-t0:.0f}s", flush=True)
q_, t_, U_ = sp.symbols("q t U")
f = {ab: sp.sympify(v, locals={"q": q_, "t": t_, "U": U_}).subs({q_: r**2, t_: s**2, U_: nu * s}) for ab, v in pickle.load(open("data/refined_kauffman/hopf_f.pkl", "rb")).items()}
def adams(e, d): return e.subs({r: r**d, s: s**d, nu: nu**d}, simultaneous=True)
L = {}
for (a, b) in f:
    g_ = sp.gcd(a, b) if a and b else (a or b)
    L[(a, b)] = sum(sp.Rational(1, d) * adams(f[(a // d, b // d)], d) for d in sp.divisors(int(g_)))
g = {}
for a in range(0, A_ + 1):
    for b in range(0, A_ + 1):
        if a + b == 0: continue
        e = L[(a, b)] - Lc[(a, b)] / 2
        g_ = sp.gcd(a, b) if a and b else (a or b)
        for d in sp.divisors(int(g_)):
            if d > 1 and d % 2 == 1:
                e -= sp.Rational(1, d) * adams(g[(a // d, b // d)], d)
        g[(a, b)] = sp.cancel(sp.together(e))
pickle.dump({ab: str(v) for ab, v in g.items()}, open("data/refined_kauffman/hopf_g_odd.pkl", "wb"))
for ab in sorted(g, key=lambda k: (k[0] + k[1], k)):
    v = g[ab]
    if v == 0:
        print(f"({ab[0]},{ab[1]}): g = 0"); continue
    num, den = sp.fraction(v)
    laurent = sp.Poly(den, r, s, nu).is_monomial
    # ν = ν~ s / r を代入して t（= s）依存性を見る
    w = sp.cancel(v.subs(nu, nut * s / r))
    t_free = not w.has(s)
    line = f"({ab[0]},{ab[1]}): Laurent={laurent}  ν~ 固定で t に依らない={t_free}"
    if t_free:
        # z_q の冪で展開：r = q^{1/2}、z_q = r − 1/r → r の Laurent 多項式を z_q と r の「偶奇」で分ける
        # 方法：w(r, ν~) を z_q の多項式に：r + 1/r = sqrt(z_q² + 4) は使わず、w を r^m で展開し係数を見る
        P = sp.Poly(sp.expand(w * (r * nut)**20), r, nut)
        coeffs = P.coeffs()
        line += f"  r,ν~ の係数: 整数={all(c.is_integer for c in coeffs)} 正{sum(1 for c in coeffs if c > 0)} 負{sum(1 for c in coeffs if c < 0)}"
        # z_q で割れる回数
        m = 0; ww = sp.expand(w)
        while True:
            quo = sp.cancel(ww / (r - 1 / r)); n_, d_ = sp.fraction(quo)
            if sp.Poly(d_, r, nut).is_monomial and quo != 0:
                ww = sp.expand(quo); m += 1
            else:
                break
        line += f"  z_q で割れる回数: {m}"
        print(line); print("     g =", sp.factor(w))
    else:
        print(line); print("     g =", sp.factor(v))
