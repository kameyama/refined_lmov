r"""K4'(3)：refined の偶数 d の引き算【仮説：AS 1202.4456 §7 の閉弦の構造 (1+q^d)/(1+t^d) を開弦に適用】
  F^{unor} = Σ_d (1/d) g(x^d) B(y^d,z^d) − Σ_d (1/2d) (1+q^d)/(1+t^d) g(x^{2d}) B(y^{2d},z^{2d})
のもとで、(2,2), (2,4), (4,2), (4,4) の g を g^{odd}（Σ_{d 奇数} の定義、hopf_g_odd.pkl）から求め、Laurent 性・整数性・ν~ 構造を見る。
  g_22 = g^odd_22 − ½ (t−q)/(1+t) g_11(x²)
  g_24 = g^odd_24 − ½ (t−q)/(1+t) g_12(x²)（(4,2) も同様）
  g_44 = g^odd_44 − ½ (t−q)/(1+t) g_22(x²) − ¼ (t²−q²)/(1+t²) g_11(x⁴)
（Schur の冪和プレシズム s_R(y²) の s_{2R} 成分の係数は +1。）"""
import sys, pickle
sys.path.insert(0, "src")
import sympy as sp
r, s, nu = sp.symbols("r s nu", positive=True)
nut = sp.Symbol("nut")
god = {ab: sp.sympify(v, locals={"r": r, "s": s, "nu": nu}) for ab, v in pickle.load(open("data/refined_kauffman/hopf_g_odd.pkl", "rb")).items()}
q, t = r**2, s**2
def ad(e, d): return e.subs({r: r**d, s: s**d, nu: nu**d}, simultaneous=True)
g = dict(god)
g[(2, 2)] = sp.cancel(god[(2, 2)] - sp.Rational(1, 2) * (t - q) / (1 + t) * ad(god[(1, 1)], 2))
for ab in ((2, 4), (4, 2)):
    g[ab] = sp.cancel(god[ab] - sp.Rational(1, 2) * (t - q) / (1 + t) * ad(god[(ab[0] // 2, ab[1] // 2)], 2))
g[(4, 4)] = sp.cancel(god[(4, 4)] - sp.Rational(1, 2) * (t - q) / (1 + t) * ad(g[(2, 2)], 2) - sp.Rational(1, 4) * (t**2 - q**2) / (1 + t**2) * ad(god[(1, 1)], 4))
for ab in ((2, 2), (2, 4), (4, 2), (4, 4)):
    v = g[ab]
    num, den = sp.fraction(v)
    laurent = sp.Poly(den, r, s, nu).is_monomial
    P = sp.Poly(sp.expand(num), r, s, nu)
    w = sp.cancel(v.subs(nu, nut * s / r))
    print(f"({ab[0]},{ab[1]}): Laurent={laurent} 分母={sp.factor(den)} 分子の係数が整数={all(c.is_integer for c in P.coeffs())} 項数={len(P.terms())} ν~固定で t 非依存={not w.has(s)}")
    if laurent:
        print("    g =", sp.factor(v) if len(P.terms()) < 40 else "(長い)")
pickle.dump({ab: str(v) for ab, v in g.items()}, open("data/refined_kauffman/hopf_g_ref.pkl", "wb"))
