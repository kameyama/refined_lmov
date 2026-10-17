r"""K4'(a)：交差帽部分を全 d（第二量子化、AS 1202.4456 §3・§7）の構造で定義する。
  F = Σ_{d≥1} (1/d) Σ_{a,b} f_{ab}(x^d) y^{ad} z^{bd}（f_{ab} は全 d の Möbius 反転、hopf_f.pkl）
  g_{ab} := f_{ab} − ½ f^{cov}_{ab}（全 d）。Σ_{d 奇数} の定義（hopf_g_odd.pkl）との関係は
  g_{ab} = g^{odd}_{ab} − ½ g^{odd}_{a/2,b/2}(x²)（a, b がともに偶数のときだけ第 2 項がある）。
予測（AS の M 理論の軌道体：probe（t-brane）に沿う z₂ → −z₂ の射影）：
  g_{ab} = ½ E_{ab}/(1+t) + C_{ab}、E：Ĩ で固定される oriented 状態の符号付き数（(偶,偶) だけ）、C：z₂ = 0 に凍結した unoriented 状態。
  したがって 2(1+t) g_{ab} は整数 Laurent、(偶,偶) 以外は g 自身が整数 Laurent、(4,4) に (1+t²) は出ない。
分離に依らない量 E_{ab}(t = −1) = [2(1+t) g_{ab}]_{t=−1} を ν~ = ν (q/t)^{1/2} で表示する。
出力：data/refined_kauffman/hopf_g_alld.pkl。"""
import pickle
import sympy as sp

r, s, nu = sp.symbols("r s nu", positive=True)
nut = sp.Symbol("nut")
t = s**2
god = {ab: sp.sympify(v, locals={"r": r, "s": s, "nu": nu}) for ab, v in pickle.load(open("data/refined_kauffman/hopf_g_odd.pkl", "rb")).items()}


def adams(e, d):
    return e.subs({r: r**d, s: s**d, nu: nu**d}, simultaneous=True)


g = {}
for (a, b), v in god.items():
    e = v
    if a % 2 == 0 and b % 2 == 0:
        e = e - sp.Rational(1, 2) * adams(god[(a // 2, b // 2)], 2)
    g[(a, b)] = sp.cancel(sp.together(e))
pickle.dump({ab: str(v) for ab, v in g.items()}, open("data/refined_kauffman/hopf_g_alld.pkl", "wb"))

for ab in sorted(g, key=lambda k: (k[0] + k[1], k)):
    v = g[ab]
    if v == 0:
        print(f"{ab}: g = 0")
        continue
    den = sp.fraction(v)[1]
    R = sp.cancel(v * 2 * (1 + t))
    Rn, Rd = sp.fraction(R)
    P = sp.Poly(sp.expand(Rn), r, s, nu)
    print(f"{ab}: 分母 {sp.factor(den)} | 2(1+t)g: Laurent={sp.Poly(Rd, r, s, nu).is_monomial} 整数={all(c.is_integer for c in P.coeffs())} 項数={len(P.terms())}")

z, zp = r - 1 / r, r + 1 / r
E = {}
for ab in ((2, 2), (2, 4), (4, 4)):
    R = sp.cancel(g[ab] * 2 * (1 + t))
    E[ab] = sp.factor(sp.cancel(sp.cancel(R.subs(nu, nut * s / r)).subs(s, sp.I)))
    print(f"E_{ab}(t=−1) [ν~] = {E[ab]}")
E22 = z * zp**2 * (nut - 1 / nut) * (nut / r - r / nut)
print("E_22(−1) = z (z^+)² (ν~−ν~⁻¹)(q^{−1/2}ν~ − q^{1/2}ν~⁻¹):", sp.simplify(E[(2, 2)] - E22) == 0)
print("E_24(−1)/E_22(−1) =", sp.factor(sp.cancel(E[(2, 4)] / E[(2, 2)])), "（= χ_1 の x² 版：q^{−1}ν~² + q ν~^{−2}）")
print("E_44(−1)/E_22(−1) =", sp.factor(sp.cancel(E[(4, 4)] / E[(2, 2)])))
