r"""K4'(a)：Ĩ = 入れ替え の表現論的な検査（成分 1 に probe 無限個、成分 2 は 1 変数、|R| ≤ 4、b ≤ 2）。
被覆 f^{cov}_{ρ,ρ',k,l}（kauffman_cover.cover_f_multi）から ½f^{cov}_{R,b} = ½ Σ_{ρ,ρ',k+l=b} N^R_{ρρ'} f^{cov}_{ρρ'kl}（Y' = Y、z1 = z2 = z）、
SO の f_{R,b}（K2'(b)、hopf_f_multi.pkl、Schur 基底）との差 g_{R,b} = f_{R,b} − ½f^{cov}_{R,b}。
Ĩ で固定される状態は (ρ,ρ) × (b/2,b/2) のセクターにあり、Ĩ は ρ⊗ρ の入れ替えとして作用するので
  E_{R,b} = Σ_ρ c_{2;ρ}^R ε_{ρ,b}（p_2∘s_ρ = Σ_R c_{2;ρ}^R s_R、Adams）。
予言：E_[2] = ε_(1)、E_[1,1] = −ε_(1)；E_[4] = ε_(2)、E_[3,1] = −ε_(2)、E_[2,2] = ε_(2) + ε_(1,1)、E_[2,1,1] = −ε_(1,1)、E_[1^4] = ε_(1,1)。
  ⇒ g_[2]+g_[1,1]、g_[4]+g_[3,1]、g_[2,1,1]+g_[1^4]、g_[2,2]−g_[4]−g_[1^4] は (1+t) の分母を持たない。ε_(2),1 = E_{42}（1 変数）、ε_(1),1 = E_{22}。
出力：data/refined_kauffman/hopf_g_multi.pkl。"""
import sys, time, pickle
sys.path.insert(0, "src")
import sympy as sp
from rcsknots.kauffman_cover import cover_f_multi, qs_, ss_, nus_
from rcsknots.characters import chi
from rcsknots.partitions import partitions, z as z_lam

r, s, nu = sp.symbols("r s nu", positive=True)
nut = sp.Symbol("nut")
t = s**2
A_, B = 4, 2
t0 = time.time()
fc = {k: v.subs({qs_: r**2, ss_: s, nus_: nu}) for k, v in cover_f_multi(A_, B).items()}
print(f"cover_f_multi({A_},{B}): {time.time()-t0:.0f}s", flush=True)


def lr(r1, r2):
    """s_{r1} s_{r2} = Σ_R N^R s_R（冪和の指標で計算）。"""
    n = sum(r1) + sum(r2)
    out = {}
    for R in partitions(n):
        c = 0
        for m1 in partitions(sum(r1)):
            for m2 in partitions(sum(r2)):
                mu = tuple(sorted(m1 + m2, reverse=True))
                c += sp.Rational(chi(r1, m1), z_lam(m1)) * sp.Rational(chi(r2, m2), z_lam(m2)) * chi(R, mu)
        if c:
            out[R] = c
    return out


half_cov = {}
for (r1, r2, k, l), v in fc.items():
    for R, c in (lr(r1, r2).items() if sum(r1) + sum(r2) else [((), 1)]):
        key = (R, k + l)
        half_cov[key] = half_cov.get(key, 0) + c * v / 2
q_, t_, U_ = sp.symbols("q t U")
fSO = {k: sp.sympify(v, locals={"q": q_, "t": t_, "U": U_}).subs({q_: r**2, t_: s**2, U_: nu * s})
       for k, v in pickle.load(open("data/refined_kauffman/hopf_f_multi.pkl", "rb"))["schur"].items()}
g = {}
for key, v in fSO.items():
    if key in half_cov or True:
        g[key] = sp.cancel(v - half_cov.get(key, 0))
pickle.dump({k: str(v) for k, v in g.items()}, open("data/refined_kauffman/hopf_g_multi.pkl", "wb"))


def den_has_1pt(e):
    d = sp.fraction(sp.cancel(e))[1]
    return sp.rem(sp.Poly(d, s), sp.Poly(1 + s**2, s)).is_zero and d != 0


def laurent(e):
    n, d = sp.fraction(sp.cancel(e))
    return sp.Poly(d, r, s, nu).is_monomial


print("各 g_{R,b} の分母（単項式以外の因子）:")
for key in sorted(g, key=lambda k: (sum(k[0]), k[1], k[0])):
    v = g[key]
    if v == 0:
        continue
    d = sp.factor(sp.fraction(v)[1])
    nonmono = [f for f in sp.Mul.make_args(d) if not sp.Poly(f, r, s, nu).is_monomial]
    print(f"  {key}: {nonmono if nonmono else '整数 Laurent' if all(c.is_integer for c in sp.Poly(sp.expand(sp.fraction(v)[0]), r, s, nu).coeffs()) else '分母は単項式（非整数係数）'}")
G = lambda R, b: g.get((tuple(R), b), 0)
print("Adams の予言（(1+t) の分母が消えるか）:")
for name, combo in (("g_[2]+g_[1,1] (b=2)", G((2,), 2) + G((1, 1), 2)),
                    ("g_[4]+g_[3,1] (b=2)", G((4,), 2) + G((3, 1), 2)),
                    ("g_[2,1,1]+g_[1^4] (b=2)", G((2, 1, 1), 2) + G((1, 1, 1, 1), 2)),
                    ("g_[2,2]−g_[4]−g_[1^4] (b=2)", G((2, 2), 2) - G((4,), 2) - G((1, 1, 1, 1), 2))):
    print(f"  {name}: Laurent = {laurent(combo)}")
E = {k: sp.sympify(v, locals={"r": r, "s": s, "nu": nu}) for k, v in pickle.load(open("data/refined_kauffman/hopf_E_fixed.pkl", "rb")).items()}
E22, E24 = E[(2, 2)], E[(2, 4)]
for R, eps, name in (((2,), E22, "E_[2],2 = E_22"), ((4,), E24, "E_[4],2 = E_42 = E_24")):
    C = sp.cancel(G(R, 2) - eps / (2 * (1 + t)))
    print(f"  {name}: g − ½E/(1+t) が Laurent = {laurent(C)}")

# ε_(1,1),(1)：E_[1^4],2 の固定セクター ((1,1),(1,1)) × (1,1) の状態の部分和として決める
sys.path.insert(0, "scripts")
from rk_k4_fixed_sector import solve_generic
P11 = (1 - t) * fc[((1, 1), (1, 1), 1, 1)]
g1111 = G((1, 1, 1, 1), 2)
g11_1_x2 = G((1, 1), 1).subs({r: r**2, s: s**2, nu: nu**2}, simultaneous=True)
Eq = sp.cancel((-(1 + t) * g11_1_x2).subs(s, r))          # t = q：ε = −(1+q) ĝ_{[1,1],1}(x²)
res = solve_generic(P11, 2 * (1 + t) * g1111, Eq)
print("ε_(1,1),(1)：解あり =", res["feasible"], " 状態の単項式", res["n"], " 一意でない変数", len(res.get("free", [])))
if res["feasible"]:
    e11 = res["E"]
    print("  P_(1,1) (ν~) =", sp.factor(sp.cancel(P11.subs(nu, nut * s / r))))
    print("  ε_(1,1) (ν~) =", sp.factor(sp.cancel(e11.subs(nu, nut * s / r))))
    print("  ε_(1,1) (ν) =", sp.factor(e11))
    print("  ε_(2)   (ν) =", sp.factor(E24))
    print("  検算 g_[2,2] − ½(ε_(2)+ε_(1,1))/(1+t) が Laurent:", laurent(sp.cancel(G((2, 2), 2) - (E24 + e11) / (2 * (1 + t)))))
    pickle.dump({"eps_11_1": str(e11)}, open("data/refined_kauffman/hopf_eps_multi.pkl", "wb"))
