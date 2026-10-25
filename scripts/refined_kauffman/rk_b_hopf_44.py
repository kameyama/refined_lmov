r"""拡張版の S3（|R| = |S| = 4）：g_{R,S}（scripts/rk_b_hopf_ext.py 4 4 の保存）から
  (1) 検算：g_{R,S} = g_{S,R}、奇数のセクターの Laurent 性、1 変数（g_{(a),(b)}）が hopf_g_alld（別の経路）と一致、
  (2) Adams：E_{R,S} = Σ_{ρ,σ ⊢ 2} c_{2;ρ}^R c_{2;σ}^S ε_{ρσ} の t = −1 での値を 25 組の留数から求め（4 つの未知数、4 組で解いて残り 21 組で検査）、
      ε_{(2),(2)}(−1) が 1 変数の E_44(−1) に一致するか、
  (3) t = q：ε_{ρσ}(t = q) = 2(1+q)(g − g^{odd})|_{t=q}（g^{odd} は Σ_{d 奇数} の定義、Mariño の ĝ）。"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rcsknots.partitions import partitions
from rk_b_classify import r, s
from rk_k4_crosscap import nu
from rk_b_hopf_ext import adams2

t = s**2
nut = sp.Symbol("nut", positive=True)
loc = {"r": r, "s": s, "nu": nu}
mono = lambda e: sp.Poly(sp.fraction(sp.cancel(e))[1], r, s, nu).is_monomial


def at_m1(e):
    return sp.cancel(sp.cancel(sp.cancel(e).subs(nu, nut * s / r)).subs(s, sp.I))


if __name__ == "__main__":
    g = {k: sp.sympify(v, locals=loc) for k, v in pickle.load(open("data/refined_kauffman/hopf_g_two_multi_44.pkl", "rb")).items()}
    print("対称性 g_{R,S} = g_{S,R}:", all(sp.cancel(g[(R, S_)] - g[(S_, R)]) == 0 for (R, S_) in g if (S_, R) in g), flush=True)
    print("奇数のセクターで Laurent でないもの:", [k for k, v in g.items() if (sum(k[0]) % 2 or sum(k[1]) % 2) and not mono(v)], flush=True)
    g42 = {k: sp.sympify(v, locals=loc) for k, v in pickle.load(open("data/refined_kauffman/hopf_g_two_multi_42.pkl", "rb")).items()}
    print("(4,2) の保存と一致:", all(sp.cancel(g[k] - v) == 0 for k, v in g42.items()), flush=True)
    galld = {k: sp.sympify(v, locals=loc) for k, v in pickle.load(open("data/refined_kauffman/hopf_g_alld.pkl", "rb")).items()}
    ok1 = all(sp.cancel(g.get(((a,) if a else (), (b,) if b else ()), 0) - v) == 0 for (a, b), v in galld.items())
    print("1 変数の g_{(a),(b)} が hopf_g_alld（別の経路、a,b ≤ 4）と一致:", ok1, flush=True)
    # (2) Adams：t = −1 の値
    P4, P2 = partitions(4), partitions(2)
    R2 = {(R, S_): at_m1(2 * (1 + t) * g[(R, S_)]) for R in P4 for S_ in P4}
    eps_sym = {(a, b): sp.Symbol(f"e_{''.join(map(str, a))}_{''.join(map(str, b))}") for a in P2 for b in P2}
    eqs = {(R, S_): sp.expand(sum(adams2(rho).get(R, 0) * adams2(sg).get(S_, 0) * eps_sym[(rho, sg)] for rho in P2 for sg in P2) - R2[(R, S_)])
           for R in P4 for S_ in P4}
    pick = [((4,), (4,)), ((4,), (1, 1, 1, 1)), ((1, 1, 1, 1), (4,)), ((1, 1, 1, 1), (1, 1, 1, 1))]
    sol = sp.solve([eqs[k] for k in pick], list(eps_sym.values()), dict=True)[0]
    bad = [k for k, e in eqs.items() if sp.cancel(e.subs(sol)) != 0]
    print(f"Adams（t = −1）：4 組で解いた ε で残りの 21 組が成り立たないもの: {bad}", flush=True)
    E1 = {k: sp.sympify(v, locals=loc) for k, v in pickle.load(open("data/refined_kauffman/hopf_E_fixed.pkl", "rb")).items()}
    print("ε_{(2),(2)}(−1) == 1 変数の E_44(−1):", sp.cancel(sol[eps_sym[((2,), (2,))]] - at_m1(E1[(4, 4)])) == 0, flush=True)
    print("ε_{(2),(1,1)}(−1) == ε_{(1,1),(2)}(−1):", sp.cancel(sol[eps_sym[((2,), (1, 1))]] - sol[eps_sym[((1, 1), (2,))]]) == 0, flush=True)
    for k, v in eps_sym.items():
        print(f"  ε_{k}(−1) = {sp.factor(sol[v])}", flush=True)
    pickle.dump({str(k): str(sol[v]) for k, v in eps_sym.items()}, open("data/refined_kauffman/hopf_eps44_m1.pkl", "wb"))
