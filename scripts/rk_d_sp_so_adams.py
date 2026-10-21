r"""SO と Sp の 2 箱の p₂ 方向の関係：Adams 項 ψ²[g_(1)] を除けば符号規則が成り立つか（progress.md 続き 54）。
  仮説：[1] で g^{Sp}_1(𝐚) = −g^{SO}_1(−𝐚) なので ψ²[g^{Sp}_1] = −ψ²[g^{SO}_1]（ψ² は 𝐚² しか見ない）。
  g_R（または p₂ 成分 D）が「Adams 項 α ψ²[g_1] ＋ 残り」で、残りが Sp(𝐚) = ε·SO(−𝐚) なら、
  X_ε = Sp(𝐚) − ε SO(−𝐚) = −α(1 − ε) ψ²[g^{SO}_1]（ε = −1 のとき −2α ψ²[g^{SO}_1]）。X_ε / ψ²[g^{SO}_1] を因数分解して見る。
変数 (r, s, ν)。ψ²：r → r², s → s², ν → ν²（x = νs → x²）。"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp

r_, s_, u_ = sp.symbols("r s u", positive=True)
nu = sp.Symbol("nu", positive=True)
loc = {"r": r_, "s": s_, "u": u_}
OUT = "data/refined_kauffman/daha"


def load():
    sp_d = pickle.load(open(f"{OUT}/largeN_2box_sp_m1.pkl", "rb"))
    so_d = pickle.load(open(f"{OUT}/largeN_2box_m1.pkl", "rb"))
    g_sp = {k: sp.cancel(sp.sympify(v, locals=loc).subs(u_, nu / s_)) for k, v in sp_d["g"].items()}
    g_so = {k: sp.cancel(sp.sympify(v, locals=loc).subs(u_, nu * s_)) for k, v in so_d["g"].items()}
    return g_sp, g_so


def psi2(e):
    return e.subs({r_: r_**2, s_: s_**2, nu: nu**2}, simultaneous=True)


def neg(e):
    return e.subs(nu, -nu)


if __name__ == "__main__":
    g_sp, g_so = load()
    P2 = sp.cancel(psi2(g_so["(1,)"]))
    print("ψ²[g^SO_1] =", sp.factor(P2), flush=True)
    zq, zqp, t = r_ - 1 / r_, r_ + 1 / r_, s_**2
    cands = {
        "g_(2)": (g_sp["(2,)"], g_so["(2,)"]), "g_(1,1)": (g_sp["(1, 1)"], g_so["(1, 1)"]),
        "g_(2) vs SO_(1,1)": (g_sp["(2,)"], g_so["(1, 1)"]), "g_(1,1) vs SO_(2)": (g_sp["(1, 1)"], g_so["(2,)"]),
        "p₂ = (g_(2)−g_(1,1))/z_q⁺": (sp.cancel((g_sp["(2,)"] - g_sp["(1, 1)"]) / zqp), sp.cancel((g_so["(2,)"] - g_so["(1, 1)"]) / zqp)),
        "p₁² = (g_(2)+g_(1,1))/z_q": (sp.cancel((g_sp["(2,)"] + g_sp["(1, 1)"]) / zq), sp.cancel((g_so["(2,)"] + g_so["(1, 1)"]) / zq)),
    }
    for lab, (a, b) in cands.items():
        for eps in (+1, -1):
            X = sp.cancel(a - eps * neg(b))
            if X == 0:
                print(f"{lab}, ε = {eps:+d}: X = 0", flush=True); continue
            ratio = sp.factor(sp.cancel(X / P2))
            n_, d_ = sp.fraction(ratio)
            nt = len(sp.Add.make_args(sp.expand(n_)))
            print(f"{lab}, ε = {eps:+d}: X/ψ²[g_1] = {ratio if nt <= 12 else f'(分子 {nt} 項)'}", flush=True)
