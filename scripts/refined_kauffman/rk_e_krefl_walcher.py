r"""Krefl–Walcher（1007.0263）の関係「β = 2 の理論は β = 1 の理論の orientifold」の結び目版の検討（progress.md 続き 61）。
  問い：unrefined の SO(2N) の [1]（Kauffman、t = q）は、refined U(N') の [1]（DGR の superpolynomial）を β = t_CS/q_CS の冪で
  特殊化したもの（β = 1/2：𝐭 = −𝐪、β = 2：𝐭 = −𝐪^{−1/2}）と、変数の単項式写像 a = λ^k q^j のもとで一致するか。
  Kauffman 側は保存データ largeN_2box_m{1,2}.pkl の f_(1)（f = −W_SO）を t = q にしたもの。unknot で割った簡約版も比べる。"""
import sys, pickle, itertools
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp

r_, s_, u_ = sp.symbols("r s u", positive=True)
lam = sp.Symbol("lambda", positive=True)       # λ = ν = U t^{-1/2}
qN = sp.Symbol("q", positive=True)             # NRZ の q = q_CS^{1/2} = r
A_, Q_, T_ = sp.symbols("A Q T")
loc = {"r": r_, "s": s_, "u": u_}
DGR = {1: A_**2 * Q_**-2 + A_**4 * T_**3 + A_**2 * Q_**2 * T_**2,
       2: A_**2 * Q_**-4 + A_**2 * T_**2 + A_**2 * Q_**4 * T_**4 + A_**4 * Q_**-2 * T_**3 + A_**4 * Q_**2 * T_**5}


def kauffman_one_box(m):
    d = pickle.load(open(f"data/refined_kauffman/daha/largeN_2box_m{m}.pkl", "rb"))
    f1 = sp.sympify(d["f"]["(1,)"], locals=loc).subs(u_, lam * s_)
    W = sp.cancel(-f1.subs(s_, r_))                                   # t = q、変数 r = q^{1/2} = q_NRZ、λ
    return sp.cancel(W.subs(r_, qN))


L_, x_ = sp.symbols("L x", positive=True)       # λ = L², q = x⁴ で半整数冪を消す


def is_mono(e):
    e = e.subs({lam: L_**2, qN: x_**4}, simultaneous=True)
    try:
        n_, d_ = sp.fraction(sp.cancel(e))
        return sp.Poly(n_, L_, x_).is_monomial and sp.Poly(d_, L_, x_).is_monomial
    except Exception:
        return False


if __name__ == "__main__":
    for m in (1, 2):
        W = kauffman_one_box(m)
        dim = sp.cancel((lam - 1 / lam) / (qN - 1 / qN) + 1)
        Wred = sp.cancel(W / dim)
        print(f"== T(2,{2*m+1})：Kauffman の [1]（unreduced）=", sp.factor(W), flush=True)
        print("   簡約 =", sp.expand(Wred), flush=True)
        hits = []
        for beta_lab, Tsub, Qsub in (("β=1/2 (t=q^{1/2})", -Q_, qN**sp.Rational(1, 2)), ("β=2 (t=q²)", -Q_**sp.Rational(-1, 2), qN**2), ("β=1 (t=q)", -1, qN)):
            P = DGR[m].subs(T_, Tsub).subs(Q_, Qsub)
            Pbar = sp.cancel((A_ * Tsub + 1 / A_) / (Qsub - 1 / Qsub) * P)          # unreduced：DGR の (at + a^{-1})/(q − q^{-1}) 倍
            for k, j in itertools.product((sp.Rational(1, 2), -sp.Rational(1, 2), 1, 2, -1, -2, 4, -4), (sp.Rational(x, 2) for x in range(-8, 9))):
                for mir in (False, True):
                    for lab, X, Pm in (("簡約", Wred, P), ("unreduced", W, Pbar)):
                        Pk = Pm.subs(A_, lam**k * qN**j)
                        if mir:
                            Pk = Pk.subs(qN, 1 / qN)
                        if is_mono(X / Pk):
                            hits.append((beta_lab, k, j, mir, lab, sp.factor(sp.cancel((X / Pk).subs({lam: L_**2, qN: x_**4}, simultaneous=True)))))
        print("   一致（β, k, j, 鏡像, 対象, 比）：", hits if hits else "なし", flush=True)


def crosscap_one_box(m):
    d = pickle.load(open(f"data/refined_kauffman/daha/largeN_2box_m{m}.pkl", "rb"))
    g1 = sp.sympify(d["g"]["(1,)"], locals=loc).subs(u_, lam * s_)
    return sp.cancel(sp.cancel(g1.subs(s_, r_)).subs(r_, qN))


def test_free_energy_form():
    """KW の形 F_unor(ε) = F(ε, −2ε) − ½F(ε, −ε) の開弦版：g_[1]|_{t=q} ?= P̄(a₂; β=2) − ½P̄(a₁; β=1)、a₂ = a₁²（a = t^{N'}）。"""
    for m in (1, 2):
        g = crosscap_one_box(m)
        P = DGR[m]
        Pb1 = sp.cancel(((A_ * T_ + 1 / A_) / (Q_ - 1 / Q_) * P).subs(T_, -1).subs(Q_, qN))                         # β = 1
        Pb2 = sp.cancel(((A_ * T_ + 1 / A_) / (Q_ - 1 / Q_) * P).subs(T_, -1 / qN).subs(Q_, qN**2))                 # β = 2：𝐪 = q_CS = q_N²、𝐭 = −q_N^{-1}
        hits = []
        for k, j in itertools.product((sp.Rational(1, 2), -sp.Rational(1, 2), 1, 2, -1, -2), (sp.Rational(x, 2) for x in range(-8, 9))):
            a1 = lam**k * qN**j
            for c in (sp.Rational(1, 2), 1, 2):
                for sgn in (1, -1):
                    Y = sp.cancel(sgn * (Pb2.subs(A_, a1**2) - c * Pb1.subs(A_, a1)))
                    for mir in (False, True):
                        Ym = Y.subs(qN, 1 / qN) if mir else Y
                        if is_mono(g / Ym):
                            hits.append((k, j, c, sgn, mir))
        print(f"T(2,{2*m+1})：g_[1]|_(t=q) = ±(P̄(a², β=2) − c·P̄(a, β=1)) × 単項式 の一致：", hits if hits else "なし", flush=True)
        # g の形（参考）
        print("   g_[1]|_(t=q) =", sp.factor(g), flush=True)


if __name__ == "__main__":
    test_free_energy_form()
