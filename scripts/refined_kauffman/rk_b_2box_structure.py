r"""タスク C の (b)：大きな N の 2 箱のクロスキャップ部分（trefoil、一様な規格化、全 d）の hatted ĝ = M(q^{−1})^{−1} g の構造。
  - 極：ĝ_(2) + ĝ_(1,1) は整数 Laurent、(1 + t) の極は p_2 方向だけ。奇数 d に直す（素朴な Adams）と (1 + q) の極が加わる。
  - G 展開（V ↔ 1/V、q → t^{−1}, t → q^{−1}、ν 固定）：c = 1 部分（V の偶奇＝ν の偶奇で分ける）は [1] で対称、2 箱では非対称。
  - 述べ直しの判定式（Hopf、report §8）：X = ĝ − ĝ|_{t=q}(ν → ν~)、ν~ = ν (q/t)^{1/2} が (q − t)(ν² − t)（太字で (𝐭² − 1)(𝐚² − 𝐪²)）で割り切れる。
変数は (r, s, ν) = (q^{1/2}, t^{1/2}, U t^{−1/2})。"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp

r, s, u = sp.symbols("r s u", positive=True)
nu = sp.symbols("nu", positive=True)


def ghat():
    d = pickle.load(open("data/refined_kauffman/daha/largeN_2box_m1.pkl", "rb"))
    g = {eval(k): sp.sympify(v, locals={"r": r, "s": s, "u": u}).subs(u, nu * s) for k, v in d["g"].items()}
    det = r**2 - r**-2                                            # M(q^{-1}) = [[r, −1/r], [−1/r, r]]
    h2 = sp.cancel((r * g[(2,)] + g[(1, 1)] / r) / det)
    h11 = sp.cancel((g[(2,)] / r + r * g[(1, 1)]) / det)
    return {(1,): sp.cancel(g[(1,)]), (2,): h2, (1, 1): h11, "sum": sp.cancel(h2 + h11)}


def sigma(e):
    return e.subs({r: 1 / s, s: 1 / r}, simultaneous=True)


def blocks(e):
    """分子を (V の偶奇, ν の偶奇)（分母に対する相対）で分ける。V = (rs)^{1/2}。"""
    n_, d_ = sp.fraction(sp.cancel(e))
    P = sp.Poly(sp.expand(n_), r, s, nu)
    Dm = sp.Poly(sp.expand(d_), r, s, nu).monoms()
    vp, npar = (Dm[0][0] + Dm[0][1]) % 2, Dm[0][2] % 2
    out = {}
    for (a, b, c), cf in P.terms():
        k = ((a + b - vp) % 2, (c - npar) % 2)
        out[k] = out.get(k, 0) + cf * r**a * s**b * nu**c
    return {k: v / d_ for k, v in out.items()}


def c1_symmetric(e):
    """c = 1 部分（V の偶数冪のブロック）が σ で対称か。"""
    return all(sp.cancel(v - sigma(v)) == 0 for k, v in blocks(e).items() if k[0] == 0)


def refined_only(e):
    """X = e − e|_{t=q}(ν → ν~) と、(q − t)(ν² − t) で割った商。"""
    X = sp.cancel(e - sp.cancel(e.subs(s, r)).subs(nu, nu * r / s))
    return X, sp.cancel(X / ((r**2 - s**2) * (nu**2 - s**2)))


def divisible_universal(e):
    X, _ = refined_only(e)
    n_ = sp.Poly(sp.fraction(X)[0], nu, r)
    return all(sp.rem(n_, sp.Poly(f, nu, r)).is_zero for f in (r**2 - s**2, nu**2 - s**2))


def sign_blocks(e):
    """商を太字（r = −𝐭𝐪, s = 𝐪, ν = 𝐚）で 𝐚 の冪ごとに分け、{冪: (項数, 符号)}。"""
    A, Q, T = sp.symbols("A Q T")
    _, Qt = refined_only(e)
    n_, d_ = sp.fraction(sp.cancel(Qt.subs({r: -T * Q, s: Q, nu: A}, simultaneous=True)))
    shift = sp.Poly(d_, A).degree()
    by = {}
    for (ea, _, _), cf in sp.Poly(sp.expand(n_), A, Q, T).terms():
        by.setdefault(ea - shift, []).append(cf)
    return {k: (len(v), "+" if all(c > 0 for c in v) else "−" if all(c < 0 for c in v) else "混在") for k, v in sorted(by.items())}


if __name__ == "__main__":
    G = ghat()
    for k in ((1,), (2,), (1, 1), "sum"):
        e = G[k]
        print(f"ĝ_{k}：分母 {sp.factor(sp.fraction(e)[1])}、c=1 が σ 対称 {c1_symmetric(e)}、"
              f"refined だけの部分が (q−t)(ν²−t) で割り切れる {divisible_universal(e)}", flush=True)
    for k in ((1,), "sum"):
        print(f"ĝ_{k} の商の 𝐚 の冪ごとの符号：", sign_blocks(G[k]), flush=True)
