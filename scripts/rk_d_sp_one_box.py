r"""案 D：refined Sp(2N)（C_N、u = t²）の [1] のクロスキャップ部分 g^{Sp}_[1]（progress.md 続き 50）。
  SO の較正（rk_d_hyper_fit.calibrate_SO）：W_SO = 鏡像[ℋ^C(u=1, a = −ν/s)]·dim_V^{SO}(ν)·(s/ν)^{2m}、g^{true}_[1] = W_SO − W^{cov}_1（T(2,3), T(2,5) で一意）。
  Sp：ν = t^{N+1/2}（A = t^{2N+1}）、a = −t^{N−1} = −ν/s³。W_Sp = 鏡像[ℋ^C(u=s⁴, a = −ν/s³)]·dim_V^{Sp}(ν)·mono、mono = ±r^c s^d ν^e を
  (1) t = 1 の極の打ち消し（g が Laurent）、(2) t = q（u = q²）で Kauffman の形（r → −1/r で r^x 倍を除いて不変）で決める。
  得た g^{Sp}_[1] を g^{SO}_[1] = g^{true} と比べる：C_{c=1} = ½(g^{SO} − g^{Sp})、C_{c=2} = ½(g^{SO} + g^{Sp}) と 𝐚 の偶奇の分割。"""
import sys, pickle, itertools
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp
import rk_d_hyper_fit as C

q, t, u, a = C.q, C.t, C.u, C.a
r_, s_ = C.r_, C.s_
OUT = "data/refined_kauffman/daha"


def load_H(m):
    d = pickle.load(open(f"{OUT}/cn_hyper_omega1_m{m}.pkl", "rb"))
    return sp.sympify(d["H"], locals={"q": q, "t": t, "u": u, "a": a})


def W_group(m, H, kind, nu):
    if kind == "SO":
        Hs = C.H_rs(H, 1).subs(a, -nu / s_)
    else:
        Hs = C.H_rs(H, s_**4).subs(a, -nu / s_**3)
    Hs = Hs.subs({r_: 1 / r_, s_: 1 / s_, nu: 1 / nu}, simultaneous=True)
    return sp.cancel(Hs * C.dim_V(kind, nu))


def is_laurent(e, nu):
    n_, d_ = sp.fraction(sp.cancel(e))
    return sp.Poly(d_, r_, s_, nu).is_monomial


def kauffman_like(e, nu):
    """t = q（s = r）、u = q² は H の代入で済んでいる。r → −1/r で r^x 倍を除いて不変か。"""
    ee = sp.cancel(e.subs(s_, r_))
    for x in range(-40, 41):
        f = sp.cancel(ee * r_**x)
        if sp.cancel(f - f.subs(r_, -1 / r_)) == 0:
            return x
    return None


def find_mono(m, H, nu, Wc, gso):
    W0 = W_group(m, H, "Sp", nu)
    hits = []
    for c, e in itertools.product(range(-12, 13), range(-14, 15)):
        g = sp.cancel(W0 * r_**c * nu**e - Wc)
        # t = 1 の極：s = ±1 で分母が 0 にならないか（s − 1, s + 1 の因子）
        d_ = sp.fraction(g)[1]
        if sp.Poly(d_, r_, s_, nu).is_monomial:
            hits.append((c, e))
    print(f"T(2,{2*m+1})：t = 1 の極が消える (c, e) = {hits}", flush=True)
    out = []
    for c, e in hits:
        for d in range(-30, 31):
            for sg in (1, -1):
                g = sp.cancel(sg * W0 * r_**c * s_**d * nu**e - Wc)
                x = kauffman_like(g, nu)
                if x is not None:
                    out.append((sg, c, d, e, x, g))
    print(f"   t = q で Kauffman の形になる (符号, c, d, e, r の冪) = {[o[:5] for o in out]}", flush=True)
    return out


if __name__ == "__main__":
    for m in (1, 2):
        H = load_H(m)
        gso, nu = C.g_known(m)
        Wc = C.cover_W1(m, nu)
        # SO の再現（確認）
        gso_rebuilt = sp.cancel(W_group(m, H, "SO", nu) * (s_ / nu)**(2 * m) - Wc)
        print(f"T(2,{2*m+1})：SO の再現 g^{{true}} − (W_SO − W^{{cov}}) = {sp.simplify(gso - gso_rebuilt)}", flush=True)
        res = find_mono(m, H, nu, Wc, gso)
        for sg, c, d, e, x, gsp in res:
            A_, Q_, T_ = sp.symbols("A Q T")
            bold = lambda f: sp.expand(sp.cancel(f.subs({r_: -T_ * Q_, nu: A_}).subs(s_, Q_)))
            gb_so, gb_sp = bold(gso), bold(gsp)
            odd = lambda f: sp.expand((f - f.subs(A_, -A_)) / 2)
            even = lambda f: sp.expand((f + f.subs(A_, -A_)) / 2)
            print(f"   候補 (符号 {sg}, c {c}, d {d}, e {e})：g^Sp の項数 {len(sp.Add.make_args(gb_sp))}", flush=True)
            print(f"      g^Sp(𝐚) = g^SO(−𝐚)：{sp.expand(gb_sp - gb_so.subs(A_, -A_)) == 0}；g^Sp = g^SO(1/𝐚)·?：{C.is_monomial_rs(sp.cancel(gsp / gso.subs(nu, 1/nu)), nu)}", flush=True)
            print(f"      ½(g^SO − g^Sp) = 𝐚 奇数部分 of g^SO：{sp.expand((gb_so - gb_sp)/2 - odd(gb_so)) == 0}；½(g^SO + g^Sp) = 𝐚 偶数部分：{sp.expand((gb_so + gb_sp)/2 - even(gb_so)) == 0}", flush=True)
            pickle.dump({"g_sp": str(gsp), "mono": (sg, c, d, e)}, open(f"{OUT}/sp_one_box_m{m}.pkl", "wb"))
