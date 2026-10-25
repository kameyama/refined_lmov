r"""framed unknot から読んだ Ĩ の符号の規則の検証（theory_orientifold_index.md §11.8、progress.md 続き 46）。
  framed unknot（x = 2, 4）では E = Σ_m (−1)^{a+b+k} c_m m（P = Σ c_m m、m = ν^{2a} q^b t^k、すなわち r^{2b} s^{2k} ν^{2a}）。
  単項式の指数 (i, j, l) = (2b, 2k, 2a) で符号は (−1)^{(i+j+l)/2}。半整数の指数（r, s, ν の奇数冪）があるときは (i+j+l)/2 が整数のときだけ定義する。
  検証：unknot（再現）、Hopf (2,2),(2,4)（既知の一意な E と比較）、結び目 T(2,3), T(2,5)（条件 (i)(ii) を満たすか）。"""
import sys, io, contextlib, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp

r, s, nu = sp.symbols("r s nu", positive=True)
t = s**2


def sign_rule(P, shift=(0, 0, 0)):
    """P の各単項式 r^i s^j ν^l に (−1)^{((i−i0)+(j−j0)+(l−l0))/2} を掛ける。指数の和の偶奇がそろわなければ None。"""
    n, d = sp.fraction(sp.cancel(P))
    Pd = sp.Poly(d, r, s, nu); (md, cd), = Pd.terms()
    out = 0
    par = None
    for (i, j, l), c in sp.Poly(sp.expand(n), r, s, nu).terms():
        e = (i - md[0] - shift[0]) + (j - md[1] - shift[1]) + (l - md[2] - shift[2])
        if par is None:
            par = e % 2
        elif e % 2 != par:
            return None
        out += (-1)**((e - par) // 2) * c * r**(i - md[0]) * s**(j - md[1]) * nu**(l - md[2]) / cd
    return sp.expand(out)


def check_conditions(E, g2, g1):
    """(i) (1+t) | 2(1+t)g_(2) − E、(ii) t = q で E = −(1+q)(g_(1)(q²,ν²) + 1)。"""
    R_ = sp.cancel(2 * (1 + t) * g2)
    ci = sp.Poly(sp.fraction(sp.cancel((R_ - E) / (1 + t)))[1], r, s, nu).is_monomial
    g1x2 = g1.subs({r: r**2, s: s**2, nu: nu**2}, simultaneous=True)
    Eq = sp.cancel((-(1 + r**2) * (g1x2 + 1)).subs(s, r))
    cii = sp.cancel(E.subs(s, r) - Eq) == 0
    return ci, cii


def hopf():
    import rk_k4_fixed_sector as H
    Ek = {k: sp.sympify(v, locals={"r": r, "s": s, "nu": nu}) for k, v in pickle.load(open("data/refined_kauffman/hopf_E_fixed.pkl", "rb")).items()}
    for ab in ((2, 2), (2, 4)):
        sec = (ab[0] // 2, ab[0] // 2, ab[1] // 2, ab[1] // 2)
        P = sp.cancel((1 - t) * H.f4[sec].subs({H.r: r, H.s: s, H.nu: nu}))
        hits = []
        for sh in [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]:
            E = sign_rule(P, sh)
            if E is not None and sp.cancel(E - Ek[ab]) == 0:
                hits.append(("+", sh))
            if E is not None and sp.cancel(E + Ek[ab]) == 0:
                hits.append(("−", sh))
        print(f"Hopf {ab}：符号の規則で既知の E を再現する (全体の符号, 基準の指数のずれ) = {hits}", flush=True)


def knot(m):
    import rk_b_knot_fixed_sector as K
    import rk_b_cover_split as C
    u = sp.Symbol("u", positive=True)
    g = K.load_g(m)
    g = {k: v.subs({K.r: r, K.s: s, K.nu: nu}) for k, v in g.items()}
    with contextlib.redirect_stdout(io.StringIO()):
        sec = C.cover_sectors(m)
    P = sp.cancel((1 - t) * sec["x"].to_sympy().subs(u, nu * s))
    res = []
    for sh in [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]:
        E = sign_rule(P, sh)
        if E is None:
            continue
        for sg in (1, -1):
            ci, cii = check_conditions(sg * E, g[(2,)], g[(1,)])
            res.append((sg, sh, ci, cii))
    ok = [x for x in res if x[2] and x[3]]
    print(f"T(2,{2*m+1})：(i)(ii) をともに満たす (符号, ずれ) = {ok}；(i) だけ {[x[:2] for x in res if x[2]]}、(ii) だけ {[x[:2] for x in res if x[3]]}", flush=True)


if __name__ == "__main__":
    hopf()
    for m in (1, 2):
        knot(m)


def sign_rule_general(P, w, shift=(0, 0, 0)):
    """P の各単項式 r^i s^j ν^l に (−1)^{(w·((i,j,l) − shift))/2} を掛ける（w ∈ {0,1}³）。偶奇がそろわなければ None。"""
    n, d = sp.fraction(sp.cancel(P))
    Pd = sp.Poly(d, r, s, nu); (md, cd), = Pd.terms()
    out, par = 0, None
    for (i, j, l), c in sp.Poly(sp.expand(n), r, s, nu).terms():
        e = w[0] * (i - md[0] - shift[0]) + w[1] * (j - md[1] - shift[1]) + w[2] * (l - md[2] - shift[2])
        if par is None:
            par = e % 2
        elif e % 2 != par:
            return None
        out += (-1)**((e - par) // 2) * c * r**(i - md[0]) * s**(j - md[1]) * nu**(l - md[2]) / cd
    return sp.expand(out)


def split_boundary(P, n):
    """P = (境界 z_q^{2n−1}) × モジュライ部分。捻った境界は z_q^{n−1}(z_q^+)^n。"""
    zq, zqp = r - 1 / r, r + 1 / r
    Pk = sp.cancel(P / zq**(2 * n - 1))
    assert sp.Poly(sp.fraction(Pk)[1], r, s, nu).is_monomial, "境界の因子で割り切れない"
    return Pk, zq**(n - 1) * zqp**n


def hopf_moduli():
    import rk_k4_fixed_sector as H
    Ek = {k: sp.sympify(v, locals={"r": r, "s": s, "nu": nu}) for k, v in pickle.load(open("data/refined_kauffman/hopf_E_fixed.pkl", "rb")).items()}
    for ab in ((2, 2), (2, 4)):
        sec = (ab[0] // 2, ab[0] // 2, ab[1] // 2, ab[1] // 2)
        P = sp.cancel((1 - t) * H.f4[sec].subs({H.r: r, H.s: s, H.nu: nu}))
        Pk, bnd = split_boundary(P, 2)
        hits = []
        for w in [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]:
            Ek_c = sign_rule_general(Pk, w)
            if Ek_c is None:
                continue
            E = sp.expand(bnd * Ek_c)
            for sg in (1, -1):
                if sp.cancel(sg * E - Ek[ab]) == 0:
                    hits.append((sg, w))
        print(f"Hopf {ab}：境界を §11.3 で分け、モジュライ部分の規則 (−1)^{{(w·e)/2}} で既知の E を再現する (符号, w) = {hits}", flush=True)


def knot_moduli(m):
    import rk_b_knot_fixed_sector as K
    import rk_b_cover_split as C
    u = sp.Symbol("u", positive=True)
    g = {k: v.subs({K.r: r, K.s: s, K.nu: nu}) for k, v in K.load_g(m).items()}
    with contextlib.redirect_stdout(io.StringIO()):
        sec = C.cover_sectors(m)
    P = sp.cancel((1 - t) * sec["x"].to_sympy().subs(u, nu * s))
    Pk, bnd = split_boundary(P, 1)
    res = []
    for w in [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]:
        for sh in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]:
            Ek_c = sign_rule_general(Pk, w, sh)
            if Ek_c is None:
                continue
            E = sp.expand(bnd * Ek_c)
            for sg in (1, -1):
                ci, cii = check_conditions(sg * E, g[(2,)], g[(1,)])
                res.append((sg, w, sh, ci, cii))
    print(f"T(2,{2*m+1})：(i)(ii) {[x[:3] for x in res if x[3] and x[4]]}；(i) だけ {[x[:3] for x in res if x[3]]}；(ii) だけ {[x[:3] for x in res if x[4]]}", flush=True)


def main_moduli():
    hopf_moduli()
    for m in (1, 2):
        knot_moduli(m)
