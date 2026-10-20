r"""refined の強い整数性の検査（theory_orientifold_index.md §11.4 (c)、progress.md 続き 45）。
  AgSc §4：f^{or}_{□,□}(q,Q) ≡ f^{unor}_□(q²,Q²) (mod 2)。全 d の設定では、brane–image のセクターの状態の指数 P = (1 − t)x と、
  [1] の凍結した状態 C_(1) = g_(1) について、P(q,t,ν) ≡ g_(1)(q², t², ν² w^{j})（mod 2、係数ごと）が成り立つかを w^{j} = (q/t)^{j/2} のずれ j で調べる。
  成り立てば「[1] の凍結した状態 ↔ brane–image の対」の mod 2 の対応が refined でも成り立つ。"""
import sys, io, contextlib
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp
import rk_b_knot_fixed_sector as K
import rk_b_cover_split as C

r, s, nu = K.r, K.s, K.nu
u = sp.Symbol("u", positive=True)


def coeffs_mod2(e):
    n, d = sp.fraction(sp.cancel(e))
    Pd = sp.Poly(d, r, s, nu)
    assert Pd.is_monomial, d
    (md, cd), = Pd.terms()
    out = {}
    for mo, c in sp.Poly(sp.expand(n), r, s, nu).terms():
        c = sp.Rational(c, cd)
        assert c.is_integer, c
        if int(c) % 2:
            out[tuple(a - b for a, b in zip(mo, md))] = 1
    return out


def run(m):
    g = K.load_g(m)
    with contextlib.redirect_stdout(io.StringIO()):
        sec = C.cover_sectors(m)
    x = sp.cancel(sec["x"].to_sympy().subs(u, nu * s))
    P = sp.cancel((1 - s**2) * x)
    cP = coeffs_mod2(P)
    print(f"T(2,{2*m+1})：P の奇数係数の単項式 {len(cP)} 個", flush=True)
    for j in range(-3, 4):
        g2 = g[(1,)].subs({r: r**2, s: s**2, nu: nu**2 * (r / s)**j}, simultaneous=True)
        cG = coeffs_mod2(sp.cancel(g2))
        diff = set(cP) ^ set(cG)
        print(f"   j = {j:+d}：ψ²[g_(1)] の奇数係数 {len(cG)} 個、mod 2 の差 {len(diff)} 個", flush=True)


if __name__ == "__main__":
    for m in (1, 2):
        run(m)


def run_unrefined(m):
    """t = q（AgSc §4 の unrefined の形）での同じ検査。"""
    g = K.load_g(m)
    with contextlib.redirect_stdout(io.StringIO()):
        sec = C.cover_sectors(m)
    x = sp.cancel(sec["x"].to_sympy().subs(u, nu * s))
    P = sp.cancel(((1 - s**2) * x).subs(s, r))
    cP = coeffs_mod2(P)
    g2 = sp.cancel(g[(1,)].subs(s, r).subs({r: r**2, nu: nu**2}, simultaneous=True))
    cG = coeffs_mod2(g2)
    print(f"T(2,{2*m+1})、t = q：P の奇数係数 {len(cP)} 個、ψ²[g_(1)] の奇数係数 {len(cG)} 個、mod 2 の差 {len(set(cP) ^ set(cG))} 個", flush=True)
