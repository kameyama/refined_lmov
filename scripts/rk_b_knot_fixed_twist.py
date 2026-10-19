r"""(a) 固定セクターで Ĩ をゼロモードの (−1)^F として扱う規則（progress.md 続き 36）。
  P = (1 − t)x を因数分解し、差の因子のうち (ν² − 1), (q − 1), (qν² − 1) を和の因子 (ν² + 1), (q + 1), (qν² + 1) に替え、普遍因子 (ν² − t) は残す。
  残りの多項式 L の単項式 r^a s^b ν^c には偶奇の符号 ε = ±(−1)^{α a + β b + γ c} を付ける。
  E の候補が条件 (i)（C が整数 Laurent）と (ii)（t = q で E = −(1 + q)(g_(1)(x²) + 1)）を満たすかを調べる。
使い方：python scripts/rk_b_knot_fixed_twist.py m"""
import sys, io, contextlib, itertools
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp
import rk_b_knot_fixed_sector as K
import rk_b_cover_split as C

r, s, nu = K.r, K.s, K.nu
u = sp.Symbol("u", positive=True)
t, q = s**2, r**2


def main(m):
    g = K.load_g(m)
    with contextlib.redirect_stdout(io.StringIO()):
        sec = C.cover_sectors(m)
    x = sp.cancel(sec["x"].to_sympy().subs(u, nu * s))
    P = sp.cancel((1 - t) * x)
    num, den = sp.fraction(P)
    c0, fl = sp.factor_list(num)
    print(f"T(2,{2*m+1})：P の分子の因子 {[(sp.factor(f), k) for f, k in fl if len(sp.Add.make_args(f)) <= 3]}、分母 {den}", flush=True)
    # P の分子 = c0 · s² (s² − ν²) (ν² − 1)(r² − 1)(ν² r² − 1) · L を確かめて分ける
    base = s**2 * (nu**2 - s**2) * (nu**2 - 1) * (r**2 - 1) * (nu**2 * r**2 - 1)
    L = sp.cancel(num / base)
    assert sp.fraction(L)[1].is_number, "想定した因子を含まない"
    keep = s**2 * (nu**2 - s**2)
    flipped = (nu**2 + 1) * (r**2 + 1) * (nu**2 * r**2 + 1)
    print(f"   そのまま残す因子 {sp.factor(keep)}、和に替えた因子 {sp.factor(flipped)}、L の項数 {len(sp.Add.make_args(sp.expand(L)))}", flush=True)
    R = sp.cancel(2 * (1 + t) * g[(2,)])
    g1x2 = g[(1,)].subs({r: r**2, s: s**2, nu: nu**2}, simultaneous=True)
    Eq = sp.cancel((-(1 + q) * (g1x2 + 1)).subs(s, r))
    PL = sp.Poly(sp.expand(L), r, s, nu)
    hits = []
    for (al, be, ga), sg in itertools.product(itertools.product((0, 1), repeat=3), (1, -1)):
        Lt = sum(sg * (-1)**(al * a + be * b + ga * c) * cf * r**a * s**b * nu**c for (a, b, c), cf in PL.terms())
        E = sp.cancel(keep * flipped * Lt / den)
        ci = sp.Poly(sp.fraction(sp.cancel((R - E) / (1 + t)))[1], r, s, nu).is_monomial
        cii = sp.cancel(E.subs(s, r) - Eq) == 0
        if ci or cii:
            hits.append(((al, be, ga), sg, ci, cii))
    print(f"   条件を満たす偶奇の規則 ((α,β,γ), 符号, (i), (ii))：{hits}", flush=True)
    # t = −1 での比較（L だけ）
    I = sp.I
    Rm1 = sp.factor(sp.cancel(R.subs(s, I)))
    base = sp.factor(sp.cancel((keep * flipped / den).subs(s, I)))
    print(f"   t = −1 で必要な L の値 = {sp.factor(sp.cancel(Rm1 / base))}；L(−1) = {sp.factor(sp.expand(L.subs(s, I)))}", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]))
