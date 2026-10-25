r"""K4'(a)：セクター別の被覆 f^{cov}_{(k,l),(m,n)}（kauffman_cover.cover_f_sym4）を計算し、検算（セクター和・Ĩ 対称性）して保存する。
出力：data/refined_kauffman/hopf_fcov_sector_A{A_}.pkl。使い方：python scripts/rk_k4_sector_cover.py 4"""
import sys, time, pickle; sys.path.insert(0, "src")
import sympy as sp
from rcsknots.refined_kauffman.kauffman_cover import cover_f_sym4, cover_f_sym, qs_, ss_, nus_
r, s, nu = sp.symbols("r s nu", positive=True); nut = sp.Symbol("nut")
t = s**2
t0 = time.time()
A_ = int(sys.argv[1]) if len(sys.argv) > 1 else 2   # A_ = 4 は約 6 分
f4 = cover_f_sym4(A_, A_)
print(f"cover_f_sym4({A_},{A_}): {time.time()-t0:.0f}s", flush=True)
sub = {qs_: r**2, ss_: s, nus_: nu}
f4 = {k: v.subs(sub) for k, v in f4.items()}
# 検算：Σ_{k+l=a, m+n=b} f4 = 1 変数の f^{cov}_{ab}
fc = {k: v.subs(sub) for k, v in cover_f_sym(A_, A_).items()}
for a in range(A_ + 1):
    for b in range(A_ + 1):
        if a + b == 0: continue
        tot = sum(v for (k, l, m, n), v in f4.items() if k + l == a and m + n == b)
        assert sp.cancel(tot - fc[(a, b)]) == 0, (a, b)
print("セクター和 = 1 変数の f^{cov}: ok")
# 対称性 (k,l,m,n) ↔ (l,k,n,m)
for key, v in f4.items():
    k, l, m, n = key
    assert sp.cancel(v - f4[(l, k, n, m)]) == 0, key
print("Ĩ 対称性 f(k,l,m,n) = f(l,k,n,m): ok")
pickle.dump({k: str(v) for k, v in f4.items()}, open(f"data/refined_kauffman/hopf_fcov_sector_A{A_}.pkl", "wb"))
g = {ab: sp.sympify(v, locals={"r": r, "s": s, "nu": nu}) for ab, v in pickle.load(open("data/refined_kauffman/hopf_g_alld.pkl", "rb")).items()}
def laurent_int(e):
    n, d = sp.fraction(sp.cancel(e))
    return sp.Poly(d, r, s, nu).is_monomial, all(c.is_integer for c in sp.Poly(sp.expand(n), r, s, nu).coeffs())
for (a, b), sec in (((2, 2), (1, 1, 1, 1)), ((2, 4), (1, 1, 2, 2)), ((4, 2), (2, 2, 1, 1)), ((4, 4), (2, 2, 2, 2))):
    if a > A_ or b > A_: continue
    P = sp.cancel((1 - t) * f4[sec])
    print(f"({a},{b}) 固定セクター {sec}: P = (1−t)f^cov_sec Laurent/整数 = {laurent_int(P)}, 項数 {len(sp.Poly(sp.expand(sp.fraction(P)[0]), r, s, nu).terms())}")
    E_m1 = sp.cancel(sp.cancel((2 * (1 + t) * g[(a, b)]).subs(nu, nut * s / r)).subs(s, sp.I))
    P_m1 = sp.cancel(sp.cancel(P.subs(nu, nut * s / r)).subs(s, sp.I))
    print(f"    P(t=−1) == E(t=−1): {sp.simplify(P_m1 - E_m1) == 0}")
    C = sp.cancel(g[(a, b)] - P / (2 * (1 + t)))
    print(f"    C = g − ½P/(1+t): Laurent/整数 = {laurent_int(C)}")
    Cw = sp.cancel(C.subs(nu, nut * s / r))
    print(f"    C は ν~ 固定で t 非依存: {not Cw.has(s)}")
    if laurent_int(C)[0]:
        print("    C =", sp.factor(Cw))
        print("    P (ν~) =", sp.factor(sp.cancel(P.subs(nu, nut * s / r))))
