r"""K4'(3)：Hopf の g_{ab} = f_{ab} − ½f^{cov}_{ab}（記号、a,b ≤ A）を求め、Laurent 性と交差帽展開を出す。
出力：data/refined_kauffman/hopf_g.pkl。"""
import sys, pickle, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rcsknots.refined_kauffman.kauffman_cover import cover_f_sym, qs_, ss_, nus_
from rk_k4_crosscap import to_VW, split_parity, in_G, V, W, nu, G
from rcsknots.refined_kauffman.kauffman_trefoil import r, s, U

A_ = int(sys.argv[1]) if len(sys.argv) > 1 else 3
t0 = time.time()
fcov = cover_f_sym(A_, A_)
print(f"f^cov 記号計算 {time.time()-t0:.0f}s", flush=True)
q_, t_, U_ = sp.symbols("q t U")
f = {ab: sp.sympify(v, locals={"q": q_, "t": t_, "U": U_}) for ab, v in pickle.load(open("data/refined_kauffman/hopf_f.pkl", "rb")).items()}
Q, T, Ab = sp.symbols("Q T A")   # DGR 太字 𝐪, 𝐭, 𝐚（𝐚 = ν）
g_all = {}
for a in range(0, A_ + 1):
    for b in range(0, A_ + 1):
        if a + b == 0 or (a, b) not in fcov: continue
        fo = f[(a, b)].subs({q_: r**2, t_: s**2, U_: U})                                   # (r, s, U)
        fc = fcov[(a, b)].subs({qs_: r**2, ss_: s, nus_: U / s})                            # ν = U/s
        g = sp.cancel(sp.together(fo - fc / 2))
        num, den = sp.fraction(g)
        laurent = sp.Poly(den, r, s, U).is_monomial
        g_all[(a, b)] = g
        print(f"\n=== ({a},{b}): Laurent = {laurent}", flush=True)
        if not laurent or g == 0:
            print("   g =", sp.factor(g)); continue
        e = to_VW(g); E, O = split_parity(e)
        P1 = in_G(E)
        print("  c=1 部分（G, w=W², ν）:", P1)
        if O != 0:
            for xname, X in (("z_q", V * W - 1 / (V * W)), ("z_t", V / W - W / V), ("w−1/w", W**2 - W**-2)):
                quo = sp.cancel(O / X); n_, d_ = sp.fraction(quo)
                if not sp.Poly(d_, V, W, nu).is_monomial:
                    print(f"  X = {xname}: 割り切れない"); continue
                P2 = in_G(sp.expand(quo))
                print(f"  X = {xname}: c=2 部分（G, w, ν）: {P2 if P2 is not None else '非対称（G で書けない）'}")
        gb = sp.expand(sp.cancel(g.subs({r: -T * Q, U: Ab * Q}).subs(s, Q)))
        coeffs = sp.Poly(sp.expand(gb * (Ab * Q * T)**20), Ab, Q, T).coeffs()
        print("  DGR 変数:", gb)
        print("  DGR 変数の係数の符号: 正", sum(1 for c in coeffs if c > 0), " 負", sum(1 for c in coeffs if c < 0), " 整数:", all(c.is_integer for c in coeffs))
pickle.dump({ab: str(v) for ab, v in g_all.items()}, open("data/refined_kauffman/hopf_g.pkl", "wb"))
