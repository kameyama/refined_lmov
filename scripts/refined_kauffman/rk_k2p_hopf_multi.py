r"""K2'(b)：Hopf 絡み目、成分 1 に probe 無限個（冪和を形式変数）、成分 2 に probe 1 個。
f_{R,b}（|R| ≤ A, b ≤ B）を Schur 基底と Macdonald P_R(y;t,q) 基底で求め、分母・整数性・Mariño との比較を出す。
出力：data/refined_kauffman/hopf_f_multi.pkl。
"""
import sys, pickle, time
from pathlib import Path
sys.path.insert(0, "src")
import sympy as sp
from rcsknots.refined_kauffman.kauffman_hopf import hopf_Ztilde_multi, log_coeffs_multi, plethystic_f, qs as q, ts as t, Us as U

A = int(sys.argv[1]) if len(sys.argv) > 1 else 4
B = int(sys.argv[2]) if len(sys.argv) > 2 else 2
t0 = time.time()
L = log_coeffs_multi(hopf_Ztilde_multi(A, B), A, B)
print(f"log Z~ 計算 {time.time()-t0:.0f}s", flush=True)
res = {}
for basis in ("schur", "macdonald"):
    f = plethystic_f(L, A, B, basis)
    res[basis] = {key: sp.cancel(v.as_expr()) for key, v in f.items()}
    print(f"{basis} 反転 {time.time()-t0:.0f}s", flush=True)
out = Path("data/refined_kauffman"); out.mkdir(parents=True, exist_ok=True)
pickle.dump({bs: {key: str(v) for key, v in fv.items()} for bs, fv in res.items()}, open(out / "hopf_f_multi.pkl", "wb"))

nu = sp.Symbol("nu")
zq = (q - 1) / sp.sqrt(q)
for basis, f in res.items():
    print(f"==== 基底 {basis}")
    for key in sorted(f, key=lambda k: (sum(k[0]), k[1], k[0])):
        R, b = key
        v = f[key]
        if v == 0:
            print(key, " f = 0"); continue
        num, den = sp.fraction(v)
        P = sp.Poly(sp.expand(num), q, t, U)
        print(key, " 分母 =", sp.factor(den), " 分子の係数が整数:", all(c.is_integer for c in P.coeffs()),
              " (q−1)|分子:", sp.rem(P, sp.Poly(q - 1, q, t, U)).is_zero, " 項数:", len(P.terms()))
    print("-- t = q（unrefined）、ν = U/q^{1/2}、z = q^{1/2} − q^{−1/2} で表示")
    for key in sorted(f, key=lambda k: (sum(k[0]), k[1], k[0])):
        if sum(key[0]) <= 2 and key[1] <= 1:
            v = sp.factor(sp.cancel(f[key].subs(t, q)).subs(U, nu * sp.sqrt(q)))
            print(key, ":", v)
