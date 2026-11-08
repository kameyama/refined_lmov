r"""保存済みの大きな N の 3 箱の g（largeN_3box_m{m}.pkl）を U = t², t³ で固定 N の g（T2_{2m+1}_fixedN_3box_g.pkl、
rk_k5_torus_fixedN_hatted.py）と全係数で照合する（progress.md 続き 139。組み立ての後に固定 N の g を作った場合に使う）。
使い方：python rk_c_3box_check_fixedN.py m"""
import sys, pickle
sys.path.insert(0, "src")
import sympy as sp
from rcsknots.characters import chi
from rcsknots.partitions import partitions

m = int(sys.argv[1])
r, s, u = sp.symbols("r s u", positive=True)
g = {eval(k): sp.sympify(v, locals={"r": r, "s": s, "u": u}) for k, v in pickle.load(open(f"data/refined_kauffman/daha/largeN_3box_m{m}.pkl", "rb"))["g"].items()}
fx = pickle.load(open(f"data/refined_kauffman/T2_{2 * m + 1}_fixedN_3box_g.pkl", "rb"))
for N in (2, 3):
    ok = {}
    for R_ in partitions(3):
        gR = sum(chi(R_, rho) * g[rho] for rho in partitions(3))
        ok[R_] = sp.cancel(gR.subs(u, s**(2 * N)) - sp.sympify(fx[N]["g_schur"][str(R_)], locals={"r": r, "s": s})) == 0
    print(f"T(2,{2 * m + 1})、U = t^{N}：g_R（|R| = 3）が固定 N の g と一致 {ok}", flush=True)
