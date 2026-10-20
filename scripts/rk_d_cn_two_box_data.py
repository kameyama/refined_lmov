r"""C_N の 2 箱の DAHA-Jones と P の係数を N ごとに計算して保存する（案 D の 2 箱、progress.md 続き 52）。
使い方：python scripts/rk_d_cn_two_box_data.py N1,N2,... m [sp]   （sp を付けると u = t² に特殊化）
出力：data/refined_kauffman/daha/cn_twobox_N{N}_m{m}[_sp].pkl：{"Z": {λ: str}, "P": {λ: str}, "coef": {"a1","a0","c0","l2","l11","l0"}}。"""
import sys, os, pickle, time
sys.path.insert(0, "src")
from rcsknots.daha_cn import daha_jones_two_box
from rcsknots.flint_rf3 import S

OUT = "data/refined_kauffman/daha"

if __name__ == "__main__":
    Ns = [int(x) for x in sys.argv[1].split(",")]
    m = int(sys.argv[2])
    sp_ = len(sys.argv) > 3 and sys.argv[3] == "sp"
    for N in Ns:
        p = os.path.join(OUT, f"cn_twobox_N{N}_m{m}{'_sp' if sp_ else ''}.pkl")
        if os.path.exists(p):
            continue
        t0 = time.time()
        out, coef = daha_jones_two_box(N, m, uhalf=(S**2 if sp_ else None))
        a1, a0 = coef["P2"]
        pickle.dump({"Z": {str(l): str(v[0].to_sympy()) for l, v in out.items()}, "P": {str(l): str(v[1].to_sympy()) for l, v in out.items()},
                     "coef": {"a1": str(a1.to_sympy()), "a0": str(a0.to_sympy()), "c0": str(coef["P11"].to_sympy()),
                              "l2": str(coef["eig"]["2"].to_sympy()), "l11": str(coef["eig"]["11"].to_sympy()), "l0": str(coef["eig"]["0"].to_sympy())}},
                    open(p, "wb"))
        print(f"N={N} m={m} sp={sp_}: {time.time()-t0:.1f}s", flush=True)
