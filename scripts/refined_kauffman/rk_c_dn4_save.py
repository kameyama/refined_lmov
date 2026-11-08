r"""SO(2N) の 4 箱（D_N の DAHA-Jones、trefoil、|λ| = 4）の Z, K を約分なしの版（daha_fast_lp.dn_jones）で計算し、
rk_c_fast_daha_dn.py と同じ形式で data/refined_kauffman/daha/Z_N{N}_m{m}_{λ}.pkl に保存する（保存済みは飛ばす）。
使い方：python rk_c_dn4_save.py m Nmin Nmax"""
import sys, os, time, pickle, resource
sys.path.insert(0, "src")
from rcsknots.refined_kauffman import daha_fast_lp as L
from rcsknots.refined_kauffman.daha_fast import partitions

OUT = "data/refined_kauffman/daha"

if __name__ == "__main__":
    m, Nmin, Nmax = map(int, sys.argv[1:4])
    for N in range(Nmin, Nmax + 1):
        lams = [l for l in partitions(4) if len(l) <= N]
        paths = {l: os.path.join(OUT, f"Z_N{N}_m{m}_{''.join(map(str, l))}.pkl") for l in lams}
        todo = [l for l in lams if not os.path.exists(paths[l])]
        if not todo:
            continue
        t0 = time.time()
        res = L.dn_jones(N, todo, m)
        for l in todo:
            Z, K = res[l]
            pickle.dump({"Z": str(Z.to_sympy()), "K": str(K.to_sympy()), "N": N, "m": m, "lam": l}, open(paths[l], "wb"))
        print(f"N={N}：{len(todo)} 個、{time.time()-t0:.1f}s、最大 RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**30:.2f} GB", flush=True)
    print("done", flush=True)
