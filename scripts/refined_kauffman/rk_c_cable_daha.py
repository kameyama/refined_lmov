r"""代数的なケーブル（K1 = trefoil の (2,13)-ケーブル、K3 = (2,15)-ケーブル）の D_N の DAHA-Jones を計算し、
data/refined_kauffman/daha/Z_N{N}_{cable}_{λ}.pkl（rk_c_fast_daha_dn.py と同じ形式）に保存する（保存済みは飛ばす。progress.md 続き 140）。
計算は rcsknots.refined_kauffman.daha_iter.iter_jones_dn_fast（CD14 の iterated の定義、companion の段は高速版）。
使い方：python rk_c_cable_daha.py K1 Nmin Nmax λ...（λ は "1", "2", "11" など）"""
import sys, time, pickle, os, resource
sys.path.insert(0, "src")
from rcsknots.refined_kauffman.daha_dn import DN
from rcsknots.refined_kauffman.daha_fast import Fast, _with_rho2
from rcsknots.refined_kauffman.daha_iter import iter_jones_dn_fast, CABLES

OUT = "data/refined_kauffman/daha"

if __name__ == "__main__":
    cable = sys.argv[1]
    Nmin, Nmax = int(sys.argv[2]), int(sys.argv[3])
    lams = [tuple(int(c) for c in a) for a in sys.argv[4:]]
    _, (m, k) = CABLES[cable]
    for N in range(Nmin, Nmax + 1):
        F = None
        for lam in lams:
            path = os.path.join(OUT, f"Z_N{N}_{cable}_{''.join(map(str, lam))}.pkl")
            if os.path.exists(path):
                continue
            if F is None:
                F = Fast(_with_rho2(DN(N)), m, k)
            t0 = time.time()
            Z, K = iter_jones_dn_fast(N, lam, cable, F=F)
            pickle.dump({"Z": str(Z.to_sympy()), "K": str(K.to_sympy()), "N": N, "cable": cable, "lam": lam}, open(path, "wb"))
            print(f"N={N} λ={lam}: {time.time()-t0:.1f}s、最大 RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.0f}MB", flush=True)
    print("done", flush=True)
