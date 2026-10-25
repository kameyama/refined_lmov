r"""タスク C：D_N の DAHA-Jones（src/rcsknots/daha_dn.py）で T(2m+1, 2) の Z = {γ̂(K_λ)}（γ = τ_+^m τ_−^2、非 reduced、鉛直 framing）を
N ごとに計算して data/refined_kauffman/daha/Z_N{N}_m{m}_{λ}.pkl に保存する（保存済みは飛ばす）。
使い方：python scripts/rk_c_daha_dn.py m Nmin Nmax λ...（λ は "1", "11", "2" など）。
検証（2026-10-19）：N = 2, 3 で Z = T_λ^{−6} × 鏡像（一様な規格化の data の W_λ）、N = 4, 5 の [1] で JD = 鏡像（AS (6.15)）。"""
import sys, time, pickle, os, resource
sys.path.insert(0, "src")
from rcsknots.refined_kauffman.daha_dn import daha_jones

OUT = "data/refined_kauffman/daha"

if __name__ == "__main__":
    m, Nmin, Nmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    lams = [tuple(int(c) for c in a) for a in sys.argv[4:]]
    for N in range(Nmin, Nmax + 1):
        for lam in lams:
            path = os.path.join(OUT, f"Z_N{N}_m{m}_{''.join(map(str, lam))}.pkl")
            if os.path.exists(path):
                print(f"N={N} λ={lam}: 保存済み", flush=True)
                continue
            t0 = time.time()
            Z, K = daha_jones(N, lam, m)
            pickle.dump({"Z": str(Z.to_sympy()), "K": str(K.to_sympy()), "N": N, "m": m, "lam": lam}, open(path, "wb"))
            print(f"N={N} λ={lam}: {time.time()-t0:.1f}s、最大 RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.0f}MB", flush=True)
    print("完了", flush=True)
