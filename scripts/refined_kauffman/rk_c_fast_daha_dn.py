r"""タスク C の高速版：scripts/rk_c_daha_dn.py と同じ Z = {γ̂(K_λ)}（D_N、γ = τ_+^m τ_−^2、T(2m+1, 2)、非 reduced、鉛直 framing）を
src/rcsknots/daha_fast.py で計算し、同じ形式で data/refined_kauffman/daha/Z_N{N}_m{m}_{λ}.pkl に保存する（保存済みは飛ばす）。
同じ N の λ は途中結果を共有する。
使い方：python scripts/rk_c_fast_daha_dn.py m Nmin Nmax λ...（λ は "1", "11", "2", "3", "21", "111" など）。
検証（2026-10-19）：保存済みの Z（現行版、m = 1 は N ≤ 9、m = 2 は N ≤ 6）と全係数一致（docs/refined_kauffman/daha_fast.md）。"""
import sys, time, pickle, os, resource
sys.path.insert(0, "src")
from rcsknots.refined_kauffman.daha_dn import DN
from rcsknots.refined_kauffman.daha_fast import Fast, daha_jones_fast, _with_rho2

OUT = "data/refined_kauffman/daha"

if __name__ == "__main__":
    m, Nmin, Nmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    lams = [tuple(int(c) for c in a) for a in sys.argv[4:]]
    for N in range(Nmin, Nmax + 1):
        F = None
        for lam in lams:
            path = os.path.join(OUT, f"Z_N{N}_m{m}_{''.join(map(str, lam))}.pkl")
            if os.path.exists(path):
                print(f"N={N} λ={lam}: 保存済み", flush=True)
                continue
            if F is None:
                F = Fast(_with_rho2(DN(N)), m)
            t0 = time.time()
            Z, K = daha_jones_fast(N, lam, m, F=F)
            pickle.dump({"Z": str(Z.to_sympy()), "K": str(K.to_sympy()), "N": N, "m": m, "lam": lam}, open(path, "wb"))
            print(f"N={N} λ={lam}: {time.time()-t0:.1f}s、最大 RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.0f}MB", flush=True)
    print("完了", flush=True)
