r"""タスク C ②：GL_N の随伴 composite (1;1) の T(2m+1, 2) の DAHA-Jones を高速版（src/rcsknots/daha_fast.gl_adj_fast）で計算し、
一様な規格化 W = 鏡像(Z / M)、M = q^{−n(b,b)} t^{−2n(b,ρ) + (m+1)|b|(N−1)}（n = 2m+1、b = (1,0,…,0,−1)、|b| = 0）で
data/refined_kauffman/daha/Wadj_N{N}_m{m}.pkl に保存する（scripts/rk_c_daha_cover.py と同じ形式。保存済みは飛ばす）。
単項式の規則は m = 1（rk_c_daha_cover.mono）と m = 2（普通の表現 9 例と随伴 U(3), U(5)、logs/refined_kauffman/c_gl_mono_m2.log）で確認。
使い方：python scripts/rk_c_fast_adj.py m Nmin Nmax"""
import sys, time, pickle, os
sys.path.insert(0, "src")
import sympy as sp
from rcsknots.refined_kauffman.daha_fast import gl_adj_fast

rr, ss = sp.symbols("r s", positive=True)
OUT = "data/refined_kauffman/daha"

if __name__ == "__main__":
    m, Nmin, Nmax = map(int, sys.argv[1:4])
    k = next((int(a[4:]) for a in sys.argv[4:] if a.startswith("--k=")), 2)      # γ = τ_+^m τ_−^k、T(mk+1, k)。k ≠ 2 は Wadj_N{N}_m{m}k{k}.pkl（続き 102）
    tag = f"m{m}" if k == 2 else f"m{m}k{k}"
    rs = (m * k + 1) * k                                       # T(r, s) の r·s（k = 2 では 2n、n = 2m + 1）
    for N in range(Nmin, Nmax + 1):
        path = os.path.join(OUT, f"Wadj_N{N}_{tag}.pkl")
        if os.path.exists(path):
            continue
        t0 = time.time()
        Z = gl_adj_fast(N, m, k).to_sympy()
        M = rr**(-rs * 2) * ss**(-2 * rs * (N - 1))          # T_adj^{−rs}：(b,b) = 2、(b,ρ) = N − 1
        W = sp.cancel((Z / M).subs({rr: 1 / rr, ss: 1 / ss}, simultaneous=True))
        pickle.dump({"W": str(W), "N": N, "m": m, "k": k}, open(path, "wb"))
        print(f"N={N}: {time.time()-t0:.1f}s", flush=True)
    print("完了", flush=True)
