r"""タスク C の被覆側：GL_N の DAHA-Jones（src/rcsknots/daha_dn.py の GLN）で随伴 composite (1;1) の T(2m+1, 2) の値を N ごとに求め、
一様な規格化の W^{true}_adj(N) = 鏡像(Z / M)、M = q^{−3(b,b)} t^{−2n(b,ρ)/2·…}（n = 2m+1、下の mono）を
data/refined_kauffman/daha/Wadj_N{N}_m{m}.pkl に保存する。
単項式 M の規則（T(3,2)）：M_b = q^{−3(b,b)} t^{−6(b,ρ) + 2|b|(N−1)}、ρ = (N−1, …, 0)。普通の表現 [1], [2], [1,1], [3], [2,1], [1,1,1]（N = 2〜5）と
随伴（N_A = 3, 5、固定 N の data）で確認した（2026-10-19）。"""
import sys, time, pickle, os
sys.path.insert(0, "src")
import sympy as sp
from rcsknots.refined_kauffman.daha_dn import GLN, gl_adj_gamma

rr, ss = sp.symbols("r s", positive=True)
OUT = "data/refined_kauffman/daha"


def mono(N, b, m=1):
    assert m == 1
    rho = [N - 1 - i for i in range(N)]
    bb = sum(x * x for x in b); br = sum(x * y for x, y in zip(b, rho)); nb = sum(b)
    return rr**(-6 * bb) * ss**(2 * (-6 * br + 2 * nb * (N - 1)))


if __name__ == "__main__":
    m, Nmin, Nmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    for N in range(Nmin, Nmax + 1):
        path = os.path.join(OUT, f"Wadj_N{N}_m{m}.pkl")
        if os.path.exists(path):
            continue
        t0 = time.time()
        D = GLN(N)
        Z = D.evaluate(gl_adj_gamma(D, m, 2)).to_sympy()
        b = [1] + [0] * (N - 2) + [-1]
        W = sp.cancel((Z / mono(N, b, m)).subs({rr: 1 / rr, ss: 1 / ss}, simultaneous=True))
        pickle.dump({"W": str(W), "N": N, "m": m}, open(path, "wb"))
        print(f"N={N}: {time.time()-t0:.1f}s", flush=True)
    print("完了", flush=True)
