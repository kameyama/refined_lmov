r"""C_N の DAHA-Jones JD(ω_1)（u は記号のまま）を N ごとに計算して保存する（案 D、progress.md 続き 50）。
出力：data/refined_kauffman/daha/cn_omega1_N{N}_m{m}.pkl（{"Z": str, "X1": str}、変数 r, s, u は q^{1/2}, t^{1/2}, u^{1/2}）。"""
import sys, os, pickle, time
sys.path.insert(0, "src")
from rcsknots.daha_cn import daha_jones_omega1

OUT = "data/refined_kauffman/daha"

if __name__ == "__main__":
    Ns = [int(x) for x in sys.argv[1].split(",")]
    ms = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [1, 2]
    for N in Ns:
        for m in ms:
            p = os.path.join(OUT, f"cn_omega1_N{N}_m{m}.pkl")
            if os.path.exists(p):
                continue
            t0 = time.time()
            Z, X1 = daha_jones_omega1(N, m)
            pickle.dump({"Z": str(Z.to_sympy()), "X1": str(X1.to_sympy())}, open(p, "wb"))
            print(f"N={N} m={m}: {time.time()-t0:.1f}s", flush=True)
