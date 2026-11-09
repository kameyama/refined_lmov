r"""代数的なケーブル K1・K3 の 2 箱の DAHA-Jones を、companion の段を約分なしの版（daha_iter.iter_jones_dn_lp・iter_adj_gl_lp）で計算して保存する
（progress.md 続き 141）。companion は K1・K3 とも trefoil なので、基底の値 V を N ごとに保存して共有する。
  dn Nmin Nmax [K…]：D_N の (2), (1,1) の Z, K を Z_N{N}_{K}_{λ}.pkl（rk_c_cable_daha.py と同じ形式）に。基底の値は Vdn_N{N}_m1.pkl と
               Vdn_N{N}_m1_{ν}.pkl（下の dnV）から読み、無いものは計算して Vdn_N{N}_m1.pkl に足す。K を並べるとそのケーブルだけ。
  dnV N ν…：基底の値 {γ̂(H_ν)} を ν ごとに計算して Vdn_N{N}_m1_{ν}.pkl に保存する（大きな N でメモリを抑えるため 1 プロセス 1 個）。
  gl Nmin Nmax：GL_N の随伴 (1;1) を W = 鏡像(Z / M)、M = T_adj^{−F} = q^{−F} t^{−F(N−1)}（F = p₂q₂、【仮説】トーラス結び目の rs を F に替えたもの。
               普通の表現での規則の確認は rk_c_cable_gl_rule.py）で Wadj_N{N}_{K}.pkl（rk_c_fast_adj.py と同じ形式）に。
               基底の値は |μ| = |ν| = 2 を保存済みの V22_N{N}_m1.pkl（rk_c_composite22_run.py）から読み、残りを計算する。
保存済みは飛ばす。最後に「done」。"""
import sys, os, time, pickle, resource
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rcsknots.refined_kauffman.flint_rf import RF, CTX, from_sympy
from rcsknots.refined_kauffman.daha_iter import iter_jones_dn_lp, iter_adj_gl_lp
import rk_torus_tag as TG

OUT = "data/refined_kauffman/daha"
CABLES = ("K1", "K3")
LAMS = ((2,), (1, 1))
rr, ss = sp.symbols("r s", positive=True)


def rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30


def ser(x):
    return ({tuple(e): int(c) for e, c in x.n.to_dict().items()}, {tuple(e): int(c) for e, c in x.d.to_dict().items()})


def deser(p):
    return RF(CTX.from_dict(p[0]), CTX.from_dict(p[1]), normalize=False)


def load_V(path):
    return {k: deser(v) for k, v in pickle.load(open(path, "rb")).items()} if os.path.exists(path) else {}


def vpath_nu(N, nu):
    return os.path.join(OUT, f"Vdn_N{N}_m1_{''.join(map(str, nu)) or '0'}.pkl")


def run_dnV(N, nus):
    from rcsknots.refined_kauffman.daha_iter import dn_basis_lp
    for nu in nus:
        p = vpath_nu(N, nu)
        if os.path.exists(p):
            continue
        t0 = time.time()
        V = dn_basis_lp(N, [nu], lowmem=True)
        pickle.dump(ser(V[nu]), open(p, "wb"))
        print(f"N={N} ν={nu}: {time.time()-t0:.1f}s、最大 RSS {rss():.2f} GB", flush=True)
    print("done", flush=True)


def run_dn(Nmin, Nmax, cables=CABLES):
    for N in range(Nmin, Nmax + 1):
        paths = {(c, l): os.path.join(OUT, f"Z_N{N}_{c}_{''.join(map(str, l))}.pkl") for c in cables for l in LAMS if len(l) <= N}
        if all(os.path.exists(p) for p in paths.values()):
            continue
        vpath = os.path.join(OUT, f"Vdn_N{N}_m1.pkl")
        V = load_V(vpath)
        for f in os.listdir(OUT):
            if f.startswith(f"Vdn_N{N}_m1_"):
                key = f[len(f"Vdn_N{N}_m1_"):-4]
                V[() if key == "0" else tuple(int(c) for c in key)] = deser(pickle.load(open(os.path.join(OUT, f), "rb")))
        n0 = len(V)
        for (c, l), p in paths.items():
            if os.path.exists(p):
                continue
            t0 = time.time()
            Z, K = iter_jones_dn_lp(N, l, c, V)
            if len(V) > n0:
                pickle.dump({k: ser(v) for k, v in V.items()}, open(vpath, "wb"))
                n0 = len(V)
            pickle.dump({"Z": str(Z.to_sympy()), "K": str(K.to_sympy()), "N": N, "cable": c, "lam": l}, open(p, "wb"))
            print(f"N={N} {c} λ={l}: {time.time()-t0:.1f}s、最大 RSS {rss():.2f} GB", flush=True)
    print("done", flush=True)


def run_gl(Nmin, Nmax):
    for N in range(Nmin, Nmax + 1):
        paths = {c: os.path.join(OUT, f"Wadj_N{N}_{c}.pkl") for c in CABLES}
        if all(os.path.exists(p) for p in paths.values()):
            continue
        V = {}
        v22 = os.path.join(OUT, f"V22_N{N}_m1.pkl")
        if os.path.exists(v22):
            V = {k: from_sympy(sp.sympify(v, locals={"r": rr, "s": ss})) for k, v in pickle.load(open(v22, "rb")).items()}
        for c, p in paths.items():
            if os.path.exists(p):
                continue
            t0 = time.time()
            Z = iter_adj_gl_lp(N, c, V).to_sympy()
            F = TG.rs(c)
            M = rr**(-2 * F) * ss**(-2 * F * (N - 1))               # T_adj^{−F}：(b,b) = 2、(b,ρ) = N − 1
            W = sp.cancel((Z / M).subs({rr: 1 / rr, ss: 1 / ss}, simultaneous=True))
            pickle.dump({"W": str(W), "N": N, "cable": c}, open(p, "wb"))
            print(f"N_A={N} {c} 随伴: {time.time()-t0:.1f}s、最大 RSS {rss():.2f} GB", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "dnV":
        run_dnV(int(sys.argv[2]), [tuple(int(c) for c in a) if a != "0" else () for a in sys.argv[3:]])
    elif mode == "dn":
        run_dn(int(sys.argv[2]), int(sys.argv[3]), tuple(sys.argv[4:]) or CABLES)
    else:
        run_gl(int(sys.argv[2]), int(sys.argv[3]))
