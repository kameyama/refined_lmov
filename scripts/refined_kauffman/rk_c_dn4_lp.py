r"""SO(2N) の 4 箱（D_N の DAHA-Jones、trefoil）の計算時間とメモリの測定（progress.md 続き 70）。
  旧版：daha_fast.daha_jones_fast（係数 RF、λ の間で Fast を共有）。新版：daha_fast_lp.dn_jones（係数 LP、最後の一段は評価の汎関数）。
使い方：
  python rk_c_dn4_lp.py verify N   … 新版と旧版の (Z, K) を |λ| ≤ 4 の全 λ で全係数比較し、両方の時間と最大 RSS
  python rk_c_dn4_lp.py time N [old|new]   … 片方だけの時間と最大 RSS（|λ| = 4 の 5 個と、それに要る |ν| ≤ 4 の基底の値）"""
import sys, time, resource
sys.path.insert(0, "src")
from rcsknots.refined_kauffman import daha_fast_lp as L
from rcsknots.refined_kauffman.daha_fast import Fast, _with_rho2, daha_jones_fast, partitions
from rcsknots.refined_kauffman.daha_dn import DN

LAM4 = [l for l in partitions(4)]


def rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30


def old(N, lams, m=1):
    F = Fast(_with_rho2(DN(N)), m)
    return {lam: daha_jones_fast(N, lam, m, F=F) for lam in lams}


def new(N, lams, m=1):
    return L.dn_jones(N, lams, m)


if __name__ == "__main__":
    mode, N = sys.argv[1], int(sys.argv[2])
    lams = [l for l in LAM4 if len(l) <= N]
    if mode == "verify":
        allp = [l for n in range(1, 5) for l in partitions(n) if len(l) <= N]
        t0 = time.time(); a = new(N, allp); t1 = time.time()
        print(f"N={N} 新版（|λ| ≤ 4 の {len(allp)} 個）：{t1-t0:.1f}s、最大 RSS {rss():.2f} GB", flush=True)
        b = old(N, allp); t2 = time.time()
        ok = all((a[l][0] - b[l][0]).is_zero() and (a[l][1] - b[l][1]).is_zero() for l in allp)
        print(f"N={N} 旧版：{t2-t1:.1f}s、最大 RSS {rss():.2f} GB；Z, K の全係数一致 {ok}", flush=True)
    else:
        which = sys.argv[3] if len(sys.argv) > 3 else "new"
        t0 = time.time()
        (new if which == "new" else old)(N, lams)
        print(f"N={N} {which}（|λ| = 4 の {len(lams)} 個）：{time.time()-t0:.1f}s、最大 RSS {rss():.2f} GB", flush=True)
