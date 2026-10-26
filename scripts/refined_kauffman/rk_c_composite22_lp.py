r"""composite (α;β)（|α| = |β| = 2）の生の値を、Laurent 多項式の係数（daha_fast_lp、約分なし）と評価の汎関数で計算する（progress.md 続き 66）。
  基底の値 V[μ; ν]（μ, ν ∈ {(2), (1,1)}）を 4 つ計算し、(2;2), (2;11), (11;11) を RF の一次結合で作る。
使い方：
  python rk_c_composite22_lp.py verify N   … 旧版（rk_c_fast_composite22.comp_raw2、RF）と全係数で比べ、両方の時間と最大 RSS を出す
  python rk_c_composite22_lp.py time N [final]   … 新版だけ（final = functional（既定）か poly）"""
import sys, time, resource
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
from rcsknots.refined_kauffman import daha_fast_lp as L

P2 = [(2,), (1, 1)]
PAIRS = (((2,), (2,)), ((2,), (1, 1)), ((1, 1), (1, 1)))


def rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30


def new(N, final="functional"):
    t0 = time.time()
    V = L.basis_values(N, P2, P2, final=final)
    t1 = time.time()
    vals = {ab: L.composite_raw(N, *ab, V) for ab in PAIRS}
    t2 = time.time()
    print(f"N={N} 新版（final={final}）：基底の値 {t1-t0:.1f}s、一次結合 {t2-t1:.1f}s、計 {t2-t0:.1f}s、最大 RSS {rss():.2f} GB", flush=True)
    return vals


if __name__ == "__main__":
    mode, N = sys.argv[1], int(sys.argv[2])
    if mode == "verify":
        import rk_c_fast_composite22 as C
        vals = new(N)
        t0 = time.time(); cache = {}
        ok = [(C.comp_raw2(N, a, b, cache=cache) - vals[(a, b)]).is_zero() for a, b in PAIRS]
        print(f"N={N} 旧版：計 {time.time()-t0:.1f}s、最大 RSS {rss():.2f} GB；全係数の一致 {ok}", flush=True)
    else:
        new(N, sys.argv[3] if len(sys.argv) > 3 else "functional")
