r"""拡張版の S1：D_N の Koornwinder 多項式 K_λ（|λ| ≤ 4、パラメータ (±1, ±q^{1/2})）を S_μ = Π S_{μ_i} の多項式として、
U = t^N を形式変数のまま評価の対称性で求める（src/rcsknots/kauffman_hopf2.koornwinder_general）。
λ ごとに data/refined_kauffman/koornwinder_DN/K_<λ>.pkl に保存し、再開時は読み込む。行列式は 8 並列。
完了の目印：ログの最後の行「保存 完了」。"""
import sys, time
sys.path.insert(0, "src")
from rcsknots.kauffman_hopf2 import koornwinder_general, koornwinder_small, _eval_S, _order

if __name__ == "__main__":
    t0 = time.time()
    Kf, Kr, cache = koornwinder_general(4, log=lambda m: print(f"{time.time()-t0:.0f}s {m}", flush=True),
                                        store="data/refined_kauffman/koornwinder_DN", workers=8)
    ks = koornwinder_small()
    print("最小版と一致:", Kf[(1, 1)].get((), 0) == ks["b"][1], Kf[(2,)].get((), 0) == ks["c"][2], flush=True)
    order = _order(4)
    bad = tot = 0
    for i, lam in enumerate(order):
        for nu in order[i + 1:]:
            tot += 1
            bad += (_eval_S(Kf[lam], nu, cache) * Kr[nu] != _eval_S(Kf[nu], lam, cache) * Kr[lam])
    print(f"評価の対称性（全 {tot} 組）: 不一致 {bad}", flush=True)
    print("保存 完了", flush=True)
