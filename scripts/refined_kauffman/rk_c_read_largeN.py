r"""タスク C ①：一般の N の簡約した不変量 R_λ(U) = Σ_j c_j(q, t) U^j を、いちばん大きな N の値から読み取る（progress.md 続き 32）。
c_j の s の次数の幅は 2N より小さいので、N が大きいと U^j = s^{2Nj} の項が s の次数の上で分離する。分離した塊ごとに c_j を読み、
読み取りに使っていないすべての N の値（R_true、scripts/rk_c_interp.py）と全係数で照合する。
使い方：python scripts/rk_c_read_largeN.py m λ...（例：1 3 21 111、2 1 2 11）。結果は data/refined_kauffman/daha/R_U_m{m}_{λ}.pkl。"""
import sys, pickle, os
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_c_interp as I

rr, ss = I.rr, I.ss


def terms(e):
    n, d = sp.fraction(sp.cancel(e))
    (dm, dc), = sp.Poly(d, rr, ss).terms()
    return {(a - dm[0], b - dm[1]): sp.Rational(c) / dc for (a, b), c in sp.Poly(sp.expand(n), rr, ss).terms()}


def read(m, lam):
    """窓 [j_min, j_max] を s の次数の最小・最大の N に対する傾きから決め（大きな N で 2N·j + 定数）、
    Newton の差分商で補間し（rk_c_interp.interpolate）、補間に使わない N で照合する。"""
    from rcsknots.refined_kauffman.flint_rf import from_sympy
    avail = [N for N in range(max(2, len(lam)), 30) if I.load(N, m, lam) is not None]          # D_N で ℓ(λ) ≤ N
    rng = {N: I.s_range(I.R_true(N, m, lam)) for N in avail[-3:]}
    a, b, c3 = avail[-3:]
    jmin = (rng[c3][0] - rng[b][0]) // 2
    jmax = (rng[c3][1] - rng[b][1]) // 2
    assert (rng[b][0] - rng[a][0]) // 2 == jmin and (rng[b][1] - rng[a][1]) // 2 == jmax, (rng, "傾きが一定でない")
    n = jmax - jmin + 1
    assert len(avail) >= n + 1, ("ノードが足りない", n, avail)
    nodes = avail[-n:]                                 # 大きな N をノードに
    c = I.interpolate(m, lam, nodes, jmin, jmax)
    ok = {N: (I.predict(c, N) - from_sympy(I.R_true(N, m, lam))).is_zero() for N in avail if N not in nodes}
    widths = {j: I.s_range(v.to_sympy()) if not v.is_zero() else None for j, v in c.items()}
    return {j: v.to_sympy() for j, v in c.items()}, nodes, widths, ok


if __name__ == "__main__":
    m = int(sys.argv[1])
    for lam in [tuple(int(x) for x in a) for a in sys.argv[2:]]:
        c, nodes, widths, ok = read(m, lam)
        print(f"m={m} λ={lam}：ノード N={nodes} で U^{min(c)}〜U^{max(c)} の {len(c)} 項を補間、c_j の s の幅 {widths}；他の N での照合 {ok}", flush=True)
        pickle.dump({j: str(v) for j, v in c.items()}, open(os.path.join(I.OUT, f"R_U_m{m}_{''.join(map(str, lam))}.pkl"), "wb"))
