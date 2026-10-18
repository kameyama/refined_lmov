r"""トーラス結び目の固定 N（N = 2, 3）の 3 箱を、一様な規格化（data × (t/q)^{δ|μ|/2}）で調べる。使い方：python scripts/rk_k5_torus_fixedN_hatted.py T2_5 2
  (1) 一様な規格化の比（共役・行列式のずらし）、(2) framing の走査：各表現の単項式を (t/q)^{β/2} に替えて Z（|Y| ≤ 3）が Laurent になる β
  （Z は各 W について 1 次なので、その表現の係数だけを足し引きする）、(3) hatted（SO の f̂、被覆の f̂ の対角と M⊗M、ĝ）の Laurent 性と整数性、
  (4) g_R の r → −r（q^{1/2} → −q^{1/2}）での不変性と ĝ_R(−r) = (−1)^{|R|−1} ĝ_R(r)。計算の本体は rk_k5_trefoil_3box.py。"""
import sys, time, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
from fractions import Fraction
import sympy as sp
import rk_k5_trefoil_2box as T
import rk_k5_trefoil_3box as X
W, ZERO, ONE = X.W_, X.ZERO, X.ONE
rr = sp.Symbol("r", positive=True)


def zdiff(a, b):
    return {k: a.get(k, ZERO) - b.get(k, ZERO) for k in set(a) | set(b)}


def scan_rep(build, base, key, size, width=8):
    """build(dict) → Zp。key の単項式を W^β（β = δ|μ| ± width）にしたとき Z が Laurent になる β。"""
    zero = dict(base); zero[key] = ZERO
    one = dict(base); one[key] = ONE
    Z0, Bk = build(zero), None
    Bk = zdiff(build(one), Z0)
    if all(v.is_laurent() for v in Bk.values()):
        return "極なし"
    c = X.DELTA * size
    return [b for b in range(c - width, c + width + 1) if all((Z0.get(k, ZERO) + W**b * Bk.get(k, ZERO)).is_laurent() for k in set(Z0) | set(Bk))]


def parity_even(e):
    s_ = e.to_sympy()
    return sp.cancel(s_ - s_.subs(rr, -rr)) == 0


if __name__ == "__main__":
    KNOT, DELTA = sys.argv[1], int(sys.argv[2])
    T.KNOT, X.DELTA = KNOT, DELTA
    out = {}
    pairs = [p for p in X.UNIFORM_PAIRS]
    ratios = {}
    for n, a, b in pairs:
        try:
            ratios[f"U({n}) {a}/{b}"] = X.uniform_ratio(n, a, b) == ONE
        except FileNotFoundError:
            pass
    print(f"{KNOT}（δ = {DELTA}）一様な規格化の比が 1: {ratios}", flush=True)
    for N in (2, 3):
        t0 = time.time()
        NA = 2 * N - 1
        fr = X.so_uniform(N)
        size = (lambda rep: sum(rep)) if N == 3 else (lambda rep: rep[0] + rep[1])
        so_scan = {rep: scan_rep(lambda f: X.so_Z(N, f), fr, rep, size(rep)) for rep in fr if size(rep) > 0}
        print(f"N={N} SO の framing の走査（予言は β = δ|μ|）: {so_scan}", flush=True)
        ords = [(1,), (1, 1), (2,), (3,), (2, 1), (1, 1, 1)]
        mo = {a: W**(DELTA * sum(a)) for a in ords}
        c_scan = {("普通", a): scan_rep(lambda m: X.cover_Z(N, {}, uniform=True, mo=m), mo, a, sum(a)) for a in ords}
        mc = {al: W**(DELTA * sum(X.comp_shift(al, NA))) for al in ((1,), (2,), (1, 1))}
        c_scan.update({("composite", al): scan_rep(lambda m: X.cover_Z(N, m, uniform=True), mc, al, sum(X.comp_shift(al, NA))) for al in mc})
        print(f"N={N} 被覆の framing の走査（予言は β = δ|μ|、composite は移した分割）: {c_scan}", flush=True)
        F = X.plog(X.so_Z(N, fr))
        Fc = X.plog(X.cover_Z(N, {}, uniform=True))
        G = {k: F.get(k, ZERO) - X.conv(Fraction(1, 2)) * Fc.get(k, ZERO) for k in set(F) | set(Fc)}
        F2c = X.plog2(X.cover_Z2(N, {}, uniform=True))
        lau = {"SO f̂": all(all(X.hatted_by_rho(F, k).values()) for k in (2, 3)),
               "被覆 f̂（対角）": all(all(X.hatted_by_rho(Fc, k).values()) for k in (2, 3)),
               "被覆 f̂（M⊗M）": all(all(X.hatted2_by_rho(F2c, *ij).values()) for ij in ((2, 0), (1, 1), (3, 0), (2, 1))),
               "ĝ": all(all(X.hatted_by_rho(G, k).values()) for k in (2, 3))}
        integ = {"f（SO）整数": all(v.is_integer_laurent() for k in (1, 2, 3) for v in X.to_schur(F, k).values()),
                 "SO f̂_R 整数": all(v.is_integer_laurent() for k in (2, 3) for v in X.hatted(F, k).values()),
                 "被覆 f̂_R 整数": all(v.is_integer_laurent() for k in (2, 3) for v in X.hatted(Fc, k).values()),
                 "g_R 整数（|R| = 1, 3）": all(v.is_integer_laurent() for k in (1, 3) for v in X.to_schur(G, k).values()),
                 "2g_R 整数（|R| = 2）": all((2 * v).is_integer_laurent() for v in X.to_schur(G, 2).values()),
                 "g_(2)+g_(1,1) 整数": (X.to_schur(G, 2)[(2,)] + X.to_schur(G, 2)[(1, 1)]).is_integer_laurent(),
                 "ĝ_R 整数（|R| = 3）": all(v.is_integer_laurent() for v in X.hatted(G, 3).values()),
                 "2ĝ_R 整数（|R| = 2）": all((2 * v).is_integer_laurent() for v in X.hatted(G, 2).values())}
        par = {"g_R が r → −r で不変": all(parity_even(v) for k in (1, 2, 3) for v in X.to_schur(G, k).values()),
               "ĝ_R(−r) = (−1)^(|R|−1) ĝ_R(r)": all(sp.cancel(v.to_sympy() - (-1)**(k - 1) * v.to_sympy().subs(rr, -rr)) == 0 for k in (2, 3) for v in X.hatted(G, k).values())}
        print(f"N={N} hatted の Laurent 性: {lau}", flush=True)
        print(f"N={N} 整数性: {integ}", flush=True)
        print(f"N={N} 偶奇: {par}", flush=True)
        out[N] = {"g_schur": {str(R_): str(v.to_sympy()) for k in (1, 2, 3) for R_, v in X.to_schur(G, k).items()},
                  "ghat_schur": {str(R_): str(v.to_sympy()) for k in (2, 3) for R_, v in X.hatted(G, k).items()}}
        print(f"  N={N}: {time.time() - t0:.0f}s", flush=True)
    pickle.dump(out, open(f"data/refined_kauffman/{KNOT}_fixedN_3box_g.pkl", "wb"))
