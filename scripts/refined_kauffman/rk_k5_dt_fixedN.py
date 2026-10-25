r"""二重ツイスト結び目 K_{m,n} の固定 N = 2（SO(4) = SU(2)²）の SO 側だけの検査（3 箱まで、SU(2) の [1], [2], [3] だけで組める）。
使い方：python scripts/rk_k5_dt_fixedN.py 名前 m n（data は scripts/rk_k5_41_data.py で作る）。
  (1) [1] の superpolynomial（KNTZ、太字 (𝐚, 𝐪, 𝐭)）の δ 次数 𝐚 + 𝐪/2 − 𝐭 の集合（薄い結び目なら一つ）。
  (2) data の行列式表現 R_[1,1]（U(2)）= (q/t)^c の c、U(2) の [2,2]・U(3) の [1,1,1] でも同じ c か。一様な規格化は data × (t/q)^{c|μ|/2}。
  (3) 一様な規格化での framing の走査（予言 β = c|μ|）、f の整数性、SO の f̂ = M^{-1}(q^{-1}) f の Laurent 性と整数性、偶奇。"""
import sys, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman"); sys.path.insert(1, "scripts/framed_knot")
from fractions import Fraction
import sympy as sp
import rk_k5_trefoil_2box as T
import rk_k5_trefoil_3box as X
from rk_k5_torus_fixedN_hatted import scan_rep
from x1_check import triply, bold_dict
W, ZERO, ONE = X.W_, X.ZERO, X.ONE
rr, ss = sp.symbols("r s", positive=True)


def delta_grading(m, n):
    P = bold_dict(triply(1, 1, m, n))
    return sorted({Fraction(a) + Fraction(q, 2) - Fraction(t) for (a, q, t) in P})


def offset_c():
    out = {}
    for n_, lab in ((2, "11"), (2, "22"), (3, "111")):
        ratio = sp.simplify(T.Rdata(lab, n_))
        # ratio = (r/s)^{2 c |μ|/2·...}：(q/t)^{c n_ k/2} の形かを見る（k は列の数）
        k = int(lab[0])
        e = sp.Poly(sp.expand(ratio * (s_ := sp.Symbol("s", positive=True))**0), rr, ss) if False else None
        num, den = sp.fraction(sp.cancel(ratio))
        ok = sp.Poly(num, rr, ss).is_monomial and sp.Poly(den, rr, ss).is_monomial
        out[f"U({n_}) [{lab}]"] = (sp.factor(ratio), ok)
    return out


def so_checks(c):
    X.DELTA = c
    X.RP.cache_clear()
    N = 2
    fr = X.so_uniform(N)
    size = lambda rep: rep[0] + rep[1]
    scan = {rep: scan_rep(lambda f: X.so_Z(N, f), fr, rep, size(rep)) for rep in fr if size(rep) > 0}
    F = X.plog(X.so_Z(N, fr))
    par = lambda sub: all(sp.cancel(v.to_sympy() - v.to_sympy().subs(sub, simultaneous=True)) == 0 for k in (1, 2, 3) for v in X.to_schur(F, k).values())
    return {"framing の走査（予言 β = c|μ|）": scan,
            "f 整数（|R| ≤ 3）": all(v.is_integer_laurent() for k in (1, 2, 3) for v in X.to_schur(F, k).values()),
            "SO f̂ Laurent（全方向）": all(all(X.hatted_by_rho(F, k).values()) for k in (2, 3)),
            "SO f̂_R 整数": all(v.is_integer_laurent() for k in (2, 3) for v in X.hatted(F, k).values()),
            "f が r → −r で不変": par({rr: -rr}), "f が (r,s) → (−r,−s) で不変": par({rr: -rr, ss: -ss})}


if __name__ == "__main__":
    knot, m, n = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    T.KNOT = knot
    t0 = time.time()
    print(f"{knot} = K_{{{m},{n}}}：[1] の δ 次数 {delta_grading(m, n)}", flush=True)
    oc = offset_c()
    print(f"{knot}：data の行列式表現の値 {oc}", flush=True)
    c = sp.Rational(sp.degree(sp.fraction(sp.cancel(T.Rdata('11', 2)))[0], rr) - sp.degree(sp.fraction(sp.cancel(T.Rdata('11', 2)))[1], rr), 2)
    print(f"{knot}：c = {c}（R_[1,1](U(2)) の r の次数の半分）", flush=True)
    if c.q == 1:
        for k, v in so_checks(int(c)).items():
            print(f"{knot} N=2 SO 側：{k}: {v}", flush=True)
    print(f"  ({time.time()-t0:.0f}s)", flush=True)
