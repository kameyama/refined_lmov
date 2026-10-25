r"""4_1 の固定 N = 2（SO(4) = SU(2)²、被覆 U(3)）。入力は scripts/rk_k5_41_data.py の data（4_1 では data が既に行列式・共役で不変なので、
一様な規格化は data のまま：δ = 0 として扱う）。
  SO 側だけの検査（3 箱まで、SU(2) の [1], [2], [3] だけで組める）：framing の走査、f の整数性、SO の f̂ = M^{-1}(q^{-1}) f の Laurent 性と整数性、偶奇。
  被覆の 2 箱（随伴 [2,1] は KN 1703.05408 の候補で、確立していない）：随伴の framing の走査、f̂^{cov}、g・ĝ の Laurent 性・整数性・偶奇。
  被覆の 3 箱は [3,1]（composite (2;1)）が無いのでできない。N = 3 は SU(4) の [2,1,1], [3,2,1] が無いのでできない。"""
import sys, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
from fractions import Fraction
import sympy as sp
import rk_k5_trefoil_2box as T
T.KNOT = "4_1"
import rk_k5_trefoil_3box as X
X.DELTA = 0
from rk_k5_torus_fixedN_hatted import scan_rep, parity_even
W, ZERO, ONE = X.W_, X.ZERO, X.ONE
rr = sp.Symbol("r", positive=True)
N = 2


def so_part():
    fr = X.so_uniform(N)
    size = lambda rep: rep[0] + rep[1]
    scan = {rep: scan_rep(lambda f: X.so_Z(N, f), fr, rep, size(rep)) for rep in fr if size(rep) > 0}
    F = X.plog(X.so_Z(N, fr))
    res = {"framing の走査（予言 β = 0）": scan,
           "f 整数（|R| ≤ 3）": all(v.is_integer_laurent() for k in (1, 2, 3) for v in X.to_schur(F, k).values()),
           "SO f̂ Laurent（全方向）": {k: X.hatted_by_rho(F, k) for k in (2, 3)},
           "SO f̂_R 整数": all(v.is_integer_laurent() for k in (2, 3) for v in X.hatted(F, k).values()),
           "f_R が r → −r で不変": [parity_even(v) for k in (1, 2, 3) for v in X.to_schur(F, k).values()]}
    return F, res


def cover_part(F):
    NA = 2 * N - 1
    Fc = X.plog(X.cover_Z(N, {}, kmax=2, uniform=True), kmax=2)
    adj_scan = scan_rep(lambda m: X.cover_Z(N, m, kmax=2, uniform=True), {(1,): ONE}, (1,), 0)
    G = {k: F.get(k, ZERO) - X.conv(Fraction(1, 2)) * Fc.get(k, ZERO) for k in set(F) | set(Fc) if sum(k) <= 2}
    res = {"随伴（KN の候補）の framing の走査（予言 β = 0）": adj_scan,
           "f^cov Laurent（|R| ≤ 2）": all(v.is_laurent() for k in (1, 2) for v in X.to_schur(Fc, k).values()),
           "被覆 f̂ Laurent（|R| = 2）": X.hatted_by_rho(Fc, 2),
           "g Laurent・整数・2g 整数（|R| = 2）": [(v.is_laurent(), v.is_integer_laurent(), (2 * v).is_integer_laurent()) for v in X.to_schur(G, 2).values()],
           "g_(2)+g_(1,1) 整数": (X.to_schur(G, 2)[(2,)] + X.to_schur(G, 2)[(1, 1)]).is_integer_laurent(),
           "ĝ Laurent（|R| = 2）": X.hatted_by_rho(G, 2),
           "g_R が r → −r で不変（|R| ≤ 2）": [parity_even(v) for k in (1, 2) for v in X.to_schur(G, k).values()]}
    return res


if __name__ == "__main__":
    t0 = time.time()
    F, so = so_part()
    for k, v in so.items():
        print(f"4_1 N=2 SO 側：{k}: {v}", flush=True)
    print(f"  ({time.time()-t0:.0f}s)", flush=True)
    for k, v in cover_part(F).items():
        print(f"4_1 N=2 被覆（[2,1] は KN の候補）：{k}: {v}", flush=True)
    print(f"  ({time.time()-t0:.0f}s)", flush=True)
