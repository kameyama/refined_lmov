r"""trefoil の 4 箱の被覆の自由エネルギー f^{cov}（大きな N、rk_c_4box_largeN_g.py）を、固定の U(N_A) での独立な計算と照合する（progress.md 続き 139）。
固定 N_A（既定 5、SO(6) の U = t³ に対応、A = t^{N_A}）：
  ⟨P_a(x) P_b(x^{−1})⟩ を GL(N_A) の Laurent 多項式として作り、既約 P_ν（ν は優勢な重み）に分解する（daha_fast.solve_symmetric、
  P_ν = e_N^{−c} P_{ν + c^N}、c = max(0, −ν_N)）。⟨P_ν⟩ は KN の data（GN の再構成。DAHA・Pieri・補間を通らない）で
  W_ν = (t/q)^{|sh|/2} R_sh P_sh(t^ρ)、sh = ν + c^N（ν_N ≥ 0 は ν のまま。rk_k5_trefoil_3box.cover_parts の uniform と同じ規則）。
  Z~ = Σ_{λ,μ} (−1)^{|λ|+|μ|} ⟨P_λ P_μ^∨⟩ P_{λ'}(Y; t, q) P_{μ'}(Y; t, q)、f = Σ_d μ(d)/d ψ_d(log Z~)（次数 ≤ 4）。
KN の data が無い形（U(5) の 10 箱など。GN の再構成は 1 件 1 時間を超える）は、GL(N_A) の DAHA の固定 N の値（rk_c_fast_composite.mono の
  一様な規格化）で代える。両者は KN の data のある 7 つの形（[2,1]〜[2,2,2,2]）の U(5) で比 1（progress 続き 139）。
  したがってこの照合は、Pieri の係数・A の補間・大きな N の組み立てを独立に確かめるもので、DAHA そのものは共有する部分がある。
使い方：python rk_c_4box_cover_fixedN.py [N_A]"""
import sys, os, pickle, time
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman"); sys.path.insert(1, "scripts/framed_knot")
import sympy as sp
from rcsknots.refined_kauffman.flint_rf import RF, R, S, ONE, conv, CTX
from rcsknots.refined_kauffman.daha_fast import solve_symmetric
from rcsknots.refined_kauffman.daha_dn import GLN, gl_P_poly
from rcsknots.macdonald import macdonald_P
from rcsknots.partitions import partitions
import rk_k5_trefoil_2box as T

ZERO = conv(0)
KMAX = 4
PARTS = [()] + [lam for n in range(1, KMAX + 1) for lam in partitions(n)]
W_ = S / R


def inv(f):
    return {tuple(-x for x in a): v for a, v in f.items()}


def mul(f, g):
    out = {}
    for a, x in f.items():
        for b, y in g.items():
            k = tuple(p + w for p, w in zip(a, b))
            out[k] = out.get(k, ZERO) + x * y
    return {k: v for k, v in out.items() if not v.is_zero()}


def weights(N, total, top, bottom, high=None):
    """優勢な重み ν（ν_1 ≥ … ≥ ν_N、Σν = total、ν_1 ≤ top、ν_N ≥ −bottom）。high を与えると支配順序で ν ≤ high のものだけ
    （P_a(x) P_b(x^{−1}) の成分は最高の重み (a, 0, …, 0, −b の逆順) 以下）。"""
    out = []

    def rec(prefix, rest_n, hi):
        if rest_n == 0:
            if sum(prefix) == total and (high is None or all(sum(prefix[:i]) <= sum(high[:i]) for i in range(1, N + 1))):
                out.append(tuple(prefix))
            return
        for x in range(hi, -bottom - 1, -1):
            rec(prefix + [x], rest_n - 1, x)
    rec([], N, top)
    return out


_PW = {}


def P_weight(D, nu):
    if (D.N, nu) in _PW:
        return _PW[(D.N, nu)]
    c = max(0, -nu[-1])
    lam = tuple(x + c for x in nu)
    P = gl_P_poly(D, tuple(x for x in lam if x > 0))
    _PW[(D.N, nu)] = {tuple(x - c for x in a): v for a, v in P.items()}
    return _PW[(D.N, nu)]


def shape(nu):
    c = max(0, -nu[-1])
    return tuple(x + c for x in nu if x + c > 0)


SOURCE = {}


def R_daha(sh, N):
    """鏡像(Z_sh / (K_sh M_sh))（GL(N) の DAHA、一様な規格化）= (t/q)^{|sh|/2} R_sh(KN)。"""
    import rk_c_fast_composite as C
    from rcsknots.refined_kauffman.daha_fast import gl_jones_fast
    from rcsknots.refined_kauffman.flint_rf import from_sympy
    from rcsknots.macdonald import principal_specialization
    rr, ss = C.rr, C.ss
    K = from_sympy(sp.cancel(principal_specialization(sh, ss**(-2 * N), rr**(-2), ss**(-2))))
    b = list(sh) + [0] * (N - len(sh))
    e = (gl_jones_fast(N, sh, 1) / (K * C.mono(N, b, 1))).to_sympy()
    return from_sympy(sp.cancel(e.subs({rr: 1 / rr, ss: 1 / ss}, simultaneous=True)))


def W_shape(sh, NA):
    if not sh:
        return ONE
    lab = "".join(map(str, sh))
    from rcsknots.refined_kauffman.flint_rf import from_sympy
    if os.path.exists(f"data/refined_kauffman/torus/T2_3_{lab}.json") or os.path.exists(f"data/framed_knot/torus/T2_3_{lab}.json"):
        SOURCE[sh] = "KN"
        return W_**sum(sh) * from_sympy(T.Rdata(lab, NA)) * from_sympy(T.P_A(sh, NA))
    t0 = time.time()
    SOURCE[sh] = "DAHA"
    cp = f"data/cache/R_daha_GL{NA}_T2_3_{lab}.pkl"          # DAHA の値のキャッシュ（10 箱は 1 件 15〜45 分）
    if os.path.exists(cp):
        from rcsknots.refined_kauffman import flint_rf as FR
        nd, dd = pickle.load(open(cp, "rb"))
        Rv = FR.RF(CTX.from_dict(nd), CTX.from_dict(dd))
    else:
        Rv = R_daha(sh, NA)
        os.makedirs("data/cache", exist_ok=True)
        pickle.dump(({tuple(e): int(c) for e, c in Rv.n.to_dict().items()}, {tuple(e): int(c) for e, c in Rv.d.to_dict().items()}), open(cp, "wb"))
    v = Rv * from_sympy(T.P_A(sh, NA))
    print(f"   [{lab}] は GL({NA}) の DAHA で評価（{time.time()-t0:.0f}s）", flush=True)
    return v


def pair_values(NA):
    global T0
    T0 = time.time()
    D = GLN(NA)
    P = {lam: gl_P_poly(D, lam) for lam in PARTS}
    dec = {}
    # 荷電共役 ⟨P_a P_b^∨⟩ = ⟨P_b P_a^∨⟩（大きな N の組み立てと同じ仮定。(1;2) と (2;1) で確かめる）を使い、
    # 重みの持ち上げが小さい側（|a| ≥ |b|、|a| = |b| = 2 では (11;2) の代わりに (2;11)）だけを分解する（U(5) の (0;4) は 16 箱になる）
    todo = [(a, b) for a in PARTS for b in PARTS if sum(a) + sum(b) <= KMAX and len(a) <= NA and len(b) <= NA
            and (sum(a) > sum(b) or (sum(a) == sum(b) and a >= b))] + [((1,), (2,))]
    for a, b in todo:
        lhs = mul(P[a], inv(P[b]))
        high = tuple(list(a) + [0] * (NA - len(a) - len(b)) + [-x for x in reversed(b)])
        ws = weights(NA, sum(a) - sum(b), a[0] if a else 0, b[0] if b else 0, high)
        basis = {nu: P_weight(D, nu) for nu in ws}
        x = solve_symmetric(D, lhs, basis)
        # 検算：分解の和が左辺に戻る
        back = {}
        for nu, c in x.items():
            if not c.is_zero():
                for k, v in basis[nu].items():
                    back[k] = back.get(k, ZERO) + c * v
        assert all((back.get(k, ZERO) - lhs.get(k, ZERO)).is_zero() for k in set(back) | set(lhs)), (a, b)
        dec[(a, b)] = {nu: c for nu, c in x.items() if not c.is_zero()}
    shapes = sorted({shape(nu) for d in dec.values() for nu in d}, key=lambda s: (sum(s), s))
    print(f"   分解 {time.time() - T0:.0f}s、形 {len(shapes)} 個", flush=True)
    Wsh = {sh: W_shape(sh, NA) for sh in shapes}
    pair = {k: sum((c * Wsh[shape(nu)] for nu, c in d.items()), ZERO) for k, d in dec.items()}
    ok = (pair[((1,), (2,))] - pair[((2,), (1,))]).is_zero()
    print(f"   荷電共役 ⟨P_1 P_2^∨⟩ = ⟨P_2 P_1^∨⟩：{ok}", flush=True)
    assert ok
    for (a, b) in list(pair):
        pair.setdefault((b, a), pair[(a, b)])
    return pair, dec


def padd(A, B, c=ONE):
    out = dict(A)
    for k, v in B.items():
        out[k] = out.get(k, ZERO) + c * v
    return out


def pmul(A, B):
    out = {}
    for a, x in A.items():
        for b, y in B.items():
            if sum(a) + sum(b) <= KMAX:
                k = tuple(sorted(a + b, reverse=True))
                out[k] = out.get(k, ZERO) + x * y
    return out


def cover_Z(pair):
    PY = {(): {(): ONE}}
    for n in range(1, KMAX + 1):
        PY.update(macdonald_P(n, S**2, R**2, conv))
    T_ = lambda lam: tuple(sum(1 for x in lam if x > j) for j in range(lam[0])) if lam else ()
    Zp = {(): ONE}
    for lam in PARTS:
        for mu in PARTS:
            n = sum(lam) + sum(mu)
            if n == 0 or n > KMAX or (lam, mu) not in pair:
                continue
            Zp = padd(Zp, pmul(PY[T_(lam)], PY[T_(mu)]), conv((-1)**n) * pair[(lam, mu)])
    return Zp


def plog(Zp):
    X = {k: v for k, v in Zp.items() if k}
    Lg, Xn = {}, dict(X)
    for n in range(1, KMAX + 1):
        Lg = padd(Lg, Xn, conv(Fraction((-1)**(n + 1), n)))
        Xn = pmul(Xn, X)
    f = {}
    for d, mu in ((1, 1), (2, -1), (3, -1)):
        for k, v in Lg.items():
            if d * sum(k) <= KMAX:
                kk = tuple(sorted((d * x for x in k), reverse=True))
                f[kk] = f.get(kk, ZERO) + conv(Fraction(mu, d)) * v.adams(d)
    return f


def from_largeN(p, N):
    """flint_rf3 の (分子, 分母) の辞書 → u = s^{2N} を代入した flint_rf。"""
    def cv(d):
        out = {}
        for e, c in d.items():                          # u = s^{2N} で指数が重なる項は足す
            k = (e[0], e[1] + 2 * N * e[2])
            out[k] = out.get(k, 0) + c
        return CTX.from_dict({k: c for k, c in out.items() if c})
    n, d = cv(p[0]), cv(p[1])
    from rcsknots.refined_kauffman import flint_rf as FR
    return FR.RF(n, d)


def main(NA=5):
    assert NA % 2 == 1
    N = (NA + 1) // 2                           # SO(2N)、U = t^N、A = t^{2N−1}
    t0 = time.time()
    pair, dec = pair_values(NA)
    print(f"U({NA}) の ⟨P_λ P_μ^∨⟩：{len(pair)} 組 {time.time()-t0:.0f}s、KN の data {sum(v == 'KN' for v in SOURCE.values())} 形・DAHA {sum(v == 'DAHA' for v in SOURCE.values())} 形", flush=True)
    for k in (((2,), (2,)), ((1, 1), (1, 1)), ((2, 1), (1,)), ((1, 1, 1), (1,))):
        if k in dec:
            print(f"   分解 {k}：{sorted(dec[k])}", flush=True)
    Fc = plog(cover_Z(pair))
    print(f"被覆の f（固定 N）{time.time()-t0:.0f}s", flush=True)
    big = pickle.load(open("data/refined_kauffman/daha/largeN_4box_m1.pkl", "rb"))["fcov"]
    res = {}
    for k in sorted(Fc, key=lambda k: (sum(k), k)):
        v = from_largeN(big[str(k)], N)
        res[k] = (v - Fc[k]).is_zero()
    print(f"U = t^{N}（U({NA})）で大きな N の f^cov と一致：{res}", flush=True)
    return res


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
