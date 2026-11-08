r"""被覆の composite (α;β)（|α| = |β| = 2）の T(2m+1, 2) の DAHA-Jones の生の値を高速版で計算し、時間とメモリを測る（progress.md 続き 65、見積もり）。
  rk_c_fast_composite.comp_raw の (α;1) を (α;β) に広げたもの：
    X^{−1} 側（β = −ε_1 の軌道、Fast(sign = −1)）で P_β(X^{−1}) を H_ν^{(−)} = Π p_{ν_i}(X^{−1}) で展開し、g = Σ_ν y_ν γ̂'(H_ν^{(−)}(Y^{−1}))·1 を作る。
    X 側（Fast(sign = +1)）の連鎖の出発点を g にして P_α(X) を H_ν で展開し、{γ̂(P_α(X) P_β(X^{−1}))} を得る。
    |α| = |β| なので総次数 c = 0 で、ρ = (N−1, …, 0) の評価への換算因子は 1。
  検算：(α;β) = ((1);(1)) は gl_adj_fast + κ、((2);(1)) と ((1,1);(1)) は rk_c_fast_composite.comp_raw に一致すること。
  既約の composite への分解（P_α(X)P_β(X^{−1}) = Σ [α′;β′]）と規格化は本計算のときに加える（見積もりには不要）。
使い方：python rk_c_fast_composite22.py check | python rk_c_fast_composite22.py time N [m]（1 プロセスで (2;2), (2;11), (11;11) を順に計算し、各時間と最大 RSS を出す）"""
import sys, time, resource
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
from rcsknots.refined_kauffman.daha_fast import Fast, _with_rho2, solve_symmetric, h_basis, gl_adj_fast
from rcsknots.refined_kauffman.daha_dn import GLN, gl_P_poly
from rcsknots.refined_kauffman.flint_rf import R, S, conv


def inv_poly(f):
    return {tuple(-x for x in a): v for a, v in f.items()}


def comp_raw2(N, alpha, beta, m=1, k=2, cache=None):
    """{γ̂(P_α(X) P_β(X^{−1}))}（ρ = (N−1, …, 0) の評価）。cache（dict）を渡すと同じ N の Fast と X^{−1} 側の g を共有する。"""
    if cache is None:
        cache = {}
    if "D" not in cache:
        D = _with_rho2(GLN(N))
        cache.update(D=D, Fm=Fast(D, m, k, -1), g={})
    D, Fm = cache["D"], cache["Fm"]
    if beta not in cache["g"]:
        y = solve_symmetric(D, inv_poly(gl_P_poly(D, beta)), h_basis(Fm, [sum(beta)]))
        g = {}
        for nu, c in y.items():
            if c.is_zero():
                continue
            for a, v in Fm.get_sym(nu).items():
                g[a] = g.get(a, conv(0)) + c * v
        cache["g"][beta] = {a: v for a, v in g.items() if not v.is_zero()}
    Fp = Fast(D, m, k, 1)
    Fp.sym = {(): cache["g"][beta]}
    x = solve_symmetric(D, gl_P_poly(D, alpha), h_basis(Fp, [sum(alpha)]))
    v = sum((c * Fp.value(nu) for nu, c in x.items() if not c.is_zero()), conv(0))
    c = (sum(alpha) - sum(beta)) * (m * k + 1)
    return v * S**(-c * (N - 1))


def maxrss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30     # macOS は bytes


def check():
    import rk_c_fast_composite as C
    q, t = R**2, S**2
    for N in (3, 4):
        kappa = (1 - t**N) * (1 - q) / ((1 - t) * (1 - q * t**(N - 1)))
        ok = (comp_raw2(N, (1,), (1,)) - (gl_adj_fast(N, 1) + kappa)).is_zero()
        print(f"N={N}：(1;1) = gl_adj_fast + κ：{ok}", flush=True)
        for alpha in ((2,), (1, 1)):
            ok = (comp_raw2(N, alpha, (1,)) - C.comp_raw(N, alpha)[0]).is_zero()
            print(f"N={N}：({''.join(map(str, alpha))};1) = comp_raw：{ok}", flush=True)


def timing(N, m=1):
    cache = {}
    tot = time.time()
    for alpha, beta in (((2,), (2,)), ((2,), (1, 1)), ((1, 1), (1, 1))):
        t0 = time.time()
        v = comp_raw2(N, alpha, beta, m, cache=cache)
        nterm = len(v.to_sympy().as_numer_denom()[0].as_ordered_terms()) if N <= 5 else -1
        print(f"N={N} m={m} ({''.join(map(str, alpha))};{''.join(map(str, beta))})：{time.time()-t0:.1f}s、最大 RSS {maxrss_gb():.2f} GB"
              + (f"、分子の項数 {nterm}" if nterm >= 0 else ""), flush=True)
    print(f"N={N} m={m} 合計 {time.time()-tot:.1f}s、最大 RSS {maxrss_gb():.2f} GB", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "check":
        check()
    else:
        timing(int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 1)
