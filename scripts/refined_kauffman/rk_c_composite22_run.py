r"""被覆の composite [α;β]（|α| = |β| = 2）の一般の N の簡約した値（progress.md 続き 67）。
  1. 生の値：基底の値 V[μ; ν]（src/rcsknots/refined_kauffman/daha_fast_lp.py）から Z_raw = {γ̂(P_α(X) P_β(X^{−1}))}。
  2. 既約への分解：P_α(X) P_β(X^{−1}) = P_b + c₁ P_adj + c₀（b = (α, 0, …, 0, −β の逆順)、P_adj = p_1(X) p_1(X^{−1}) − κ）。
     係数は二つの評価から解く【導出】：
       ev(f) = f(x)、x = (t^{−(N−1)}, …, t^{−1}, 1)（DAHA の評価 X_b ↦ t^{−(b, ρ)}、ρ = (N−1, …, 0)）。
       Macdonald の作用素 D = Σ_i Π_{j≠i}(t x_i − x_j)/(x_i − x_j) T_{q, x_i}（固有値 e_b = Σ_i q^{b_i} t^{N−i}）は、点 x では i = N の項だけが残り
       ev(D f) = [N]_t f(x′)、x′ = (t^{−(N−1)}, …, t^{−1}, q)。
       よって f(x) = K_b + c₁ K_adj + c₀、[N]_t f(x′) = e_b K_b + c₁ e_adj K_adj + c₀ [N]_t（K = ev(P)）。
     P_α(x) は P_α の冪和展開（macdonald_P）と p_d(x) = Σ x_i^d から、K_b は主特殊化の積公式（P_b = det^{−k} P_{b+k}）から。
     小さな N では多項式の分解（solve_symmetric）と照合する（check）。
  3. 規格化：Z_b = Z_raw − c₁ Z_adj − c₀、R_b(N) = 鏡像(Z_b / (K_b M_b))（M_b は rk_c_fast_composite.mono、(α;1) と同じ規則）。
使い方：
  python rk_c_composite22_run.py check            … 分解の係数と固有値の公式を N = 3, 4 で多項式と照合
  python rk_c_composite22_run.py run Nmin Nmax [m]  … R_b(N) を計算して保存（保存済みは飛ばす。最後に「done」。m は既定 1、T(2,5) は 2）
  python rk_c_composite22_run.py window [m]         … 保存済みの R_b(N) から U の冪の範囲（x = s^N の冪）を傾きで決め、補間して残りの N で照合"""
import sys, os, time, pickle, resource
from fractions import Fraction
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rcsknots.refined_kauffman import daha_fast_lp as L
from rcsknots.refined_kauffman.flint_rf import R, S, ONE, conv, from_sympy
from rcsknots.macdonald import macdonald_P, boxes, arm, leg

OUT = "data/refined_kauffman/daha"
P2 = [(2,), (1, 1)]
PAIRS = (((2,), (2,)), ((2,), (1, 1)), ((1, 1), (1, 1)))
tag = lambda a, b: f"{''.join(map(str, a))}_{''.join(map(str, b))}"
q, t = R**2, S**2
rr, ss = sp.symbols("r s", positive=True)


def rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30


def bvec(N, alpha, beta):
    return list(alpha) + [0] * (N - len(alpha) - len(beta)) + [-x for x in reversed(beta)]


def ps_rf(lam, N):
    """P_λ(1, t^{−1}, …, t^{−(N−1)}; q, t) = 主特殊化（a = t^{−N}, q^{−1}, t^{−1}）を RF で。"""
    a, qi, ti = S**(-2 * N), 1 / q, 1 / t
    out = ONE
    for i, j in boxes(lam):
        out = out * (ti**i - a * qi**j) / (1 - qi**arm(lam, i, j) * ti**(leg(lam, i, j) + 1))
    return out


def K_of(N, b):
    """ev(P_b)（b は降順の整数列）。P_b = det^{−k} P_{b+k}、ev(det) = t^{−N(N−1)/2}。"""
    k = -min(b)
    lam = tuple(x + k for x in b if x + k > 0)
    return S**(k * N * (N - 1)) * ps_rf(lam, N)


def e_of(N, b):
    return sum((q**b[i] * t**(N - 1 - i) for i in range(N)), conv(0))


_MP = {}


def P_at(lam, xs):
    """P_λ(x)（xs は RF の点、冪和展開から）。"""
    n = sum(lam)
    if n not in _MP:
        _MP[n] = macdonald_P(n, q, t, conv)
    tot = conv(0)
    for rho, c in _MP[n][tuple(lam)].items():
        v = ONE
        for d in rho:
            v = v * sum((x**d for x in xs), conv(0))
        tot = tot + c * v
    return tot


def decomposition(N, alpha, beta):
    """(c₁, c₀, K_b, K_adj) を二つの評価の式から。"""
    x = [S**(-2 * (N - 1 - i)) for i in range(N)]
    xp = x[:-1] + [q]
    f = lambda pt: P_at(alpha, pt) * P_at(beta, [1 / y for y in pt])
    kappa = (1 - t**N) * (1 - q) / ((1 - t) * (1 - q * t**(N - 1)))
    adj = lambda pt: sum(pt, conv(0)) * sum((1 / y for y in pt), conv(0)) - kappa
    b = bvec(N, alpha, beta)
    badj = [1] + [0] * (N - 2) + [-1]
    Kb, Kadj = K_of(N, b), adj(x)
    eN = e_of(N, [0] * N)
    eb, eadj = e_of(N, b), e_of(N, badj)
    # 検算：ev(D P_adj) = [N]_t P_adj(x′) = e_adj K_adj（固有値の式と点 x′ の規約の確認）
    assert (eN * adj(xp) - eadj * Kadj).is_zero(), "ev(D P_adj) の式が合わない"
    r1 = f(x) - Kb                      # c₁ K_adj + c₀
    r2 = eN * f(xp) - eb * Kb           # c₁ e_adj K_adj + c₀ e_N
    det = Kadj * eN - eadj * Kadj
    c1 = (r1 * eN - r2) / det
    c0 = r1 - c1 * Kadj
    return c1, c0, Kb, Kadj


def check():
    from rcsknots.refined_kauffman.daha_dn import GLN, gl_P_poly, gl_adj_poly
    from rcsknots.refined_kauffman.daha_fast import solve_symmetric, poly_mul
    for N in (3, 4):
        D = GLN(N)
        inv = lambda f: {tuple(-x for x in a): v for a, v in f.items()}
        for alpha, beta in PAIRS:
            if len(alpha) + len(beta) > N:
                continue
            t0 = time.time()
            b = bvec(N, alpha, beta)
            k = -min(b)
            Pb = {tuple(x - k for x in a): v for a, v in gl_P_poly(D, tuple(x + k for x in b if x + k > 0)).items()}
            lhs = poly_mul(gl_P_poly(D, alpha), inv(gl_P_poly(D, beta)))
            sol = solve_symmetric(D, lhs, {"b": Pb, "adj": gl_adj_poly(D), "0": {tuple([0] * N): ONE}})
            c1, c0, Kb, Kadj = decomposition(N, alpha, beta)
            ok = [sol["b"] == conv(1), (sol["adj"] - c1).is_zero(), (sol["0"] - c0).is_zero(),
                  (D.evaluate(Pb) - Kb).is_zero(), (e_of(N, b) * Kb - eval_D(D, Pb)).is_zero()]
            print(f"N={N} ({tag(alpha, beta)})：[P_b の係数 1, c₁, c₀, K_b, ev(D P_b) = e_b K_b] {ok}、c₀ = 0 {c0.is_zero()}（{time.time()-t0:.1f}s）", flush=True)


def eval_D(D, f):
    """ev(D f) を多項式から直接：[N]_t f(x′)（ここでは x′ の評価を多項式で行う）。"""
    N = D.N
    xp = [S**(-2 * (N - 1 - i)) for i in range(N - 1)] + [q]
    tot = conv(0)
    for a, v in f.items():
        m = ONE
        for xi, ai in zip(xp, a):
            m = m * xi**ai
        tot = tot + v * m
    return e_of(N, [0] * N) * tot


def reduced(N, V, alpha, beta, Zadj, m=1):
    from rk_c_fast_composite import mono
    c1, c0, Kb, Kadj = decomposition(N, alpha, beta)
    Zraw = L.composite_raw(N, alpha, beta, V)
    Zb = Zraw - c1 * Zadj - c0
    e = (Zb / (Kb * mono(N, bvec(N, alpha, beta), m))).to_sympy()
    return sp.cancel(e.subs({rr: 1 / rr, ss: 1 / ss}, simultaneous=True)), c0.is_zero()


def run(Nmin, Nmax, m=1):
    from rcsknots.refined_kauffman.daha_fast import gl_adj_fast
    for N in range(Nmin, Nmax + 1):
        paths = {ab: os.path.join(OUT, f"Rcomp22_N{N}_m{m}_{tag(*ab)}.pkl") for ab in PAIRS}
        if all(os.path.exists(p) for p in paths.values()):
            continue
        t0 = time.time()
        vpath = os.path.join(OUT, f"V22_N{N}_m{m}.pkl")
        if os.path.exists(vpath):
            V = {k: from_sympy(sp.sympify(v, locals={"r": rr, "s": ss})) for k, v in pickle.load(open(vpath, "rb")).items()}
        else:
            V = L.basis_values(N, P2, P2, m=m, final="functional")
            pickle.dump({k: str(v.to_sympy()) for k, v in V.items()}, open(vpath, "wb"))
        t1 = time.time()
        Zadj = gl_adj_fast(N, m)
        for ab, p in paths.items():
            if os.path.exists(p):
                continue
            Rb, c0zero = reduced(N, V, *ab, Zadj, m)
            pickle.dump({"R": str(Rb), "N": N, "m": m, "alpha": ab[0], "beta": ab[1]}, open(p, "wb"))
        print(f"N={N}：基底の値 {t1-t0:.1f}s、分解と規格化 {time.time()-t1:.1f}s、最大 RSS {rss():.2f} GB", flush=True)
    print("done", flush=True)


def window(m=1):
    import rk_c_interp as I
    for ab in PAIRS:
        avail = [N for N in range(2, 40) if os.path.exists(os.path.join(OUT, f"Rcomp22_N{N}_m{m}_{tag(*ab)}.pkl"))]
        vals = {N: sp.sympify(pickle.load(open(os.path.join(OUT, f"Rcomp22_N{N}_m{m}_{tag(*ab)}.pkl"), "rb"))["R"], locals={"r": rr, "s": ss})
                for N in avail}
        rng = {N: I.s_range(vals[N]) for N in avail}
        slopes = [(rng[avail[i + 1]][0] - rng[avail[i]][0], rng[avail[i + 1]][1] - rng[avail[i]][1]) for i in range(len(avail) - 1)]
        print(f"({tag(*ab)})：N = {avail}、s の次数の範囲 {[rng[N] for N in avail]}、隣の N との差 {slopes}", flush=True)
        jmin, jmax = slopes[-1]
        if not (len(slopes) >= 2 and slopes[-1] == slopes[-2]):
            print("   傾きがまだ一定でない", flush=True)
            continue
        # x = s^N の冪の窓。偶数冪だけなら U = x² = t^N の冪で補間する（項の数が半分になる）
        for unit in ((2, 1) if jmin % 2 == 0 and jmax % 2 == 0 else (1,)):
            lo, hi = jmin // unit, jmax // unit
            n = hi - lo + 1
            if len(avail) < n + 1:
                print(f"   unit={unit}：窓 {lo}〜{hi}（{n} 項）の補間にはノードが足りない（{len(avail)} 個）", flush=True)
                continue
            nodes = avail[-n - 2:-2] if len(avail) >= n + 3 else avail[-n:]      # 大きい方の 2 つを補間の外側の照合に残す
            c = I.interpolate_generic({N: vals[N] for N in nodes}, unit, lo, hi)
            ok = {N: (I.predict_generic(c, unit, N) - from_sympy(vals[N])).is_zero() for N in avail if N not in nodes}
            var = "U = t^N" if unit == 2 else "x = s^N"
            print(f"   {var} の冪 {lo}〜{hi}（{n} 項）、ノード {nodes}、照合 {ok}", flush=True)
            if all(ok.values()):
                pickle.dump({"unit": unit, "c": {j: str(v.to_sympy()) for j, v in c.items()}},
                            open(os.path.join(OUT, f"Rcomp22_A_m{m}_{tag(*ab)}.pkl"), "wb"))
                break


if __name__ == "__main__":
    if sys.argv[1] == "check":
        check()
    elif sys.argv[1] == "run":
        run(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]) if len(sys.argv) > 4 else 1)
    else:
        window(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
