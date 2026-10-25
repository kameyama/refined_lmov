r"""rk_d_hyper_fit.fit_hyper の flint 版（有理関数体 rcsknots.refined_kauffman.flint_rf3 での Gauss 消去。sympy 版は遅すぎる）。
  JD(N)·M(N) = ℋ(a = −t^{N−1})、M(N) = q^α t^{β(N−1)} u^γ。ℋ = Σ_{j≤r} h_j a^j を r + 1 個の N で解き（h_j は q,t,u の Laurent 多項式。a⁰ の係数は 1 + qt、Cherednik (super-c-2-3)）、残りの N で照合する。"""
import sys, pickle, itertools
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rcsknots.refined_kauffman.flint_rf3 import RF, R, S, U, conv, from_sympy

OUT = "data/refined_kauffman/daha"
r_, s_, u_ = sp.symbols("r s u", positive=True)
q, t, u, a = sp.symbols("q t u a")


def jd_rf(N, m):
    d = pickle.load(open(f"{OUT}/cn_omega1_N{N}_m{m}.pkl", "rb"))
    loc = {"r": r_, "s": s_, "u": u_}
    return from_sympy(sp.sympify(d["Z"], locals=loc)) / from_sympy(sp.sympify(d["X1"], locals=loc))


def solve_rf(A, b):
    """Gauss 消去（A は RF の正方行列、b はベクトル）。"""
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = next(i for i in range(c, n) if not M[i][c].is_zero())
        M[c], M[p] = M[p], M[c]
        inv = conv(1) / M[c][c]
        M[c] = [x * inv for x in M[c]]
        for i in range(n):
            if i != c and not M[i][c].is_zero():
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    return [M[i][n] for i in range(n)]


def fit(m, al, be, ga, Ns=range(2, 8), check=(8,)):
    r = 2 * m + 1
    vals = {N: jd_rf(N, m) for N in list(Ns) + list(check)}
    M = lambda N: R**(2 * al) * S**(2 * be * (N - 1)) * U**(2 * ga)
    nodes = list(Ns)[: r + 1]
    aN = [-(S**(2 * (N - 1))) for N in nodes]
    A = [[aN[i]**j for j in range(r + 1)] for i in range(r + 1)]
    h = solve_rf(A, [vals[N] * M(N) for N in nodes])
    if not all(x.is_laurent() for x in h):
        return None
    for N in check:
        val = sum((h[j] * (-(S**(2 * (N - 1))))**j for j in range(r + 1)), conv(0))
        if not (val - vals[N] * M(N)).is_zero():
            return None
    return h


if __name__ == "__main__":
    for m in (1, 2):
        found = None
        for al, be, ga in [(2 * m + 1, 4 * m, m + 1)] + [(al, be, ga) for al in range(0, 13) for be in range(0, 13) for ga in range(0, 7)]:
            h = fit(m, al, be, ga)
            if h is not None:
                H = sum(sp.sympify(str(x.to_sympy()), locals={"r": r_, "s": s_, "u": u_}).subs({r_: sp.sqrt(q), s_: sp.sqrt(t), u_: sp.sqrt(u)}) * a**j for j, x in enumerate(h))
                H = sp.expand(sp.cancel(H))
                print(f"T({2*m+1},2)：M(N) = q^{al} t^{{{be}(N−1)}} u^{ga}、ℋ^C = {H}", flush=True)
                pickle.dump({"H": str(H), "M": (al, be, ga)}, open(f"{OUT}/cn_hyper_omega1_m{m}.pkl", "wb"))
                found = (al, be, ga)
                break
        if found is None:
            print(f"T({2*m+1},2)：当てはまらない", flush=True)
