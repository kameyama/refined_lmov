r"""3 箱のクロスキャップ用：GL_N の composite (α;1)（|α| = 2）の T(2m+1, 2) の DAHA-Jones を高速版で計算する（progress.md 続き 33）。
  生の値 Zraw = {γ̂(P_α(X) p_1(X^{−1}))}：X^{−1} 側（β = −ε_1、符号 −1）の対称化で g = γ̂'(p_1(Y))·1 を作り、X 側（β = ε_1）の連鎖の出発点にする
  （src/rcsknots/daha_fast.gl_adj_fast と同じ考え方を P_α(X) に広げたもの）。総次数 c = (|α| − 1)(mk + 1) で ρ = (N−1, …, 0) の評価に直す。
  既約への分解 P_α(X) p_1(X^{−1}) = P_b(X) + d_α P_(1)(X)（b = α − ε_N）：Z_b = Zraw − d_α Z_(1)。
  一様な規格化 W_b = 鏡像(Z_b / M_b)、M_b = q^{−n(b,b)} t^{−2n(b,ρ) + (m+1)|b|(N−1)}（n = 2m + 1）。
  d_α と K_b = ev(P_b) は小さな N では多項式から直接、一般の N では主特殊化の閉じた形から求め、両者を照合する。"""
import sys, os, pickle, time
sys.path.insert(0, "src")
import sympy as sp
from rcsknots.daha_fast import Fast, _with_rho2, solve_symmetric, h_basis, poly_mul, evaluate
from rcsknots.daha_dn import GLN, gl_P_poly
from rcsknots.flint_rf import R, S, ONE, conv

OUT = "data/refined_kauffman/daha"
rr, ss = sp.symbols("r s", positive=True)


def comp_raw(N, alpha, m=1, k=2):
    D = _with_rho2(GLN(N))
    Fm = Fast(D, m, k, -1)
    g = Fm.get_sym((1,))                                   # γ̂'(p_1(Y))·1
    Fp = Fast(D, m, k, 1)
    Fp.sym = {(): g}
    x = solve_symmetric(D, gl_P_poly(D, alpha), h_basis(Fp, [sum(alpha)]))
    v = sum((c * Fp.value(nu) for nu, c in x.items() if not c.is_zero()), conv(0))
    c = (sum(alpha) - 1) * (m * k + 1)
    return v * S**(-c * (N - 1)), D


def ev_rho(D, f):
    """ρ = (N−1, …, 0) での評価（daha_dn.DN.evaluate と同じ規約）。"""
    return D.evaluate(f)


def det_inv(D, f):
    return {tuple(x - 1 for x in a): v for a, v in f.items()}


def decomposition_small(N, alpha):
    """小さな N で P_α(X) p_1(X^{−1}) = P_b + d P_(1) を多項式から解き、(d, K_b, K_1) を返す（K は ρ = (N−1,…,0) の評価）。"""
    D = GLN(N)
    lhs = poly_mul(gl_P_poly(D, alpha), {tuple(-1 if j == i else 0 for j in range(N)): ONE for i in range(N)})
    bp = tuple(a + 1 for a in alpha) + (1,) * (N - len(alpha) - 1)       # b + (1^N)、b = α − ε_N
    Pb = det_inv(D, gl_P_poly(D, bp))
    P1 = gl_P_poly(D, (1,))
    x = solve_symmetric(D, lhs, {"b": Pb, "1": P1})
    assert x["b"] == conv(1)
    return x["1"], D.evaluate(Pb), D.evaluate(P1)


def mono(N, b, m):
    n = 2 * m + 1
    rho = [N - 1 - i for i in range(N)]
    bb = sum(x * x for x in b); br = sum(x * y for x, y in zip(b, rho)); nb = sum(b)
    return R**(-2 * n * bb) * S**(2 * (-2 * n * br + (m + 1) * nb * (N - 1)))


def comp_reduced(N, alpha, m=1, check=False):
    """一様な規格化の簡約した値 R_b(N) = 鏡像(Z_b / (K_b M_b))（sympy）。d_α, K_b は主特殊化の閉じた形から（ev は乗法的）：
    ev(P_λ) = P_λ(1, t^{−1}, …; q, t) = PS(λ; a = t^{−N}, q^{−1}, t^{−1})、K_b = t^{N(N−1)/2} ev(P_{b+1})、d_α = (ev(P_α) ev(p_1(X^{−1})) − K_b)/ev(P_1)。"""
    from rcsknots.macdonald import principal_specialization
    from rcsknots.flint_rf import from_sympy
    from rcsknots.daha_fast import gl_jones_fast
    q_, t_ = rr**2, ss**2
    ps = lambda lam: from_sympy(sp.cancel(principal_specialization(lam, t_**(-N), 1 / q_, 1 / t_)))
    bp = tuple(a + 1 for a in alpha) + (1,) * (N - len(alpha) - 1)
    Kb = from_sympy(t_**(N * (N - 1) // 2)) * ps(bp)
    p1inv = sum((S**(2 * j) for j in range(N)), conv(0))
    d = (ps(alpha) * p1inv - Kb) / ps((1,))
    if check:
        d0, Kb0, _ = decomposition_small(N, alpha)
        assert (d - d0).is_zero() and (Kb - Kb0).is_zero(), "主特殊化の閉じた形が多項式の分解と合わない"
    Zraw, _ = comp_raw(N, alpha, m)
    Zb = Zraw - d * gl_jones_fast(N, (1,), m)
    b = list(alpha) + [0] * (N - len(alpha) - 1) + [-1]
    e = (Zb / (Kb * mono(N, b, m))).to_sympy()
    return sp.cancel(e.subs({rr: 1 / rr, ss: 1 / ss}, simultaneous=True))


def run(m, Nmin, Nmax, alphas=((2,), (1, 1))):
    for N in range(Nmin, Nmax + 1):
        for alpha in alphas:
            if N < len(alpha) + 1:
                continue
            path = os.path.join(OUT, f"Rcomp_N{N}_m{m}_{''.join(map(str, alpha))}.pkl")
            if os.path.exists(path):
                continue
            t0 = time.time()
            R_ = comp_reduced(N, alpha, m, check=N <= 5)
            pickle.dump({"R": str(R_), "N": N, "m": m, "alpha": alpha}, open(path, "wb"))
            print(f"N={N} α={alpha}: {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        run(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
        sys.exit(0)
    # 検証：U(3) の (2;1) ↔ [3,1]、(11;1) ↔ [2,2]（KN の data、一様な規格化）
    sys.path.insert(0, "scripts")
    import rk_k5_trefoil_2box as T
    from rcsknots.flint_rf import from_sympy
    from rcsknots.daha_fast import gl_jones_fast
    mir = lambda e: e.subs({rr: 1 / rr, ss: 1 / ss}, simultaneous=True)
    for N in (3, 5):
        for alpha, shc in (((2,), {3: "31", 5: "3111"}), ((1, 1), {3: "22", 5: "2211"})):
            t0 = time.time()
            Zraw, D = comp_raw(N, alpha)
            d, Kb, K1 = decomposition_small(N, alpha)
            Z1 = gl_jones_fast(N, (1,), 1)
            Zb = Zraw - d * Z1
            b = list(alpha) + [0] * (N - len(alpha) - 1) + [-1]
            sh = shc[N]
            W = (ss / rr)**(sum(int(c) for c in sh)) * T.Rdata(sh, N)            # 一様な規格化（δ = 1、移した分割の箱の数）
            rat = sp.factor(sp.cancel((Zb / Kb / mono(N, b, 1)).to_sympy() / mir(W)))
            print(f"U({N}) の ({''.join(map(str, alpha))};1) ↔ [{sh}]：Z_b/(K_b·M_b) / 鏡像(W) = {rat}（{time.time()-t0:.1f}s）", flush=True)
