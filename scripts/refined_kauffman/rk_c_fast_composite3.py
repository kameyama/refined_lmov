r"""4 箱のクロスキャップ用：GL_N の composite (α;1̄)（|α| = 3）の T(2m+1, 2) の DAHA-Jones（progress.md 続き 139）。
  生の値 Zraw = {γ̂(P_α(X) p_1(X^{−1}))} は rk_c_fast_composite.comp_raw（|α| によらない）。
  既約への分解 P_α(X) p_1(X^{−1}) = P_b(X) + Σ_μ c_{α,μ} P_μ(X)（b = α − ε_N、μ = α − ε_i は分割、i < ℓ(α)）：
    p_1(X^{−1}) = e_{N−1}(X)/e_N(X) と縦帯の Pieri P_α e_{N−1} = Σ_λ ψ'_{λ/α} P_λ から、c_{α,μ} = ψ'_{λ/α}、λ = α + (1^N) − ε_i。
    ψ' は kauffman_cover.pieri_vertical（数値で特定済み）。A = t^N の閉じた形は psi_A（列 0 の脚 N − i − 1 だけが N に依る）。
  照合：主特殊化（ev は乗法的）で ev(P_α) ev(p_1(X^{−1})) = K_b + Σ_μ c_{α,μ} ev(P_μ)、小さな N では多項式の分解と比較。
  一様な規格化 W_b = 鏡像(Z_b / M_b)（rk_c_fast_composite.mono）。保存：data/refined_kauffman/daha/Rcomp_N{N}_m{m}_{α}.pkl（同じ形式）。
使い方：python rk_c_fast_composite3.py m Nmin Nmax [α ...]（α は "3", "21", "111"）"""
import sys, os, pickle, time, resource
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_c_fast_composite as C
from rcsknots.refined_kauffman.daha_fast import solve_symmetric, poly_mul, gl_jones_fast
from rcsknots.refined_kauffman.daha_dn import GLN, gl_P_poly
from rcsknots.refined_kauffman.flint_rf import R, S, ONE, conv, from_sympy
from rcsknots.refined_kauffman.kauffman_cover import pieri_vertical
from rcsknots.partitions import conjugate

rr, ss = C.rr, C.ss
OUT = C.OUT


def lower(alpha):
    """[(i, μ)]：μ = α − ε_i が分割になる i < ℓ(α)。"""
    out = []
    for i in range(len(alpha)):
        if i == len(alpha) - 1 or alpha[i] > alpha[i + 1]:
            mu = tuple(x for x in (alpha[:i] + (alpha[i] - 1,) + alpha[i + 1:]) if x > 0)
            out.append((i, mu))
    return out


def strip_lam(alpha, i, N):
    a = list(alpha) + [0] * (N - len(alpha))
    return tuple(a[r] + (0 if r == i else 1) for r in range(N))


def c_pieri(alpha, i, N):
    """c_{α,μ}(N) = ψ'_{λ/α}（sympy、q = r², t = s²）。"""
    return sp.cancel(pieri_vertical(strip_lam(alpha, i, N), alpha, rr**2, ss**2))


def psi_A(alpha, i, q, t, A):
    """c_{α,μ} を A = t^N の閉じた形で（q, t, A は任意の環の元）。
    ψ'_{λ/α} = Π_{j < α_i, j ∈ cols} b_λ(i, j)/b_α(i, j)、cols = {α_r : r ≠ i} ∪ {0}（ℓ(α) < N）。
    b(i, j) = (1 − q^a t^{l+1})/(1 − q^{a+1} t^l)、a = α_i − j − 1。l_λ は j ≥ 1 で N に依らず、j = 0 で t^{l_λ} = A t^{−i−1}。"""
    cols = {alpha[r] for r in range(len(alpha)) if r != i} | {0}
    aT = conjugate(alpha)
    out = 1
    for j in range(alpha[i]):
        if j not in cols:
            continue
        a = alpha[i] - j - 1
        la = aT[j] - i - 1
        if j == 0:
            tl = A * t**(-i - 1)
        else:
            lamT_j = sum(1 for r in range(len(alpha)) if (alpha[r] + (0 if r == i else 1)) > j)
            tl = t**(lamT_j - i - 1)
        out = out * ((1 - q**a * tl * t) / (1 - q**(a + 1) * tl)) / ((1 - q**a * t**(la + 1)) / (1 - q**(a + 1) * t**la))
    return out


def decomposition_small(N, alpha):
    """小さな N で P_α(X) p_1(X^{−1}) = P_b + Σ_μ c_μ P_μ を多項式から解く。{μ: c_μ} と K_b を返す。"""
    D = GLN(N)
    lhs = poly_mul(gl_P_poly(D, alpha), {tuple(-1 if j == i else 0 for j in range(N)): ONE for i in range(N)})
    bp = tuple(a + 1 for a in alpha) + (1,) * (N - len(alpha) - 1)
    Pb = C.det_inv(D, gl_P_poly(D, bp))
    basis = {"b": Pb}
    for i, mu in lower(alpha):
        basis[mu] = gl_P_poly(D, mu)
    x = solve_symmetric(D, lhs, basis)
    assert x["b"] == conv(1)
    return {mu: x[mu] for _, mu in lower(alpha)}, D.evaluate(Pb)


def comp_reduced(N, alpha, m=1, check=False):
    from rcsknots.macdonald import principal_specialization
    q_, t_ = rr**2, ss**2
    ps = lambda lam: from_sympy(sp.cancel(principal_specialization(lam, t_**(-N), 1 / q_, 1 / t_)))
    bp = tuple(a + 1 for a in alpha) + (1,) * (N - len(alpha) - 1)
    Kb = from_sympy(t_**(N * (N - 1) // 2)) * ps(bp)
    p1inv = sum((S**(2 * j) for j in range(N)), conv(0))
    cs = {mu: from_sympy(c_pieri(alpha, i, N)) for i, mu in lower(alpha)}
    # 主特殊化の恒等式
    lhs = ps(alpha) * p1inv
    rhs = Kb + sum((cs[mu] * ps(mu) for mu in cs), conv(0))
    assert (lhs - rhs).is_zero(), ("主特殊化の恒等式が合わない", N, alpha)
    # A の閉じた形
    for i, mu in lower(alpha):
        cA = from_sympy(sp.cancel(sp.sympify(psi_A(alpha, i, q_, t_, sp.Symbol("A"))).subs(sp.Symbol("A"), t_**N)))
        assert (cA - cs[mu]).is_zero(), ("psi_A が合わない", N, alpha, mu)
    if check:
        c0, Kb0 = decomposition_small(N, alpha)
        assert all((c0[mu] - cs[mu]).is_zero() for mu in cs) and (Kb - Kb0).is_zero(), "Pieri の係数が多項式の分解と合わない"
    Zraw, _ = C.comp_raw(N, alpha, m)
    Zb = Zraw
    for mu, c in cs.items():
        Zb = Zb - c * gl_jones_fast(N, mu, m)
    b = list(alpha) + [0] * (N - len(alpha) - 1) + [-1]
    e = (Zb / (Kb * C.mono(N, b, m))).to_sympy()
    return sp.cancel(e.subs({rr: 1 / rr, ss: 1 / ss}, simultaneous=True))


def run(m, Nmin, Nmax, alphas=((3,), (2, 1), (1, 1, 1))):
    for N in range(Nmin, Nmax + 1):
        for alpha in alphas:
            if N < len(alpha) + 1:
                continue
            path = os.path.join(OUT, f"Rcomp_N{N}_m{m}_{''.join(map(str, alpha))}.pkl")
            if os.path.exists(path):
                continue
            t0 = time.time()
            R_ = comp_reduced(N, alpha, m, check=N <= 4)
            pickle.dump({"R": str(R_), "N": N, "m": m, "alpha": alpha}, open(path, "wb"))
            print(f"N={N} α={alpha}: {time.time()-t0:.1f}s、最大 RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.0f}MB", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    m, Nmin, Nmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    al = [tuple(int(c) for c in a) for a in sys.argv[4:]] or [(3,), (2, 1), (1, 1, 1)]
    run(m, Nmin, Nmax, al)
