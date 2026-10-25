r"""タスク C：D_N の DAHA-Jones の Z(N)（scripts/rk_c_daha_dn.py の保存値）から、一様な規格化の reduced 不変量
  R_λ(N) = T_λ(N)^6 × 鏡像(Z(N)) / K_λ(t^ρ)、T_λ = q^{(λ,λ)/2} t^{(λ,ρ)}（鏡像：r → 1/r, s → 1/s）
を求め、U = t^N の Laurent 多項式 R_λ(U) = Σ_j c_j(q, t) U^j として補間する（Vandermonde、ノードは N の値）。
補間に使わない N の値で予言を確かめる。"""
import sys, pickle, os
sys.path.insert(0, "src")
import sympy as sp

rr, ss, U = sp.symbols("r s U", positive=True)
OUT = "data/refined_kauffman/daha"


def load(N, m, lam):
    path = os.path.join(OUT, f"Z_N{N}_m{m}_{''.join(map(str, lam))}.pkl")
    if not os.path.exists(path):
        return None
    d = pickle.load(open(path, "rb"))
    return sp.sympify(d["Z"], locals={"r": rr, "s": ss}), sp.sympify(d["K"], locals={"r": rr, "s": ss})


def R_true(N, m, lam):
    got = load(N, m, lam)
    if got is None:
        return None
    Z, K = got
    lamp = list(lam) + [0] * (N - len(lam))
    ll = sum(x * x for x in lamp)
    lr = sum(x * (N - 1 - i) for i, x in enumerate(lamp))
    Tl = rr**ll * ss**(2 * lr)
    n = 2 * m + 1                       # T(n, 2)：framing の冪は 2n
    mir = Z.subs({rr: 1 / rr, ss: 1 / ss}, simultaneous=True)
    return sp.cancel(Tl**(-2 * n) * mir / K)          # W = T_λ^{−2n} × 鏡像(Z)（N = 2, 3 の data で確認）


def s_range(e):
    n, d = sp.fraction(sp.cancel(e))
    (dm, _), = sp.Poly(d, rr, ss).terms()
    degs = [mon[1] - dm[1] for mon in sp.Poly(sp.expand(n), rr, ss).monoms()]
    return min(degs), max(degs)


def interpolate(m, lam, nodes, jmin, jmax):
    """nodes の N で R(N) = Σ_{j=jmin}^{jmax} c_j U^j（U = s^{2N}）を解く。R/U^{jmin} を x = U の多項式として Newton の差分商で補間
    （有理関数は rcsknots.refined_kauffman.flint_rf.RF）。{j: c_j（sympy）} を返す。"""
    from rcsknots.refined_kauffman.flint_rf import from_sympy, S as Sf, ONE, conv
    js = list(range(jmin, jmax + 1))
    assert len(nodes) == len(js)
    xs = [Sf**(2 * N) for N in nodes]
    ys = [from_sympy(R_true(N, m, lam)) / xs[k]**jmin for k, N in enumerate(nodes)]
    n = len(xs)
    coef = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) / (xs[i] - xs[i - j])
    # Newton 形 → 単項式の係数
    poly = [conv(0)] * n
    basis = [ONE]                                   # Π_{i<j}(x − x_i) の係数
    for j in range(n):
        for d, b_ in enumerate(basis):
            poly[d] = poly[d] + coef[j] * b_
        nb = [conv(0)] * (len(basis) + 1)
        for d, b_ in enumerate(basis):
            nb[d + 1] = nb[d + 1] + b_
            nb[d] = nb[d] - xs[j] * b_
        basis = nb
    return {jmin + d: poly[d] for d in range(n)}


def predict(c, N):
    from rcsknots.refined_kauffman.flint_rf import S as Sf, conv
    return sum((v * Sf**(2 * N * j) for j, v in c.items()), conv(0))


if __name__ == "__main__":
    m = int(sys.argv[1])
    for lam in [tuple(int(x) for x in a) for a in sys.argv[2:]]:
        avail = [N for N in range(2, 20) if load(N, m, lam) is not None]
        print(f"λ={lam}、m={m}：N = {avail}", flush=True)
        for N in avail:
            R = R_true(N, m, lam)
            lau = sp.Poly(sp.fraction(sp.cancel(R))[1], rr, ss).is_monomial
            print(f"   N={N}: R_true は Laurent {lau}、s の次数の範囲 {s_range(R) if lau else '-'}", flush=True)


def interpolate_generic(Rvals, unit, jmin, jmax):
    """Rvals = {N: R(N)（sympy）}。R(N) = Σ_{j=jmin}^{jmax} c_j x^j、x = s^{unit·N} を Newton の差分商で解く（ノードは Rvals の全部）。"""
    from rcsknots.refined_kauffman.flint_rf import from_sympy, S as Sf, ONE, conv
    nodes = sorted(Rvals)
    assert len(nodes) == jmax - jmin + 1
    xs = [Sf**(unit * N) for N in nodes]
    ys = [from_sympy(Rvals[N]) / xs[k]**jmin for k, N in enumerate(nodes)]
    n = len(xs)
    coef = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) / (xs[i] - xs[i - j])
    poly = [conv(0)] * n
    basis = [ONE]
    for j in range(n):
        for d, b_ in enumerate(basis):
            poly[d] = poly[d] + coef[j] * b_
        nb = [conv(0)] * (len(basis) + 1)
        for d, b_ in enumerate(basis):
            nb[d + 1] = nb[d + 1] + b_
            nb[d] = nb[d] - xs[j] * b_
        basis = nb
    return {jmin + d: poly[d] for d in range(n)}


def predict_generic(c, unit, N):
    from rcsknots.refined_kauffman.flint_rf import S as Sf, conv
    return sum((v * Sf**(unit * N * j) for j, v in c.items()), conv(0))


def R_adj(N, m=1):
    path = os.path.join(OUT, f"Wadj_N{N}_m{m}.pkl")
    if not os.path.exists(path):
        return None
    W = sp.sympify(pickle.load(open(path, "rb"))["W"], locals={"r": rr, "s": ss})
    A2 = ss**N                                    # A^{1/2} = t^{N/2}
    qn = (A2 - 1 / A2) / (ss - 1 / ss)
    q_, t_ = rr**2, ss**2
    kappa = (1 - t_**N) * (1 - q_) / ((1 - t_) * (1 - q_ * t_**(N - 1)))
    return sp.cancel(W / (qn**2 - kappa))
