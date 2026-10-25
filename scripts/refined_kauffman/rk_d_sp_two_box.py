r"""案 D の 2 箱：refined Sp(2N)（C_N、u = t²）の trefoil の 2 箱のクロスキャップ部分（progress.md 続き 52）。
  入力：data/refined_kauffman/daha/cn_twobox_N{N}_m1.pkl（N = 2..9、u は記号。rk_d_cn_two_box_data.py）。
  手順：
  (1) R_λ(N) = T_λ^{−2n}·鏡像(Z_λ)/P_λ(q^{−ρ_k})（T_λ = q^{(λ,λ)/2} t^{(λ,ρ_sht)} u^{(λ,ρ_lng)/2}、u = t²）を x = t^N の Laurent 多項式として補間
      （窓は s の次数の傾き、ノードは大きな N、残りで照合）。u = 1 なら D_N の R_U_m1 データに一致するはず（検算）。
  (2) 核の展開係数：P_(2) = (a1/2)S_1² + (1 − a1/2)S_2 + a0′、P_(1,1) = ½S_1² − ½S_2 + b0（a0′ = a0 − a1 N、b0 = c0 − N）と P_λ(q^{−ρ_k}) を
      x の有理関数として当てはめる（分子・分母の次数の窓を走査、8 点のうち 1 点で照合）。
  (3) f^{Sp} = so_f(W, co)（SO と同じ形、W_λ = R_λ P_λ(q^{−ρ_k}) mono_λ）、g = f − ½f^{cov}（被覆は U(2N+1)、A = ν²、ν = t^{N+1/2} = x s）。
      mono_λ は [1] の較正（rk_d_sp_one_box）と g の Laurent 性で決める。
  (4) g^{Sp}_R（R = (2), (1,1)）を g^{SO}_R（largeN_2box_m1.pkl）と比べる：g^{Sp}_R(𝐚) = −g^{SO}_R(−𝐚) か、E^{Sp} と E^{SO} の関係。"""
import sys, os, pickle, itertools
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from fractions import Fraction
from rcsknots.refined_kauffman.flint_rf3 import RF, R, S, U, conv, from_sympy

OUT = "data/refined_kauffman/daha"
r_, s_, u_ = sp.symbols("r s u", positive=True)
nu = sp.Symbol("nu", positive=True)
LAMS = ((1,), (2,), (1, 1))
loc = {"r": r_, "s": s_, "u": u_}


def load(N, m=1):
    d = pickle.load(open(os.path.join(OUT, f"cn_twobox_N{N}_m{m}.pkl"), "rb"))
    Z = {lam: sp.sympify(d["Z"][str(lam)], locals=loc) for lam in LAMS}
    P = {lam: sp.sympify(d["P"][str(lam)], locals=loc) for lam in LAMS}
    coef = {k: sp.sympify(v, locals=loc) for k, v in d["coef"].items()}
    return Z, P, coef


def T_lam(lam, N, usub):
    """framing 固有値 T_λ = q^{(λ,λ)/2} t^{(λ,ρ_sht)} u^{(λ,ρ_lng)/2}（sympy、u = usub）。"""
    lamp = list(lam) + [0] * (N - len(lam))
    ll = sum(x * x for x in lamp)
    lr = sum(x * (N - 1 - i) for i, x in enumerate(lamp))
    lg = sum(lamp)
    return r_**ll * s_**(2 * lr) * sp.sqrt(usub)**lg


def R_of_N(N, lam, usub, m=1):
    """R_λ(N) = T_λ^{−2n}·鏡像(Z_λ)/P_λ(q^{−ρ_k})、鏡像：r, s → 1/r, 1/s（u は usub を代入してから）。sympy（r, s）。"""
    Z, P, _ = load(N, m)
    n = 2 * m + 1
    uh = sp.sqrt(usub)                                            # pkl の u は u^{1/2}
    Zs = sp.cancel(Z[lam].subs(u_, uh)); Ps = sp.cancel(P[lam].subs(u_, uh))
    mir = Zs.subs({r_: 1 / r_, s_: 1 / s_}, simultaneous=True)
    return sp.cancel(T_lam(lam, N, usub)**(-2 * n) * mir / Ps)


def s_range(e):
    n_, d_ = sp.fraction(sp.cancel(e))
    (dm, _), = sp.Poly(d_, r_, s_).terms()
    degs = [mon[1] - dm[1] for mon in sp.Poly(sp.expand(n_), r_, s_).monoms()]
    return min(degs), max(degs)


def newton_interp(xs, ys):
    """x の多項式 Σ_d poly[d] x^d を Newton の差分商で（RF の体の上で）。"""
    n = len(xs)
    coef = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) / (xs[i] - xs[i - j])
    poly = [conv(0)] * n
    basis = [conv(1)]
    for j in range(n):
        for d, b_ in enumerate(basis):
            poly[d] = poly[d] + coef[j] * b_
        nb = [conv(0)] * (len(basis) + 1)
        for d, b_ in enumerate(basis):
            nb[d + 1] = nb[d + 1] + b_
            nb[d] = nb[d] - xs[j] * b_
        basis = nb
    return poly


def interp_laurent(vals, avail):
    """vals[N]（sympy、r, s）を x = s^{2N} の Laurent 多項式として補間。窓は最後の 3 点の s の次数の傾きから。
    返り値：{j: RF（r, s）}（x^j の係数）、ノード、照合結果。"""
    rng = {N: s_range(vals[N]) for N in avail[-3:]}
    a, b, c = avail[-3:]
    jmin = (rng[c][0] - rng[b][0]) // 2; jmax = (rng[c][1] - rng[b][1]) // 2
    assert (rng[b][0] - rng[a][0]) // 2 == jmin and (rng[b][1] - rng[a][1]) // 2 == jmax, ("傾きが一定でない", rng)
    n = jmax - jmin + 1
    assert len(avail) >= n + 1, ("ノードが足りない", n, avail)
    nodes = avail[-n:]
    xs = [S**(2 * N) for N in nodes]
    ys = [from_sympy(vals[N]) / xs[k]**jmin for k, N in enumerate(nodes)]
    poly = newton_interp(xs, ys)
    c_ = {jmin + d: poly[d] for d in range(n)}
    pred = lambda N: sum((v * S**(2 * N * j) for j, v in c_.items()), conv(0))
    ok = {N: (pred(N) - from_sympy(vals[N])).is_zero() for N in avail if N not in nodes}
    return c_, nodes, ok


def rational_fit(vals, avail, max_num=4, max_den=3):
    """vals[N]（sympy、r, s）を x = s^{2N} の有理関数 N(x)/D(x)（D は monic、定数項 1 に正規化）として当てはめる。
    分子の窓 [n0, n1]、分母の次数 d1 を走査し、最後の 1 点で照合。返り値：(分子 {j: RF}, 分母 {j: RF})。"""
    pts = avail[:-1]; check = avail[-1]
    for d1 in range(0, max_den + 1):
        for n0 in range(-max_num, max_num + 1):
            for n1 in range(n0, max_num + 1):
                nn = n1 - n0 + 1
                if nn + d1 > len(pts):
                    continue
                # 未知数：num_j (j = n0..n1)、den_j (j = 1..d1)（den_0 = 1）。方程式：Σ num_j x^j − f(x) Σ_{j≥1} den_j x^j = f(x)
                rows, rhs = [], []
                for N in pts[: nn + d1]:
                    x = S**(2 * N); f = from_sympy(vals[N])
                    rows.append([x**j for j in range(n0, n1 + 1)] + [-f * x**j for j in range(1, d1 + 1)])
                    rhs.append(f)
                try:
                    sol = _solve(rows, rhs)
                except ZeroDivisionError:
                    continue
                if sol is None:
                    continue
                num = {n0 + k: sol[k] for k in range(nn)}
                den = {0: conv(1)}; den.update({k + 1: sol[nn + k] for k in range(d1)})
                good = True
                for N in pts[nn + d1:] + [check]:
                    x = S**(2 * N); f = from_sympy(vals[N])
                    val = sum((v * x**j for j, v in num.items()), conv(0)) / sum((v * x**j for j, v in den.items()), conv(0))
                    if not (val - f).is_zero():
                        good = False; break
                if good:
                    return num, den
    return None


def _solve(A, b):
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    m_ = len(M[0]) - 1
    for c in range(m_):
        p = next((i for i in range(c, n) if not M[i][c].is_zero()), None)
        if p is None:
            return None
        M[c], M[p] = M[p], M[c]
        inv = conv(1) / M[c][c]
        M[c] = [x * inv for x in M[c]]
        for i in range(n):
            if i != c and not M[i][c].is_zero():
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    return [M[i][m_] for i in range(m_)]


def rf_eval(num, den, x):
    return sum((v * x**j for j, v in num.items()), conv(0)) / sum((v * x**j for j, v in den.items()), conv(0))


def to_nu(e):
    """flint_rf3（R, S, U）の式で U を x = t^N とみなし、ν = x s に置き換えて sympy に（r, s, ν）。"""
    return sp.cancel(e.to_sympy().subs(u_, nu / s_))


if __name__ == "__main__":
    avail = [N for N in range(2, 20) if os.path.exists(os.path.join(OUT, f"cn_twobox_N{N}_m1.pkl"))]
    print("利用できる N:", avail, flush=True)
    for usub, lab in ((1, "u = 1（SO）"), (s_**4, "u = t²（Sp）")):
        for lam in LAMS:
            vals = {N: R_of_N(N, lam, usub) for N in avail}
            c_, nodes, ok = interp_laurent(vals, avail)
            print(f"{lab} λ={lam}：x の冪 {min(c_)}〜{max(c_)}、ノード {nodes}、照合 {ok}", flush=True)
            pickle.dump({j: str(v.to_sympy()) for j, v in c_.items()}, open(os.path.join(OUT, f"R_x_{'so' if usub == 1 else 'sp'}_m1_{''.join(map(str, lam))}.pkl"), "wb"))


# ---- (2)〜(4) ----
def coef_functions(avail, usub, m=1):
    """a1(x)、a0′(x) = a0 − a1 N、b0(x) = c0 − N、P_λ(q^{−ρ_k})(x) を x の有理関数として当てはめる。"""
    vals = {"a1": {}, "a0p": {}, "b0": {}}
    Pv = {lam: {} for lam in LAMS}
    for N in avail:
        Z, P, coef = load(N, m)
        uh = sp.sqrt(usub)
        a1 = sp.cancel(coef["a1"].subs(u_, uh)); a0 = sp.cancel(coef["a0"].subs(u_, uh)); c0 = sp.cancel(coef["c0"].subs(u_, uh))
        vals["a1"][N] = a1; vals["a0p"][N] = sp.cancel(a0 - a1 * N); vals["b0"][N] = sp.cancel(c0 - N)
        for lam in LAMS:
            Pv[lam][N] = sp.cancel(P[lam].subs(u_, uh))
    out = {}
    for k, v in vals.items():
        fit = rational_fit(v, avail)
        assert fit is not None, ("当てはまらない", k)
        out[k] = fit
    for lam in LAMS:
        fit = rational_fit(Pv[lam], avail)
        assert fit is not None, ("当てはまらない", lam)
        out[("P", lam)] = fit
    return out


def assemble(m, Rx, cf, mono):
    """f^{Sp}_R（flint_rf3、U = x = t^N）。Rx：{λ: {j: RF}}、cf：coef_functions の出力、mono：{λ: RF}。"""
    import rk_c_2box_largeN as L
    x = U
    ev = lambda fit: rf_eval(fit[0], fit[1], x)
    a1, a0p, b0 = ev(cf["a1"]), ev(cf["a0p"]), ev(cf["b0"])
    W = {}
    for lam in LAMS:
        Rl = sum((v * x**j for j, v in Rx[lam].items()), conv(0))
        W[lam] = Rl * ev(cf[("P", lam)]) * mono[lam]
    co = {"S1sq": {(2,): conv(1), (1, 1): conv(2) - a1, (1,): conv(0), (): -a0p + (a1 - conv(2)) * b0},
          "S2": {(2,): conv(1), (1, 1): -a1, (1,): conv(0), (): -a0p + a1 * b0}}
    return L.so_f(W, co), W, co


def cover_in_x(m):
    """被覆（rk_c_2box_largeN_m.cover_f、変数 U_D = ν s）を U = x = t^N（ν = x s、U_D = x s²）に置き換える。"""
    import rk_c_2box_largeN_m as Lm
    fc = Lm.cover_f(m, f"T2_{2*m+1}")
    out = {}
    for k, v in fc.items():
        e = v.to_sympy().subs(u_, u_ * s_**2)
        out[k] = from_sympy(sp.cancel(e))
    return out


def bold_nu(e):
    """flint_rf3（U = x、ν = x s）→ 太字 (A, Q, T)：r = −TQ、s = Q、ν = A。"""
    A_, Q_, T_ = sp.symbols("A Q T")
    ex = e.to_sympy().subs(u_, nu / s_)
    return sp.expand(sp.cancel(ex.subs({r_: -T_ * Q_, nu: A_}).subs(s_, Q_))), (A_, Q_, T_)


def bold_so(e_str):
    """largeN_2box_m1.pkl の g（r, s, u = U_D = ν s）→ 太字。"""
    A_, Q_, T_ = sp.symbols("A Q T")
    ex = sp.sympify(e_str, locals=loc).subs(u_, nu * s_)
    return sp.expand(sp.cancel(ex.subs({r_: -T_ * Q_, nu: A_}).subs(s_, Q_)))


def main2(m=1):
    avail = [N for N in range(2, 20) if os.path.exists(os.path.join(OUT, f"cn_twobox_N{N}_m{m}.pkl"))]
    usub = s_**4
    Rx = {}
    for lam in LAMS:
        c_, nodes, ok = interp_laurent({N: R_of_N(N, lam, usub, m) for N in avail}, avail)
        assert all(ok.values()), (lam, ok)
        Rx[lam] = c_
    cf = coef_functions2(avail, usub, m)
    print("係数の当てはめ：", {k: (sorted(v[0]), sorted(v[1])) for k, v in cf.items()}, flush=True)
    fc = cover_in_x(m)
    gso = pickle.load(open(os.path.join(OUT, f"largeN_2box_m{m}.pkl"), "rb"))["g"]
    gso_b = {k: bold_so(v) for k, v in gso.items()}
    A_, Q_, T_ = sp.symbols("A Q T")
    half = conv(Fraction(1, 2))
    # [1] の単項式（rk_d_sp_one_box の較正）：W_Sp,(1) = R_(1) P_(1) mono_(1)。mono_(1) は main で求めた値を使う（引数で上書き可）
    results = {}
    for e1 in (0, 1, 2, -1, -2):
        for e2 in (0, 1, 2, -1, -2):
            mono = {(1,): S**(2 * e1) * U**0, (2,): S**(2 * e2), (1, 1): S**(2 * e2)}
            f, W, co = assemble(m, Rx, cf, mono)
            g = {k: f[k] - half * fc[k] for k in f}
            lau = {k: g[k].is_laurent() for k in g}
            if all(lau.values()):
                results[(e1, e2)] = g
                gb = {k: bold_nu(g[k])[0] for k in g}
                rel = {k: sp.expand(gb[k] + gso_b[str(k)].subs(A_, -A_)) == 0 for k in g}
                print(f"mono: s^{2*e1}（[1]）, s^{2*e2}（2 箱）：g は Laurent。g^Sp_R(𝐚) = −g^SO_R(−𝐚)：{rel}", flush=True)
    return results


def P_eval_formula(N, lam, usub):
    """Cherednik (pebeb)：P_b(q^{−ρ_k}) = q^{(ρ_k,b)} Π_{α>0} Π_{j=0}^{(α∨,b)−1} (1 − q_α^j t_α X_α(q^{ρ_k}))/(1 − q_α^j X_α(q^{ρ_k}))。
    C_N：短い根 ε_i ± ε_j（q_α = q、t_α = t、α∨ = α）、長い根 2ε_i（q_α = q²、t_α = u、α∨ = ε_i）。X_b(q^{ρ_k}) = t^{(b,ρ_sht)} u^{(b,ρ_lng)/2}。sympy（r, s）。"""
    q2, t2 = r_**2, s_**2
    b = list(lam) + [0] * (N - len(lam))
    rho_s = [N - 1 - i for i in range(N)]
    X = lambda v: t2**sum(x * y for x, y in zip(v, rho_s)) * sp.sqrt(usub)**sum(v)       # X_v(q^{ρ_k})
    val = 1 / X(b)                                   # q^{−(ρ_k,b)}（Cherednik (pebeb) の印刷は q^{(ρ_k,b)} だが、データ（N = 2..9）に合うのはこちら）
    roots = []
    for i in range(N):
        roots.append(([2 if k == i else 0 for k in range(N)], q2**2, usub, [1 if k == i else 0 for k in range(N)]))   # 長い根、α∨ = ε_i
        for j in range(i + 1, N):
            for sg in (1, -1):
                al = [1 if k == i else (sg if k == j else 0) for k in range(N)]
                roots.append((al, q2, t2, al))
    for al, qa, ta, alv in roots:
        n_ = sum(x * y for x, y in zip(alv, b))
        Xa = X(al)
        for j in range(n_):
            val = val * (1 - qa**j * ta * Xa) / (1 - qa**j * Xa)
    return sp.cancel(val)


def coef_functions2(avail, usub, m=1, Nmax=24):
    """coef_functions の P_λ(q^{−ρ_k}) を評価公式（N ≤ Nmax）で当てはめる版。"""
    out = coef_functions_abc(avail, usub, m)
    for lam in LAMS:
        # u = 1 では長い根 2ε_N の因子が 0/0 になるので、u を記号のまま計算してから代入する
        vals = {N: sp.cancel(P_eval_formula(N, lam, u_**2).subs(u_, sp.sqrt(usub))) for N in range(2, Nmax + 1)}
        # データとの照合
        for N in avail:
            Z, P, coef = load(N, m)
            assert sp.cancel(vals[N] - P[lam].subs(u_, sp.sqrt(usub))) == 0, ("評価公式がデータと合わない", lam, N)
        fit = rational_fit(vals, list(range(2, Nmax + 1)), max_num=8, max_den=6)
        assert fit is not None, ("当てはまらない", lam)
        out[("P", lam)] = fit
    return out


def coef_functions_abc(avail, usub, m=1):
    vals = {"a1": {}, "a0p": {}, "b0": {}}
    for N in avail:
        Z, P, coef = load(N, m)
        uh = sp.sqrt(usub)
        a1 = sp.cancel(coef["a1"].subs(u_, uh)); a0 = sp.cancel(coef["a0"].subs(u_, uh)); c0 = sp.cancel(coef["c0"].subs(u_, uh))
        vals["a1"][N] = a1; vals["a0p"][N] = sp.cancel(a0 - a1 * N); vals["b0"][N] = sp.cancel(c0 - N)
    out = {}
    for k, v in vals.items():
        fit = rational_fit(v, avail)
        assert fit is not None, ("当てはまらない", k)
        out[k] = fit
    return out
