r"""(b) framed unknot（framing x）の大きな N の 2 箱のクロスキャップ部分と、固定セクターでの E の決定（progress.md 続き 45）。
  framing 0 の unknot では f_R = f^{cov}_R = 0（|R| ≥ 2）で E = C = 0 なので、「loop quiver では E が決まる」という予言は自明に成り立つ。
  自明でない版として framing x を入れる：SO 側 W_λ = T_λ^x K_λ(t^ρ)、被覆 W_μ = T_μ^x P_μ(t^ρ)、T は DAHA の framing 固有値
  （D_N：T_λ = q^{(λ,λ)/2} t^{(λ,ρ)}、GL_{N_A}：同じ式。両者は U = t^N、A = t^{2N−1} で一致：T_(1) = q^{1/2}U/t、T_(2) = q²U²/t²、T_(1,1) = qU²/t³、T_adj = qU²/t²）。
  【仮説】framed unknot の t = q の条件 (ii) は結び目と同じ形 E|_{t=q} = −(1+q)(g_(1)(q², ν²) + 1) とする。"""
import sys, io, contextlib, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp
from fractions import Fraction
from rcsknots.flint_rf3 import RF, R, S, U, from_sympy, conv
import rk_c_2box_largeN as L
import rk_c_2box_largeN_m as Lm

q, t = R * R, S * S
nu = U / S
half = conv(Fraction(1, 2))
r_, s_, u_ = sp.symbols("r s u", positive=True)
nus = sp.Symbol("nu", positive=True)


def T_so(lam):
    return {(1,): R * U / t, (2,): q * q * U * U / (t * t), (1, 1): q * U * U / (t**3)}[lam]


def run(x):
    from rcsknots.kauffman_hopf2 import koornwinder_small, K_rho
    from rcsknots.kauffman_hopf import qs, ts, Us
    sub = {qs: r_**2, ts: s_**2, Us: u_}
    cv = lambda v: from_sympy(sp.sympify(v.as_expr() if hasattr(v, "as_expr") else v).subs(sub))
    ks = koornwinder_small()
    co = {name: {lam: cv(v) for lam, v in ks[name].items()} for name in ("S1sq", "S2")}
    W = {lam: T_so(lam)**x * cv(K_rho(lam)) for lam in ((1,), (2,), (1, 1))}
    f = L.so_f(W, co)
    # 被覆
    Wc = {lam: T_so(lam)**x * L.Psym(lam) for lam in ((1,), (1, 1), (2,))}
    A = nu * nu
    qn = (nu - 1 / nu) / (S - 1 / S)
    kappa = (1 - A) * (1 - q) / ((1 - t) * (1 - q * A / t))
    Tadj = q * U * U / (t * t)
    Wadj = Tadj**x * (qn * qn - kappa)
    c = (t - q) / (1 - q * t)
    A1 = -2 * Wc[(1,)]
    A2 = 2 * Wc[(1, 1)] + Wadj + kappa
    A11 = 2 * (Wc[(2,)] + c * Wc[(1, 1)]) + Wadj + kappa
    Z1, Z11, Z2 = A1, (A2 + A11) * half, (A2 - A11) * half
    f1 = Z1
    H11 = Z11 - Z1 * Z1 * half
    H2 = Z2 - f1.adams(2) * half
    fc = {(1,): f1, (2,): H11 + H2, (1, 1): H11 - H2}
    g = {k: f[k] - half * fc[k] for k in f}
    xsec = Wadj + kappa - Wc[(1,)] * Wc[(1,)]                     # brane–image のセクター
    out = {}
    for k, v in g.items():
        e = sp.cancel(v.to_sympy())
        out[k] = e
        print(f"framing {x}：g_{k} の分母 {sp.factor(sp.fraction(e)[1])}、項数 {len(sp.Add.make_args(sp.fraction(e)[0]))}", flush=True)
    print(f"   g_(2)+g_(1,1) 整数 Laurent {(g[(2,)] + g[(1,1)]).is_integer_laurent()}、(1+t)(g_(2)−g_(1,1)) 整数 Laurent {((1 + t) * (g[(2,)] - g[(1,1)])).is_integer_laurent()}", flush=True)
    return out, sp.cancel(xsec.to_sympy())


def solve_E(x, variant="A", quiet=False):
    """variant：t = q の条件 (ii) の形。A：E = −(1+q)(g_(1)(q²,ν²) + 1)（結び目の形）、B：−(1+q) g_(1)(q²,ν²)、C：−(1+q)(g_(1)(q²,ν²) − 1)。"""
    import rk_k4_fixed_sector as H
    if quiet:
        with contextlib.redirect_stdout(io.StringIO()):
            g, xsec = run(x)
    else:
        g, xsec = run(x)
    r, s, nu_ = H.r, H.s, H.nu
    sub = lambda e: sp.cancel(e.subs({r_: r, s_: s}).subs(u_, nu_ * s))
    g = {k: sub(v) for k, v in g.items()}
    P = sp.cancel((1 - s**2) * sub(xsec))
    R_ = sp.cancel(2 * (1 + s**2) * g[(2,)])
    g1x2 = g[(1,)].subs({r: r**2, s: s**2, nu_: nu_**2}, simultaneous=True)
    shift = {"A": 1, "B": 0, "C": -1}[variant]
    Eq = sp.cancel((-(1 + r**2) * (g1x2 + shift)).subs(s, r))
    print(f"   P の分子の因子 {sp.factor(sp.fraction(P)[0])}", flush=True)
    try:
        res = H.solve_generic(P, R_, Eq)
    except AssertionError as ex:
        print(f"   (ii) の形 {variant}：{ex}（偶奇が合わないので不可）", flush=True)
        return None
    print(f"   (ii) の形 {variant}：固定セクターの単項式 {res['n']} 個、(i)+(ii) の解あり = {res['feasible']}" + (f"、一意でない変数 {len(res['free'])} 個" if res["feasible"] else ""), flush=True)
    if res["feasible"]:
        print(f"   E = {sp.factor(res['E'])}", flush=True)
    return res


if __name__ == "__main__":
    for x in (1, 2, 3, 4):
        for variant in ("A", "B", "C"):
            solve_E(x, variant, quiet=(variant != "A"))
