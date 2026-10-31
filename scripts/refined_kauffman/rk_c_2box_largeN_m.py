r"""タスク C ②：大きな N の 2 箱のクロスキャップ部分 g = f − ½f^{cov} を T(2,2m+1) に一般化する（rk_c_2box_largeN.py は m = 1 の版、progress.md 続き 32）。
  SO：W_λ(U) = R_λ(U) K_λ(t^ρ)、R_λ(U) は data/refined_kauffman/daha/R_U_m{m}_{λ}.pkl（scripts/rk_c_read_largeN.py）。
  被覆 U(N_A)：普通の表現は W_μ = (t/q)^{δ|μ|/2} R_μ(KN の data) P_μ(t^ρ)（δ = m）、随伴は R_adj(A^{1/2}) を GL の DAHA（scripts/rk_c_fast_adj.py）から
  A^{1/2} = t^{N_A/2} の Laurent 多項式として補間し（傾きで窓を決め、補間に使わない N_A で照合）、data/refined_kauffman/daha/Radj_A_m{m}.pkl に保存。
  検算：U = t², t³ で固定 N の g（data/refined_kauffman/T2_{2m+1}_fixedN_3box_g.pkl）と全係数比較。
使い方：python scripts/rk_c_2box_largeN_m.py m [--k=K]（K ≠ 2 は T(mK+1, K)、タグ m{m}k{K}、δ は rk_torus_tag.DELTA。続き 102）"""
import sys, pickle, os
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from fractions import Fraction
from rcsknots.refined_kauffman.flint_rf3 import RF, R, S, U, from_sympy, conv
import rk_c_2box_largeN as L
import rk_c_interp as I
import rk_torus_tag as TG

r_, s_, u_ = sp.symbols("r s u", positive=True)
q, t = R * R, S * S
nu = U / S
half = conv(Fraction(1, 2))
OUT = "data/refined_kauffman/daha"


def so_parts(m, k=2):
    from rcsknots.refined_kauffman.kauffman_hopf2 import koornwinder_small, K_rho
    from rcsknots.refined_kauffman.kauffman_hopf import qs, ts, Us
    sub = {qs: r_**2, ts: s_**2, Us: u_}
    cv = lambda v: from_sympy(sp.sympify(v.as_expr() if hasattr(v, "as_expr") else v).subs(sub))
    ks = koornwinder_small()
    co = {name: {lam: cv(v) for lam, v in ks[name].items()} for name in ("S1sq", "S2")}
    W = {}
    for lam in ((1,), (2,), (1, 1)):
        c = pickle.load(open(os.path.join(OUT, f"R_U_{TG.tag(m, k)}_{''.join(map(str, lam))}.pkl"), "rb"))
        Rf = sum((from_sympy(sp.sympify(v, locals={"r": r_, "s": s_})) * U**j for j, v in c.items()), conv(0))
        W[lam] = Rf * cv(K_rho(lam))
    return W, co


def adj_interp(m):
    path = os.path.join(OUT, f"Radj_A_m{m}.pkl")
    avail = [N for N in range(2, 40) if os.path.exists(os.path.join(OUT, f"Wadj_N{N}_m{m}.pkl"))]
    vals = {N: I.R_adj(N, m) for N in avail}
    a, b, c3 = avail[-3:]
    rng = {N: I.s_range(vals[N]) for N in (a, b, c3)}
    jmin, jmax = rng[c3][0] - rng[b][0], rng[c3][1] - rng[b][1]          # x = s^N の冪 j で s の次数は N·j + 定数
    assert rng[b][0] - rng[a][0] == jmin and rng[b][1] - rng[a][1] == jmax, rng
    n = jmax - jmin + 1
    nodes = avail[-n:]
    c = I.interpolate_generic({N: vals[N] for N in nodes}, 1, jmin, jmax)
    from rcsknots.refined_kauffman.flint_rf import from_sympy as fs2
    ok = {N: (I.predict_generic(c, 1, N) - fs2(vals[N])).is_zero() for N in avail if N not in nodes}
    pickle.dump({j: str(v.to_sympy()) for j, v in c.items()}, open(path, "wb"))
    print(f"随伴 m={m}：A^{{1/2}} の冪 {jmin}〜{jmax}（ノード N_A = {nodes}）、照合 {ok}", flush=True)
    return c


def cover_f(m, knot, k=2):
    w = (S / R)**TG.delta(m, k)                          # (t/q)^{δ/2}、δ = m（k = 2）
    Wc = {lam: w**sum(lam) * L.Rdata_nu("".join(map(str, lam)), knot) * L.Psym(lam) for lam in ((1,), (1, 1), (2,))}
    radj = pickle.load(open(os.path.join(OUT, f"Radj_A_{TG.tag(m, k)}.pkl"), "rb"))
    Radj = sum((from_sympy(sp.sympify(v, locals={"r": r_, "s": s_})) * nu**j for j, v in radj.items()), conv(0))
    A = nu * nu
    qn = (nu - 1 / nu) / (S - 1 / S)
    kappa = (1 - A) * (1 - q) / ((1 - t) * (1 - q * A / t))
    Wadj = Radj * (qn * qn - kappa)
    c = (t - q) / (1 - q * t)
    A1 = -2 * Wc[(1,)]
    A2 = 2 * Wc[(1, 1)] + Wadj + kappa
    A11 = 2 * (Wc[(2,)] + c * Wc[(1, 1)]) + Wadj + kappa
    Z1, Z11, Z2 = A1, (A2 + A11) * half, (A2 - A11) * half
    f1 = Z1
    H11 = Z11 - Z1 * Z1 * half
    H2 = Z2 - f1.adams(2) * half
    return {(1,): f1, (2,): H11 + H2, (1, 1): H11 - H2}


def main(m, k=2):
    knot, tg = TG.knot_name(m, k), TG.tag(m, k)
    if os.path.exists(os.path.join(OUT, f"Radj_A_{tg}.pkl")):
        print(f"Radj_A_{tg}.pkl があるので随伴の補間は飛ばす（rk_c_bounded_save.py --adj で作った場合）", flush=True)
    else:
        assert k == 2
        adj_interp(m)
    W, co = so_parts(m, k)
    f = L.so_f(W, co)
    fc = cover_f(m, knot, k)
    g = {k: f[k] - half * fc[k] for k in f}
    for k, v in g.items():
        d = sp.factor(sp.fraction(sp.cancel(v.to_sympy()))[1])
        print(f"g_{k}：分母 {d}", flush=True)
    print("g_(2) + g_(1,1) は整数 Laurent", (g[(2,)] + g[(1, 1)]).is_integer_laurent(), "、(1+t)(g_(2) − g_(1,1)) は整数 Laurent",
          ((1 + t) * (g[(2,)] - g[(1, 1)])).is_integer_laurent(), flush=True)
    fxp = f"data/refined_kauffman/{knot}_fixedN_3box_g.pkl"
    fx = pickle.load(open(fxp, "rb")) if os.path.exists(fxp) else {}
    if not fx:
        print(f"固定 N の g（{fxp}）が無いので U = t², t³ の照合は飛ばす", flush=True)
    for N in (2, 3):
        if N not in fx:
            continue
        ok = all(sp.cancel(g[k].to_sympy().subs(u_, s_**(2 * N)) - sp.sympify(fx[N]["g_schur"][str(k)], locals={"r": r_, "s": s_})) == 0 for k in ((1,), (2,), (1, 1)))
        print(f"U = t^{N} で固定 N の g と一致：{ok}", flush=True)
    pickle.dump({"f": {str(k): str(v.to_sympy()) for k, v in f.items()}, "fcov": {str(k): str(v.to_sympy()) for k, v in fc.items()},
                 "g": {str(k): str(v.to_sympy()) for k, v in g.items()}}, open(os.path.join(OUT, f"largeN_2box_{tg}.pkl"), "wb"))
    return f, fc, g


if __name__ == "__main__":
    main(int(sys.argv[1]), next((int(a[4:]) for a in sys.argv[2:] if a.startswith("--k=")), 2))
