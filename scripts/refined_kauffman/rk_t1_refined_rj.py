"""T1 の refined 版（続き 108）：随伴の refined composite Rosso–Jones の仮説
  W_adj = (枠) [Σ_{μ,ν ⊢ r} c_μ c_ν T_μ^{s/r} T_ν^{s/r} dim^{ref}[μ;ν] + c₀(q,t)]
（c_μ：p_r = Σ c_μ P_μ(q,t) の係数、T_μ = q^{n(μ′)} t^{−n(μ)}、dim^{ref}：kauffman_refined.dim_composite、A = ν² = (u/s)²）
を DAHA の W_adj = R_adj·K_adj（K_adj = q_n² − κ）と比べる。枠の単項式は最上行で合わせ、差が ν によらない定数 c₀ になるかを見る。
flint の有理関数（r = q^{1/2}、s = t^{1/2}、u = νs）で計算する。T^{s/r} の冪が整数になる 2 本撚り（r = 2）が対象。
実行：PYTHONPATH=src:scripts/refined_kauffman python scripts/refined_kauffman/rk_t1_refined_rj.py"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rcsknots.refined_kauffman.flint_rf3 import R, S, U, from_sympy, conv
from rcsknots.macdonald import macdonald_P
from rcsknots.partitions import partitions, n_stat, conjugate
from rcsknots.refined_kauffman import kauffman_refined as KR
import rk_f_frozen_basis as FB
import rk_torus_tag as TG

q, t, A = KR.q, KR.t, KR.A
r_, s_, u_ = sp.symbols("r s u", positive=True)
nu = U / S
qq, tt = R * R, S * S
OUT = "data/refined_kauffman/daha"
SUB = {q: r_**2, t: s_**2, A: (u_ / s_)**2}


def adams_coeffs(r):
    """p_r = Σ_μ c_μ P_μ（q, t の有理式）。"""
    cv = lambda c: sp.Rational(c.numerator, c.denominator)
    P = macdonald_P(r, q, t, cv); parts = list(P)
    cs = sp.symbols(f"c0:{len(parts)}")
    eqs = [sum(c * P[lam].get(n_, 0) for c, lam in zip(cs, parts)) - (1 if n_ == (r,) else 0) for n_ in partitions(r)]
    sol = sp.solve(eqs, cs, dict=True)[0]
    return {lam: sp.factor(sol[c]) for c, lam in zip(cs, parts)}


def rf(e):
    return from_sympy(sp.sympify(e).subs(SUB))


def dim_ref(al, be):
    if not al and not be:
        return conv(1)
    return rf(KR.dim_composite(list(al), list(be)))


def Tpow(mu, s, r):
    """T_μ^{s/r} = q^{n(μ′)s/r} t^{−n(μ)s/r}。"""
    eq_, et_ = sp.Rational(2 * n_stat(conjugate(mu)) * s, r), sp.Rational(-2 * n_stat(mu) * s, r)
    assert eq_.q == 1 and et_.q == 1, "分数冪（r ≠ 2）は未対応"
    return R**int(eq_) * S**int(et_)


def hypothesis(r, s, image="same"):
    c = adams_coeffs(r)
    tot = conv(0)
    for mu, cm in c.items():
        for nv, cn in c.items():
            cnu = cn if image == "same" else cn.subs({q: t, t: q}, simultaneous=True)
            tot = tot + rf(cm * cnu) * Tpow(mu, s, r) * Tpow(nv, s, r) * dim_ref(mu, nv)
    return tot, c


def daha_W(m, k=2):
    radj = pickle.load(open(f"{OUT}/Radj_A_{TG.tag(m, k)}.pkl", "rb"))
    Rr = sum((from_sympy(sp.sympify(v, locals={"r": r_, "s": s_})) * nu**int(j) for j, v in radj.items()), conv(0))
    qn = (nu - 1 / nu) / (S - 1 / S)
    kappa = (1 - nu * nu) * (1 - qq) / ((1 - tt) * (1 - qq * nu * nu / tt))
    return Rr * (qn * qn - kappa)


def rows(e):
    X = sp.cancel(e.to_sympy().subs(u_, FB.nu * s_))
    return FB.nu_blocks(X)


def _lead_u(poly):
    """flint の多項式の u の最高次：(u の次数, その係数の多項式（u を含まない）)。"""
    from rcsknots.refined_kauffman.flint_rf3 import CTX
    d = poly.to_dict()
    top = max(e[2] for e in d)
    return top, CTX.from_dict({(e[0], e[1], 0): c for e, c in d.items() if e[2] == top})


def leading(e):
    """u → ∞ の主要項：(係数の RF（r, s）, u の次数)。"""
    from rcsknots.refined_kauffman.flint_rf3 import RF
    tn, cn = _lead_u(e.n); td, cd = _lead_u(e.d)
    return RF(cn, cd), tn - td


def has_u(e):
    return any(k[2] != 0 for k in e.n.to_dict()) or any(k[2] != 0 for k in e.d.to_dict())


if __name__ == "__main__":
    for m, (r, s) in ((1, (2, 3)), (2, (2, 5)), (3, (2, 7))):
        print(f"== T({r},{s})", flush=True)
        W = daha_W(m)
        lw, dw = leading(W)
        for image in ("same", "swap"):
            H, c = hypothesis(r, s, image)
            if image == "same":
                print("  c_μ：", {k_: str(v) for k_, v in c.items()}, flush=True)
            lh, dh = leading(H)
            ratio = lw / lh
            print(f"  仮説（image={image}）：主要項の比 W/H = {sp.factor(ratio.to_sympy())} × u^{dw - dh}", flush=True)
            D = W - ratio * U**(dw - dh) * H
            print(f"    W − 枠·H は u に依らない：{not has_u(D)}", flush=True)
            if not has_u(D):
                print(f"    c₀ = {sp.factor(D.to_sympy())}", flush=True)
            else:
                ld, dd = leading(D)
                print(f"    残りの主要項：{sp.factor(ld.to_sympy())} × u^{dd}（W の主要項は u^{dw}）"[:300], flush=True)
