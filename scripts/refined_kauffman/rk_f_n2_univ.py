r"""n = 2（trefoil）の N_λ を式から読む（theory §12.25、progress.md 続き 91）：結び目のデータ（R_b, R_adj, R_c, R_λ）を 1 にした
「unknot の参照」の N^{univ}_λ を、枠 w = (t/q)^{1/2}（trefoil と同じ、枠付き unknot）と w = 1（枠なし unknot）で計算する。
  univ：枠付き unknot の N^{univ}_λ の有理式の分母の因数分解と、ν^{-k} の形式級数のブロック（分母＝塔、分子）。
  w1  ：枠なし unknot で x₁₁ = 0、N^{univ}_λ = 0 になることと、n = 1 の unknot（w = 1 で 0、枠付きで非 Laurent）の確認。
使い方：python rk_f_n2_univ.py univ | w1"""
import sys, time, inspect, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_f_cover_fixed22 as CF
import rk_f_twisted22_lam as L2
import rk_c_2box_largeN as L
from rcsknots.refined_kauffman.flint_rf3 import conv

r_, s_, u_ = sp.symbols("r s u", positive=True)
NU, y = sp.symbols("nu y")
OUT = "data/refined_kauffman/daha"


def set_unknot(w_one=False):
    """結び目のデータを 1 にする（w_one なら枠 w も 1 に）。"""
    CF.load_poly_A = lambda path, unit: CF.ONE
    L.Rdata_nu = lambda c, knot="T2_3": conv(1)
    if w_one:
        src = inspect.getsource(CF.pair_values).replace("w = S / R", "w = ONE")
        ns = dict(CF.__dict__); exec(src, ns); CF.pair_values = ns["pair_values"]


def hatted_N():
    """rk_f_twisted22_lam.cover と同じ hatted の P̂_λλ を N_λ = P̂/(s(t − ν²)(ν² − 1)(ν²q − 1)) にして返す（flint RF）。"""
    P, x11, Lg = CF.fixed_sector()
    t, zq = CF.t, CF.R - 1 / CF.R
    phi = {(1, 1): zq, (2,): CF.R + 1 / CF.R}
    a = {(rho, sig): Lg.get((rho, sig), CF.ZERO) for rho in ((1, 1), (2,)) for sig in ((1, 1), (2,))}
    a[((2,), (2,))] = a[((2,), (2,))] - CF.half * x11.adams(2)
    ah = {k: (1 - t) * v / (zq * phi[k[0]] * phi[k[1]]) for k, v in a.items()}
    Mh = L2.schur_split(ah)
    den = CF.S * (t - CF.nu**2) * (CF.nu**2 - 1) * (CF.nu**2 * CF.q - 1)
    return {"2": Mh[("2", "2")] / den, "11": Mh[("11", "11")] / den}, x11


def series_y(ex, nmax):
    """ex（y の有理式、係数は r, s）の y = 0 での冪級数（多項式の割り算）。{冪: 係数}。"""
    num, dd_ = sp.fraction(sp.cancel(sp.together(ex)))
    Pn, Pd = sp.Poly(sp.expand(num), y), sp.Poly(sp.expand(dd_), y)
    nd = {k[0]: v for k, v in zip(Pn.monoms(), Pn.coeffs())}; dd = {k[0]: v for k, v in zip(Pd.monoms(), Pd.coeffs())}
    d0 = min(dd); dd = {k - d0: v for k, v in dd.items()}
    c = {}
    for k in range(nmax + d0 + 1):
        v = nd.get(k, 0) - sum(dd[j] * c[k - j] for j in dd if j >= 1 and k - j in c)
        c[k] = sp.cancel(v / dd[0])
    return {k - d0: c[k] for k in c}


def univ():
    set_unknot(False)
    Ns, _ = hatted_N()
    out = {}
    for lab, N in Ns.items():
        print(f"λ = {lab}：枠付き unknot の N^univ は Laurent {N.is_laurent()}", flush=True)
        ex = sp.cancel(N.to_sympy().subs(u_, NU * s_))
        num, dd = sp.fraction(ex)
        print(f"   分母：{sp.factor(dd)}", flush=True)
        Pn = sp.Poly(sp.expand(num), NU, r_, s_)
        print(f"   分子：ν の次数 {min(m[0] for m in Pn.monoms())}..{max(m[0] for m in Pn.monoms())}、項数 {len(Pn.monoms())}、Σ|c| {sum(abs(c) for c in Pn.coeffs())}", flush=True)
        out[lab] = {"num": str(num), "den": str(dd)}
        t0 = time.time()
        ser = series_y(N.to_sympy().subs(u_, s_ / y), 26)
        print(f"   ν^(-k) の形式級数の係数（{time.time() - t0:.0f} 秒）：", flush=True)
        for k in sorted(ser):
            if k < 0 or ser[k] == 0:
                continue
            c = sp.factor(ser[k]); n_, d_ = sp.fraction(sp.cancel(c))
            print(f"     ν^-{k}：分母 {sp.factor(d_)}、分子 {str(sp.factor(n_))[:200]}", flush=True)
    pickle.dump(out, open(f"{OUT}/n2_N_univ_m1.pkl", "wb"))
    print(f"保存：{OUT}/n2_N_univ_m1.pkl", flush=True)


def w1():
    RR, SS = L.R, L.S; qq, tt = RR * RR, SS * SS; nu_ = L.nu; A = nu_ * nu_
    qn = (nu_ - 1 / nu_) / (SS - 1 / SS)
    kappa = (1 - A) * (1 - qq) / ((1 - tt) * (1 - qq * A / tt))
    D = SS**2 * (tt - nu_**2) * (nu_**2 - 1) * (qq - 1) * (nu_**2 * qq - 1)
    for name, w in (("w = 1", conv(1)), ("w = (t/q)^{1/2}", SS / RR)):
        W1 = w * L.Psym((1,))
        N1 = (1 - tt) * (qn * qn - W1 * W1) / D                  # R_adj = 1：K_adj + κ = q_n²
        print(f"n = 1 の unknot（{name}）：N は Laurent {N1.is_laurent()}；{sp.factor(N1.to_sympy().subs(u_, NU * s_))}", flush=True)
    set_unknot(True)
    Ns, x11 = hatted_N()
    print(f"n = 2 の unknot（w = 1）：x₁₁ = 0 ? {x11.is_zero()}", flush=True)
    for lab, N in Ns.items():
        print(f"  λ = {lab}：N^univ(w = 1) = 0 ? {N.is_zero()}、Laurent {N.is_laurent()}", flush=True)


if __name__ == "__main__":
    {"univ": univ, "w1": w1}[sys.argv[1]]()
