r"""Hopf 絡み目：凍結した状態 C_22（電荷 ν~^β）と、固定セクターの指数 P_44（電荷 ν~^{2β}、被覆の Ĩ 不変な曲線）の関係（theory §12）。
  予想：凍結した状態の被覆は Ĩ 不変な曲線で、両端のブロックでは P_44|_{2β} = 単項式 × (C_22|_β / z_q)²（Krefl–Walcher 型の平方根）。
  入力：hopf_E_fixed.pkl（E_22、一意）、hopf_g_alld.pkl、hopf_fcov_sector_A4.pkl（固定セクター (2,2,2,2)）。"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp

r, s, nu = sp.symbols("r s nu", positive=True)
nut = sp.Symbol("nut", positive=True)
loc = {"r": r, "s": s, "nu": nu}
t = s**2; zq = r - 1 / r
E = {k: sp.sympify(v, locals=loc) for k, v in pickle.load(open("data/refined_kauffman/hopf_E_fixed.pkl", "rb")).items()}
g = {k: sp.sympify(v, locals=loc) for k, v in pickle.load(open("data/refined_kauffman/hopf_g_alld.pkl", "rb")).items()}
f4 = {k: sp.sympify(v, locals=loc) for k, v in pickle.load(open("data/refined_kauffman/hopf_fcov_sector_A4.pkl", "rb")).items()}


def blocks(X):
    X = sp.cancel(X.subs(nu, nut * s / r)); n_, d_ = sp.fraction(X)
    P = sp.Poly(sp.expand(n_), nut); dn = sp.Poly(d_, nut).monoms()[0][0]
    return {int(k[0] - dn): sp.cancel(c / (d_ / nut**dn)) for k, c in zip(P.monoms(), P.coeffs())}


if __name__ == "__main__":
    C22 = sp.cancel(g[(2, 2)] - E[(2, 2)] / (2 * (1 + t)))
    bC = blocks(C22 / zq)                                   # 商の境界モード 1 個を除く
    bP = blocks((1 - t) * f4[(2, 2, 2, 2)] / zq**3)         # 被覆の境界モード 3 個を除く
    print("C_22/z_q のブロック（ν~^β）:", {b: sp.factor(v) for b, v in sorted(bC.items())})
    for b in sorted(bC):
        if 2 * b not in bP:
            print(f"  β={b}: P_44 に ν~^{2*b} のブロックなし"); continue
        ratio = sp.cancel(bP[2 * b] / bC[b]**2)
        mono = sp.Poly(sp.fraction(ratio)[0], r, s).is_monomial and sp.Poly(sp.fraction(ratio)[1], r, s).is_monomial
        print(f"  β={b}: P_44|_{{2β}} / (C_22|_β)² = {ratio if mono else '単項式でない'}", flush=True)
    # (1,1) ↔ (2,2)：C_11 = g_11（E = 0）
    bC1 = blocks(g[(1, 1)] / zq); bP2 = blocks((1 - t) * f4[(1, 1, 1, 1)] / zq**3)
    print("C_11/z_q:", {b: sp.factor(v) for b, v in sorted(bC1.items())}, " P_22/z_q³:", {b: sp.factor(v) for b, v in sorted(bP2.items())})
    for b in sorted(bC1):
        if 2 * b in bP2:
            print(f"  β={b}: P_22|_{{2β}} / (C_11|_β)² = {sp.cancel(bP2[2 * b] / bC1[b]**2)}")
