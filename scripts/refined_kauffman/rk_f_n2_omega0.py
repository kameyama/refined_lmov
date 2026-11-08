"""単一の brane の 2 箱のキュムラントを Schur の λ に射影した Ω^{(0)}_λ(ν)（続き 119、theory §12.49）。
Z_b(Y) = 1 + Σ_λ W_λ P_λ(Y)（W_λ = R_λ Psym(λ)、trefoil の実際の行）を冪和基底に直し、log と Adams（p₂ の成分から ½ψ²(p₁ の成分)を引く）で
単一粒子の自由エネルギー f_ρ（ρ = (1), (1,1), (2)）を作り、Ω_λ = Σ_ρ χ_λ(ρ) f_ρ/φ_ρ とする。ν の最上の係数を ν^{−12} の X_λ・Ŷ_λ と比べる。
実行：PYTHONPATH=src:scripts/refined_kauffman python scripts/refined_kauffman/rk_f_n2_omega0.py"""
import sys
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_f_n2_rows_sym as RS
from rk_f_n2_rows_sym import r, s, q, t, nu, Psym, macdonald_P
import rk_f_n2_grows as G


def omega0():
    val = RS.actual_rows()
    R = {(1,): val[RS.a2] * nu**-2 + val[RS.a4] * nu**-4,
         (2,): val[RS.b4] * nu**-4 + val[RS.b6] * nu**-6 + val[RS.b8] * nu**-8,
         (1, 1): val[RS.c4] * nu**-4 + val[RS.c6] * nu**-6 + val[RS.c8] * nu**-8}
    W = {lam: sp.cancel(R[lam] * Psym(lam)) for lam in R}
    PY = {}
    for n in (1, 2):
        for lam, v in macdonald_P(n, t, q, lambda c: sp.Rational(c.numerator, c.denominator)).items():
            PY[lam] = {a_: sp.sympify(x) for a_, x in v.items()}
    T_ = lambda lam: tuple(sum(1 for x in lam if x > j) for j in range(lam[0])) if lam else ()
    # 冪和基底の係数：Z = 1 + Σ_λ W_λ P_{T(λ)}（rows_sym と同じく転置と符号 (−1)^{|λ|}）
    Zp = {}
    for lam, w in W.items():
        sg = (-1)**sum(lam)
        for a_, x in PY[T_(lam)].items():
            Zp[a_] = Zp.get(a_, 0) + sg * w * x
    psi2 = lambda e: sp.cancel(e.subs({r: r**2, s: s**2, nu: nu**2}, simultaneous=True))
    f = {(1,): Zp[(1,)], (1, 1): sp.cancel(Zp[(1, 1)] - Zp[(1,)]**2 / 2), (2,): sp.cancel(Zp[(2,)] - psi2(Zp[(1,)]) / 2)}
    zq = r - 1 / r
    phi = {(1, 1): zq, (2,): r + 1 / r}
    CHI = {"2": {(1, 1): 1, (2,): 1}, "11": {(1, 1): 1, (2,): -1}}
    Om = {lab: sp.cancel(sum(CHI[lab][rho] * f[rho] / phi[rho] for rho in ((1, 1), (2,)))) for lab in ("2", "11")}
    return Om, f, Zp


def top(e, nmax=16):
    """ν^{−1} の級数の先頭の係数 {k: c}。"""
    ser = RS.rat_series(e, nmax)
    return {k: v.as_expr() for k, v in ser.c.items()}


if __name__ == "__main__" and "--series" not in sys.argv:
    RS.SYMS = [sp.Symbol("_dummy")]; RS.ZERO = RS.P(0)
    Om, f, Zp = omega0()
    X = {lab: sp.cancel(G.X_of(lab).subs(RS.actual_rows())) for lab in ("2", "11")}
    for lab in ("11", "2"):
        ser = top(Om[lab])
        ks = sorted(ser)
        print(f"λ = {lab}：Ω^(0) の ν の冪 {ks[:6]}", flush=True)
        for k in ks[:3]:
            c = sp.factor(ser[k])
            print(f"   ν^-{k}：{c}", flush=True)
        k0 = ks[0]
        for name, Y in (("X_λ", X[lab]), ("Ŷ_λ", G.Y_LAM[lab])):
            print(f"   最上の係数/{name} = {sp.factor(sp.cancel(ser[k0] / Y))}", flush=True)
        # 2 乗の最上と、(s³/r)X² の比
        print(f"   (最上)²/((s³/r)X²) = {sp.factor(sp.cancel(ser[k0]**2 / (s**3 / r * X[lab]**2)))}", flush=True)
    for rho, v in f.items():
        print(f"f_{rho} の先頭：{ {k: sp.factor(c) for k, c in list(top(v).items())[:2]} }", flush=True)


def series_test():
    """Ω^{(0)} の準位を Ŷ の規格化に直し（Ω_k′）、σ 共変な部分 W₁₄・W₁₆ が Ŷ(ν)²F(ν) の形か（W₁₄ − 2ŶΩ₂′ が Ŷ で割れるか）を試す。"""
    RS.SYMS = [sp.Symbol("_dummy")]; RS.ZERO = RS.P(0)
    Om, f, Zp = omega0()
    val = RS.actual_rows(); valrf = G.rows_rf(val)
    D20 = G.blocks("data/refined_kauffman/daha/n2_rows_sym_m1_k20.pkl")
    for lab in ("11", "2"):
        Y = G.Y_LAM[lab]
        ser = top(Om[lab], 20)
        fac = sp.cancel(ser[4] / Y)                  # Ω^(0)_4 = fac · Ŷ
        Omn = {k: sp.cancel(ser[k] / fac) for k in sorted(ser)}   # 規格化：Omn[4] = Ŷ
        print(f"λ = {lab}：fac = {sp.factor(fac)}；Ω′_k = Ω^(0)_k/fac：", flush=True)
        for k in (2, 4, 6, 8):
            print(f"   Ω′_{k} = {sp.factor(Omn[k])}", flush=True)
        W = {}
        for K in (12, 14, 16):
            Wk, Ck = G.split_sigma(lab, K, valrf, D20)
            W[K] = G.to_sym(Wk)
        print(f"   W₁₂ − Ŷ² = {sp.cancel(W[12] - Y**2)}", flush=True)
        for sg in (1, -1):
            rem = sp.cancel(W[14] - 2 * sg * Y * Omn[6])
            qtnt = sp.cancel(rem / Y)
            lau = len(sp.Poly(sp.fraction(qtnt)[1], r, s).terms()) == 1
            print(f"   σ = {sg:+d}：W₁₄ − 2σŶΩ′₆ は Ŷ で割れる：{lau}" + (f"；f₁ := 商/Ŷ = {sp.factor(sp.cancel(qtnt / Y))}（Laurent：{len(sp.Poly(sp.fraction(sp.cancel(qtnt / Y))[1], r, s).terms()) == 1}）" if lau else f"；余りの因子 {sp.factor(sp.fraction(qtnt)[1])}"), flush=True)
        # 2 乗の級数そのもの：Ω′(ν)² の ν^{-8}（Ŷ²）, ν^{-10}（2ŶΩ′₆）, ν^{-12} と W の比
        Osq = {8: Y**2, 10: sp.expand(2 * Y * Omn[6]), 12: sp.expand(2 * Y * Omn[8] + Omn[6]**2)}
        for K, k2 in ((14, 10), (16, 12)):
            d = sp.cancel(W[K] - Osq[k2])
            print(f"   W_{K} − [Ω′²]_{k2} = {sp.factor(d)}"[:400], flush=True)
            print(f"      /Ŷ Laurent：{len(sp.Poly(sp.fraction(sp.cancel(d / Y))[1], r, s).terms()) == 1}；/(q−t) Laurent：{len(sp.Poly(sp.fraction(sp.cancel(d / (q - t)))[1], r, s).terms()) == 1}；/(q−1)(t−1) Laurent：{len(sp.Poly(sp.fraction(sp.cancel(d / ((q - 1) * (t - 1))))[1], r, s).terms()) == 1}", flush=True)


if __name__ == "__main__" and "--series" in sys.argv:
    series_test()
