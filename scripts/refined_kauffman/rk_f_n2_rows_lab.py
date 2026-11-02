"""brane と image の行を別の記号にした n = 2 の記号パイプライン（続き 116、theory §12.46）。
rk_f_n2_rows_sym.build と同じ組み立てで、W_λ の brane 側は a2b, b4b, …、image 側は a2i, b4i, …、composite [α;β̄] は R_α^b R_β^i + (q−1)(t−1)G_{αβ}
（G は brane–image の連結部分で共通の記号）、[α;1̄] の縮約 d_α W₁ は brane の W₁、[1;ᾱ] のは image の W₁。
入れ替え Ĩ（brane ↔ image）のトレース：G を含まない単項式は brane の内容と image の内容が一致するときだけ残り ψ² になる。
実行（サンドボックス外）：PYTHONPATH=src:scripts/refined_kauffman python scripts/refined_kauffman/rk_f_n2_rows_lab.py [kmax]"""
import sys, pickle, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_f_n2_rows_sym as RS
from rk_f_n2_rows_sym import r, s, q, t, nu, Ser, rat_series, Psym, dimc, decomposition_flint, macdonald_P, P2, LAMS

BR = {k: sp.Symbol(k + "b") for k in ("a2", "a4", "b4", "b6", "b8", "c4", "c6", "c8")}
IM = {k: sp.Symbol(k + "i") for k in ("a2", "a4", "b4", "b6", "b8", "c4", "c6", "c8")}
GS = [RS.gadj, RS.g21, RS.g111, RS.g22, RS.g211, RS.g1111]
SYMS_LAB = list(BR.values()) + list(IM.values()) + [x for d in GS for x in d.values()]


def build_lab(kmax=14, verbose=True):
    RS.SYMS = SYMS_LAB
    RS.ZERO = RS.P(0)
    RS.KMAX = kmax
    P = RS.P
    t0 = time.time()
    qn = (nu - 1 / nu) / (s - 1 / s)
    kappa_e = (1 - RS.A) * (1 - q) / ((1 - t) * (1 - q * RS.A / t))
    Rb = {(1,): Ser({2: P(BR["a2"]), 4: P(BR["a4"])}), (2,): Ser({4: P(BR["b4"]), 6: P(BR["b6"]), 8: P(BR["b8"])}), (1, 1): Ser({4: P(BR["c4"]), 6: P(BR["c6"]), 8: P(BR["c8"])})}
    Ri = {(1,): Ser({2: P(IM["a2"]), 4: P(IM["a4"])}), (2,): Ser({4: P(IM["b4"]), 6: P(IM["b6"]), 8: P(IM["b8"])}), (1, 1): Ser({4: P(IM["c4"]), 6: P(IM["c6"]), 8: P(IM["c8"])})}
    qt1 = (q - 1) * (t - 1)
    Gadj = Ser({j: P(v) for j, v in RS.gadj.items()})
    Gc = {(2,): Ser({j: P(v) for j, v in RS.g21.items()}), (1, 1): Ser({j: P(v) for j, v in RS.g111.items()})}
    Gb = {((2,), (2,)): Ser({j: P(v) for j, v in RS.g22.items()}), ((2,), (1, 1)): Ser({j: P(v) for j, v in RS.g211.items()}), ((1, 1), (1, 1)): Ser({j: P(v) for j, v in RS.g1111.items()})}
    K = {lam: rat_series(Psym(lam)) for lam in Rb}
    kappa = rat_series(kappa_e)
    Kadj = rat_series(qn * qn - kappa_e)
    Wb = {lam: Rb[lam] * K[lam] for lam in Rb}
    Wi = {lam: Ri[lam] * K[lam] for lam in Ri}
    Wadj = (Rb[(1,)] * Ri[(1,)] + Gadj * qt1) * Kadj
    pair = {((), ()): Ser.const(1), ((1,), (1,)): Wadj + kappa}
    for lam in Rb:
        pair[(lam, ())] = Wb[lam]
        pair[((), lam)] = Wi[lam]
    for alpha in P2:
        dc = dimc(alpha, (1,))
        Sdc = rat_series(dc)
        d = (K[alpha] * K[(1,)] - Sdc) * rat_series(1 / Psym((1,)))
        pair[(alpha, (1,))] = (Rb[alpha] * Ri[(1,)] + Gc[alpha] * qt1) * Sdc + d * Wb[(1,)]      # 縮約 α/1 は brane の状態
        pair[((1,), alpha)] = (Rb[(1,)] * Ri[alpha] + Gc[alpha] * qt1) * Sdc + d * Wi[(1,)]      # 縮約は image の状態
    for alpha in P2:
        for beta in P2:
            key = (alpha, beta) if (alpha, beta) != ((1, 1), (2,)) else ((2,), (1, 1))
            c1, c0, Kb = decomposition_flint(*key)
            pair[(alpha, beta)] = (Rb[alpha] * Ri[beta] + Gb[key] * qt1) * rat_series(Kb) + rat_series(c1) * Wadj + rat_series(c0)
    if verbose:
        print(f"  核と対の値：{time.time() - t0:.0f} 秒", flush=True)
    PY = {(): {(): 1}}
    for n in (1, 2):
        for lam, v in macdonald_P(n, t, q, lambda c: sp.Rational(c.numerator, c.denominator)).items():
            PY[lam] = v
    T_ = lambda lam: tuple(sum(1 for x in lam if x > j) for j in range(lam[0])) if lam else ()
    Z = {}
    for (lam, mu), val in pair.items():
        sg = (-1)**(sum(lam) + sum(mu))
        for a_, x in PY[T_(lam)].items():
            for b_, y in PY[T_(mu)].items():
                Z[(a_, b_)] = Z.get((a_, b_), Ser()) + val * (sg * x * y)
    def mul2(A_, B_):
        out = {}
        for (a1, a2_), x in A_.items():
            for (b1, b2), y in B_.items():
                if sum(a1) + sum(b1) <= 2 and sum(a2_) + sum(b2) <= 2:
                    k = (tuple(sorted(a1 + b1, reverse=True)), tuple(sorted(a2_ + b2, reverse=True)))
                    out[k] = out.get(k, Ser()) + x * y
        return out
    X = {k: v for k, v in Z.items() if k != ((), ())}
    Lg, Xn = {}, dict(X)
    for n in range(1, 5):
        for k, v in Xn.items():
            Lg[k] = Lg.get(k, Ser()) + v * sp.Rational((-1)**(n + 1), n)
        if n < 4:
            Xn = mul2(Xn, X)
    if verbose:
        print(f"  log Z：{time.time() - t0:.0f} 秒", flush=True)
    x11 = Lg[((1,), (1,))]
    rs = ((1, 1), (2,))
    a = {(rho, sig): Lg.get((rho, sig), Ser()) for rho in rs for sig in rs}
    zq = r - 1 / r
    phi = {(1, 1): zq, (2,): r + 1 / r}
    CHI = {"2": {(1, 1): 1, (2,): 1}, "11": {(1, 1): 1, (2,): -1}}
    den = s * (t - nu**2) * (nu**2 - 1) * (nu**2 * q - 1)
    inv_den = rat_series(1 / den)
    N = {}
    for lab in LAMS:
        M = Ser()
        for (rho, sig), v in a.items():
            M = M + v * ((1 - t) * CHI[lab][rho] * CHI[lab][sig] / (zq * phi[rho] * phi[sig]))
        N[lab] = M * inv_den
    adams_coef = {lab: sp.Rational(-1, 2) * (1 - t) * CHI[lab][(2,)]**2 / (zq * phi[(2,)]**2) for lab in LAMS}
    if verbose:
        print(f"  N_λ（Adams の項を除く）：{time.time() - t0:.0f} 秒", flush=True)
    return N, x11, adams_coef, inv_den


if __name__ == "__main__":
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    t0 = time.time()
    N, x11, adams_coef, inv_den = build_lab(kmax)
    out = {"N": {(lab, k): str(v.as_expr()) for lab in LAMS for k, v in N[lab].c.items()},
           "x11": {k: str(v.as_expr()) for k, v in x11.c.items()},
           "adams_coef": {lab: str(v) for lab, v in adams_coef.items()},
           "inv_den": {k: str(v.as_expr()) for k, v in inv_den.c.items()}}
    f = f"data/refined_kauffman/daha/n2_rows_lab_m1_k{kmax}.pkl"
    pickle.dump(out, open(f, "wb"))
    print(f"保存：{f}（{time.time() - t0:.0f} 秒）", flush=True)
