"""第 4 期 V4 の残り：ツイスト結び目 K_p のニュートン形の係数 c_k を、GNSSS の入力の q-Pochhammer の和から直接書き下す。

入力（scripts/framed_knot/t3_pretzel.py の quadruple、s = 1）は p > 0 で
  P_r(K_p) = Σ_{j=0}^{r} B^{(r)}_j M_j^p,   M_0 = 1, M_1 = X = 𝐚²𝐭², M_2 = Y = 𝐚⁴𝐪⁴𝐭⁸
  B^{(r)}_j = (−𝐭)^{−r} Σ_{i=j}^{r} ϖ(1, i) V(i, j) W^{(r)}_i
（ϖ(p, i) = Σ_j V(i, j) M_j^p、V, W は q-Pochhammer と q 二項係数の積、ϖ(1, i) は単項式）。p < 0 は前因子 (−t_r t_c)^{−r·sgn} だけが変わる。
\\hat f_ρ = Σ_k π_{ρ,k} R_k（R = (P_[2], P_[1,1], P_[1]², P_[1](x²))、π_{ρ,k} は 2 箱の公式の有理関数）なので
  \\hat f_ρ(K_p) = ε Σ_j A_j m_j^p,   A_j = Σ_k π_{ρ,k} B^{R_k}_j    （ε = 1（p > 0）、𝐭²（p < 0）、reports/REPORT_phase4.md §5.1）
ニュートン形（並び μ = (1, X, X², Y, Ȳ)^{±1}、n0 = 1）の係数は
  c_k = ε Σ_{j ≥ k} A_j μ_j Π_{i<k} (μ_j − μ_i) = Σ_{j≥k} F_j / Π_{k≤i≠j}(μ_j − μ_i)   （F_j は residues() の留数）。
main() は、この閉じた式が直接計算と一致すること、c_4 = F_Ȳ、F_Ȳ = σ(F_Y)、c_3 の q 整数（sl_2 指標）による表示、
F_{X²} の Sym² 型の式、c_0 = \\hat f(K_{±1})・c_1 = \\hat f(K_{±2}) − \\hat f(K_{±1})、c_0..c_3 の sl_2 展開の符号を確かめる（logs/framed_knot/v4_ck.log）。
使い方：python scripts/framed_knot/v4_ck.py
"""
import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from v4_closed import RHOS, hat_direct, load_inputs, newton_coeffs, roots_for  # noqa: E402
from rcsknots.bps import C, C1, M_matrix, Q as QS, T as TS, f_formulas  # noqa: E402
from rcsknots.partitions import arm, boxes, leg  # noqa: E402

a, q, t = sp.symbols("a q t")                       # 太字 (𝐚, 𝐪, 𝐭)
Qq, tr, tc = sp.symbols("Q t_r t_c")
X, Y, Yb = a**2 * t**2, a**4 * q**4 * t**8, a**4 * q**-4 * t**4
ROOTS = {"1": sp.Integer(1), "X": X, "X2": X**2, "Y": Y, "Yb": Yb}
ORDER = ["1", "X", "X2", "Y", "Yb"]
SUBS = {Qq: q, tr: 1 / q, tc: t * q}


def qp(x, k, base):
    out = sp.Integer(1)
    for i in range(k):
        out *= 1 - x * base**i
    return out


def qbinom(n, k, base):
    return sp.cancel(qp(base, n, base) / (qp(base, k, base) * qp(base, n - k, base)))


def V(i, j):
    """ϖ(p, i) = Σ_j V(i, j) M_j^p（t3_pretzel.varpi の p によらない部分）。"""
    return a**i * tc**sp.Rational(i * (i - 1), 2) * (-1)**j * tc**(j * (j - 1)) \
        * (1 - a**2 * tr**2 * tc**(4 * j)) / qp(a**2 * tr**2 * tc**(2 * j), i + 1, tc**2) * qbinom(i, j, tc**2)


def M(j):
    return (a * tr * tc)**(2 * j) * tc**(2 * j * (j - 1))


def varpi1(i):
    return sum(V(i, j) * M(j) for j in range(i + 1))


def W(i, r):
    return (-a**2 * tc**(2 * (r - 1)))**(-i) * qp(-a**2 * Qq**-2 * tr * tc, i, tc**2) \
        * qp(-a**2 * Qq**2 * tr**3 * tc**(1 + 2 * r), i, tc**2) * qbinom(r, i, tc**2)


def blocks(r):
    """{j: B^{(r)}_j}（太字、p > 0 の前因子 (−t_r t_c)^{−r} 込み）。"""
    return {j: sp.factor(sp.cancel((sum(varpi1(i) * V(i, j) * W(i, r) for i in range(j, r + 1)) * (-tr * tc)**(-r)).subs(SUBS)))
            for j in range(r + 1)}


def parts():
    """{ρ: [π_{ρ,k}]}（太字）。\\hat f_ρ = Σ_k π_{ρ,k} R_k、R = (P_[2], P_[1,1], P_[1]², P_[1](x²))（いずれも reduced、太字）。"""
    A = sp.Symbol("A")

    def unknot(lam, d=1):
        k = sum(lam)
        out = TS**(d * k) * A**(-d * k)
        for i, j in boxes(lam):
            out *= (TS**(2 * d * i) - A**(2 * d) * QS**(2 * d * j)) / (1 - QS**(2 * d * arm(lam, i, j)) * TS**(2 * d * (leg(lam, i, j) + 1)))
        return out

    F = f_formulas(2)
    hat = M_matrix(2).inv() * sp.Matrix([F[(2,)][0], F[(1, 1)][0]])
    syms = [C((2,)), C((1, 1)), C((1,)), C1(2)]
    us = [unknot((2,)), unknot((1, 1)), unknot((1,)), unknot((1,), 2)]
    out = {}
    for rho, h in zip(RHOS, hat):
        num, den = sp.fraction(sp.together(h))
        P = sp.Poly(sp.expand(num), *syms)
        ps = [0, 0, 0, 0]
        for mon, c in P.terms():
            k = {(1, 0, 0, 0): 0, (0, 1, 0, 0): 1, (0, 0, 2, 0): 2, (0, 0, 0, 1): 3}[mon]
            ps[k] = c * (us[k] if k != 2 else us[2]**2) / den
        # 小文字 → 太字：A² = −𝐚²𝐭, Q = −𝐪𝐭, T = 𝐪（A は偶数冪でだけ現れる）
        Ab = sp.Symbol("Ab")
        out[rho] = [sp.factor(sp.cancel(sp.expand(sp.cancel(x).subs(A, sp.sqrt(Ab))).subs({Ab: -a**2 * t, QS: -q * t, TS: q})))
                    for x in ps]
    return out


def input_blocks():
    """{入力名: {根: ブロック}}（太字、p > 0）。"""
    b1, b2 = blocks(1), blocks(2)
    mir = lambda e: sp.factor(sp.cancel(e.subs({q: 1 / (q * t)}, simultaneous=True)))   # noqa: E731  (e_𝐪, e_𝐭) ↦ (−e_𝐪, e_𝐭 − e_𝐪)
    adams = lambda e: sp.factor(sp.cancel(e.subs({a: a**2, q: q**2, t: -t**2}, simultaneous=True)))  # noqa: E731  P_[1](x²)
    return {
        "P2": {"1": b2[0], "X": b2[1], "Y": b2[2]},
        "P11": {"1": mir(b2[0]), "X": mir(b2[1]), "Yb": mir(b2[2])},
        "P1sq": {"1": sp.factor(b1[0]**2), "X": sp.factor(2 * b1[0] * b1[1]), "X2": sp.factor(b1[1]**2)},
        "P1d": {"1": adams(b1[0]), "X2": adams(b1[1])},
    }


def amplitudes(PI, IB):
    """{ρ: {根: A_j}}。"""
    names = ["P2", "P11", "P1sq", "P1d"]
    return {rho: {r: sp.factor(sp.cancel(sum(PI[rho][k] * IB[n].get(r, 0) for k, n in enumerate(names)))) for r in ORDER} for rho in RHOS}


def ck(Aj, sign):
    eps = 1 if sign == 1 else t**2
    mu = {r: ROOTS[r]**sign for r in ORDER}
    return [sp.factor(sp.cancel(eps * sum(Aj[ORDER[j]] * mu[ORDER[j]] * sp.prod([mu[ORDER[j]] - mu[ORDER[i]] for i in range(k)])
                                          for j in range(k, 5)))) for k in range(5)]


def to_dict(e):
    e = sp.cancel(e)
    n, d = sp.fraction(e)
    dp = sp.Poly(d, a, q, t)
    assert len(dp.terms()) == 1, "単項式でない分母"
    (de, dc), = dp.terms()
    return {(i - de[0], j - de[1], k - de[2]): sp.Rational(c, dc) for (i, j, k), c in sp.Poly(sp.expand(n), a, q, t).terms()}


def residues(Aj, sign):
    """F_j = ε A_j μ_j Π_{i≠j}(μ_j − μ_i)（生成関数 Σ_{n≥1} s_n z^{−n} = N(z)/Π(z − μ_i) の分子 N の節点での値）。
    c_k = N[μ_k, …, μ_4] = Σ_{j≥k} F_j / Π_{k≤i≠j}(μ_j − μ_i)（差分商）。"""
    eps = 1 if sign == 1 else t**2
    mu = {r: ROOTS[r]**sign for r in ORDER}
    return {j: sp.factor(sp.cancel(eps * Aj[j] * mu[j] * sp.prod([mu[j] - mu[i] for i in ORDER if i != j]))) for j in ORDER}


# ---- (a, v, τ) 座標：v = 𝐪²𝐭（σ は v ↦ 1/v、τ = 𝐭 は固定）。単項式 𝐚^α 𝐪^β 𝐭^γ = 𝐚^α v^{β/2} τ^{γ−β/2} ----------------

v, tau = sp.symbols("v tau")


def to_avt(e):
    d = to_dict(e)
    assert all(be % 2 == 0 for (_, be, _) in d)
    return sum(c * a**al * v**(be // 2) * tau**(ga - be // 2) for (al, be, ga), c in d.items())


def avt_terms(e):
    """(a, v, τ) の Laurent 多項式 → {(α, n, γ): c}。"""
    e = sp.cancel(e)
    n_, d_ = sp.fraction(e)
    (de, dc), = sp.Poly(d_, a, v, tau).terms()
    return {(i - de[0], j - de[1], k - de[2]): sp.Rational(c, dc) for (i, j, k), c in sp.Poly(sp.expand(n_), a, v, tau).terms()}


def qint(n):
    """[n]_v = (v^n − v^{−n})/(v − v^{−1})（sl_2 の n 次元表現の指標。n > 0 で非負係数、n < 0 で −[|n|]_v）。"""
    return sp.cancel((v**n - v**-n) / (v - 1 / v))


S_Y = (1 + a**2 * tau) * (1 + a**2 * tau**3) * (1 + a**2 * tau**2 * v) * (1 + a**2 * tau**2 / v)   # σ 不変


def c3_structure(F, sign, c3):
    """F_Y = −κ_F S_Y (v + v^{−1}) h(v)（h は正係数）と書き、c_3 = (F_Y − σF_Y)/(μ_Y − μ_Ȳ) = −(κ_F/λ) S_Y Σ_n h_n [n]_v
    （μ_Y − μ_Ȳ = λ (v² − v^{−2})）を確かめる。全ての n が同符号なら c_3 は一様な符号。"""
    FY = to_avt(F["Y"])
    h = sp.cancel(-FY / (S_Y * (v + 1 / v)))
    hd = avt_terms(h)
    lam = sp.cancel((to_avt(ROOTS["Y"]**sign) - to_avt(ROOTS["Yb"]**sign)) / (v**2 - v**-2))
    rec = -S_Y / lam * sum(c * a**al * tau**ga * qint(n) for (al, n, ga), c in hd.items())
    ok = sp.expand(sp.cancel(rec - to_avt(c3))) == 0
    return {"h の係数は全て正": all(c > 0 for c in hd.values()), "h の v の指数": sorted({n for (_, n, _) in hd}),
            "λ": lam, "c_3 の再構成が一致": ok}


def sym2_check(rho, F, PI):
    """F_{X²} = (X²−Y)(X²−Ȳ) X [π_{1²}(X+1) x² + π_{1d}(X−1) x^{(2)}]·(p<0 の単項式)、x = B^{[1]}_X X(X−1)（[1] の上端の留数）。
    π の分母 (1−𝐪⁴)(1−𝐪⁴𝐭⁴) は (X²−Ȳ)(X²−Y) と打ち消し、Z = (1+𝐚²𝐭)x とおくと
    F_{X²} ∝ [(1+X)(1+𝐪²)(1+𝐪²𝐭²) Z² ± (1−X)(1−𝐪²)(1−𝐪²𝐭²) Z^{(2)}]/2（Sym²/Λ² 型。単項式ごとに正）。ここでは p > 0 の恒等式を確かめる。"""
    b1 = blocks(1)
    x = sp.cancel(b1[1] * X * (X - 1))
    xd = sp.cancel(x.subs({a: a**2, q: q**2, t: -t**2}, simultaneous=True))
    Fx2 = (X**2 - Y) * (X**2 - Yb) * X * (PI[rho][2] * (X + 1) * x**2 + PI[rho][3] * (X - 1) * xd)
    return sp.cancel(Fx2 - F["X2"]) == 0, sp.factor(x)


def sl2_positive(c):
    """σ 不変な c を (𝐚, τ) の行ごとに v の sl_2 指標 [k]_v で展開したときの係数 f_k = e_{k−1} − e_{k+1} が全て c の符号と同じか。"""
    rows = {}
    for (al, n, ga), x in avt_terms(to_avt(c)).items():
        rows.setdefault((al, ga), {})[n] = x
    sym = all(r.get(-n, 0) == x for r in rows.values() for n, x in r.items())
    fs = [r.get(k - 1, 0) - r.get(k + 1, 0) for r in rows.values() for k in {abs(n) + 1 for n in r}]
    return sym, sum(1 for f in fs if f > 0), sum(1 for f in fs if f < 0)


def main():
    import time
    t0 = time.time()
    PI = parts()
    for rho in RHOS:
        print(f"π_{rho}：", [str(x) for x in PI[rho]], flush=True)
    IB = input_blocks()
    for n, bl in IB.items():
        for r, e in bl.items():
            print(f"入力 {n} のブロック {r}: {e}", flush=True)
    print("ϖ(1, i) =", [sp.factor(sp.cancel(varpi1(i).subs(SUBS))) for i in range(3)], "（単項式）", flush=True)
    print("V(2,2) =", sp.factor(sp.cancel(V(2, 2).subs(SUBS))), "  W(2,2) =", sp.factor(sp.cancel(W(2, 2).subs(SUBS))), flush=True)
    AJ = amplitudes(PI, IB)
    for rho in RHOS:
        for r in ORDER:
            print(f"A_{r}({rho}) = {AJ[rho][r]}", flush=True)
    cache = load_inputs()
    for sign in (1, -1):
        H = {n: hat_direct(sign * n, cache) for n in range(1, 6)}
        for rho in RHOS:
            cs = ck(AJ[rho], sign)
            direct = newton_coeffs([H[n][rho] for n in range(1, 6)], roots_for(sign))
            print(f"\n=== p の符号 {sign} ρ={rho} ===", flush=True)
            print("閉じた式の c_k が直接計算のニュートン係数と全係数一致:", all(to_dict(c) == direct[k] for k, c in enumerate(cs)), flush=True)
            print("c_0 = \\hat f(K_{±1})、c_1 = \\hat f(K_{±2}) − \\hat f(K_{±1}):",
                  to_dict(cs[0]) == H[1][rho],
                  to_dict(cs[1]) == {k: H[2][rho].get(k, 0) - H[1][rho].get(k, 0) for k in set(H[1][rho]) | set(H[2][rho]) if H[2][rho].get(k, 0) != H[1][rho].get(k, 0)}, flush=True)
            F = residues(AJ[rho], sign)
            for j in ORDER:
                print(f"  F_{j} = {F[j]}", flush=True)
            sig = lambda e: sp.cancel(e.subs({q: 1 / (q * t)}, simultaneous=True))   # noqa: E731
            print("  F_Ȳ = σ(F_Y):", sp.cancel(sig(F["Y"]) - F["Yb"]) == 0, "  c_4 = F_Ȳ:", sp.cancel(cs[4] - F["Yb"]) == 0, flush=True)
            print("  c_3 の構造:", c3_structure(F, sign, cs[3]), flush=True)
            if sign == 1:
                ok, x = sym2_check(rho, F, PI)
                print("  F_{X²} の Sym² 型の式が一致:", ok, "  x =", x, flush=True)
            for k in range(4):
                print(f"  c_{k} の sl_2 展開（σ 不変, 正の係数の数, 負の係数の数）:", sl2_positive(cs[k]), flush=True)
    import resource
    print(f"{time.time() - t0:.0f} 秒、最大メモリ {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20:.0f} MB", flush=True)


if __name__ == "__main__":
    main()
