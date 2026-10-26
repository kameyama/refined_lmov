r"""被覆 U(N_A) の brane 2・image 2 の固定セクター（大きな N、trefoil）と、凍結した状態との β ↔ 2β の対応（theory §12.6、progress.md 続き 68）。
  被覆の生成関数（brane の probe Y、image の probe Y′、双対 Cauchy）：
    Z(Y, Y′) = Σ_{λ,μ} (−1)^{|λ|+|μ|} ⟨P_λ(x) P_μ(x^{−1})⟩ P_{λ′}(Y; t, q) P_{μ′}(Y′; t, q)。
  ⟨P_λ⟩ = W_λ（一様な規格化）、⟨P_1 P_1^∨⟩ = W_adj + κ、⟨P_α P_1^∨⟩ = W_{[α;1]} + d_α W_1（rk_c_3box_largeN と同じ）、
  ⟨P_α P_β^∨⟩ = W_{[α;β]} + c₁ W_adj + c₀（|α| = |β| = 2、W_{[α;β]} = R_{[α;β]}(A) K_{[α;β]}、R は Rcomp22_A_m1_{α}_{β}.pkl、
  K は kauffman_refined.dim_composite、c₁, c₀ は rk_c_composite22_run.decomposition の式を A = t^{N_A} の閉じた形にしたもの）。
  固定セクター：P = (1 − t) [log Z]_{p_1(Y)² p_1(Y′)²}（Adams の項はこの成分に入らないので f^{cov} の成分と同じ）。
  凍結した状態：trefoil の 2 箱 p₁² 方向の hatted (g_(2) + g_(1,1))/z_q の ν 奇（c = 1）のブロック（rk_f_frozen_basis.knot_two_box）。
変数は (r, s, u) = (q^{1/2}, t^{1/2}, U)、ν = U/s、A = ν²（rcsknots.refined_kauffman.flint_rf3）。"""
import sys, os, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from fractions import Fraction
from rcsknots.refined_kauffman.flint_rf3 import RF, R, S, U, from_sympy, conv
from rcsknots.macdonald import macdonald_P
import rk_c_2box_largeN as L

r_, s_, u_ = sp.symbols("r s u", positive=True)
q, t = R * R, S * S
nu = U / S
A = nu * nu
OUT = "data/refined_kauffman/daha"
ZERO, ONE = conv(0), conv(1)
half = conv(Fraction(1, 2))
P2 = ((2,), (1, 1))


def sym_rf(e):
    """sympy（r, s, A または q, t, A）→ flint_rf3（A = ν²）。"""
    qs, ts, As = sp.symbols("q t A", positive=True)
    return from_sympy(sp.cancel(sp.sympify(e).subs({qs: r_**2, ts: s_**2, As: (u_ / s_)**2})))


def load_poly_A(path, unit):
    """{j: c_j} を Σ c_j x^j（unit = 1：x = A^{1/2} = ν、unit = 2：x = A = ν²）に。"""
    d = pickle.load(open(path, "rb"))
    c = d["c"] if "c" in d else d
    return sum((from_sympy(sp.sympify(v, locals={"r": r_, "s": s_})) * nu**(unit * int(j)) for j, v in c.items()), ZERO)


# ---- 分解の係数（A の閉じた形） -------------------------------------------------------------------
def psum(d, prime, inverse):
    """p_d(x)（x = (t^{−(N−1)}, …, 1)、prime なら最後の座標を q に、inverse なら 1/x）。A = t^N。"""
    if not inverse:
        base = (1 - A**(-d)) / (1 - t**(-d))
        return base + (q**d - 1 if prime else ZERO)
    base = (1 - A**d) / (1 - t**d)
    return base + (q**(-d) - 1 if prime else ZERO)


_MP = {}


def P_pt(lam, prime, inverse):
    n = sum(lam)
    if n not in _MP:
        _MP[n] = macdonald_P(n, q, t, conv)
    tot = ZERO
    for rho, c in _MP[n][tuple(lam)].items():
        v = ONE
        for d in rho:
            v = v * psum(d, prime, inverse)
        tot = tot + c * v
    return tot


def e_val(alpha, beta):
    """e_b = Σ q^{b_i} t^{N−i}（b = (α, 0, …, 0, −β の逆順)）。"""
    eN = (A - 1) / (t - 1)
    e = eN
    for i, a in enumerate(alpha):                     # 位置 i+1：t^{N−1−i} = A t^{−1−i}
        e = e + (q**a - 1) * A * t**(-1 - i)
    for k, b in enumerate(beta):                      # 位置 N−k：t^{k}
        e = e + (q**(-b) - 1) * t**k
    return e


def decomposition_A(alpha, beta):
    from rcsknots.refined_kauffman.kauffman_refined import dim_composite
    Kb = sym_rf(dim_composite(alpha, beta))
    qn = (nu - 1 / nu) / (S - 1 / S)
    kappa = (1 - A) * (1 - q) / ((1 - t) * (1 - q * A / t))
    Kadj = qn * qn - kappa
    eN = (A - 1) / (t - 1)
    f = P_pt(alpha, False, False) * P_pt(beta, False, True)
    fp = P_pt(alpha, True, False) * P_pt(beta, True, True)
    eb, eadj = e_val(alpha, beta), e_val((1,), (1,))
    r1 = f - Kb
    r2 = eN * fp - eb * Kb
    det = Kadj * eN - eadj * Kadj
    c1 = (r1 * eN - r2) / det
    c0 = r1 - c1 * Kadj
    return c1, c0, Kb, Kadj


def check_decomposition():
    """A の閉じた形の (c₁, c₀, K_b) が、N ごとの式（rk_c_composite22_run.decomposition）と A = t^N で一致するか（N = 5, 7）。"""
    import rk_c_composite22_run as C
    for alpha, beta in ((P2[0], P2[0]), (P2[0], P2[1]), (P2[1], P2[1])):
        cA = decomposition_A(alpha, beta)
        for N in (5, 7):
            cN = C.decomposition(N, alpha, beta)
            sub = lambda e: sp.cancel(e.to_sympy().subs(u_, s_**(N + 1)))       # A = U²/t = t^N ⇔ U = s^{N+1}
            ok = [sp.cancel(sub(a) - b.to_sympy()) == 0 for a, b in zip(cA[:3], cN[:3])]
            print(f"[{''.join(map(str, alpha))};{''.join(map(str, beta))}] N={N}：c₁, c₀, K_b が N ごとの式と一致 {ok}", flush=True)


# ---- ⟨P_λ P_μ^∨⟩ ----------------------------------------------------------------------------------
def pair_values():
    from rcsknots.refined_kauffman.kauffman_refined import dim_composite
    w = S / R
    Wo = {lam: w**sum(lam) * L.Rdata_nu("".join(map(str, lam))) * L.Psym(lam) for lam in ((1,), (2,), (1, 1))}
    Radj = load_poly_A(os.path.join(OUT, "Radj_A_m1.pkl"), 1)
    qn = (nu - 1 / nu) / (S - 1 / S)
    kappa = (1 - A) * (1 - q) / ((1 - t) * (1 - q * A / t))
    Wadj = Radj * (qn * qn - kappa)
    pair = {((), ()): ONE, ((1,), (1,)): Wadj + kappa}
    for lam, v in Wo.items():
        pair[(lam, ())] = pair[((), lam)] = v
    for alpha in P2:                                       # [α;1]
        Rc = load_poly_A(os.path.join(OUT, f"Rcomp_A_m1_{''.join(map(str, alpha))}.pkl"), 1)
        dimc = sym_rf(dim_composite(alpha, (1,)))
        d = (L.Psym(alpha) * L.Psym((1,)) - dimc) / L.Psym((1,))
        pair[(alpha, (1,))] = pair[((1,), alpha)] = Rc * dimc + d * Wo[(1,)]
    for alpha in P2:                                       # [α;β]、|α| = |β| = 2
        for beta in P2:
            key = (alpha, beta) if (alpha, beta) != ((1, 1), (2,)) else ((2,), (1, 1))
            tag = f"{''.join(map(str, key[0]))}_{''.join(map(str, key[1]))}"
            Rb = load_poly_A(os.path.join(OUT, f"Rcomp22_A_m1_{tag}.pkl"), 2)
            c1, c0, Kb, _ = decomposition_A(*key)
            pair[(alpha, beta)] = Rb * Kb + c1 * Wadj + c0
    return pair


# ---- 二変数の冪和（キーは (Y の分割, Y′ の分割)） ------------------------------------------------------
def kadd(acc, k, v):
    acc[k] = acc.get(k, ZERO) + v


def mul2(A_, B_, cap=(2, 2)):
    out = {}
    for (a1, a2), x in A_.items():
        for (b1, b2), y in B_.items():
            if sum(a1) + sum(b1) <= cap[0] and sum(a2) + sum(b2) <= cap[1]:
                kadd(out, (tuple(sorted(a1 + b1, reverse=True)), tuple(sorted(a2 + b2, reverse=True))), x * y)
    return out


def cover_Z():
    pair = pair_values()
    PY = {(): {(): ONE}}
    for n in (1, 2):
        for lam, v in macdonald_P(n, t, q, conv).items():           # P_λ(Y; t, q) の冪和表示
            PY[lam] = v
    T_ = lambda lam: tuple(sum(1 for x in lam if x > j) for j in range(lam[0])) if lam else ()
    Z = {}
    for (lam, mu), val in pair.items():
        sg = conv((-1)**(sum(lam) + sum(mu)))
        for a, x in PY[T_(lam)].items():
            for b, y in PY[T_(mu)].items():
                kadd(Z, (a, b), sg * val * x * y)
    return Z


def log2(Z):
    X = {k: v for k, v in Z.items() if k != ((), ())}
    out, Xn = {}, dict(X)
    for n in range(1, 5):
        for k, v in Xn.items():
            kadd(out, k, conv(Fraction((-1)**(n + 1), n)) * v)
        Xn = mul2(Xn, X)
    return out


def fixed_sector():
    Z = cover_Z()
    Lg = log2(Z)
    x11 = Lg[((1,), (1,))]
    P = (1 - t) * Lg[((1, 1), (1, 1))]
    return P, x11, Lg


def to_rs_nu(e):
    nus = sp.Symbol("nu", positive=True)
    return sp.cancel(e.to_sympy().subs(u_, nus * s_)), nus


def analyze():
    """保存済みの P と凍結した状態 C の比較（theory §12.6）。"""
    import rk_f_frozen_basis as FB
    from rk_k4_crosscap import W, G, in_G
    r, s, nu_ = FB.r, FB.s, FB.nu
    d = pickle.load(open(os.path.join(OUT, "largeN_cover_fixed22_m1.pkl"), "rb"))
    P = sp.sympify(d["P"], locals={"r": r, "s": s, "nu": nu_})
    x11 = sp.sympify(d["x11"], locals={"r": r, "s": s, "nu": nu_})
    zq_ = r - 1 / r
    cases = (("[1] ↔ brane–image (1),(1)", FB.knot_one_box(1), sp.cancel((1 - s**2) * x11 / zq_)),
             ("2 箱 p₁² ↔ brane 2・image 2", FB.knot_two_box(1), sp.cancel(4 * P / zq_**3)))
    for lab, C, Pk in cases:
        bC, bP = FB.nu_blocks(C), FB.nu_blocks(Pk)
        print(f"== {lab}", flush=True)
        for b, blk in sorted(bP.items()):
            Gx = in_G(FB.to_vw(blk / s))
            info = "G 展開できない（σ 非対称）"
            if Gx is not None:
                cs = sp.Poly(sp.expand(Gx * W**400), G, W).coeffs()
                info = f"σ 対称、G の次数 {sp.Poly(Gx, G).degree()}、符号 {FB.sign_of(cs)}"
            sq = ""
            if b % 2 == 0 and b // 2 in bC:
                rat = sp.cancel(blk / bC[b // 2]**2)
                n_, d_ = sp.fraction(rat)
                if sp.Poly(n_, r, s).is_monomial and sp.Poly(d_, r, s).is_monomial:
                    sq = f"、P_{{2β}} = {sp.factor(rat)} × C_β²"
            print(f"  P の ν^{b}：{info}{sq}", flush=True)
    low = FB.nu_blocks(sp.cancel(4 * P / zq_**3))[-24]
    pred = -s**20 * r**-41 * sum(r**(2 * max(i, j)) * s**(2 * min(i, j)) for i in range(12) for j in range(12))
    print("最低のブロック ν^{−24} = −s^20 r^{−41} Σ_{i,j=0}^{11} q^{max(i,j)} t^{min(i,j)}:", sp.expand(sp.cancel(low - pred)) == 0, flush=True)


if __name__ == "__main__" and sys.argv[1:] == ["analyze"]:
    analyze()
elif __name__ == "__main__":
    import rk_f_frozen_basis as FB
    check_decomposition()
    P, x11, Lg = fixed_sector()
    # 検算：(1,1) の成分は brane–image の x = W_adj + κ − W_1²（rk_b_cover_split と同じ）
    import rk_b_cover_split as CS
    print("brane–image (1),(1) の成分が rk_b_cover_split の x と一致:", (x11 - CS.cover_sectors()["x"]).is_zero(), flush=True)
    Ps, nus = to_rs_nu(P)
    n_, d_ = sp.fraction(Ps)
    print("P の分母:", sp.factor(d_), flush=True)
    zq = r_ - 1 / r_
    for k in (0, 1, 2, 3):
        Q = sp.cancel(Ps / zq**k)
        print(f"  P / z_q^{k} は Laurent: {sp.Poly(sp.fraction(Q)[1], r_, s_, nus).is_monomial}", flush=True)
    pickle.dump({"P": str(Ps), "x11": str(to_rs_nu(x11)[0])}, open(os.path.join(OUT, "largeN_cover_fixed22_m1.pkl"), "wb"))
    print("P の因数分解:", sp.factor(Ps), flush=True)
