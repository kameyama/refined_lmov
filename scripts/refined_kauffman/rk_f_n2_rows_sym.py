"""n = 2（brane 2・image 2、trefoil）の N_λ を、[1]・2 箱・随伴・composite の行を記号にして ν^{−1} の級数で導出する（続き 110、theory §12.40）。
rk_f_cover_fixed22 のパイプライン（⟨P_αP_β^∨⟩ の composite の分解、log Z のキュムラント展開、(1−t)/(z_q φ φ)、Schur 分解、
N_λ = M̂_λλ/(s(t − ν²)(ν² − 1)(ν²q − 1))）を、入力の R を記号の行の多項式に置き換えて sympy で再現する。
  R₁ = a₂ν^{−2} + a₄ν^{−4}（w = s/r を含む鏡像の行）、R_(2) = b₄ν^{−4} + b₆ν^{−6} + b₈ν^{−8}、R_(1,1) = c₄ν^{−4} + …、
  R_adj = R₁² + (q−1)(t−1)G_adj、R_{[α;1]} = R_αR₁ + (q−1)(t−1)G_{α1}、R_{[α;β]} = R_αR_β + (q−1)(t−1)G_{αβ}（§12.36 の因子化）。
級数の係数は行の記号の多項式（係数体 Q(q,t)、r = q^{1/2}, s = t^{1/2} を使うので Q(r,s)）。
実行：PYTHONPATH=src:scripts/refined_kauffman python scripts/refined_kauffman/rk_f_n2_rows_sym.py [--kmax=20]"""
import sys, json, pickle, time
from fractions import Fraction
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rcsknots.macdonald import macdonald_P, principal_specialization
from rcsknots.refined_kauffman import kauffman_refined as KR

r, s = sp.symbols("r s", positive=True)
q, t = r**2, s**2
nu = sp.Symbol("nu", positive=True)
KMAX = 20
LAMS = ("2", "11")
P2 = ((2,), (1, 1))
# 行の記号
a2, a4 = sp.symbols("a2 a4")
b4, b6, b8 = sp.symbols("b4 b6 b8")
c4, c6, c8 = sp.symbols("c4 c6 c8")
gadj = {j: sp.Symbol(f"gadj{j}") for j in (6, 8, 10)}
g21 = {j: sp.Symbol(f"g21_{j}") for j in (8, 10, 12, 14)}
g111 = {j: sp.Symbol(f"g111_{j}") for j in (8, 10, 12, 14)}
g22 = {j: sp.Symbol(f"g22_{j}") for j in (10, 12, 14, 16, 18, 20)}
g211 = {j: sp.Symbol(f"g211_{j}") for j in (10, 12, 14, 16, 18)}
g1111 = {j: sp.Symbol(f"g1111_{j}") for j in (10, 12, 14, 16, 18, 20)}
SYMS = [a2, a4, b4, b6, b8, c4, c6, c8] + list(gadj.values()) + list(g21.values()) + list(g111.values()) + list(g22.values()) + list(g211.values()) + list(g1111.values())
DOM = sp.QQ.frac_field(r, s)


def P(e):
    return sp.Poly(sp.sympify(e), *SYMS, domain=DOM)


ZERO = P(0)


class Ser:
    """Σ_k c_k ν^{−k}（k ≤ KMAX、負の k も可）。c_k は Poly。"""
    __slots__ = ("c",)

    def __init__(self, c=None):
        self.c = {k: v for k, v in (c or {}).items() if k <= KMAX and not v.is_zero}

    @staticmethod
    def const(e):
        return Ser({0: P(e)})

    def __add__(self, o):
        o = o if isinstance(o, Ser) else Ser.const(o)
        c = dict(self.c)
        for k, v in o.c.items():
            c[k] = c.get(k, ZERO) + v
        return Ser(c)

    __radd__ = __add__

    def __neg__(self):
        return Ser({k: -v for k, v in self.c.items()})

    def __sub__(self, o):
        return self + (-(o if isinstance(o, Ser) else Ser.const(o)))

    def __mul__(self, o):
        if not isinstance(o, Ser):
            p = P(o)
            return Ser({k: v * p for k, v in self.c.items()})
        c = {}
        for k1, v1 in self.c.items():
            for k2, v2 in o.c.items():
                if k1 + k2 <= KMAX:
                    c[k1 + k2] = c.get(k1 + k2, ZERO) + v1 * v2
        return Ser(c)

    __rmul__ = __mul__

    def adams2(self):
        """r → r², s → s²、ν → ν²、記号 → 記号（ψ² は行の記号には作用させない：log の Adams 項は結び目の行を含む p₁p₁′ の 2 乗なので別扱い）。"""
        raise NotImplementedError

    def is_zero(self):
        return not self.c


def rat_series(expr, kmax=None):
    """ν の有理式（係数 r, s）→ Ser。ν = 1/w として w の級数。"""
    kmax = KMAX if kmax is None else kmax
    w = sp.Symbol("w", positive=True)
    e = sp.cancel(sp.together(sp.sympify(expr).subs(nu, 1 / w)))
    num, den = sp.fraction(e)
    Pn, Pd = sp.Poly(sp.expand(num), w), sp.Poly(sp.expand(den), w)
    nd = {m[0]: c for m, c in zip(Pn.monoms(), Pn.coeffs())}
    dd = {m[0]: c for m, c in zip(Pd.monoms(), Pd.coeffs())}
    d0 = min(dd); dd = {k - d0: v for k, v in dd.items()}
    c = {}
    for k in range(kmax + d0 + 1):
        v = nd.get(k, 0) - sum(dd[j] * c[k - j] for j in dd if j >= 1 and k - j in c)
        c[k] = sp.cancel(v / dd[0])
    return Ser({k - d0: P(v) for k, v in c.items() if k - d0 <= kmax})


# ---- 核（有理式） ----
A = nu**2
w_fr = s / r


def Psym(lam):
    k = sum(lam)
    return s**k * nu**(-k) * principal_specialization(lam, A, q, t)


def psum(d, prime, inverse):
    if not inverse:
        base = (1 - A**(-d)) / (1 - t**(-d))
        return base + (q**d - 1 if prime else 0)
    base = (1 - A**d) / (1 - t**d)
    return base + (q**(-d) - 1 if prime else 0)


_MP = {}


def P_pt(lam, prime, inverse):
    n = sum(lam)
    if n not in _MP:
        _MP[n] = macdonald_P(n, q, t, lambda c: sp.Rational(c.numerator, c.denominator))
    tot = 0
    for rho, c in _MP[n][tuple(lam)].items():
        v = 1
        for d in rho:
            v = v * psum(d, prime, inverse)
        tot = tot + c * v
    return tot


def e_val(alpha, beta):
    eN = (A - 1) / (t - 1)
    e = eN
    for i, a_ in enumerate(alpha):
        e = e + (q**a_ - 1) * A * t**(-1 - i)
    for k, b_ in enumerate(beta):
        e = e + (q**(-b_) - 1) * t**k
    return e


def dimc(alpha, beta):
    return KR.dim_composite(list(alpha), list(beta)).subs({KR.q: q, KR.t: t, KR.A: A})


def decomposition_A(alpha, beta):
    Kb = dimc(alpha, beta)
    qn = (nu - 1 / nu) / (s - 1 / s)
    kappa = (1 - A) * (1 - q) / ((1 - t) * (1 - q * A / t))
    Kadj = qn * qn - kappa
    eN = (A - 1) / (t - 1)
    f = P_pt(alpha, False, False) * P_pt(beta, False, True)
    fp = P_pt(alpha, True, False) * P_pt(beta, True, True)
    eb, eadj = e_val(alpha, beta), e_val((1,), (1,))
    r1 = f - Kb
    r2 = eN * fp - eb * Kb
    det = Kadj * eN - eadj * Kadj
    c1 = sp.cancel((r1 * eN - r2) / det)
    c0 = sp.cancel(r1 - c1 * Kadj)
    return c1, c0, Kb, Kadj


def decomposition_flint(alpha, beta):
    """rk_f_cover_fixed22.decomposition_A（flint RF）→ ν の有理式（sympy）。"""
    import rk_f_cover_fixed22 as CF
    u_ = sp.Symbol("u", positive=True)
    c1, c0, Kb, _ = CF.decomposition_A(list(alpha), list(beta))
    cv = lambda x: x.to_sympy().subs(u_, nu * s)
    return cv(c1), cv(c0), cv(Kb)


def build(kmax=KMAX, verbose=True, val=None):
    """val = {記号: 値} を与えると記号の代わりに値（r, s の式）で計算する（検証用）。"""
    global KMAX
    KMAX = kmax
    t0 = time.time()
    V = (lambda x: val.get(x, x)) if val else (lambda x: x)
    qn = (nu - 1 / nu) / (s - 1 / s)
    kappa_e = (1 - A) * (1 - q) / ((1 - t) * (1 - q * A / t))
    # 記号の行
    R = {(1,): Ser({2: P(V(a2)), 4: P(V(a4))}), (2,): Ser({4: P(V(b4)), 6: P(V(b6)), 8: P(V(b8))}), (1, 1): Ser({4: P(V(c4)), 6: P(V(c6)), 8: P(V(c8))})}
    qt1 = (q - 1) * (t - 1)
    Gadj = Ser({j: P(V(v)) for j, v in gadj.items()})
    Radj = R[(1,)] * R[(1,)] + Gadj * qt1
    Gc = {(2,): Ser({j: P(V(v)) for j, v in g21.items()}), (1, 1): Ser({j: P(V(v)) for j, v in g111.items()})}
    Gb = {((2,), (2,)): Ser({j: P(V(v)) for j, v in g22.items()}), ((2,), (1, 1)): Ser({j: P(V(v)) for j, v in g211.items()}), ((1, 1), (1, 1)): Ser({j: P(V(v)) for j, v in g1111.items()})}
    # 核の級数
    K = {lam: rat_series(Psym(lam)) for lam in R}
    kappa = rat_series(kappa_e)
    Kadj = rat_series(qn * qn - kappa_e)
    Wo = {lam: R[lam] * K[lam] for lam in R}
    Wadj = Radj * Kadj
    pair = {((), ()): Ser.const(1), ((1,), (1,)): Wadj + kappa}
    for lam, v in Wo.items():
        pair[(lam, ())] = pair[((), lam)] = v
    for alpha in P2:
        dc = dimc(alpha, (1,))
        Sdc = rat_series(dc)
        d = (K[alpha] * K[(1,)] - Sdc) * rat_series(1 / Psym((1,)))      # d_α = (P_α(t^ρ)P_1(t^ρ) − dim[α;1])/P_1(t^ρ) を級数で
        Rc = R[alpha] * R[(1,)] + Gc[alpha] * qt1
        pair[(alpha, (1,))] = pair[((1,), alpha)] = Rc * Sdc + d * Wo[(1,)]
    for alpha in P2:
        for beta in P2:
            key = (alpha, beta) if (alpha, beta) != ((1, 1), (2,)) else ((2,), (1, 1))
            c1, c0, Kb = decomposition_flint(*key)              # flint で約分してから級数化（sympy の約分は遅い）
            Rb = R[key[0]] * R[key[1]] + Gb[key] * qt1
            pair[(alpha, beta)] = Rb * rat_series(Kb) + rat_series(c1) * Wadj + rat_series(c0)
    if verbose:
        print(f"  核と対の値：{time.time() - t0:.0f} 秒", flush=True)
    # Z(Y, Y′)
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
    # log
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
    # Adams の項：½ x11.adams(2)。x11 は記号の行を含むので、ψ² を記号にも作用させる必要がある → 新しい記号 ψ²(a) を導入せず、
    # 実際の値の検証では x11 の adams を数値で引く。ここでは記号のまま「adams 項」を別に返す。
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
    # Adams の項の係数（λ 成分、符号付き）：−½ (1−t)/(z_q φ_(2) φ_(2)) × χχ × ψ²(x11) → 呼び出し側で足す
    adams_coef = {lab: sp.Rational(-1, 2) * (1 - t) * CHI[lab][(2,)]**2 / (zq * phi[(2,)]**2) for lab in LAMS}
    if verbose:
        print(f"  N_λ（Adams の項を除く）：{time.time() - t0:.0f} 秒", flush=True)
    return N, x11, adams_coef, inv_den


def actual_rows():
    """被覆データの行（鏡像、w^{|λ|} を含む）と、§12.36 の G の行（Radj_A、Rcomp_A、Rcomp22_A から）。"""
    import rk_f_comp_factor as CFp
    import pickle, os
    OUT = "data/refined_kauffman/daha"
    w = s / r
    Rd = {c: sp.expand(w**k * CFp.R_data(c, "T2_3")) for c, k in (("1", 1), ("2", 2), ("11", 2))}
    rows = {c: CFp.rows(v) for c, v in Rd.items()}
    val = {a2: rows["1"][-2], a4: rows["1"][-4], b4: rows["2"][-4], b6: rows["2"][-6], b8: rows["2"][-8],
           c4: rows["11"][-4], c6: rows["11"][-6], c8: rows["11"][-8]}
    radj = pickle.load(open(f"{OUT}/Radj_A_m1.pkl", "rb"))
    Radj = sp.expand(sum(sp.sympify(v, locals={"r": r, "s": s}) * CFp.nu**int(j) for j, v in radj.items()))
    def Gof(Rab, Ra, Rb, syms):
        G, _, _ = CFp.G_rows(Rab, Ra, Rb)
        for j, sym in syms.items():
            val[sym] = sp.expand(G.get(-j, 0))
        assert all(k in [-j for j in syms] for k in G if G[k] != 0), (syms, sorted(G))
    Gof(Radj, Rd["1"], Rd["1"], gadj)
    Gof(CFp.R_comp(f"{OUT}/Rcomp_A_m1_2.pkl", 1), Rd["2"], Rd["1"], g21)
    Gof(CFp.R_comp(f"{OUT}/Rcomp_A_m1_11.pkl", 1), Rd["11"], Rd["1"], g111)
    Gof(CFp.R_comp(f"{OUT}/Rcomp22_A_m1_2_2.pkl", 2), Rd["2"], Rd["2"], g22)
    Gof(CFp.R_comp(f"{OUT}/Rcomp22_A_m1_2_11.pkl", 2), Rd["2"], Rd["11"], g211)
    Gof(CFp.R_comp(f"{OUT}/Rcomp22_A_m1_11_11.pkl", 2), Rd["11"], Rd["11"], g1111)
    return val


def substitute(ser, val):
    return {k: sp.expand(v.as_expr().subs(val)) for k, v in ser.c.items()}


def verify(kmax=14):
    """実際の行を代入した build と、保存済みの N_λ（rk_f_n2_rterms の N）をブロックごとに比べる。"""
    import rk_f_n2_rterms as NR
    val = actual_rows()
    N, x11, adams_coef, inv_den = build(kmax, val=val)
    D = NR.load(0)
    for lab in LAMS:
        ser = NR.series_u(D["N"][(NR.GROUPS, lab)], kmax)
        x11v = {k: v.as_expr() for k, v in x11.c.items()}
        ad = Ser({2 * k: P(v.subs({r: r**2, s: s**2}, simultaneous=True)) for k, v in x11v.items()}) * adams_coef[lab] * inv_den
        tot = N[lab] + ad
        ok = {}
        for k in range(0, kmax + 1):
            mine = tot.c.get(k, ZERO).as_expr()
            ref = NR.to_sym(NR.nu_block(ser[k], k)) if k in ser and not ser[k].is_zero() else 0
            ok[k] = sp.cancel(mine - ref) == 0
        print(f"λ = {lab}：ν^-k のブロックの一致 {ok}", flush=True)


def structure(kmax=14, kshow=(8, 10, 12)):
    """記号版のブロックの単項式構造。"""
    N, x11, adams_coef, inv_den = build(kmax, verbose=False)
    for lab in LAMS:
        print(f"== λ = {lab}（Adams の項を除く）", flush=True)
        for k in kshow:
            p = N[lab].c.get(k, ZERO)
            terms = p.terms()
            print(f"  ν^-{k}：{len(terms)} 単項式", flush=True)
            for mon, c in terms[:40]:
                mstr = "·".join(f"{str(g)}^{e}" if e > 1 else str(g) for g, e in zip(SYMS, mon) if e) or "1"
                print(f"     [{mstr}] {sp.factor(sp.sympify(c).as_expr() if hasattr(c, 'as_expr') else c)}"[:200], flush=True)


if __name__ == "__main__":
    kmax = next((int(a_[7:]) for a_ in sys.argv[1:] if a_.startswith("--kmax=")), 14)
    if "--structure" in sys.argv:
        structure(kmax)
    else:
        verify(kmax)
