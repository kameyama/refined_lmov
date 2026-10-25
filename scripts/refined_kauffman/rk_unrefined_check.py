r"""refined Kauffman（第 3 の研究）K0：unrefined の Kauffman large N duality（Mariño 0904.1088）の再現。

Mariño の変数 (ν, t)、z = t − t^{−1}。
- trefoil：\hat g_[2], \hat g_[1,1] を Mariño (4.8) と全係数で比較。
- 4_1：\hat h_[2], \hat g_[2], \hat g_[1,1] を NRZ 1310.2240 付録 C と比較。
入力：colored HOMFLY（trefoil は src/rcsknots/homfly.py の Rosso–Jones、4_1 は NRZ/FGS のツイスト結び目の式）、
colored Kauffman（NRZ (6.11), (6.26) の対称表現の式の t = −1）、composite ([1],[1]) の HOMFLY（Mariño (4.4)、Hadji–Morton）。
[1,1] の Kauffman は F_{[1,1]}(ν, t) = F_{[2]}(−ν, t^{−1}) の対称性で作る（unknot で確認、trefoil の (4.8) で検証）。
"""
import sys
import sympy as sp

sys.path.insert(0, "src")
from rcsknots import homfly

nu, t = sp.symbols("nu t")
z = t - 1/t


def qpoch(x, qq, k):
    return sp.prod([1 - x*qq**j for j in range(k)]) if k > 0 else sp.Integer(1)


def qbinom(n, k, qq):
    if k < 0 or k > n:
        return sp.Integer(0)
    num = sp.prod([1 - qq**(n-j) for j in range(k)]) if k > 0 else sp.Integer(1)
    den = sp.prod([1 - qq**(j+1) for j in range(k)]) if k > 0 else sp.Integer(1)
    return sp.cancel(num/den)


# ---- unknot（Mariño の規約）----
p = {n: (nu**n - nu**(-n))/(t**n - t**(-n)) for n in (1, 2)}
s1 = p[1]; s2 = (p[1]**2 + p[2])/2; s11 = (p[1]**2 - p[2])/2
H_unk = {(1,): s1, (2,): s2, (1, 1): s11}
G_unk = {(1,): s1 + 1, (2,): s2 + s1, (1, 1): s11 + s1}
H_adj_unk = p[1]**2 - 1


# ---- colored Kauffman（NRZ、対称表現、t = −1）----
def nrz_trefoil(r, lam, qq, tt):
    """NRZ (triple-Kauffman-trefoil)：正規化された [r] 彩色の trefoil の Kauffman 超多項式（t = −1 で Kauffman）。"""
    tot = 0
    for k in range(r+1):
        for j in range(k+1):
            for i in range(r-k+1):
                tot += (lam**(i-k+3*r) * qq**(3*k - 2*j*(1+r) + r*(2*r-3) + i*(2*j+2*r-1)) * tt**(2*(i-j+r))
                        * qbinom(r, k, qq**2) * qbinom(k, j, qq**2) * qbinom(r-k, i, qq**2)
                        * qpoch(-lam**2*tt**3*qq**(2*r-2), qq**2, j) * qpoch(-qq**2*tt, qq**2, r-k) * qpoch(-lam*tt/qq, qq**2, i))
    return tot


def nrz_fig8(r, lam, qq, tt):
    """NRZ (triple-Kauffman-fig8)。"""
    tot = 0
    for k in range(r+1):
        for j in range(k+1):
            for i in range(r-k+1):
                tot += (lam**(2*(k-r-j)+i) * qq**(4*j - 2*k**2 - 2*k*(1+i-j-3*r) + (i-2*r)*(2*r+2*j-1)) * tt**(2*(2*k-2*r+i-j))
                        * qbinom(r, k, qq**2) * qbinom(r-k, i, qq**2) * qbinom(k, j, qq**2)
                        * qpoch(-lam**2*qq**(2*(r-k-i)-4)*tt, qq**2, j) * qpoch(-lam**2*qq**(-2+4*r-2*k)*tt**3, qq**2, j)
                        * qpoch(-qq**2*tt, qq**2, r-k) * qpoch(-lam**2*qq**(-2+2*r)*tt**3, qq**2, r-k)
                        * qpoch(-lam*tt/qq, qq**2, r-k-i) * qpoch(-lam*qq**(1+2*(r-k+j))*tt**3, qq**2, r-k-i))
    return tot


# ---- colored HOMFLY of twist knots（FGS/NRZ の式、ノートブック supertwistneg）----
def supertwist(n, pp, aa, qq, tt, sign):
    """sign=+1: supertwistpos（(−t)^{1−n} 付き）、sign=−1: supertwistneg。n = r + 1。"""
    tot = 0
    for k in range(n):
        inner = 0
        for l in range(k+1):
            inner += ((-1)**l * aa**(pp*l) * tt**(2*pp*l) * qq**(sp.Rational(2*pp+1, 2)*l*(l-1)) * (1 - aa*tt**2*qq**(2*l-1))
                      / qpoch(aa*tt**2*qq**(l-1), qq, k+1) * qbinom(k, l, qq))
        tot += (qq**k / qpoch(qq, qq, k) * qpoch(-aa*tt/qq, qq, k) * qpoch(qq**(1-n), qq, k) * qpoch(-aa*tt**3*qq**(n-1), qq, k)) * inner
    if sign > 0:
        tot *= (-tt)**(1-n)
    return tot


def simp(e):
    return sp.factor(sp.cancel(sp.together(e)))


def hat2(g2, g11):
    return simp((t*g2 + g11/t)/(t**2 - t**(-2))), simp((g2/t + t*g11)/(t**2 - t**(-2)))


def reformulated(H, Hadj, G):
    """2 箱までの g_R, h_R（Mariño (3.14)、NRZ (C.3)）。H, G は unreduced。"""
    g1 = G[(1,)] - H[(1,)]
    g2 = G[(2,)] - G[(1,)]**2/2 - H[(2,)] + H[(1,)]**2 - Hadj/2
    g11 = G[(1, 1)] - G[(1,)]**2/2 - H[(1, 1)] + H[(1,)]**2 - Hadj/2
    H1sq = H[(1,)].subs({nu: nu**2, t: t**2})
    h2 = 2*H[(2,)] + Hadj - 2*H[(1,)]**2 - H1sq
    h11 = 2*H[(1, 1)] + Hadj - 2*H[(1,)]**2 + H1sq
    return simp(g1), simp(g2), simp(g11), simp(h2), simp(h11)


def in_ring(expr, kind):
    r"""kind='g': Z[z, ν^{±1}]（\hat g）、'h': z^{-1} Z[z², ν^{±1}]（\hat h）。z = t − 1/t での展開で判定。"""
    zz = sp.symbols("zz")
    e = sp.cancel(sp.together(expr))
    num, den = sp.fraction(e)
    if kind == "h":
        e = sp.cancel(e * z)
    # t を z で表す：t を根とする二次式 t² − z t − 1 = 0。e が z の多項式かどうかは、t → 1/t で不変（z 不変）かつ…
    # 単純に：e(t) が t ↔ −1/t で不変（z ↔ z）なら e は z の有理式。さらに分母が t の冪だけなら Laurent。
    inv = sp.simplify(e - e.subs(t, -1/t))
    num, den = sp.fraction(sp.cancel(sp.together(e)))
    den_ok = den.is_polynomial(nu, t) and all(sp.degree(den, s) >= 0 for s in (nu, t)) and sp.factor_list(den)[1] == [] or all(f[0] in (nu, t) for f in sp.factor_list(den)[1])
    return inv == 0 and den_ok


def z_expansion(expr):
    r"""\hat g を Σ c_{k,Q} z^k ν^Q に展開（辞書）。t = (z + sqrt(z²+4))/2 の代わりに、t の Laurent 多項式を z の多項式に書き直す。"""
    zz = sp.symbols("zz")
    e = sp.expand(sp.cancel(sp.together(expr)))
    # e は t の Laurent 多項式で t ↔ −1/t 不変 → 各 ν 冪の係数は z の多項式。ν の係数ごとに、t^k + (−1)^k t^{−k} を z で表す。
    out = {}
    e = sp.expand(e)
    for Qpow in range(-30, 31):
        c = e.coeff(nu, Qpow)
        if c == 0:
            continue
        c = sp.expand(c)
        poly = sp.Integer(0)
        # Chebyshev 型：t^k + (−1)^k t^{−k} = P_k(z)
        P = {0: sp.Integer(2), 1: zz}
        for k in range(2, 40):
            P[k] = sp.expand(zz*P[k-1] + P[k-2])  # t^k+(-1)^k t^-k = z (t^{k-1}+(-1)^{k-1}t^{-(k-1)}) + (t^{k-2}+...)
        rem = c
        for k in range(39, -1, -1):
            ck = rem.coeff(t, k)
            if ck != 0:
                if k == 0:
                    poly += ck
                    rem = sp.expand(rem - ck)
                else:
                    poly += ck*P[k]
                    rem = sp.expand(rem - ck*(t**k + (-1)**k*t**(-k)))
        assert rem == 0, (Qpow, rem)
        for k, ck in sp.Poly(poly, zz).as_dict().items():
            out[(k[0], Qpow)] = int(ck)
    return out


if __name__ == "__main__":
    # ===== trefoil =====
    # HOMFLY（package の Rosso–Jones、KN 変数 a = ν², Q = t）。符号は H_[1] を Mariño の P で固定。
    P31_marino = 2*nu**2 - nu**4 - z**2*nu**2          # Mariño の本文の値
    P31_std = 2*nu**2 - nu**4 + z**2*nu**2             # 量子群の標準（a = ν², q^{1/2} = t）
    sign, sub = -1, {homfly.a: nu**2, homfly.Q: t}
    h1 = sp.simplify(homfly.torus_reduced((1,), 2, 3, sign).subs(sub))
    print("package H_[1](3_1) =", sp.expand(h1), "| == std:", sp.simplify(h1 - P31_std) == 0, "| == Mariño 本文:", sp.simplify(h1 - P31_marino) == 0)
    H = {lam: simp(homfly.torus_reduced(lam, 2, 3, sign).subs(sub) * H_unk[lam]) for lam in ((1,), (2,), (1, 1))}
    Hadj = simp(H_adj_unk * (4*nu**4 - 4*nu**6 + nu**8 + z**2*(4*nu**4 - 7*nu**6 + 2*nu**8 + nu**10) + z**4*(nu**4 - 2*nu**6 + nu**8)))
    F31 = 2*nu**2 - nu**4 + z*(-nu**3 + nu**5) + z**2*(nu**2 - nu**4)
    print("NRZ r=1 (t=-1) == Mariño F(3_1):", sp.simplify(nrz_trefoil(1, nu, t, sp.Integer(-1)) - F31) == 0, "|", sp.expand(nrz_trefoil(1, nu, t, sp.Integer(-1))))
    G = {(1,): simp(F31 * G_unk[(1,)]),
         (2,): simp(nrz_trefoil(2, nu, t, sp.Integer(-1)) * G_unk[(2,)]),
         (1, 1): simp(nrz_trefoil(2, -nu, 1/t, sp.Integer(-1)) * G_unk[(1, 1)])}
    g1, g2, g11, h2, h11 = reformulated(H, Hadj, G)
    gh2, gh11 = hat2(g2, g11)
    hh2, hh11 = hat2(h2, h11)
    # Mariño (4.8)
    M2 = (-21*nu**11 + 79*nu**9 - 111*nu**7 + 69*nu**5 - 16*nu**3 + 21*z*(nu**12 - 3*nu**10 + 3*nu**8 - nu**6)
          + z**2*(-70*nu**11 + 251*nu**9 - 307*nu**7 + 146*nu**5 - 20*nu**3) + 7*z**3*(10*nu**12 - 33*nu**10 + 33*nu**8 - 10*nu**6)
          - 2*z**4*(42*nu**11 - 165*nu**9 + 183*nu**7 - 64*nu**5 + 4*nu**3) + 14*z**5*(6*nu**12 - 23*nu**10 + 23*nu**8 - 6*nu**6)
          + z**6*(-45*nu**11 + 220*nu**9 - 230*nu**7 + 56*nu**5 - nu**3) + 3*z**7*(15*nu**12 - 73*nu**10 + 73*nu**8 - 15*nu**6)
          + z**8*(-11*nu**11 + 78*nu**9 - 79*nu**7 + 12*nu**5) + z**9*(-11*nu**11 + 78*nu**9 - 79*nu**7 + 12*nu**5)
          + z**10*(-nu**11 + 14*nu**9 - 14*nu**7 + nu**5) + z**11*(nu**12 - 14*nu**10 + 14*nu**8 - nu**6)
          + z**12*(nu**9 - nu**7) + z**13*(nu**8 - nu**10))
    M11 = (-15*nu**11 + 53*nu**9 - 69*nu**7 + 39*nu**5 - 8*nu**3 + 15*z*(nu**12 - 3*nu**10 + 3*nu**8 - nu**6)
           + z**2*(-35*nu**11 + 126*nu**9 - 146*nu**7 + 61*nu**5 - 6*nu**3) + 5*z**3*(7*nu**12 - 24*nu**10 + 24*nu**8 - 7*nu**6)
           + z**4*(-28*nu**11 + 120*nu**9 - 128*nu**7 + 37*nu**5 - nu**3) + 7*z**5*(4*nu**12 - 17*nu**10 + 17*nu**8 - 4*nu**6)
           + z**6*(-9*nu**11 + 55*nu**9 - 56*nu**7 + 10*nu**5) + z**7*(9*nu**12 - 55*nu**10 + 55*nu**8 - 9*nu**6)
           + z**8*(-nu**11 + 12*nu**9 - 12*nu**7 + nu**5) + z**9*(nu**12 - 12*nu**10 + 12*nu**8 - nu**6)
           + z**10*(nu**9 - nu**7) + z**11*(nu**8 - nu**10))
    print("trefoil  g_[1] =", g1)
    print("trefoil  \\hat g_[2]  == Mariño (4.8):", sp.simplify(gh2 - M2) == 0)
    d = z_expansion(sp.expand(gh2 - M2)); print("   差の z 展開 {(z冪, ν冪): 係数}:", d)
    print("   計算した \\hat g_[2] の z^8, z^9 の係数:", {k: v for k, v in z_expansion(gh2).items() if k[0] in (8, 9)})
    print("trefoil  \\hat g_[1,1]== Mariño (4.8):", sp.simplify(gh11 - M11) == 0)
    print("trefoil  \\hat h_[2]  =", hh2)
    print("trefoil  \\hat h_[1,1]=", hh11)
    print("trefoil  z^-1 \\hat h_[2] in Z[z^2,nu]:", in_ring(hh2, "h"), " \\hat g_[2] in Z[z,nu]:", in_ring(gh2, "g"))

    # ===== 4_1 =====
    # HOMFLY：supertwistneg（a → ν², q → t², t → −1、n = r+1、p = −1）× unknot
    H41 = {}
    for lam, r in (((1,), 1), ((2,), 2)):
        H41[lam] = simp(supertwist(r+1, -1, nu**2, t**2, sp.Integer(-1), -1) * H_unk[lam])
    H41[(1, 1)] = simp(supertwist(3, -1, nu**2, t**-2, sp.Integer(-1), -1) * H_unk[(1, 1)])
    adj41 = (3 - 2*z**2 - 6*z**4 - 2*z**6 + (-2 - z**2 + 2*z**4 + z**6)*(nu**-2 + nu**2) + (1 + 2*z**2 + z**4)*(nu**-4 + nu**4))
    Hadj41 = simp(adj41 * H_adj_unk)
    G41 = {(1,): simp(nrz_fig8(1, nu, t, sp.Integer(-1)) * G_unk[(1,)]),
           (2,): simp(nrz_fig8(2, nu, t, sp.Integer(-1)) * G_unk[(2,)]),
           (1, 1): simp(nrz_fig8(2, -nu, 1/t, sp.Integer(-1)) * G_unk[(1, 1)])}
    g1, g2, g11, h2, h11 = reformulated(H41, Hadj41, G41)
    gh2, gh11 = hat2(g2, g11)
    hh2, hh11 = hat2(h2, h11)
    # NRZ (C.10)（λ = ν、q = t）
    lam, qq = nu, t
    N2 = -1/(qq**11*lam**6)*(qq-lam)*(-1+lam)*(1+lam)*(qq+lam)*(-1+qq*lam)*(1+qq*lam)*(
        qq**2+qq**6-qq**12-qq**16-qq**3*lam-qq**7*lam-qq**9*lam-qq**11*lam-qq**15*lam-qq**6*lam**2+qq**8*lam**2
        -qq**10*lam**2+qq**12*lam**2+qq**5*lam**3+qq**7*lam**3+qq**11*lam**3+qq**13*lam**3+qq**4*lam**4-qq**8*lam**4+qq**10*lam**4
        -qq**14*lam**4-qq*lam**5-qq**5*lam**5-qq**7*lam**5-qq**9*lam**5-qq**11*lam**5-qq**13*lam**5-qq**17*lam**5-lam**6
        -qq**4*lam**6-qq**6*lam**6+qq**12*lam**6+qq**14*lam**6+qq**18*lam**6)
    N11 = 1/(qq**11*lam**6)*(qq-lam)*(-1+lam)*(1+lam)*(qq+lam)*(-1+qq*lam)*(1+qq*lam)*(
        -1-qq**4-qq**6+qq**12+qq**14+qq**18+qq*lam+qq**5*lam+qq**7*lam+qq**9*lam+qq**11*lam+qq**13*lam
        +qq**17*lam+qq**4*lam**2-qq**8*lam**2+qq**10*lam**2-qq**14*lam**2-qq**5*lam**3-qq**7*lam**3-qq**11*lam**3-qq**13*lam**3-qq**6*lam**4
        +qq**8*lam**4-qq**10*lam**4+qq**12*lam**4+qq**3*lam**5+qq**7*lam**5+qq**9*lam**5+qq**11*lam**5+qq**15*lam**5+qq**2*lam**6
        +qq**6*lam**6-qq**12*lam**6-qq**16*lam**6)
    NH2 = 1/(lam**6*qq**7*(qq-1/qq))*(-lam**4+2*lam**6-lam**8-2*lam**6*qq**2+4*lam**8*qq**2-2*lam**10*qq**2+qq**4-4*lam**2*qq**4
        +4*lam**4*qq**4+2*lam**8*qq**4-6*lam**10*qq**4+3*lam**12*qq**4+2*qq**6-2*lam**2*qq**6-5*lam**6*qq**6+6*lam**8*qq**6
        -2*lam**10*qq**6+qq**8-4*lam**2*qq**8+4*lam**4*qq**8+2*lam**8*qq**8-6*lam**10*qq**8+3*lam**12*qq**8-2*lam**6*qq**10
        +4*lam**8*qq**10-2*lam**10*qq**10-lam**4*qq**12+2*lam**6*qq**12-lam**8*qq**12)
    print("4_1  g_[1] =", g1)
    print("4_1  \\hat g_[2]  == NRZ (C.10):", sp.simplify(gh2 - N2) == 0)
    print("4_1  \\hat g_[1,1]== NRZ (C.10):", sp.simplify(gh11 - N11) == 0)
    print("4_1  \\hat h_[2]  == NRZ (C.9) :", sp.simplify(hh2 - NH2) == 0)
    exp2 = z_expansion(gh2)
    neg = {k: v for k, v in exp2.items()}
    print("4_1  \\hat g_[2] の z 展開（(z 冪, ν 冪): 係数）の一部:", sorted(exp2.items())[:8])
