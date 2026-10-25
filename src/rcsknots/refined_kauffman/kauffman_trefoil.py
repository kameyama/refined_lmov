r"""K3'：trefoil [1] の refined Kauffman（AS 1202.4456 §6）と、refined の h/g 分解。

変数：r = q^{1/2}, s = t^{1/2}, U = t^N。AS の太字 𝐪 = s, 𝐭 = −r/s, 𝐚 = U/s（= ν、Mariño の ν² = Q = U²/t）。
- Z_V(T(2,2m+1)) は AS §6 の一般式（unknot で正規化）。m = 1 で (6.15) に一致【計算】。
- t = q で (6.15) は Mariño 0904.1088 の F_{3_1}(z,ν) の鏡像（ν → 1/ν、z → −z）× ν^{−6}（framing 6 単位）【計算】。
- refined の h/g 分解【仮説】：f_1 = −W_SO、½h^{ref}_[1] = 𝒫(𝐚_U,𝐪,𝐭) (ν − ν^{−1})/(𝐪 − 𝐪^{−1})、𝐚_U² = A (t/q)^{1/2}、A = ν²
  （unknot の f_{10} = −1 − (ν−ν^{−1})/(t^{1/2}−t^{−1/2}) から A = ν² が決まる）、𝒫 は DGR の超多項式の鏡像。
  g^{ref}_[1] = W_SO − ½h^{ref} が Laurent 多項式になるための (6.15) への単項式は r^{−5} s^{d} U^{−6}（d 奇数）に限られ、
  d = 11（= T_[1]^{−6} (q/t)^{1/2}、T_[1] = q^{1/2} t^{N−1}）で g が最短（11 項）【計算】。
"""
import sympy as sp

r, s, U = sp.symbols("r s U", positive=True)
bq, bt, ba = s, -r / s, U / s
nu = ba


def Z_V_T2(m, bq=bq, bt=bt, ba=ba):
    """AS §6 の T(2,2m+1) の一般式（[1]、unknot で正規化）。"""
    inner = sum(bq**(-4 * i) * bt**(-2 * i) for i in range(0, m + 1))
    for j in range(1, 2 * m):
        c1 = sum((bq**4 * bt**2)**(-i + 1) for i in range(1, sp.ceiling(sp.Rational(2 * m + 1 - j, 2)) + 1))
        c2 = sum((bq**4 * bt**2)**(-i + 1) for i in range(1, sp.floor(sp.Rational(2 * m + 1 - j, 2)) + 1))
        inner += ba**(-j) * bq**(-j) * bt**(-2 * j) * (c1 - bq**(-2) * c2)
    return ba**2 * bt**2 + ba**(2 * m + 2) * bq**(2 * m) * bt**(4 * m + 2) * (1 - ba**(-2) * bt**(-2)) * inner


Z615 = (bt**6 * bq**2 * ba**4 + bt**4 * bq**(-2) * ba**4 + bt**4 * bq * ba**3 - bt**4 * bq**(-1) * ba**3
        - bt**4 * bq**2 * ba**2 + bt**2 * ba**2 - bt**2 * bq**(-2) * ba**2 - bt**2 * bq * ba + bt**2 * bq**(-1) * ba)
unknot_SO = (ba - 1 / ba) / (bq - 1 / bq) + 1


def marino_F31(nu_, z_):
    """Mariño 0904.1088 の Kauffman 多項式 F_{3_1}(z, ν)。"""
    return 2 * nu_**2 - nu_**4 + z_ * (-nu_**3 + nu_**5) + z_**2 * (nu_**2 - nu_**4)


def P_DGR(A2, Q, T):
    """DGR の超多項式 𝒫(T_{2,3}) = 𝐚²𝐪⁻² + 𝐚²𝐪²𝐭² + 𝐚⁴𝐭³（𝐚² を A2 として渡す）。"""
    return A2 * Q**-2 + A2 * Q**2 * T**2 + A2**2 * T**3


def h_half_ref(A=nu**2, mirror=True):
    """½h^{ref}_[1](3_1) = 𝒫(𝐚_U, 𝐪, 𝐭) (√A − 1/√A)/(𝐪 − 1/𝐪)、𝐚_U² = A (t/q)^{1/2}。"""
    A2U = A * s / r
    P = P_DGR(1 / A2U, 1 / bq, 1 / bt) if mirror else P_DGR(A2U, bq, bt)
    return P * (sp.sqrt(A) - 1 / sp.sqrt(A)) / (bq - 1 / bq)


def W_SO(mono=r**-5 * s**11 * U**-6):
    """標準 framing の unnormalized refined SO trefoil [1]（仮説の単項式 mono を掛ける）。"""
    return Z615 * mono * unknot_SO


def g_ref(mono=r**-5 * s**11 * U**-6):
    return sp.cancel(W_SO(mono) - h_half_ref())


def is_laurent(expr):
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    return sp.Poly(den, r, s, U).is_monomial
