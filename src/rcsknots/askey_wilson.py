r"""Askey–Wilson 多項式（1 変数 Koornwinder）を三項漸化式で厳密に計算する。

Koekoek–Swarttouw (3.1.4)：x = (y + 1/y)/2、基底 Q、パラメータ (a,b,c,d)：
  2x p_n(x) = A_n p_{n+1} + B_n p_n + C_n p_{n−1}、p_0 = 1、
  A_n = (1 − abcd Q^{n−1})(1 − ab Q^n)(1 − ac Q^n)(1 − ad Q^n) / (a (1 − abcd Q^{2n−1})(1 − abcd Q^{2n})),
  C_n = a (1 − Q^n)(1 − bc Q^{n−1})(1 − bd Q^{n−1})(1 − cd Q^{n−1}) / ((1 − abcd Q^{2n−2})(1 − abcd Q^{2n−1})),
  B_n = a + 1/a − A_n − C_n.
KS の漸化式は p̃_n = p_n/(a^{−n}(ab,ac,ad;Q)_n) に対するものなので、monic の K_{(n)}(y) = p̃_n · a^{−n}(ab,ac,ad;Q)_n/(abcd Q^{n−1};Q)_n。
多項式は {j: 係数}（y^j、j = −n..n）の辞書で表す。
"""
from fractions import Fraction


def _qpoch(x, Q, n):
    r = 1
    for i in range(n):
        r *= 1 - x * Q**i
    return r


def askey_wilson_monic(nmax, Q, a, b, c, d):
    """K_{(n)}(y) for n = 0..nmax を辞書 {j: coeff} のリストで返す（係数は Fraction か sympy の数）。"""
    abcd = a * b * c * d
    p = [{0: Q**0}]  # p_0（Q と同じ型の 1）
    p_prev = {}
    out = []
    for n in range(0, nmax + 1):
        # monic 化
        # KS (3.1.4) の漸化式は p̃_n = p_n / (a^{−n} (ab, ac, ad; Q)_n) に対するもの。
        # p̃_n の (y + 1/y)^n の係数は (abcd Q^{n−1}; Q)_n / (a^{−n} (ab, ac, ad; Q)_n)。
        lead = (_qpoch(abcd * Q**(n - 1), Q, n) / (a**(-n) * _qpoch(a*b, Q, n) * _qpoch(a*c, Q, n) * _qpoch(a*d, Q, n))) if n > 0 else 1
        out.append({j: v / lead for j, v in p[n].items()})
        if n == nmax:
            break
        An = ((1 - abcd * Q**(n - 1)) * (1 - a * b * Q**n) * (1 - a * c * Q**n) * (1 - a * d * Q**n)
              / (a * (1 - abcd * Q**(2 * n - 1)) * (1 - abcd * Q**(2 * n))))
        Cn = (a * (1 - Q**n) * (1 - b * c * Q**(n - 1)) * (1 - b * d * Q**(n - 1)) * (1 - c * d * Q**(n - 1))
              / ((1 - abcd * Q**(2 * n - 2)) * (1 - abcd * Q**(2 * n - 1))))
        Bn = a + 1 / a - An - Cn
        # p_{n+1} = ((y + 1/y) p_n − B_n p_n − C_n p_{n−1}) / A_n
        new = {}
        for j, v in p[n].items():
            new[j + 1] = new.get(j + 1, 0) + v
            new[j - 1] = new.get(j - 1, 0) + v
            new[j] = new.get(j, 0) - Bn * v
        for j, v in (p[n - 1] if n >= 1 else {}).items():
            new[j] = new.get(j, 0) - Cn * v
        p.append({j: v / An for j, v in new.items() if v != 0})
    return out


def top_coeff_closed(j, Un, q, t):
    r"""probe 側の Askey–Wilson 多項式 K_{(n)}(y; t, q; 1, −1, q^{1/2}, −q^{1/2})（基底 t、monic）の
    y^{n−2j} の係数の閉じた形（U_n = t^n の有理式）：
      c_{2j}(n) = t^j (q/t; t²)_j / (t²; t²)_j × (1 − U_n)(1 + U_n t^{−2j}) Π_{i=1}^{j−1}(1 − U_n² t^{−2i}) / Π_{i=1}^{j}(1 − q U_n² t^{−2i−1}).
    【予想】有理数 (q,t) での有理補間（n ≤ 40、j ≤ 4）と記号計算（n ≤ 9、j ≤ 4）で確認。
    n → ∞（U_n → 0）で (q y²; t²)_∞/(t y²; t²)_∞ の係数（Rains–Warnaar の極限）に一致する。"""
    if j == 0:
        return Un**0
    lim = t**j * _qpoch(q / t, t**2, j) / _qpoch(t**2, t**2, j)
    num = (1 - Un) * (1 + Un * t**(-2 * j))
    for i in range(1, j):
        num *= 1 - Un**2 * t**(-2 * i)
    den = Un**0
    for i in range(1, j + 1):
        den *= 1 - q * Un**2 * t**(-2 * i - 1)
    return lim * num / den
