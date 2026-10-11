"""Nawata–Oblomkov（arXiv:1510.01795、§2.3、TeX knot-review.tex の 1084〜1093 行）の 6_2, 6_3 の [r] の四重次数の閉公式（予想式）を、
t_c 次数の太字に直す：P(𝐚, 𝐪, 𝐭) = P̃(a = 𝐚, Q = 𝐪, t_r = 𝐪⁻¹, t_c = 𝐭𝐪)（scripts/t3_pretzel.py と同じ規則）。
D10 の予言の照合に使う（docs/dgr_lmov/research_plan.md）。
"""
import sympy as sp

a, Q, tr, tc = sp.symbols("a Q t_r t_c")
q, t = sp.symbols("q t")


def poch(x, n, base):
    out = sp.Integer(1)
    for i in range(n):
        out *= 1 - x * base ** i
    return out


def qbin(n, k, base):
    return sp.cancel(poch(base, n, base) / (poch(base, k, base) * poch(base, n - k, base)))


def quad(knot, r):
    b = tc ** 2
    tot = 0
    for k in range(r + 1):
        for j in range(k + 1):
            for i in range(j + 1):
                if knot == "6_2":
                    pref = (-1) ** (r - k) * Q ** (2 * (i + j - k)) * tr ** (i + j - k - r) * tc ** (i * i - j * j + 2 * j * k - 2 * k * r + k - r)
                else:
                    pref = a ** (-2 * k) * Q ** (2 * (i - 2 * j + k)) * tr ** (i - 2 * j - k) * tc ** (i * i + k * (-2 * j + k - 2 * r))
                tot += pref * qbin(r, k, b) * qbin(k, j, b) * qbin(j, i, b) \
                    * poch(-a ** 2 * Q ** -2 * tr * tc, j, b) * poch(-a ** 2 * Q ** -2 * tr * tc ** (2 * i + 1), k - j, b) \
                    * poch(-a ** 2 * Q ** 2 * tr ** 3 * tc ** (2 * r + 1), k, b)
    return sp.expand(sp.cancel(tot))


def tc_bold(knot, r):
    e = sp.expand(quad(knot, r).subs({Q: q, tr: 1 / q, tc: t * q}, simultaneous=True) * (a * q * t) ** 200)
    return {(i - 200, j - 200, k - 200): int(c) for (i, j, k), c in sp.Poly(e, a, q, t).terms()}
