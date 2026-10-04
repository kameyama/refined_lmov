"""Reference formulas transcribed from the ancillary notebook of 1703.05408 (refined-LMOV-*.nb).

torusknotsym[n,p] / torusknotanti[n,p]: unreduced invariants of T(2,2p+1) colored by [n] / [1^n],
with the symmetric unknot normalization a^{-n/2} t^{n/2} (a;q)_n/(t;q)_n.
Arguments are field elements; `sqrt_*` give the chosen square roots q^{1/2}, t^{1/2}, a^{1/2}.
"""


def qp(x, y, n):
    """QPochhammer[x, y, n] = Π_{k=0}^{n-1} (1 - x y^k)."""
    out = 1
    for k in range(n):
        out = out * (1 - x * y ** k)
    return out


def _pw(x, h, e2):
    """x^{e2/2} given h = x^{1/2}."""
    return h ** e2


def torusknotsym(n, p, a, q, t, ha, hq, ht):
    pre = _pw(a, ha, -n) * _pw(t, ht, n) * qp(a, q, n) / qp(t, q, n) * (a * ht / hq) ** n
    s = 0
    for l in range(n + 1):
        e_q2 = -3 * l - l * l * (1 + 2 * p) + n * (n - p + 2 * n * p)   # exponent of q is e_q2/2
        e_t2 = n * (-3 + p) - 2 * l * p                                 # exponent of t is e_t2/2
        num = (-1) ** ((l + 3 * n + 4 * l * p) % 2) * a ** (n * p) * _pw(q, hq, e_q2) * _pw(t, ht, e_t2) \
            * (-1 + q ** (2 * l) * t) * qp(1 / a, 1 / q, l + n) * qp(1 / q, 1 / q, n) \
            * qp(1 / t, 1 / q, l) * qp(t / a, 1 / q, n - l)
        den = (-1 + t) * qp(1 / a, 1 / q, n) * qp(1 / q, 1 / q, l) * qp(1 / q, 1 / q, n - l) \
            * qp(1 / (q * t), 1 / q, l + n)
        s = s + num / den
    return pre * s


def torusknotanti(n, p, a, q, t, ha, hq, ht):
    pre = _pw(a, ha, -n) * _pw(t, ht, n * n) * qp(a, 1 / t, n) / qp(t, t, n) * (t / q) ** (p * n)
    s = 0
    for l in range(n + 1):
        e_q2 = n * (-2 + p) + 2 * l * (1 + p)
        e_t2 = -3 * l + l * l * (1 + 2 * p) - n * (1 + n + p + 2 * n * p)
        num = (-1) ** ((-l + n * (3 + 2 * p)) % 2) * a ** (n * p) * _pw(q, hq, e_q2) * _pw(t, ht, e_t2) \
            * (-1 + q * t ** (2 * l)) * qp(a, 1 / t, l + n) * qp(1 / q, 1 / t, l) \
            * qp(a * q, 1 / t, n - l) * qp(1 / t, 1 / t, n)
        den = (-1 + q) * qp(a, 1 / t, n) * qp(1 / t, 1 / t, l) * qp(1 / t, 1 / t, n - l) \
            * qp(1 / (q * t), 1 / t, l + n)
        s = s + num / den
    return pre * s
