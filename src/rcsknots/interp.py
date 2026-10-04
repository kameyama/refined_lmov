"""Polynomial interpolation over F_p (vectorized with numpy, p < 2^31)."""
import numpy as np


def _inv(x, p):
    return pow(int(x) % p, p - 2, p)


def newton_to_monomial(xs, vals, p):
    """Coefficients c (shape (N, ...)) with Σ_k c_k x^k = vals at xs (vals shape (N, ...))."""
    xs = [int(v) % p for v in xs]
    N = len(xs)
    d = np.array(vals, dtype=np.int64) % p
    # divided differences, in place
    for j in range(1, N):
        for i in range(N - 1, j - 1, -1):
            inv = _inv(xs[i] - xs[i - j], p)
            d[i] = (d[i] - d[i - 1]) % p * inv % p
    # Newton form -> monomial form (Horner)
    c = np.zeros_like(d)
    for k in range(N - 1, -1, -1):
        # c(x) <- c(x)·(x − x_k) + d_k
        shifted = np.zeros_like(c)
        shifted[1:] = c[:-1]
        c = (shifted - xs[k] * c) % p
        c[0] = (c[0] + d[k]) % p
    return c


def symmetric_lift(c, p):
    c = np.asarray(c, dtype=np.int64) % p
    return np.where(c > p // 2, c - p, c)
