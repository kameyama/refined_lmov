"""Polynomial interpolation over F_p (vectorized with numpy, p < 2^31)."""
import numpy as np


def _inv(x, p):
    return pow(int(x) % p, p - 2, p)


def _powmod_vec(a, e, p):
    result = np.ones_like(a)
    base = a % p
    while e:
        if e & 1:
            result = result * base % p
        base = base * base % p
        e >>= 1
    return result


def newton_to_monomial_direct(xs, vals, p):
    """Coefficients c (shape (N, ...)) with Σ_k c_k x^k = vals at xs (vals shape (N, ...))."""
    xs = [int(v) % p for v in xs]
    N = len(xs)
    d = np.array(vals, dtype=np.int64) % p
    # divided differences: level j uses the level j−1 values d[i], d[i−1] (vectorized over i)
    xa = np.array(xs, dtype=np.int64)
    shape = (-1,) + (1,) * (d.ndim - 1)
    for j in range(1, N):
        invs = _powmod_vec((xa[j:] - xa[:N - j]) % p, p - 2, p).reshape(shape)
        d[j:] = (d[j:] - d[j - 1:N - 1]) % p * invs % p
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


_VINV = {}


def _matmul_mod(A, B, p):
    """A @ B mod p for int64 arrays with entries in [0, p), p < 2^31.

    16 ビットずつに分けた 4 つの積を float64 の行列積（BLAS）で計算する。各積の成分は < 2^32·K（K = A.shape[1]）で、
    K < 2^21 なら 2^53 未満に収まり float64 で厳密。
    """
    assert A.shape[1] < 2**21
    a0, a1 = (A & 0xFFFF).astype(np.float64), (A >> 16).astype(np.float64)
    b0, b1 = (B & 0xFFFF).astype(np.float64), (B >> 16).astype(np.float64)

    def mm(x, y):
        return (x @ y).astype(np.int64) % p

    mid = (mm(a0, b1) + mm(a1, b0)) % p
    return (mm(a0, b0) + mid * 65536 % p + mm(a1, b1) * (65536 * 65536 % p) % p) % p


def newton_to_monomial(xs, vals, p):
    """Coefficients c (shape (N, ...)) with Σ_k c_k x^k = vals at xs, as (inverse Vandermonde) @ vals mod p.

    The inverse Vandermonde matrix is computed once per point set with newton_to_monomial_direct and cached.
    """
    xs = [int(v) % p for v in xs]
    N = len(xs)
    key = (p, tuple(xs))
    if key not in _VINV:
        if len(_VINV) > 8:
            _VINV.clear()
        _VINV[key] = newton_to_monomial_direct(xs, np.eye(N, dtype=np.int64), p)
    d = np.array(vals, dtype=np.int64) % p
    return _matmul_mod(_VINV[key], d.reshape(N, -1), p).reshape(d.shape)
