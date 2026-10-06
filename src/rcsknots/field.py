"""Elementwise arithmetic in F_p on numpy arrays: one entry per evaluation point.

`V` behaves like a field element (+, -, *, /, ** with integer exponents), so the
generic code in macdonald.py etc. runs on many evaluation points at once.
"""
from fractions import Fraction

import numpy as np

P31 = 2**31 - 1          # primes below 2^31 keep products inside int64
P31B = 2**31 - 19
P31C = 2147483587          # third prime, used only for verification


class V:
    __slots__ = ("a", "p", "_iv")      # _iv：逆元のキャッシュ（同じオブジェクトの負のべき・割り算で再利用）
    __array_priority__ = 1000

    def __init__(self, a, p):
        self.a = np.asarray(a, dtype=np.int64) % p
        self.p = p
        self._iv = None

    # -- helpers ------------------------------------------------------------
    def _coerce(self, other):
        if isinstance(other, V):
            return other.a
        if isinstance(other, Fraction):
            return np.int64(other.numerator % self.p) * _inv_scalar(other.denominator, self.p) % self.p
        return np.int64(int(other) % self.p)

    def _new(self, a):
        out = V.__new__(V)
        out.a, out.p, out._iv = a, self.p, None
        return out

    # -- arithmetic -----------------------------------------------------------
    def __add__(self, o):
        return self._new((self.a + self._coerce(o)) % self.p)

    __radd__ = __add__

    def __sub__(self, o):
        return self._new((self.a - self._coerce(o)) % self.p)

    def __rsub__(self, o):
        return self._new((self._coerce(o) - self.a) % self.p)

    def __neg__(self):
        return self._new((-self.a) % self.p)

    def __mul__(self, o):
        return self._new(self.a * self._coerce(o) % self.p)

    __rmul__ = __mul__

    def inv(self):
        if self._iv is None:
            if np.any(self.a == 0):
                raise ZeroDivisionError("division by zero at some evaluation point")
            self._iv = self._new(_powmod(self.a, self.p - 2, self.p))
            self._iv._iv = self
        return self._iv

    def __truediv__(self, o):
        o = o if isinstance(o, V) else self._new(np.broadcast_to(self._coerce(o), self.a.shape))
        return self * o.inv()

    def __rtruediv__(self, o):
        return self.inv() * o

    def __pow__(self, e):
        e = int(e)
        if e < 0:
            return self.inv() ** (-e)
        return self._new(_powmod(self.a, e, self.p))

    def __eq__(self, o):
        return bool(np.all(self.a == self._coerce(o)))

    def __ne__(self, o):
        return not self.__eq__(o)

    def __repr__(self):
        return f"V({self.a[:4]}..., p={self.p})"


def _powmod(a, e, p):
    result = np.ones_like(a)
    base = a % p
    while e:
        if e & 1:
            result = result * base % p
        base = base * base % p
        e >>= 1
    return result


def _inv_scalar(d, p):
    return np.int64(pow(int(d) % p, p - 2, p))


def make_conv(p, n):
    """Fraction -> constant V of length n."""
    def conv(c):
        return V(np.full(n, int(c.numerator % p) * pow(int(c.denominator) % p, p - 2, p) % p), p)
    return conv
