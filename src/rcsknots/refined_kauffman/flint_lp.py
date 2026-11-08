r"""変数 (r, s) の Laurent 多項式を、python-flint の fmpz_mpoly と単項式のずらし (e_r, e_s) で表す小さいクラス（値は n · r^{e_r} s^{e_s}）。
flint_rf.RF との違いは約分（gcd）をしないこと。DAHA の作用素（T_i、π_1、X、評価）の係数はすべて Laurent 多項式なので、
有理関数の係数（Macdonald 多項式の展開係数など）を最後の一次結合だけに回せば、途中の計算はこのクラスで閉じる。
表示は一意でない（n が r, s で割り切れてもよい）。0 の判定は n.is_zero()、等号は差の 0 判定で行う。"""
import flint

from .flint_rf import CTX, RF, _ONE

_MONO = {}


def _mono(a, b):
    """r^a s^b（a, b ≥ 0）の fmpz_mpoly（キャッシュする）。"""
    key = (a, b)
    m = _MONO.get(key)
    if m is None:
        m = _MONO[key] = CTX.from_dict({(a, b): 1})
    return m


class LP:
    __slots__ = ("n", "er", "es")

    def __init__(self, n, er=0, es=0):
        self.n, self.er, self.es = n, er, es

    @staticmethod
    def coerce(x):
        if isinstance(x, LP):
            return x
        if isinstance(x, int):
            return LP(CTX.from_dict({(0, 0): x}) if x else CTX.from_dict({}))
        raise TypeError(type(x))

    def __add__(a, b):
        b = LP.coerce(b)
        if a.er == b.er and a.es == b.es:
            return LP(a.n + b.n, a.er, a.es)
        er, es = min(a.er, b.er), min(a.es, b.es)
        na = a.n if (a.er == er and a.es == es) else a.n * _mono(a.er - er, a.es - es)
        nb = b.n if (b.er == er and b.es == es) else b.n * _mono(b.er - er, b.es - es)
        return LP(na + nb, er, es)

    __radd__ = __add__

    def __neg__(a):
        return LP(-a.n, a.er, a.es)

    def __sub__(a, b):
        return a + (-LP.coerce(b))

    def __rsub__(a, b):
        return LP.coerce(b) + (-a)

    def __mul__(a, b):
        if isinstance(b, int):
            return LP(a.n * b, a.er, a.es)
        return LP(a.n * b.n, a.er + b.er, a.es + b.es)

    __rmul__ = __mul__

    def _unit(self):
        """単項式 ±r^i s^j なら (符号, i, j)、そうでなければ None。"""
        d = self.n.to_dict()
        if len(d) != 1:
            return None
        (i, j), c = next(iter(d.items()))
        if c not in (1, -1):
            return None
        return int(c), i + self.er, j + self.es

    def __pow__(a, k):
        if k >= 0:
            return LP(a.n ** k, a.er * k, a.es * k)
        u = a._unit()
        if u is None:
            raise ValueError("負冪は ±単項式だけ")
        c, i, j = u
        return LP(CTX.from_dict({(0, 0): c ** (-k)}), i * k, j * k)          # (±1)^k = (±1)^{−k}

    def __truediv__(a, b):
        u = LP.coerce(b)._unit()
        if u is None:
            raise ValueError("割り算は ±単項式だけ")
        c, i, j = u
        return LP(a.n * c, a.er - i, a.es - j)

    def __rtruediv__(a, b):
        return LP.coerce(b) * (a ** -1)

    def is_zero(self):
        return self.n.is_zero()

    def __eq__(a, b):
        return (a - LP.coerce(b)).is_zero()

    def to_rf(self):
        """RF（flint_rf）に直す。"""
        num, den = self.n, _ONE
        if self.er >= 0 and self.es >= 0:
            num = num * _mono(self.er, self.es)
        else:
            num = num * _mono(max(self.er, 0), max(self.es, 0))
            den = _mono(max(-self.er, 0), max(-self.es, 0))
        return RF(num, den)

    def __repr__(self):
        return f"LP({self.n})·r^{self.er}·s^{self.es}"


R = LP(_mono(1, 0))
S = LP(_mono(0, 1))
ONE = LP(_ONE)
ZERO = LP(CTX.from_dict({}))
