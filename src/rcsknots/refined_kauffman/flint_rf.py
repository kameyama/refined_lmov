r"""変数 (r, s) の有理関数を python-flint の fmpz_mpoly の分子・分母で表す小さいクラス。
sympy の疎多項式の gcd は発見的 gcd が失敗すると止まる（次数 6 の Macdonald P で起きた）ので、gcd は flint に任せる。
分母は先頭係数が正になるように正規化する。Laurent 多項式 ⇔ 分母が 1 項。整数係数の Laurent 多項式 ⇔ 分母が係数 1 の 1 項。"""
from fractions import Fraction

import flint
import sympy as sp

CTX = flint.fmpz_mpoly_ctx.get(("r", "s"), "lex")
_R, _S = CTX.gens()
_ONE = CTX.from_dict({(0, 0): 1})


def _const(k):
    return CTX.from_dict({(0, 0): int(k)}) if k else CTX.from_dict({})


class RF:
    __slots__ = ("n", "d")

    def __init__(self, n, d=None, normalize=True):
        if d is None:
            d = _ONE
        if normalize:
            if n.is_zero():
                d = _ONE
            else:
                g = n.gcd(d)
                if g != _ONE:
                    n, d = n / g, d / g
                if d.coeffs()[0] < 0:
                    n, d = -n, -d
        self.n, self.d = n, d

    @staticmethod
    def coerce(x):
        if isinstance(x, RF):
            return x
        if isinstance(x, Fraction):
            return RF(_const(x.numerator), _const(x.denominator))
        if isinstance(x, int):
            return RF(_const(x), _ONE, normalize=False)
        if isinstance(x, sp.Basic):
            return from_sympy(x)
        raise TypeError(type(x))

    def __add__(a, b):
        b = RF.coerce(b)
        if a.d == b.d:
            return RF(a.n + b.n, a.d)
        return RF(a.n * b.d + b.n * a.d, a.d * b.d)

    __radd__ = __add__

    def __neg__(a):
        return RF(-a.n, a.d, normalize=False)

    def __sub__(a, b):
        return a + (-RF.coerce(b))

    def __rsub__(a, b):
        return RF.coerce(b) + (-a)

    def __mul__(a, b):
        b = RF.coerce(b)
        return RF(a.n * b.n, a.d * b.d)

    __rmul__ = __mul__

    def __truediv__(a, b):
        b = RF.coerce(b)
        if b.n.is_zero():
            raise ZeroDivisionError
        return RF(a.n * b.d, a.d * b.n)

    def __rtruediv__(a, b):
        return RF.coerce(b) / a

    def __pow__(a, k):
        if k < 0:
            return RF(a.d ** (-k), a.n ** (-k))
        return RF(a.n ** k, a.d ** k, normalize=False)

    def __eq__(a, b):
        b = RF.coerce(b)
        return a.n * b.d == b.n * a.d

    def __hash__(self):
        return hash((str(self.n), str(self.d)))

    def is_zero(self):
        return self.n.is_zero()

    def is_laurent(self):
        return len(self.d.coeffs()) == 1

    def is_integer_laurent(self):
        return len(self.d.coeffs()) == 1 and self.d.coeffs()[0] == 1

    def adams(self, k):
        """r → r^k、s → s^k。"""
        f = lambda p: CTX.from_dict({(e[0] * k, e[1] * k): c for e, c in p.to_dict().items()}) if not p.is_zero() else p
        return RF(f(self.n), f(self.d), normalize=False)

    def num_divisible(self, x):
        """分子が多項式 x（RF の分子として与える）で割り切れるか。"""
        x = RF.coerce(x).n
        g = self.n.gcd(x)
        return g == x or g == -x

    def to_sympy(self):
        r, s = sp.symbols("r s", positive=True)
        cv = lambda p: sp.Add(*[int(c) * r**e[0] * s**e[1] for e, c in p.to_dict().items()])   # 一度に足す（項ごとの sum は大きな式で数百倍遅い）
        return cv(self.n) / cv(self.d)

    def __repr__(self):
        return f"RF({self.n})/({self.d})"


def from_sympy(e):
    """sympy の有理式（変数は r, s の名前の記号）→ RF。負冪は分母へ。"""
    num, den = sp.fraction(sp.together(sp.sympify(e)))
    syms = sorted(sp.sympify(e).free_symbols, key=str)
    names = {str(x) for x in syms}
    if not names <= {"r", "s"}:
        raise ValueError(names)
    rs = {str(x): x for x in syms}
    gens = [rs.get("r", sp.Symbol("r")), rs.get("s", sp.Symbol("s"))]

    def cv(p):
        P = sp.Poly(sp.expand(p), *gens, domain="QQ")
        den_ = 1
        for c in P.coeffs():
            den_ = den_ * int(sp.Rational(c).q) // __import__("math").gcd(den_, int(sp.Rational(c).q))
        return CTX.from_dict({m: int(sp.Rational(c) * den_) for m, c in P.terms()}), den_

    n, dn = cv(num)
    d, dd = cv(den)
    return RF(n * dd, d * dn)


R = RF(_R, _ONE, normalize=False)
S = RF(_S, _ONE, normalize=False)
ONE = RF(_ONE, _ONE, normalize=False)


def conv(c):
    return RF.coerce(Fraction(c))
