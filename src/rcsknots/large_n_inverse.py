"""refined large N の式（KN の (refined-GT1), (refined-GT2)）を逆に解き、\\hat f_ρ から rCS_λ を求める（段階 D8、docs/research_plan_dgr_lmov.md）。

評価点 (A, Q, T) = (a, q, t)^{1/2}（有理数）で数値として解く。|λ| ≤ 3。
  GT1：f^q = M(t) \\hat f、f^q_μ は bps.f_formulas の記号式（\\bar{rCS}_λ = C[λ] と Adams の像 C1[d] の多項式）。
  GT2：f^{\\bar t} = M(q^{−1}) \\hat f、同じく f_formulas の第 2 成分。
各次数 n で f_μ（|μ| = n）は C[λ]（|λ| = n）について一次なので、低い次数を代入した一次方程式を解く。
C1[d] = \\bar{rCS}_{[1]}(a^d, q^d, t^d) は、\\hat f_{[1]} を点 (A^d, Q^d, T^d) で評価して次数 1 の式から求める。
reduced は rCS_λ = C[λ] / unknot_λ（KN のノートブックの対称な正規化、conventions.md §2）。
"""
from functools import lru_cache
from fractions import Fraction as Fr

import sympy as sp

from .bps import C, C1, M_matrix, Q, T, f_formulas
from .partitions import arm, boxes, leg, partitions

Asym = sp.Symbol("A")


def unknot(lam, A, Qv, Tv):
    """\\bar{rCS}_λ(○)（対称な正規化）を数値で。"""
    k = sum(lam)
    out = Fr(Tv) ** k / Fr(A) ** k
    for i, j in boxes(lam):
        out *= (Fr(Tv) ** (2 * i) - Fr(A) ** 2 * Fr(Qv) ** (2 * j)) / (1 - Fr(Qv) ** (2 * arm(lam, i, j)) * Fr(Tv) ** (2 * (leg(lam, i, j) + 1)))
    return out


def ev_poly(poly, A, Qv, Tv):
    """{(e_A, e_Q, e_T): c} を点で評価。"""
    return sum(Fr(int(c)) * Fr(A) ** e[0] * Fr(Qv) ** e[1] * Fr(Tv) ** e[2] for e, c in poly.items())


@lru_cache(maxsize=None)
def _formulas(nmax):
    return f_formulas(nmax)


def _num(expr, Qv, Tv):
    return expr.subs({Q: sp.Rational(Qv.numerator, Qv.denominator), T: sp.Rational(Tv.numerator, Tv.denominator)})


def rcs_from_hat(hat, A, Qv, Tv, which="q", nmax=3):
    """hat = {ρ: \\hat f_ρ の多項式（小文字 (A,Q,T)、a 次数のずらしの前）}。点での \\bar{rCS}_λ を GT1（which="q"）か GT2（"t"）から解く。

    返り値 {λ: \\bar{rCS}_λ の値}（unreduced）。
    """
    F = _formulas(nmax)
    comp = 0 if which == "q" else 1
    Qv, Tv, A = Fr(Qv), Fr(Tv), Fr(A)
    # C1[d]：点 (A^d, Q^d, T^d) で次数 1 の式を解く
    c1 = {}
    for d in range(2, nmax + 1):
        h1 = ev_poly(hat[(1,)], A ** d, Qv ** d, Tv ** d)
        m1 = sp.nsimplify(M_matrix(1, T if comp == 0 else 1 / Q)[0, 0]).subs({Q: sp.Rational(Qv.numerator ** d, Qv.denominator ** d), T: sp.Rational(Tv.numerator ** d, Tv.denominator ** d)})
        f1 = m1 * sp.Rational(h1.numerator, h1.denominator)
        e1 = F[(1,)][comp].subs({Q: sp.Rational(Qv.numerator ** d, Qv.denominator ** d), T: sp.Rational(Tv.numerator ** d, Tv.denominator ** d)})
        sol = sp.solve(sp.Eq(e1, f1), C((1,)))
        c1[C1(d)] = sol[0]
    known = dict(c1)
    for n in range(1, nmax + 1):
        parts = partitions(n)
        Mm = M_matrix(n, T if comp == 0 else 1 / Q)
        h = sp.Matrix([sp.Rational(*(lambda v: (v.numerator, v.denominator))(ev_poly(hat[r], A, Qv, Tv))) for r in parts])
        fvals = _num(Mm, Qv, Tv) * h
        eqs = [sp.Eq(_num(F[mu][comp], Qv, Tv).subs(known), fvals[i]) for i, mu in enumerate(parts)]
        unknowns = [C(lam) for lam in parts]
        sol = sp.solve(eqs, unknowns, dict=True)
        if len(sol) != 1:
            raise ArithmeticError(f"degree {n}: {len(sol)} solutions")
        known.update(sol[0])
    return {lam: Fr(str(known[C(lam)])) for n in range(1, nmax + 1) for lam in partitions(n)}


def reduced(unred, A, Qv, Tv):
    return {lam: v / unknot(lam, A, Qv, Tv) for lam, v in unred.items()}
