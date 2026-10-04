"""Stage 1: Macdonald foundation."""
import random
from fractions import Fraction

import pytest
import sympy as sp
from flint import nmod

from partitions import partitions, z, dominates
from characters import chi
from macdonald import macdonald_P, macdonald_norm, principal_specialization, pairing_weight, _p_in_m

PR = 2**61 - 1
q, t = sp.symbols("q t")


def nconv(c):
    return nmod(c.numerator, PR) / nmod(c.denominator, PR)


def rand_qt():
    rng = random.Random(1)
    return nmod(rng.randrange(2, PR), PR), nmod(rng.randrange(2, PR), PR)


@pytest.mark.parametrize("n", range(1, 9))
def test_schur_limit(n):
    """q = t: P_λ = s_λ = Σ_ν χ_λ(ν)/z_ν p_ν."""
    x = nmod(12345, PR)
    P = macdonald_P(n, x, x, nconv)
    for lam in partitions(n):
        for nu in partitions(n):
            assert P[lam].get(nu, 0) == nconv(Fraction(chi(lam, nu), z(nu)))


@pytest.mark.parametrize("n", range(1, 10))
def test_orthogonality_and_norm(n):
    qq, tt = rand_qt()
    P = macdonald_P(n, qq, tt, nconv)
    w = {nu: pairing_weight(nu, qq, tt, nconv) for nu in partitions(n)}
    for lam in partitions(n):
        for mu in partitions(n):
            ip = sum((P[lam].get(nu, 0) * P[mu].get(nu, 0) * w[nu] for nu in partitions(n)), nmod(0, PR))
            expected = macdonald_norm(lam, qq, tt) if lam == mu else 0
            assert ip == expected


@pytest.mark.parametrize("n", range(1, 9))
def test_triangularity(n):
    """Converting back to monomials: P_λ = m_λ + (dominated terms)."""
    qq, tt = rand_qt()
    P = macdonald_P(n, qq, tt, nconv)
    R = _p_in_m(n)
    for lam in partitions(n):
        coeffs = {mu: sum((c * R[nu][mu] for nu, c in P[lam].items()), nmod(0, PR)) for mu in partitions(n)}
        assert coeffs[lam] == 1
        for mu, c in coeffs.items():
            if c != 0:
                assert dominates(lam, mu)


def test_explicit_low_degree():
    """Appendix A of 1703.05408 (P in power sums)."""
    P = {n: macdonald_P(n, q, t, lambda c: sp.Rational(c.numerator, c.denominator)) for n in (2, 3)}
    expected = {
        (2,): {(1, 1): (1 - t) * (1 + q) / (1 - t * q) / 2, (2,): (1 + t) * (1 - q) / (1 - t * q) / 2},
        (1, 1): {(1, 1): sp.Rational(1, 2), (2,): -sp.Rational(1, 2)},
        (3,): {(1, 1, 1): (1 + q) * (1 - q**3) * (1 - t) ** 2 / ((1 - q) * (1 - t * q) * (1 - t * q**2)) / 6,
               (2, 1): (1 - q) * (1 - t**2) * (1 - q**3) / ((1 - q) * (1 - t * q) * (1 - t * q**2)) / 2,
               (3,): (1 - q) * (1 - q**2) * (1 - t**3) / ((1 - t) * (1 - t * q) * (1 - t * q**2)) / 3},
        (2, 1): {(1, 1, 1): (1 - t) * (2 * q * t + q + t + 2) / (1 - q * t**2) / 6,
                 (2, 1): (1 + t) * (t - q) / (1 - q * t**2) / 2,
                 (3,): -(1 - q) * (1 - t**3) / ((1 - t) * (1 - q * t**2)) / 3},
        (1, 1, 1): {(1, 1, 1): sp.Rational(1, 6), (2, 1): -sp.Rational(1, 2), (3,): sp.Rational(1, 3)},
    }
    for lam, exp in expected.items():
        got = P[sum(lam)][lam]
        for nu, val in exp.items():
            assert sp.simplify(got.get(nu, 0) - val) == 0


@pytest.mark.parametrize("N", [2, 3, 4])
def test_principal_specialization(N):
    """Closed formula vs direct substitution x = (1, t, ..., t^{N-1})."""
    qq, tt = rand_qt()
    a = tt ** N
    for n in range(1, 7):
        P = macdonald_P(n, qq, tt, nconv)
        for lam in partitions(n):
            direct = nmod(0, PR)
            for nu, c in P[lam].items():
                term = c
                for k in nu:
                    term *= sum((tt ** (k * i) for i in range(N)), nmod(0, PR))
                direct += term
            assert direct == principal_specialization(lam, a, qq, tt)
