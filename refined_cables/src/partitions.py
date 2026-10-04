"""Partitions and their combinatorial statistics."""
from functools import lru_cache
from math import factorial
from collections import Counter


@lru_cache(maxsize=None)
def partitions(n, max_part=None):
    """All partitions of n as tuples, in decreasing lexicographic order."""
    if max_part is None:
        max_part = n
    if n == 0:
        return ((),)
    out = []
    for k in range(min(n, max_part), 0, -1):
        for rest in partitions(n - k, k):
            out.append((k,) + rest)
    return tuple(out)


def conjugate(lam):
    if not lam:
        return ()
    return tuple(sum(1 for x in lam if x > j) for j in range(lam[0]))


def boxes(lam):
    """Boxes (i, j), 0-based row i and column j."""
    return [(i, j) for i, row in enumerate(lam) for j in range(row)]


def arm(lam, i, j):
    return lam[i] - j - 1


def leg(lam, i, j):
    return conjugate(lam)[j] - i - 1


def n_stat(lam):
    """n(λ) = Σ (i-1) λ_i."""
    return sum(i * x for i, x in enumerate(lam))


def z(lam):
    """z_λ = Π i^{m_i} m_i!."""
    out = 1
    for part, mult in Counter(lam).items():
        out *= part ** mult * factorial(mult)
    return out


def dominates(lam, mu):
    """λ ≥ μ in dominance order (same size assumed)."""
    s1 = s2 = 0
    for k in range(max(len(lam), len(mu))):
        s1 += lam[k] if k < len(lam) else 0
        s2 += mu[k] if k < len(mu) else 0
        if s1 < s2:
            return False
    return True


def scale(lam, m):
    """mλ = (mλ_1, mλ_2, ...): the partition of p_λ[p_k -> p_{mk}]."""
    return tuple(m * x for x in lam)
