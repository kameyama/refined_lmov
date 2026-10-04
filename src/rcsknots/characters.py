"""Characters of symmetric groups (Murnaghan–Nakayama rule)."""
from functools import lru_cache


@lru_cache(maxsize=None)
def chi(lam, nu):
    """χ_λ evaluated on the conjugacy class of cycle type ν."""
    if sum(lam) != sum(nu):
        raise ValueError("size mismatch")
    if not nu:
        return 1
    k, rest = nu[0], nu[1:]
    # beta-numbers of λ: λ_i + (ℓ - i), strictly decreasing
    ell = len(lam)
    beta = [lam[i] + ell - 1 - i for i in range(ell)]
    beads = set(beta)
    total = 0
    for b in beta:
        if b - k >= 0 and (b - k) not in beads:
            sign = (-1) ** sum(1 for c in beta if b - k < c < b)
            new = sorted((beads - {b}) | {b - k}, reverse=True)
            mu = tuple(x for x in (new[i] - (ell - 1 - i) for i in range(ell)) if x > 0)
            total += sign * chi(mu, rest)
    return total
