"""Modified Macdonald polynomials H̃_μ(q,t) in the power-sum basis, and change of basis.

GN Thm 4.x: H̃_μ(q,t) = t^{n(μ)} φ[J_μ(q,t^{-1})],  J_μ = h_μ P_μ,  h_μ(q,t) = Π (1 − q^a t^{l+1}),
φ: p_k ↦ p_k / (1 − t^{−k}).
"""
from fractions import Fraction

from .partitions import partitions, boxes, arm, leg, n_stat
from .macdonald import macdonald_P


def h(lam, q, t):
    out = 1
    for i, j in boxes(lam):
        out = out * (1 - q ** arm(lam, i, j) * t ** (leg(lam, i, j) + 1))
    return out


def phi_factor(nu, t):
    """φ(p_ν) = p_ν · Π 1/(1 − t^{−ν_i})."""
    out = 1
    for k in nu:
        out = out / (1 - t ** (-k))
    return out


def modified_H(n, q, t, conv):
    """{μ: {ν: coefficient of p_ν in H̃_μ(q,t)}} for μ ⊢ n."""
    tinv = 1 / t
    P = macdonald_P(n, q, tinv, conv)
    out = {}
    for mu in partitions(n):
        c = t ** n_stat(mu) * h(mu, q, tinv)
        out[mu] = {nu: c * x * phi_factor(nu, t) for nu, x in P[mu].items()}
    return out


def solve_in_basis(target, basis, zero):
    """Coefficients x with target = Σ_μ x_μ basis[μ] (all vectors are dicts over the p-basis). Gaussian elimination."""
    keys = list(basis)
    coords = sorted({nu for v in basis.values() for nu in v} | set(target))
    # rows = coordinates, columns = basis elements, augmented with target
    rows = [[basis[mu].get(nu, zero) for mu in keys] + [target.get(nu, zero)] for nu in coords]
    ncol = len(keys)
    pivot_row = 0
    piv_cols = []
    for col in range(ncol):
        sel = None
        for r in range(pivot_row, len(rows)):
            if rows[r][col] != 0:
                sel = r
                break
        if sel is None:
            raise ValueError("singular basis")
        rows[pivot_row], rows[sel] = rows[sel], rows[pivot_row]
        inv = 1 / rows[pivot_row][col]
        rows[pivot_row] = [x * inv for x in rows[pivot_row]]
        for r in range(len(rows)):
            if r != pivot_row and rows[r][col] != 0:
                f = rows[r][col]
                rows[r] = [a - f * b for a, b in zip(rows[r], rows[pivot_row])]
        piv_cols.append(col)
        pivot_row += 1
    return {keys[c]: rows[i][ncol] for i, c in enumerate(piv_cols)}


def multiply_p(f, g):
    """Product of two symmetric functions in the p-basis."""
    out = {}
    for a, x in f.items():
        for b, y in g.items():
            key = tuple(sorted(a + b, reverse=True))
            out[key] = out[key] + x * y if key in out else x * y
    return out
