"""Compare our normalized refined invariants with DAHA superpolynomials (CD14) up to a monomial.

DAHA variables: (q, t) are those of 1703.05408 and a_DAHA = −a_KN (docs/conventions.md §1.2).
"""
import sympy as sp

q, t, a = sp.symbols("q t a")


def daha_to_kn(poly_expr):
    """DAHA polynomial (sympy in q, t, a) -> {(e_A, e_Q, e_T): c} in KN's half-variables, with a_DAHA = −a_KN."""
    e = sp.expand(poly_expr.subs(a, -a))
    out = {}
    for term, c in e.as_coefficients_dict().items():
        pw = term.as_powers_dict()
        key = (2 * int(pw.get(a, 0)), 2 * int(pw.get(q, 0)), 2 * int(pw.get(t, 0)))
        out[key] = out.get(key, 0) + int(c)
    return {k: v for k, v in out.items() if v}


def monomial_ratio(ours, theirs):
    """(sign, shift) with ours = sign · A^{s0} Q^{s1} T^{s2} · theirs exactly, or None."""
    top_o, top_t = max(ours), max(theirs)
    shift = tuple(x - y for x, y in zip(top_o, top_t))
    sign = 1 if ours[top_o] * theirs[top_t] > 0 else -1
    moved = {tuple(k[i] + shift[i] for i in range(3)): sign * v for k, v in theirs.items()}
    return (sign, shift) if moved == ours else None
