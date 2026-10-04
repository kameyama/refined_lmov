"""Laurent polynomials in A = a^{1/2}, Q = q^{1/2}, T = t^{1/2} with integer coefficients.

Stored as {(eA, eQ, eT): coeff}.
"""
import sympy as sp

A, Q, T = sp.symbols("A Q T")


def parse(text):
    expr = sp.expand(sp.sympify(text, locals={"A": A, "Q": Q, "T": T}))
    out = {}
    for term, coeff in expr.as_coefficients_dict().items():
        powers = term.as_powers_dict()
        key = (int(powers.get(A, 0)), int(powers.get(Q, 0)), int(powers.get(T, 0)))
        out[key] = out.get(key, 0) + int(coeff)
    return {k: v for k, v in out.items() if v}


def evaluate(poly, b, x, y):
    """Value at A = b, Q = x, T = y (field elements)."""
    s = 0
    for (ea, eq, et), c in poly.items():
        s = s + c * b ** ea * x ** eq * y ** et
    return s


def load_table(path):
    """{(knot, color): poly} from a TAB-separated reference file."""
    out = {}
    with open(path) as f:
        for line in f:
            knot, color, text = line.rstrip("\n").split("\t")
            out[(knot, color)] = parse(text)
    return out
