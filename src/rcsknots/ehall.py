"""Elliptic Hall algebra operators on the K-theory of Hilbert schemes (Gorsky–Negut 1304.3328, §4).

Basis: fixed points I_λ ↔ modified Macdonald H̃_λ(q,t).  Box weight χ(□) = q^{col} t^{row} (0-based).
Matrix elements of P̃_{kn,km} (n > 0) are GN eq. (fix3):

  (I_μ|P̃_{kn,km}|I_λ) = γ^{kn}/[k] · g_λ/g_μ · Σ_{SYT μ = λ + □_1 + … + □_{kn}} [Σ_{j<k} (qt)^j Π_{l≤j} χ_{n(k-l)+1}/χ_{n(k-l)}]
        · Π_i χ_i^{S_{m/n}(i)} (qtχ_i − 1) / Π_{i<kn} (1 − qt χ_{i+1}/χ_i) · Π_{i<j} ω^{-1}(χ_j/χ_i) · Π_{i, □∈λ} ω^{-1}(χ(□)/χ_i)

with ω(x) = (x−1)(x−qt)/((x−q)(x−t)), γ = (q−1)(t−1)/(qt(qt−1)), [k] = (q^k−1)(t^k−1)/((q−1)(t−1)),
g_λ = Π (1 − q^a t^{−l−1})(1 − q^{−a−1} t^l), S_{m/n}(i) = ⌊im/n⌋ − ⌊(i−1)m/n⌋.

In GN's SYT convention the labels decrease to the right and up, so □_{kn} is added first and □_1 last
(λ + □_i + … + □_{kn} is a partition for every i).  This orientation is checked by the tests
(P̃_{k,0} must be multiplication by φ(p_k) = p_k/(1 − t^{−k})).
"""
from functools import lru_cache

from .partitions import arm, leg, boxes


def addable(shape):
    """Boxes (row, col) that can be added to `shape`."""
    out = []
    for r in range(len(shape) + 1):
        c = shape[r] if r < len(shape) else 0
        if r == 0 or shape[r - 1] > c:
            out.append((r, c))
    return out


def add_box(shape, box):
    r, _ = box
    s = list(shape)
    if r == len(s):
        s.append(1)
    else:
        s[r] += 1
    return tuple(s)


class EHA:
    """Matrix elements of P̃_{kn,km} at the evaluation points (q, t) (field elements, e.g. `field.V`)."""

    def __init__(self, q, t, one):
        self.q, self.t, self.one = q, t, one
        self.qt = q * t
        self.gamma = (q - 1) * (t - 1) / (self.qt * (self.qt - 1))
        self._chi, self._g, self._omega_inv = {}, {}, {}

    def chi(self, box):
        if box not in self._chi:
            r, c = box
            self._chi[box] = self.q ** c * self.t ** r * self.one
        return self._chi[box]

    def g(self, lam):
        if lam not in self._g:
            out = self.one
            for i, j in boxes(lam):
                a, l = arm(lam, i, j), leg(lam, i, j)
                out = out * (1 - self.q ** a * self.t ** (-l - 1)) * (1 - self.q ** (-a - 1) * self.t ** l)
            self._g[lam] = out
        return self._g[lam]

    def omega_inv(self, old, new):
        """ω^{-1}(χ(old)/χ(new))."""
        key = (old, new)
        if key not in self._omega_inv:
            x = self.chi(old) / self.chi(new)
            q, t, qt = self.q, self.t, self.qt
            self._omega_inv[key] = (x - q) * (x - t) / ((x - 1) * (x - qt))
        return self._omega_inv[key]

    def bracket_k(self, k):
        return (self.q ** k - 1) * (self.t ** k - 1) / ((self.q - 1) * (self.t - 1))

    def apply(self, k, n, m, vec):
        """P̃_{kn,km} applied to vec = {λ: coefficient of I_λ}. Returns {μ: coefficient of I_μ}."""
        L = k * n
        S = [None] + [(i * m) // n - ((i - 1) * m) // n for i in range(1, L + 1)]
        r_steps = {n * (k - l): l for l in range(1, k)}      # label i at which r_l = χ_{i+1}/χ_i is completed
        pref = self.gamma ** L / self.bracket_k(k)
        out = {}
        for lam, coeff in vec.items():
            # states: (shape, last box) -> [A, B];  labels are assigned from L down to 1
            states = {(lam, None): [self.one, self.one]}
            for label in range(L, 0, -1):
                new_states = {}
                for (shape, last), (A, B) in states.items():
                    for box in addable(shape):
                        chi_b = self.chi(box)
                        w = chi_b ** S[label] * (self.qt * chi_b - 1)
                        for old in boxes_of(shape):
                            w = w * self.omega_inv(old, box)
                        if last is not None:            # pair (label, label+1): 1/(1 − qt χ_{label+1}/χ_label)
                            w = w / (1 - self.qt * self.chi(last) / chi_b)
                        if label in r_steps:
                            Bn = w * self.qt * (self.chi(last) / chi_b) * B
                            An = w * A + Bn
                        else:
                            An, Bn = w * A, w * B
                        key = (add_box(shape, box), box)
                        if key in new_states:
                            new_states[key][0] = new_states[key][0] + An
                            new_states[key][1] = new_states[key][1] + Bn
                        else:
                            new_states[key] = [An, Bn]
                states = new_states
            for (mu, _), (A, _) in states.items():
                term = coeff * pref * self.g(lam) / self.g(mu) * A
                out[mu] = out[mu] + term if mu in out else term
        return out


@lru_cache(maxsize=None)
def boxes_of(shape):
    return tuple((r, c) for r, row in enumerate(shape) for c in range(row))
