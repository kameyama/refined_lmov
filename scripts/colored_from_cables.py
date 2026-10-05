"""2 本ケーブルの HOMFLY-PT からケーブル公式を逆に解き、[2], [1,1] 彩色 HOMFLY-PT（reduced、framing 0）を求める（ステージ T4）。

入力：data/reference/two_cable_homfly.txt（scripts/two_cable_homfly.sage、s = 1, 3）。
W_□(C_{2,s}(K)) = θ_□^{2s}·Σ_μ c^{(2)}_{□μ} θ_μ^{−s/2} W_μ(K)（src/rcsknots/homfly.py の cable_unreduced、sign = −1、ステージ 4 で検証済み）。
変数は KN の (a, Q = q^{1/2}) で、Sage の値は tests/test_stage4_cable_homfly.py の sage_to_kn で変換する（DGR の向き）。
使い方：python scripts/colored_from_cables.py 3_1 4_1 9_42
"""
import sys

import sympy as sp

from rcsknots.homfly import Q, a, adams_coeffs, theta, unknot

L, M, Abold = sp.symbols("L M A")


def sage_to_kn(expr):
    e = sp.expand(expr.subs({L: sp.I * Abold, M: -sp.I * (Q - 1 / Q)}, simultaneous=True))
    e = e.subs({Abold: 1 / Abold, Q: 1 / Q}, simultaneous=True)
    return sp.expand(e.subs(Abold, sp.sqrt(a)))


def load(path="data/reference/two_cable_homfly.txt"):
    out = {}
    for line in open(path):
        name, s, h = line.rstrip("\n").split("\t")
        out[(name, int(s))] = sage_to_kn(sp.sympify(h, locals={"L": L, "M": M}))
    return out


def colored(name, data=None):
    """{(2,): P_[2], (1,1): P_[1,1]}（reduced、framing 0、KN の (a, Q)）。"""
    data = data or load()
    c = adams_coeffs((1,), 2)
    mus = [(2,), (1, 1)]
    W = sp.symbols("W2 W11")
    eqs = []
    for s in (1, 3):
        lhs = data[(name, s)] * unknot((1,))
        rhs = theta((1,), 2 * s) * sum(c[mu] * theta(mu, sp.Rational(-s, 2)) * w for mu, w in zip(mus, W))
        eqs.append(sp.Eq(lhs, rhs))
    sol = sp.solve(eqs, W, dict=True)[0]
    return {mu: sp.factor(sp.simplify(sol[w] / unknot(mu))) for mu, w in zip(mus, W)}


if __name__ == "__main__":
    data = load()
    for name in sys.argv[1:]:
        for mu, p in colored(name, data).items():
            print(name, mu, sp.expand(p))
