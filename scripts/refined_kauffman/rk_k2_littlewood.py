r"""K2：refined SO(2N) unknot の母関数 Z(v) = Σ_λ K_λ w_λ B_λ(v) が積の形（refined Littlewood 恒等式）をもつか。

積の形 Π_i f(v_i) Π_{i<j} h(v_i v_j) なら log Z = Σ_n α_n p_n + Σ_n β_n (p_n² − p_{2n})/2 なので、
3 次の部分に p_1³, p_1 p_2 の項が無く、p_3 だけになる。これを重み w_λ と基底 B_λ の候補ごとに調べる。
K_λ = dim^{D_N}_λ（AS）、U = t^N。
"""
import sys
sys.path.insert(0, "src")
import sympy as sp
from fractions import Fraction
from rcsknots.refined_kauffman.kauffman_refined import dim_DN, metric_DN, q, t, U
from rcsknots.macdonald import macdonald_P, macdonald_norm
from rcsknots.partitions import partitions, conjugate

conv = lambda c: sp.Rational(c.numerator, c.denominator)
P = {}
for n in range(0, 4):
    P.update(macdonald_P(n, q, t, conv))
p = {n: sp.Symbol(f"p{n}") for n in (1, 2, 3)}
eps = sp.Symbol("eps")


def basis(lam, kind):
    f = P[lam]
    if kind == "schur":
        f = macdonald_P(sum(lam), q, q, conv)[lam]
        return sum(c*sp.prod([p[k] for k in nu]) for nu, c in f.items())
    if kind == "P":       # P_λ(v; q, t)
        return sum(c*sp.prod([p[k] for k in nu]) for nu, c in f.items())
    if kind == "PT":      # P_{λ^T}(−v; t, q)
        f = P[conjugate(lam)]
        e = sum(c*sp.prod([p[k] for k in nu]) for nu, c in f.items())
        e = e.subs({q: t, t: q}, simultaneous=True)
        return e.subs({p[k]: (-1)**k * p[k] for k in p}, simultaneous=True)
    raise ValueError(kind)


def bA(lam):
    """type A の b_λ = 1/<P_λ,P_λ>（安定）。"""
    return 1/macdonald_norm(lam, q, t)


def weight(lam, kind):
    if kind == "1":
        return sp.Integer(1)
    if kind == "bA":
        return bA(lam)
    if kind == "bA_qt":
        return bA(lam)*(q/t)**sp.Rational(sum(lam), 2)
    if kind == "1/g":
        return 1/metric_DN(lam)
    if kind == "bA/g":
        return bA(lam)/metric_DN(lam)
    raise ValueError(kind)


def log_coeffs(wkind, bkind, Kfun=dim_DN, maxdeg=3, point=None):
    """point = {q: .., t: .., U: ..}（有理数）を与えると数値で計算する。"""
    Z = sp.Integer(1)
    for n in range(1, maxdeg+1):
        for lam in partitions(n):
            term = Kfun(lam)*weight(lam, wkind)*basis(lam, bkind)
            if point:
                term = term.subs(point)
            Z += sp.expand(term)*eps**n
    Z = sp.expand(Z)
    # log の級数：log(1 + X) = X − X²/2 + X³/3（X は eps の 1 次以上）
    X = Z - 1
    L = 0
    for k in range(1, maxdeg+1):
        L += (-1)**(k+1) * X**k / k
    L = sp.expand(L)
    out = {}
    for n in range(1, maxdeg+1):
        cn = sp.expand(L.coeff(eps, n))
        for nu in partitions(n):
            mono = sp.prod([p[k] for k in nu])
            c = sp.Poly(cn, *p.values()).coeff_monomial(mono)
            out[nu] = sp.factor(sp.cancel(c)) if not point else sp.nsimplify(c)
    return out


if __name__ == "__main__":
    pts = [{q: sp.Rational(3, 7), t: sp.Rational(5, 11), U: sp.Rational(2, 13)}, {q: sp.Rational(7, 5), t: sp.Rational(3, 2), U: sp.Rational(9, 4)}]
    for wkind in ("1", "bA", "bA_qt", "1/g", "bA/g"):
        for bkind in ("schur", "P", "PT"):
            res = []
            for pt in pts:
                out = log_coeffs(wkind, bkind, point=pt)
                res.append((out[(1, 1, 1)] == 0, out[(2, 1)] == 0, out[(1, 1)]))
            print(f"w = {wkind:6s} basis = {bkind:5s}: p1^3 = 0? {[r[0] for r in res]}  p1 p2 = 0? {[r[1] for r in res]}  (p1^2 係数: {res[0][2]})")
