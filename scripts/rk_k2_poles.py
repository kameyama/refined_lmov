r"""K2：unknot の 2 次（と 3 次）で、Kauffman 側 F_G と composite 側 F_R の Q（U）の極が
定数 c（q, t の有理式）の組み合わせ F_G − c F_R で打ち消せるかを調べる。

Kauffman 側の unknot：AS の D_N（Koornwinder ±1, ±q^{1/2}）か LRW 型（±1, ±t^{1/2}）。
composite 側：U(M) の Macdonald 次元、A ∈ {Q, Qt, Q/t, Qt², U, U²}（Q = U²/t）。
基底は P_λ(v;q,t)、重み w_λ ∈ {1, m_λ = (q/t)^{|λ|/2} b_λ}。
有理点 (q, t) で U を記号のまま計算し、各成分（p_1², p_2）の U の極と留数を比べる。
"""
import sys
sys.path.insert(0, "src")
import sympy as sp
from rcsknots.kauffman_refined import dim_DN, dim_composite, koornwinder_eval, q, t, U, A
from rcsknots.macdonald import macdonald_P, macdonald_norm
from rcsknots.partitions import partitions

conv = lambda c: sp.Rational(c.numerator, c.denominator)
P = {}
for n in range(0, 3):
    P.update(macdonald_P(n, q, t, conv))
p1, p2 = sp.symbols("p1 p2")
pv = {1: p1, 2: p2}


def Pbasis(lam):
    return sum(c*sp.prod([pv[k] for k in nu]) for nu, c in P[lam].items())


def weight(lam, kind):
    if kind == "1":
        return sp.Integer(1)
    if kind == "m":
        return (q/t)**sp.Rational(sum(lam), 2)/macdonald_norm(lam, q, t)
    raise ValueError


def F2(Z1, Z2):
    """Z = 1 + Z1 + Z2（次数 1, 2）の log の次数 2 の部分。"""
    return sp.expand(Z2 - Z1**2/2)


def components(expr):
    e = sp.expand(expr)
    return {"p1^2": e.coeff(p1, 2).subs(p2, 0), "p2": e.coeff(p2, 1).subs(p1, 0)}


def poles(expr, point):
    """U の有理式として、点 point（q, t の値）での分母の因子（U の多項式）を返す。"""
    e = sp.factor(sp.cancel(sp.together(expr.subs(point))))
    num, den = sp.fraction(e)
    fl = sp.factor_list(den, U)
    return {sp.Poly(f, U).as_expr(): m for f, m in fl[1] if sp.Poly(f, U).degree() > 0}


if __name__ == "__main__":
    sq, st = sp.sqrt(q), sp.sqrt(t)
    Kside = {"DN": lambda lam: koornwinder_eval(lam, (1, -1, sq, -sq)),
             "LRW": lambda lam: koornwinder_eval(lam, (1, -1, st, -st))}
    Q = U**2/t
    Achoices = {"Q": Q, "Qt": Q*t, "Q/t": Q/t, "Qt^2": Q*t**2, "U": U, "U^2": U**2}
    point = {q: sp.Rational(3, 7), t: sp.Rational(5, 11)}
    for wkind in ("1", "m"):
        for kname, K in Kside.items():
            Z1 = K((1,))*weight((1,), wkind)*Pbasis((1,))
            Z2 = sum(K(lam)*weight(lam, wkind)*Pbasis(lam) for lam in partitions(2))
            FG = components(F2(Z1, Z2))
            print(f"=== w = {wkind}, Kauffman = {kname}")
            for comp, e in FG.items():
                print(f"   F_G[{comp}] の U 極: {poles(e, point)}")
            for aname, Aval in Achoices.items():
                H = lambda lam, mu: dim_composite(lam, mu).subs(A, Aval)
                R1 = (H((1,), ()) + H((), (1,)))*weight((1,), wkind)*Pbasis((1,))
                R2 = sum((H(lam, ()) + H((), lam))*weight(lam, wkind)*Pbasis(lam) for lam in partitions(2))
                R2 += H((1,), (1,))*weight((1,), wkind)**2*Pbasis((1,))**2
                FR = components(F2(R1, R2))
                print(f"   A = {aname}:")
                for comp in FG:
                    print(f"      F_R[{comp}] の U 極: {poles(FR[comp], point)}")
