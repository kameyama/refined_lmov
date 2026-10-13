"""open_problems D1：BGLX（arXiv:1404.4616）の Q_{m,n} と GN（1304.3328）式 (fix3) の演算子の対応を、演算子として（1 以外にも作用させて）確かめる。

記号（世界は GN の (q, τ)、gn.py では τ = t_w = 1/t）：
  F_{n,m}：(fix3) を H̃ 基底の係数として文字どおりに読んだ演算子（GN Thm fixed の I_λ = H̃_λ）。
  O_{n,m} = ∇^{-1} F_{n,m} ∇：gn.py の読み方（I_μ = H̃_μ/∇_μ）。gn.py の op は O/c_k（k = 1 では c_1 = 1 − 1/τ）。
照合（有理数の点で厳密、入力は次数 ≤ 3 の H̃_λ すべて）：
  P1  O_{1,0} = e_1 の掛け算（GN の Remark :889「(fix3) は Pieri 則」）。F_{1,0} は e_1 の掛け算ではない。
  P2  ∇ O_{n,m} ∇^{-1} = O_{n,m+n}（(fix3) で S_{m/n}(i) が 1 増えると各項が Π χ_i = T_μ/T_λ 倍になることから）。
  P3  正の次数の擬空三角形（GN (commutation) :631 の向き、m1 n2 − m2 n1 = 1）で [O_{n1,m1}, O_{n2,m2}] = −M O_{n1+n2,m1+m2}（M = (1−q)(1−τ)）。
  P4  Q_{m,n} = (−1)^n O_{n,m}（演算子として）。
使い方：PYTHONPATH=src:scripts/framed_knot:scripts/dgr_lmov python scripts/framed_knot/d1_ops.py（ログ logs/framed_knot/d1_ops.log）。
"""
from fractions import Fraction as Fr

from rcsknots.ehall import EHA
from rcsknots.modified import modified_H, solve_in_basis
from rcsknots.partitions import partitions, n_stat, conjugate
from d1_bglx import Qop, add, mul

Q0, T0 = Fr(3, 7), Fr(5, 11)


class World:
    def __init__(self, q, t):
        self.q, self.t = q, t
        self.E = EHA(q, t, Fr(1))
        self.Hs = {}

    def H(self, d):
        if d not in self.Hs:
            self.Hs[d] = modified_H(d, self.q, self.t, lambda x: x) if d else {(): {(): Fr(1)}}
        return self.Hs[d]

    def T(self, mu):
        return self.q ** n_stat(conjugate(mu)) * self.t ** n_stat(mu)

    def to_H(self, f):
        """p 基底の斉次式 → H̃ 基底の係数。"""
        if not f:
            return {}
        d = sum(next(iter(f)))
        return solve_in_basis(f, self.H(d), Fr(0)) if d else {(): f[()]}

    def to_p(self, v):
        out = {}
        for mu, c in v.items():
            out = add(out, {r: c * x for r, x in self.H(sum(mu))[mu].items()})
        return out

    def F(self, n, m, v):
        return self.E.apply(1, n, m, v)

    def O(self, n, m, v):
        w = self.E.apply(1, n, m, {mu: c * self.T(mu) for mu, c in v.items()})
        return {mu: c / self.T(mu) for mu, c in w.items()}

    def nabla(self, v, s=1):
        return {mu: c * self.T(mu) ** s for mu, c in v.items()}


def clean(v):
    return {k: c for k, c in v.items() if c != 0}


def main():
    log = open("logs/framed_knot/d1_ops.log", "w")

    def out(s):
        print(s, flush=True)
        log.write(s + "\n")

    W = World(Q0, T0)
    M = (1 - Q0) * (1 - T0)
    basis = [lam for d in range(0, 4) for lam in (partitions(d) if d else [()])]
    e1 = {(1,): Fr(1)}

    p1o = p1f = True
    for lam in basis:
        v = {lam: Fr(1)}
        target = W.to_H(mul(W.to_p(v), e1))
        p1o &= clean(W.O(1, 0, v)) == clean(target)
        p1f &= clean(W.F(1, 0, v)) == clean(target)
    out(f"P1：O_(1,0) = e_1 の掛け算 {p1o}、F_(1,0) = e_1 の掛け算 {p1f}")

    p2 = True
    for n, m in [(1, 0), (1, 1), (2, 1), (2, -1), (3, 1), (3, -2), (2, 3)]:
        for lam in basis[:5]:
            v = {lam: Fr(1)}
            p2 &= clean(W.nabla(W.O(n, m, W.nabla(v, -1)))) == clean(W.O(n, m + n, v))
    out(f"P2：∇ O_(n,m) ∇^(-1) = O_(n,m+n)：{p2}")

    p3 = []
    for x, y in [((1, 1), (1, 0)), ((1, 2), (1, 1)), ((2, 1), (1, 0)), ((1, 0), (1, -1)), ((1, 1), (2, 1)),
                 ((2, 3), (1, 1)), ((1, -1), (2, -3)), ((1, 2), (2, 3))]:
        (n1, m1), (n2, m2) = x, y
        assert m1 * n2 - m2 * n1 == 1
        ok = True
        for lam in basis[:5]:
            v = {lam: Fr(1)}
            lhs = add(W.O(n1, m1, W.O(n2, m2, v)), W.O(n2, m2, W.O(n1, m1, v)), -1)
            rhs = {k: -M * c for k, c in W.O(n1 + n2, m1 + m2, v).items()}
            ok &= clean(add(lhs, rhs, -1)) == {}
        p3.append((x, y, ok))
    out("P3：[O_x, O_y] = −M O_(x+y)：" + "、".join(f"{x},{y} {ok}" for x, y, ok in p3))

    p4 = []
    for m, n in [(0, 1), (1, 1), (2, 1), (3, 1), (1, 2), (3, 2), (1, 3), (2, 3), (4, 3), (3, 4)]:
        ok = True
        for lam in basis[:4]:
            f = W.to_p({lam: Fr(1)})
            if m == 0:
                Qf = mul(f, {(1,): Fr(-1)})          # Q_{0,1} = −e_1（spl(1,n) の (0,1)）
            else:
                Qf = Qop(m, n, f, Q0, T0)
            ok &= clean(W.to_H(Qf)) == clean({k: (-1) ** n * c for k, c in W.O(n, m, {lam: Fr(1)}).items()})
        p4.append(((m, n), ok))
    out("P4：Q_(m,n) = (−1)^n O_(n,m)：" + "、".join(f"{mn} {ok}" for mn, ok in p4))


if __name__ == "__main__":
    main()
