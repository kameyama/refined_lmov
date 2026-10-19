r"""(2) E を 1 箱のクロスキャップの Adams 倍化で書く候補の走査（progress.md 続き 35）。
  候補 E = F · (ψ²[g_(1)] + 1)、ψ²[g](q, t, ν) = g(q², t², ν² (q/t)^{j/2}) と、前因子 F ∈ {−(1+q), −(1+t), −(1+q) (t/q)^{k/2}, −(q^{1/2} t^{1/2} + q^{−1/2} t^{−1/2}) (qt)^{1/2}, …}。
  t = q で AgSc の強い整数性の形 E = −(1 + q)(g_(1)(x²) + 1) に戻る。条件 (i)：C = g − E/(2(1+t)) が整数 Laurent（(1 + t) | (2(1+t)g_(2) − E)）。
使い方：python scripts/rk_b_knot_adams_E.py m"""
import sys, itertools
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp
import rk_b_knot_fixed_sector as K

r, s, nu = K.r, K.s, K.nu
t, q = s**2, r**2


def scan(m):
    g = K.load_g(m)
    R = sp.cancel(2 * (1 + t) * g[(2,)])
    F = {"−(1+q)": -(1 + q), "−(1+t)": -(1 + t), "−(q+t)": -(q + t), "−(1+qt)": -(1 + q * t),
         "−(r s + 1/(r s)) r s": -(r * s + 1 / (r * s)) * r * s, "−(1+q)·t/q": -(1 + q) * t / q, "−(1+q)·(t/q)^{1/2}": -(1 + q) * s / r,
         "−(1+t)·q/t": -(1 + t) * q / t, "−(1+q)·q/t": -(1 + q) * q / t}
    hits = []
    for (fl, f), j, sg in itertools.product(F.items(), range(-2, 3), (1, -1)):
        psi = g[(1,)].subs({r: r**2, s: s**2, nu: nu**2 * (r / s)**j}, simultaneous=True)
        E = sp.cancel(sg * f * (psi + 1))
        ok_i = sp.Poly(sp.fraction(sp.cancel((R - E) / (1 + t)))[1], r, s, nu).is_monomial
        if ok_i:
            hits.append((fl, j, sg))
    print(f"T(2,{2*m+1})：条件 (i) を満たす候補 (F, j, 符号)：{hits}", flush=True)
    return hits, g


if __name__ == "__main__":
    scan(int(sys.argv[1]))
