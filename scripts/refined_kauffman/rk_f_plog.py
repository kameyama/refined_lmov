"""n = 1 の母関数の plethystic logarithm（theory §12.31、progress.md 続き 101）。
Fock 空間の解釈（有限個の生成元）があるかを、X_{j′}・F_{j′}・上の系列 u_j = Sym^j V + (q − t) t^{m−1−j′} F_{j′}・E の上の系列について調べる。
plethystic log：G(z) = Π (1 − q^a t^b z^n)^{−g_{n,a,b}} の g を、PL = Σ_k μ(k)/k log G(z^k; q^k, t^k) で求める。

実行：PYTHONPATH=src:scripts/refined_kauffman .venv/bin/python scripts/refined_kauffman/rk_f_plog.py [m ...]
"""
import sys
import sympy as sp

q, t, z = sp.symbols("q t z")
V = 1 + t + q * t
NUM_F = 1 + (2 * t + t**2 + q * t + 2 * q * t**2) * z + q * t**3 * z**2
DEN_F = (1 - z) * (1 - q * t * z) * (1 - t**2 * z) * (1 - q**2 * t**2 * z)


def series_coeffs(G, order):
    s = sp.series(G, z, 0, order + 1).removeO()
    return [sp.expand(s.coeff(z, n)) for n in range(order + 1)]


def plog(g):
    """g_0 = 1 の z 展開の係数列から plethystic log の係数 PL_1..PL_order（q, t の多項式）。"""
    assert g[0] == 1
    order = len(g) - 1
    L = [0] * (order + 1)
    for n in range(1, order + 1):
        L[n] = sp.expand((n * g[n] - sum(k * L[k] * g[n - k] for k in range(1, n))) / n)
    return [0] + [sp.expand(sum(sp.Rational(sp.mobius(d), d) * L[n // d].subs({q: q**d, t: t**d}, simultaneous=True)
                                for d in sp.divisors(n))) for n in range(1, order + 1)]


def gen_X():
    return (1 + t * z) / ((1 - z) * (1 - t * z) * (1 - q * t * z))


def gen_F():
    return NUM_F / DEN_F


def gen_upper_N(m):
    """上の系列 u_j = Sym^j V + (q − t) t^{m−1−j′} F_{j′}（j = m + j′）の母関数。"""
    return 1 / ((1 - z) * (1 - t * z) * (1 - q * t * z)) + (q - t) * t**(m - 1) * z**m * gen_F().subs(z, z / t)


def gen_upper_E(m):
    """E の上の系列（符号 −t^{K−2j} を除き、全体に −1 を掛けて定数項 1 にしたもの）：
    Σ_cells (−1)^{b+c} m − (q − t)(−1)^{m+j′+1} t^{m−1−j′} ψ²X_{j′}。"""
    X2 = gen_X().subs({q: q**2, t: t**2}, simultaneous=True).subs(z, -z / t)
    return 1 / ((1 - z) * (1 + t * z) * (1 + q * t * z)) - (q - t) * (-1)**(m + 1) * t**(m - 1) * z**m * X2


def generators_N(m, jmax):
    """上の系列の plethystic log の z^{m+j′} の係数を (q − t) t^{m−1−j′} で割った c_{j′}（j′ ≤ jmax < m）。"""
    PL = plog(series_coeffs(gen_upper_N(m), m + jmax))
    return [sp.expand(sp.cancel(PL[m + jp] / ((q - t) * t**(m - 1 - jp)))) for jp in range(jmax + 1)]


def _show(name, PL):
    print(f"== {name}")
    for n in range(1, len(PL)):
        if PL[n] != 0:
            print(f"  z^{n}: {sp.factor(PL[n])}")


if __name__ == "__main__":
    ms = [int(a) for a in sys.argv[1:]] or [5]
    _show("X（有限：V ⊕ t の Fock に関係 t²z²）", plog(series_coeffs(gen_X(), 6)))
    _show("F", plog(series_coeffs(gen_F(), 6)))
    for m in ms:
        _show(f"N の上の系列 m={m}", plog(series_coeffs(gen_upper_N(m), m + 3)))
        for jp, c in enumerate(generators_N(m, min(m - 1, 5))):
            cs = sp.Poly(c, q, t).coeffs()
            print(f"   c_{jp}：Σ|c| = {sum(abs(v) for v in cs)}、負の係数 {sum(1 for v in cs if v < 0)}")
        _show(f"E の上の系列 m={m}", plog(series_coeffs(gen_upper_E(m), m + 3)))
