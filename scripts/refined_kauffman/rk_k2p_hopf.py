r"""K2'：t-brane（双対 Cauchy、Mimachi）型の refined orientifold OV 核での Hopf 絡み目。

ホロノミー U ∈ SO(2n)（AS の D_n = Koornwinder (±1, ±q^{1/2})）、probe の数 m = 1（変数 y, z）。
Mimachi の恒等式（役割を入れ替えた読み）：
  Σ_{k=0}^{n} (−1)^k K_{(n−k)}(y; t,q; t̄) K_{1^k}(x; q,t; t̄) = y^{−n} Π_i (1 − y x_i)(1 − y/x_i)
なので、probe 側の基底 B_k(y) := K_{(n−k)}(y;t,q;t̄)（1 変数 Askey–Wilson）とホロノミー側の標準表現 1^k が対になる。
Hopf 絡み目（AS：S_{αβ}/S_{00} = K_α(t^ρ) K_β(t^ρ q^α)）：
  Z(y,z) = Σ_{k,l} (−1)^{k+l} K_{1^k}(t^ρ) K_{1^l}(t^ρ q^{1^k}) B_k(y) B_l(z)
         = Σ_k (−1)^k K_{1^k}(t^ρ) B_k(y) · z^{−n} Π_i (1 − z x_i)(1 − z/x_i)|_{x = t^ρ q^{1^k}}.
(a) 直接の二重和と Mimachi で縮約した式が一致するか（n = 2, 3）。
(b) log Z の y z, y z², y² z, y² z² の係数を n = 2..7 で数値計算し、U = t^n の有理式に合うかを見る。
q, t は有理数、y, z は記号。
"""
import sys
sys.path.insert(0, "src")
import sympy as sp
from rcsknots.refined_kauffman.koornwinder import koornwinder
from rcsknots.refined_kauffman.kauffman_refined import koornwinder_eval, q as qs, t as ts, U as Us

y, z = sp.symbols("y z")


def setup(nval, q, t, sq):
    tbar = (1, -1, sq, -sq)
    rho = [t**(nval - i) for i in range(1, nval + 1)]
    ys = (y,)
    B = {k: koornwinder((nval - k,), 1, t, q, tbar, ys) if nval - k > 0 else sp.Integer(1) for k in range(nval + 1)}
    W = {k: koornwinder_eval((1,)*k if k else (), tbar).subs({qs: q, ts: t, Us: t**nval}) for k in range(nval + 1)}
    return tbar, rho, B, W


def Z_reduced(nval, q, t, sq):
    tbar, rho, B, W = setup(nval, q, t, sq)
    Z = 0
    for k in range(nval + 1):
        x = [rho[i] * (q if i < k else 1) for i in range(nval)]
        Z += (-1)**k * W[k] * B[k] * z**(-nval) * sp.prod([(1 - z*xi)*(1 - z/xi) for xi in x])
    return sp.expand(Z)


def Z_direct(nval, q, t, sq):
    tbar, rho, B, W = setup(nval, q, t, sq)
    xs = sp.symbols(f"x1:{nval+1}")
    Kpoly = {l: koornwinder((1,)*l if l else (), nval, q, t, tbar, xs) for l in range(nval + 1)}
    Z = 0
    for k in range(nval + 1):
        x = {xs[i]: rho[i] * (q if i < k else 1) for i in range(nval)}
        for l in range(nval + 1):
            Z += (-1)**(k + l) * W[k] * Kpoly[l].subs(x) * B[k] * B[l].subs(y, z)
    return sp.expand(Z)


def log_coeffs(Z, maxdeg=2):
    """Z を y, z の Laurent 級数（|q|<|y|,|z|<1/|q| の環）でなく、Z が多項式になるよう y^n z^n を掛けてから log を展開する。
    ここでは Z は y, z の Laurent 多項式なので、Z0 = Z(y=0,z=0) で割って log を級数展開（y, z の正冪の部分だけを見る）。"""
    e = sp.Symbol("e")
    Zs = sp.expand(Z.subs({y: e*y, z: e*z}))
    # y, z の負冪は e の負冪になる：負冪の項は別の級数展開になるので、ここでは y^n z^n を掛けて多項式にする
    return None


if __name__ == "__main__":
    q = sp.Rational(9, 49); sq = sp.Rational(3, 7); t = sp.Rational(5, 11)
    for nval in (2, 3):
        d = sp.simplify(Z_reduced(nval, q, t, sq) - Z_direct(nval, q, t, sq))
        print(f"(a) n={nval}: 直接 − 縮約 =", d)
    # (b) y^n z^n Z は y, z の多項式。その log の低次の係数（Z の y^0 z^0 の項で正規化）
    for nval in range(2, 8):
        Z = sp.expand(Z_reduced(nval, q, t, sq) * y**nval * z**nval)
        Z0 = Z.subs({y: 0, z: 0})
        P = sp.Poly(sp.expand(Z / Z0), y, z)
        c = lambda a, b: P.coeff_monomial(y**a * z**b)
        # log(1 + X) の係数：y z: c11 − c10 c01；y^2 z: c21 − c20 c01 − c11 c10 + c10^2 c01；…
        L11 = c(1, 1) - c(1, 0)*c(0, 1)
        L21 = c(2, 1) - c(2, 0)*c(0, 1) - c(1, 1)*c(1, 0) + c(1, 0)**2*c(0, 1)
        print(f"(b) n={nval}: U={t**nval}  L11={L11}  L21={L21}  c10={c(1,0)}  c01={c(0,1)}")
