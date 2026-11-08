"""T1（続き 107）：unrefined（t = q）の随伴（composite [1;1̄]）のトーラス結び目の不変量を、大きな N の Rosso–Jones 型の公式で a を記号のまま計算し、
DAHA の R_adj|_{t=q} と照合する。目的は α₀ = ρ₀²・g₁ の関係・最深行の由来を公式の項から読むこと。
  ψ^r(s_{[1;1]}) = ψ^r(s_1 s̄_1) − 1 = Σ_{μ,ν：r のフック} (−1)^{k_μ+k_ν} s_μ s̄_ν − 1、
  s_μ s̄_ν = Σ_κ Σ_{α,β} c^μ_{κα} c^ν_{κβ} s_{[α;β]}（Littlewood–Koike）、
  W_{[1;1]}(T(r,s)) = θ_{[1;1]}^{−rs} Σ_{[α;β]} c θ_{[α;β]}^{s/r} dim_q[α;β]、θ_{[α;β]} = q^{e0} a^{e1}（U(N) の Casimir κ + N|D|/2 を N の一次式に当てる）。
次元は kauffman_refined.dim_composite（対称な正規化、A = a）を t = q で使う。
実行：PYTHONPATH=src:scripts/refined_kauffman python scripts/refined_kauffman/rk_t1_composite_rj.py"""
import sys, pickle
from fractions import Fraction
from functools import lru_cache
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rcsknots.characters import chi
from rcsknots.partitions import partitions, z
from rcsknots.refined_kauffman import kauffman_refined as KR

q, a = sp.symbols("q a", positive=True)
r_, s_ = sp.symbols("r s", positive=True)


def hooks(n):
    return [((n - k,) + (1,) * k, (-1)**k) for k in range(n)]


@lru_cache(maxsize=None)
def lr(mu, kappa, alpha):
    """Littlewood–Richardson 係数 c^μ_{κα} = Σ_{σ,τ} χ_κ(σ) χ_α(τ) χ_μ(σ∪τ) / (z_σ z_τ)。"""
    if sum(kappa) + sum(alpha) != sum(mu):
        return 0
    tot = Fraction(0)
    for sg in partitions(sum(kappa)) if kappa else [()]:
        for ta in partitions(sum(alpha)) if alpha else [()]:
            rho = tuple(sorted(sg + ta, reverse=True))
            cs, ct = (chi(kappa, sg) if kappa else 1), (chi(alpha, ta) if alpha else 1)
            tot += Fraction(cs * ct * chi(mu, rho), z(sg) * z(ta)) if (sg or ta) else Fraction(cs * ct * chi(mu, rho))
    assert tot.denominator == 1
    return int(tot)


def subparts(mu):
    """μ の部分分割（κ ⊆ μ）。"""
    n = sum(mu)
    return [k for m in range(n + 1) for k in (partitions(m) if m else [()]) if len(k) <= len(mu) and all(k[i] <= mu[i] for i in range(len(k)))]


def koike(mu, nu):
    """s_μ(x) s_ν(x̄) = Σ c_{[α;β]} s_{[α;β]}。"""
    out = {}
    for kappa in subparts(mu):
        if not (len(kappa) <= len(nu) and all(kappa[i] <= nu[i] for i in range(len(kappa)))):
            continue
        for al in (partitions(sum(mu) - sum(kappa)) if sum(mu) > sum(kappa) else [()]):
            ca = lr(mu, kappa, al)
            if not ca:
                continue
            for be in (partitions(sum(nu) - sum(kappa)) if sum(nu) > sum(kappa) else [()]):
                cb = lr(nu, kappa, be)
                if cb:
                    out[(al, be)] = out.get((al, be), 0) + ca * cb
    return out


def psi_adjoint(r):
    """ψ^r(s_{[1;1]}) の composite 指標への展開 {(α, β): 係数}。"""
    out = {}
    for mu, smu in hooks(r):
        for nu, snu in hooks(r):
            for ab, c in koike(mu, nu).items():
                out[ab] = out.get(ab, 0) + smu * snu * c
    out[((), ())] = out.get(((), ()), 0) - 1
    return {k: v for k, v in out.items() if v}


def casimir_exponent(al, be):
    """θ_{[α;β]} = q^{e0} a^{e1}：U(N) = SU(N) × U(1) の Casimir
    κ(D_N) + N|D_N|/2 − |D_N|²/(2N) + n²/(2N)（n = |α| − |β| は U(1) 電荷）= e0 + e1 N。
    [1] と [∅;1] で a^{1/2}（同じ枠の因子）、随伴 [1;1] で a。"""
    n_ = sum(al) - sum(be)
    def E(N):
        D = KR.composite_diagram(list(al), list(be), N)
        kappa = sum(j - i for i, row in enumerate(D) for j in range(row))
        d = sum(D)
        return Fraction(kappa) + Fraction(N * d, 2) - Fraction(d * d, 2 * N) + Fraction(n_ * n_, 2 * N)
    N0 = len(al) + len(be) + 3
    e0, e1, e2 = E(N0), E(N0 + 1), E(N0 + 2)
    assert e2 - e1 == e1 - e0, "Casimir が N の一次式でない"
    c1 = e1 - e0
    return e0 - c1 * N0, c1


def dim_tq(al, be):
    """composite [α;β] の次元（対称な正規化）を t = q、A = a で。"""
    if not al and not be:
        return sp.Integer(1)
    d = KR.dim_composite(list(al), list(be))
    return sp.cancel(d.subs({KR.t: q, KR.q: q, KR.A: a}))


def W_torus(ps, r, s):
    """W_V(T(r,s)) = θ_V^{−rs} Σ c θ^{s/r} dim（ps = ψ^r(χ_V) の展開、V の θ は ps を与える前に別途）。"""
    tot = 0
    for (al, be), c in ps.items():
        e0, e1 = casimir_exponent(al, be)
        tot += c * q**(sp.Rational(e0) * sp.Rational(s, r)) * a**(sp.Rational(e1) * sp.Rational(s, r)) * dim_tq(al, be)
    return tot


def adjoint_rj(r, s):
    ps = psi_adjoint(r)
    e0, e1 = casimir_exponent((1,), (1,))
    W = q**(-e0 * r * s) * a**(-e1 * r * s) * W_torus(ps, r, s)
    return sp.cancel(sp.together(W / dim_tq((1,), (1,)))), ps


def fund_rj(r, s):
    ps = {(mu, ()): sg for mu, sg in hooks(r)}
    e0, e1 = casimir_exponent((1,), ())
    W = q**(-e0 * r * s) * a**(-e1 * r * s) * W_torus(ps, r, s)
    return sp.cancel(sp.together(W / dim_tq((1,), ())))


def rows_a(expr):
    """a の冪ごとの行（Laurent 多項式 in a^{1/2}）。"""
    e = sp.expand(expr)
    nu = sp.Symbol("nu", positive=True)
    e = sp.expand(sp.powdenest(e.subs(a, nu**2), force=True))
    P = sp.Poly(sp.expand(e * nu**400), nu)
    return {int(m[0]) - 400: sp.factor(c) for m, c in zip(P.monoms(), P.coeffs()) if c != 0}


if __name__ == "__main__":
    for r, s in ((2, 3), (2, 5), (3, 4)):
        print(f"== T({r},{s})")
        F = fund_rj(r, s)
        print("  基本表現（Rosso–Jones、t = q、reduced）の a の行：", rows_a(F))
        Rad, ps = adjoint_rj(r, s)
        print(f"  ψ^{r}(s_[1;1]) の composite の項：{len(ps)} 個：", ps)
        print("  随伴（Rosso–Jones、t = q、reduced）の a の行：")
        for k, v in sorted(rows_a(Rad).items(), reverse=True):
            print(f"    a^{{{k}/2}}: {v}"[:200])
