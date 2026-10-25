r"""(b)：N を固定した refined SO(2N) の trefoil（1^k、k ≤ N）で、自由エネルギー f_a（a ≤ amax）がすべて Laurent になる
framing の単項式 m_k = r^{α_k} s^{β_k} を探す（r = q^{1/2}、s = t^{1/2}）。
  W_{1^k} = Σ m · R · K（R：data/refined_kauffman/torus の reduced 不変量の鏡像（AS の向き）、K：unknot の主特殊化）。
  SO(4) = SU(2)²：W_1 = m_1 R_[1]² K_1、W_2 = m_2 R_[2] · 2P_[2]。
  SO(6) = SU(4)：W_1 = m_1 R_[1,1] K_1、W_2 = m_2 R_[2,1,1] K_{11}、W_3 = m_3 R_[2] P_[2] + m_3' R_[2,2,2] P_[2,2,2]。
方法：s を有理数 s0 に固定し、log Z~ の係数 L_b を M（= m）の多項式・r の有理式として求める。Möbius 反転の d 乗の項は
  s0 → s0^d、r → r^d、M_k → r^{dα_k} s0^{dβ_k}。f_a = Σ_i C_i(r) r^{e_i(α)} s0^{g_i(β)} と書き、分母の非単項式部分 D を法とする剰余が 0 か
  （⇔ f_a が r の Laurent 多項式）を判定する。s0 を 2 つ取って両方で成り立つ組を残す。"""
import sys, json, itertools, time
from functools import reduce
sys.path.insert(0, "src")
import sympy as sp
from sympy import mobius
from rcsknots.refined_kauffman.askey_wilson import askey_wilson_monic
from rcsknots.refined_kauffman.kauffman_knot import unknot_W, r, s
from rcsknots.macdonald import principal_specialization


def Rdata(c, Nsu):
    d = json.load(open(f"data/refined_kauffman/torus/T2_3_{c}.json"))
    A = s**Nsu
    return sp.expand(sum(v * A**(-int(k.split(",")[0])) * r**(-int(k.split(",")[1])) * s**(-int(k.split(",")[2])) for k, v in d["poly"].items()))


def P_sym(lam, Nsu):
    k = sum(lam)
    return sp.cancel(s**k / s**(Nsu * k) * principal_specialization(lam, s**(2 * Nsu), r**2, s**2))


def setup(N):
    if N == 2:
        return {1: [Rdata("1", 2)**2 * unknot_W(2, 1)], 2: [Rdata("2", 2) * 2 * P_sym((2,), 2)]}
    if N == 3:
        return {1: [Rdata("11", 4) * unknot_W(3, 1)], 2: [Rdata("211", 4) * unknot_W(3, 2)],
                3: [Rdata("2", 4) * P_sym((2,), 4), Rdata("222", 4) * P_sym((2, 2, 2), 4)]}
    raise ValueError(N)


def L_terms(N, parts, s_val, bmax):
    """{b: {指数ベクトル n: C(r)}}：L_b = Σ_n C_n(r) Π M^n（s = s_val）。有理関数体 Q(r, M) の上で計算する。"""
    from sympy import QQ
    from sympy.polys.fields import field
    nM = sum(len(v) for v in parts.values())
    K, *gens = field(",".join(["r"] + [f"M{i}" for i in range(nM)]), QQ)
    rK, Ms = gens[0], gens[1:]
    plain_r = sp.Symbol("r")
    conv = lambda e: K.from_expr(sp.cancel(e).subs(r, plain_r)) if e != 0 else K.zero
    W, i = [K.one], 0
    for k in range(1, N + 1):
        e = K.zero
        for part in parts[k]:
            e += Ms[i] * conv(part.subs(s, s_val)); i += 1
        W.append(e)
    aw = askey_wilson_monic(N, s_val**2, sp.Integer(1), sp.Integer(-1), r, -r)
    Z = {}
    for k in range(N + 1):
        n = N - k
        for j in range(n + 1):
            c = aw[n].get(n - 2 * j, 0)
            if c:
                Z[k + 2 * j] = Z.get(k + 2 * j, K.zero) + (-1)**k * W[k] * conv(c)
    X = {a: v for a, v in Z.items() if 0 < a <= bmax}
    L, P = {}, dict(X)
    for m in range(1, bmax + 1):
        for a, v in P.items():
            L[a] = L.get(a, K.zero) + K(sp.Rational((-1)**(m + 1), m)) * v
        Q = {}
        for a1, v1 in P.items():
            for a2, v2 in X.items():
                if a1 + a2 <= bmax:
                    Q[a1 + a2] = Q.get(a1 + a2, K.zero) + v1 * v2
        P = Q
    out = {}
    Rr = sp.Poly(1, plain_r)
    for b, v in L.items():
        # 分母は M を含まない（W は M に線形、Askey–Wilson の係数は M を含まない）
        num, den = v.numer, v.denom
        den_e = den.as_expr().subs(plain_r, r)
        terms = {}
        for monom, c in num.terms():
            key = tuple(monom[1:])
            terms[key] = terms.get(key, 0) + c * r**monom[0]
        out[b] = {key: sp.cancel(c / den_e) for key, c in terms.items()}
    return out, nM


def f_terms(N, parts, s0, amax):
    """{a: [(C(r), n, d)]}：f_a = Σ (μ(d)/d) C(r^d)|_{s=s0^d} · Π (r^{dα} s0^{dβ})^n。"""
    Ld = {}
    for d in range(1, amax + 1):
        Ld[d], nM = L_terms(N, parts, s0**d, amax // d)
    out = {}
    for a in range(1, amax + 1):
        terms = []
        for d in sp.divisors(a):
            for mon, C in Ld[d].get(a // d, {}).items():
                terms.append((sp.Rational(mobius(d), d) * C.subs(r, r**d), mon, d))
        out[a] = terms
    return out, nM


def _strip_r(pd):
    """pd = r^k · pd'（pd'(0) ≠ 0）の (k, pd')。"""
    k = 0
    while pd.eval(0) == 0:
        pd = sp.Poly(sp.cancel(pd.as_expr() / r), r); k += 1
    return k, pd


def prepare(terms):
    """f_a·D ≡ Σ A_i r^{e_i − k_i} s0^{g_i} (mod D) の A_i, k_i を作る（D：分母の非単項式部分の最小公倍数）。"""
    raw = []
    for C, mon, d in terms:
        num, den = sp.fraction(sp.cancel(C))
        k, pd = _strip_r(sp.Poly(den, r))
        raw.append((sp.Poly(num, r), k, pd, mon, d))
    D = reduce(lambda x, y: x.lcm(y), [x[2] for x in raw], sp.Poly(1, r))
    prepped = [((pn * sp.Poly(sp.cancel(D.as_expr() / pd.as_expr()), r)).rem(D), k, mon, d) for pn, k, pd, mon, d in raw]
    return D, prepped


def _powmod(b, e, D):
    out, base = sp.Poly(1, r), b.rem(D)
    while e:
        if e & 1:
            out = (out * base).rem(D)
        base = (base * base).rem(D)
        e >>= 1
    return out


def laurent_ok(D, prepped, alpha, beta, s0):
    """(α, β) で f_a が Laurent か（D を法として f_a·D ≡ 0）。"""
    if D.degree() == 0:
        return True
    rP = sp.Poly(r, r)
    rinv = sp.Poly(sp.invert(r, D.as_expr(), r), r)
    tot = sp.Poly(0, r)
    for A, k, mon, d in prepped:
        e = d * sum(n * a for n, a in zip(mon, alpha)) - k
        g = d * sum(n * b for n, b in zip(mon, beta))
        tot = (tot + A * _powmod(rP if e >= 0 else rinv, abs(e), D) * s0**g).rem(D)
    return tot.is_zero
