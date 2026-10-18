"""refined Kauffman K3'：trefoil [1]（AS 1202.4456 §6）と refined の h/g 分解。"""
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rcsknots.kauffman_trefoil import (Z_V_T2, Z615, r, s, U, marino_F31, g_ref, is_laurent, W_SO, h_half_ref)  # noqa: E402


def test_general_formula_reproduces_615():
    assert sp.simplify(Z_V_T2(1) - Z615) == 0


def test_unrefined_is_mirror_of_marino_with_framing():
    """t = q：AS (6.15) = F_{3_1}(ν → 1/ν, z → −z) × ν^{−6}、ν = 1/𝐚。"""
    Zu = Z615.subs(r, s)
    FM = marino_F31(s / U, 1 / s - s)
    assert sp.cancel(Zu / FM - U**6 / s**6) == 0


def test_refined_rudolph_with_identified_monomial():
    """g^{ref}_[1] = W_SO − ½h^{ref} は単項式 r^{−5} s^{11} U^{−6} で Laurent 多項式（11 項、整数係数）、
    T_[1]^{−6} = r^{−6} s^{12} U^{−6} では Laurent にならない。"""
    g = g_ref()
    assert is_laurent(g)
    num, den = sp.fraction(g)
    P = sp.Poly(num, r, s, U)
    assert len(P.terms()) == 11 and all(c.is_integer for c in P.coeffs())
    assert not is_laurent(g_ref(r**-6 * s**12 * U**-6))
    # unrefined でも Laurent
    assert is_laurent(g.subs(r, s))


def test_as_615_matches_gn_at_N2_N3():
    """SO(4) = SU(2)×SU(2)、SO(6) = SU(4)：AS (6.15) は、A 型の Gorsky–Negut（楕円 Hall 代数）で計算した
    SU(2) の [1] の 2 乗（N = 2）、SU(4) の [1,1]（N = 3）の正規化した trefoil に、単項式 q^{−1}、q^{−1}t^{−4} 倍で一致する。"""
    from fractions import Fraction as F
    from rcsknots.gn import GN
    for rr, ss in ((F(3, 7), F(5, 11)), (F(2, 3), F(7, 13))):
        qq, tt = rr**2, ss**2
        R = lambda x: sp.Rational(x.numerator, x.denominator)
        G2 = GN(qq, tt, tt**2, lambda x: x)
        v2 = (G2.torus((1,), 2, 3) / G2.torus((1,), 1, 0))**2
        assert R(v2) == Z615.subs({r: R(rr), s: R(ss), U: R(ss)**4}) * R(rr)**-2
        G4 = GN(qq, tt, tt**4, lambda x: x)
        v4 = G4.torus((1, 1), 2, 3) / G4.torus((1, 1), 1, 0)
        assert R(v4) == Z615.subs({r: R(rr), s: R(ss), U: R(ss)**6}) * R(rr)**-2 * R(ss)**-8


def test_trefoil_fixed_N_pole_cancellation():
    """(b)：N を固定した refined SO(2N) の trefoil（1^k、k ≤ N）。Z~ の N 依存の極が打ち消すには、reduced 不変量が極の上で単項式（倍）でなければならない。
    SO(4) = SU(2)²：R_[2] は q t = 1（r = ±1/s）で s^{−4}；m_2 = s^4 で f_2 は整数 Laurent、m_2 = 1 では Laurent でない。
    SO(6) = SU(4)：R_[2,1,1] は q t³ = 1（r = ±s^{−3}）で s^{−16}。R は data/refined_kauffman/torus の鏡像（AS の向き）。"""
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "scripts"))
    import rk_k5_framing_scan as Fs
    from rcsknots.kauffman_knot import free_energy, r as r_, s as s_
    for sign in (1, -1):
        assert sp.expand(Fs.Rdata("2", 2).subs(r_, sign / s_) - s_**-4) == 0
        assert sp.expand(Fs.Rdata("211", 4).subs(r_, sign * s_**-3) - s_**-16) == 0
    parts = Fs.setup(2)
    mono = lambda e: sp.Poly(sp.fraction(sp.cancel(e))[1], r_, s_).is_monomial
    for m2, expect in ((s_**4, True), (sp.Integer(1), False)):
        W = [sp.Integer(1), parts[1][0], sp.cancel(m2 * parts[2][0])]
        f = free_energy(2, W, 2)
        assert mono(f[2]) == expect


def test_torus_family_framing_monomial():
    """A：T(2,5) の [1]。g^{ref} を Laurent にし、t = q で Kauffman 型（r → −1/r 不変）にする単項式は
    r^{−8} s^{18} U^{−10} = T_[1]^{−10} (q/t)^{δ/2}（δ = 2、T_[1] = q^{1/2}t^{N−1} = r U s^{−2}）。AS 1105.5117 の SU(N) の因子と同じ形。
    t = q の g^{ref} は z の偶数冪だけ（N^{c=2} = 0）。"""
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "scripts"))
    import rk_k3_torus_family as Tf
    assert Tf.kauffman_q_power(2) == 0
    Tm = r * U / s**2
    assert sp.simplify(r**-8 * s**18 * U**-10 - Tm**-10 * (r / s)**2) == 0
    g = Tf.g_family(2, 18)
    assert is_laurent(g)
    # d の偶奇を変えると Laurent でない
    assert not is_laurent(sp.cancel(Tf.Z_V_T2(2) * r**-8 * s**19 * U**-10 * Tf.unknot_SO - Tf.h_half(2)))
    Q, Tb, A = sp.symbols("Q T A")
    gq = sp.expand(sp.cancel(g.subs({r: -Tb * Q, U: A * Q}).subs(s, Q)).subs(Tb, -1))
    assert sp.expand(gq - gq.subs(Q, -1 / Q)) == 0


def test_torus_family_refined_only_part():
    """A：g^{ref}_[1](T(2,2m+1)) の 𝐚 の奇数冪（t = q で消える refined の部分）は
    (−1)^{m+1}(𝐚²−𝐪²)(𝐭²−1) 𝐚^{−(4m+1)}𝐪^{−(2m−1)}𝐭^{−3m} Π_m、Π_m = Σ_{j<m} 𝐚^{2j}𝐪^{2(m−1−j)}𝐭^{2j} Σ_{i≤j}(𝐪⁴𝐭²)^i（係数はすべて +1）。m = 2, 3。"""
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "scripts"))
    import rk_k3_torus_family as Tf
    for m in (2, 3):
        Pi, (Q, Tb, A) = Tf.refined_only_part(m)
        guess = sum(A**(2 * j) * Q**(2 * (m - 1 - j)) * Tb**(2 * j) * sum((Q**4 * Tb**2)**i for i in range(j + 1)) for j in range(m))
        assert sp.expand(Pi - guess) == 0


def test_refined_crosscap_c1_blocks():
    """B：T(2,5) の [1] の c = 1 部分（G = z_q z_t、W = (q/t)^{1/4} の多項式）を ν^β ごとに分けると、
    unrefined の台のブロックは係数が同符号（+, −, +）、台の外のブロック（β = −8）は (𝐭² − 1) × 正の多項式。"""
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "scripts"))
    import rk_k3_torus_family as Tf
    assert Tf.c1_blocks(2) == {-10: "coherent+", -8: "tfactor+", -6: "coherent−", -4: "coherent+"}


def test_refined_crosscap_special_points():
    """B の因子の意味：T(2,5) の g^{ref}_[1] は SO(2)（𝐚 = ±𝐪、ν² = t）で (−𝐭)^m = (q/t)^{δ/2}（U(1) の自明さ）、
    𝐭 → −𝐭 で (−1)^m 倍（だから t = q で消える部分は (𝐭² − 1) を因子にもつ）。"""
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "scripts"))
    import rk_k3_torus_family as Tf
    so2, par, (Q, Tb, A) = Tf.special_points(2)
    assert so2 == {Tb**2} and par == 1


def test_trefoil_fixed_N_two_probe():
    """固定 N の 2 箱（probe 2 個、Y の Schur |R| ≤ 2）。m_(1,1) = (t/q)^{N−1}（極の上の値に t = q で 1 の条件）のとき、
    f_(1,1)（Sym²V、D_N の (2)）を Laurent にする m_(2) は (t/q)^2（scripts/rk_k5_trefoil_2box.py の走査で一意、t = q で自動的に 1）。
    N = 2 では m_(2) = m_(1,1)²（SO(4) = SU(2)² の積の形）。近くの単項式・m_(2) = 1 では Laurent でない。
    f は整数係数で、f_(1), f_(2) は 1 probe の計算（別の展開）の f_1, f_2 と一致。"""
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "scripts"))
    import rk_k5_trefoil_2box as T2
    import rk_k5_framing_scan as Fs
    from rcsknots.kauffman_knot import free_energy, r as r_, s as s_
    N = 2
    assert sp.expand(T2.M2[N] - T2.M11[N]**2) == 0
    fu = T2.free_energy_2(N, {lam: T2.unknot_K(N, lam) for lam in ((1,), (1, 1), (2,))})
    assert fu[(2,)] == 0 and fu[(1, 1)] == 0
    f = T2.free_energy_2(N, T2.trefoil_W(N, T2.M2[N]))
    for R in ((1,), (2,), (1, 1)):
        st = T2.stats(f[R])
        assert st["laurent"] and st["integer"]
    for m in (sp.Integer(1), T2.M2[N] * r_**2, T2.M2[N] / r_**2, T2.M2[N] * s_**2, T2.M2[N] / s_**2, -T2.M2[N]):
        assert not T2.laurent(T2.free_energy_2(N, T2.trefoil_W(N, m))[(1, 1)])
    parts = Fs.setup(N)
    W = [sp.Integer(1), sp.cancel(T2.M1 * parts[1][0]), sp.cancel(T2.M11[N] * parts[2][0])]
    f1p = free_energy(N, W, 2)
    assert sp.cancel(f1p[1] - f[(1,)]) == 0 and sp.cancel(f1p[2] - f[(2,)]) == 0


def test_trefoil_fixed_N_cover_and_hatted():
    """固定 N の被覆 U(3)（N = 2）とクロスキャップ部分 g = f − ½f^{cov}（scripts/rk_k5_trefoil_2box_cover.py）。
    unknot は g = (−1, 0, 0)。随伴の単項式 (t/q)^{3/2}（t = q で 1）で f^{cov} は Laurent、1 では Laurent でない。g_(1) は K3' の −g^{ref}_[1]。
    2g は整数で g_(2) + g_(1,1) は整数（半整数は p_2 方向だけ：Adams の予言 E_(2) = −E_(1,1)）。
    ĝ = M^{-1}(q^{-1}) g は m_(1) = t/q で z_q ĝ が Laurent、K3' の m_(1) = (t/q)^{1/2} では (1+q) の分母が残る。"""
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "scripts"))
    import rk_k5_trefoil_2box as T2
    import rk_k5_trefoil_2box_cover as C2
    from rcsknots.kauffman_knot import r as r_, s as s_
    from rcsknots.bps import M_matrix
    N = 2
    fu = C2.cover_f(N, unknot=True)
    gu = {R: sp.cancel(v - fu[R] / 2) for R, v in C2.so_f(N, unknot=True).items()}
    assert gu == {(1,): -1, (2,): 0, (1, 1): 0}
    assert not all(T2.laurent(v) for v in C2.cover_f(N).values())
    Minv = M_matrix(2, 1 / r_).inv()
    zq = r_ - 1 / r_
    for m1, hatted_ok in ((T2.M1_HAT, True), (T2.M1, False)):
        f, fc, g = C2.crosscap(N, m1)
        assert all(T2.laurent(v) for v in fc.values())
        assert C2.p2_divisible(fc, r_**2 + 1) and C2.p2_divisible(fc, s_**2 + 1)
        assert all(T2.stats(2 * v)["integer"] for v in g.values())
        assert T2.stats(g[(2,)] + g[(1, 1)])["integer"] and not T2.stats(g[(2,)])["integer"]
        gh = Minv * sp.Matrix([g[(2,)], g[(1, 1)]])
        assert all(T2.laurent(zq * gh[i]) for i in range(2)) == hatted_ok
        assert C2.p2_divisible(g, r_**2 + 1) == hatted_ok
    assert sp.cancel(C2.crosscap(N, T2.M1)[2][(1,)] + g_ref().subs({r: r_, s: s_, U: s_**4})) == 0



def test_flint_rf_basic():
    """rcsknots.flint_rf.RF：約分・Laurent 判定・整数判定・Adams・sympy との往復。"""
    from rcsknots.flint_rf import RF, R as R_, S as S_, from_sympy
    rr, ss = sp.symbols("r s", positive=True)
    e = from_sympy((rr**2 * ss - 1 / ss) / (1 + rr**2))
    assert not e.is_laurent() and (e * (1 + R_**2)).is_laurent()
    assert (e * (1 + R_**2)).is_integer_laurent() and not (e * (1 + R_**2) / 2).is_integer_laurent()
    assert sp.simplify(e.adams(2).to_sympy() - ((rr**4 * ss**2 - ss**-2) / (1 + rr**4))) == 0
    assert ((R_**2 - S_**2) / (R_ - S_) - (R_ + S_)).is_zero()


def test_trefoil_fixed_N_three_box():
    """固定 N = 2 の 3 箱（scripts/rk_k5_trefoil_3box.py、SO(4) = SU(2)² を型 A の経路で、被覆 U(3)）。
    型 A の経路は 2 箱で D_N の経路と一致し、unknot の f_R（|R| = 2, 3）は SO・被覆とも 0。
    被覆の (2;1), (11;1) の単項式は (t/q)^{3/2}（t = q で 1 の族で一意）。m_(1) = t/q のとき SO の [3]⊗[3], [3]⊗[1] は (t/q)^3, (t/q)^2
    （どれも (t/q)^{|μ|/2}）で、SO の f̂ = M^{-1}(q^{-1}) f（|R| ≤ 3）は整数 Laurent。K3' の m_(1) = (t/q)^{1/2} では整数 Laurent でない。
    g_R（|R| = 3）は整数 Laurent（奇数のセクターに E は無い）。被覆の f̂^{cov} は p_3 方向では Laurent、p_2 p_1・p_1³ 方向では Laurent でない。"""
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "scripts"))
    import rk_k5_trefoil_3box as X
    import rk_k5_trefoil_2box as T2
    import rk_k5_trefoil_2box_cover as C2
    from rcsknots.flint_rf import from_sympy
    N, W = 2, X.W_
    fu, fcu = X.plog(X.so_Z(N, None, unknot=True)), X.plog(X.cover_Z(N, None, unknot=True))
    assert all(v.is_zero() for k in (2, 3) for v in list(X.to_schur(fu, k).values()) + list(X.to_schur(fcu, k).values()))
    F2 = X.plog(X.so_Z(N, X.so_framing(N, from_sympy(T2.M1), (X.ZERO, X.ZERO)), kmax=2), kmax=2)
    old = C2.so_f(N)
    assert all((X.to_schur(F2, sum(R_))[R_] - from_sympy(old[R_])).is_zero() for R_ in old)
    assert X.scan_pairs(lambda ma, mb: X.cover_Z(N, {(1,): W**3, (2,): ma, (1, 1): mb}), rng=range(0, 7)) == [(3, 3)]
    assert X.so_new_framing(N, W**2) == [(6, 4)]
    for b1, ok in ((2, True), (1, False)):
        new = (W**(b1 + 4), W**(b1 + 2))
        F, Fc, G = X.crosscap3(N, W**b1, new)
        assert all(v.is_integer_laurent() for k in (2, 3) for v in X.hatted(F, k).values()) == ok
        assert all(v.is_integer_laurent() for v in X.to_schur(G, 3).values())
    assert X.hatted_by_rho(Fc, 3) == {(3,): True, (2, 1): False, (1, 1, 1): False}


def test_trefoil_uniform_normalization():
    """一様な規格化 W_μ = (t/q)^{|μ|/2} R_μ P_μ(t^ρ)（R は data の trefoil）。共役・行列式のずらしで不変（U(3), U(4) の例）。
    N = 2 でこの規格化を SO 側（SU(2)²）と被覆 U(3) の両方に使うと、composite の単項式は走査で予言どおり (t/q)^2、
    SO の f̂・被覆の f̂（対角と M⊗M）・ĝ が |R| ≤ 3 の全方向で Laurent、ĝ_R（|R| = 3）は整数。
    大きな N では W_SO = (6.15)·T_[1]^{−6}·unknot と (t/q)^{1/2}·½h^{ref} の差が Laurent で、(t/q)^{1/2} g^{ref}_[1]（K3'）に等しい。"""
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "scripts"))
    import rk_k5_trefoil_3box as X
    for n, a, b in ((3, (2, 2), (2,)), (3, (1, 1, 1), ()), (4, (2, 2, 2), (2,)), (4, (2, 1, 1, 1), (1,))):
        assert X.uniform_ratio(n, a, b) == X.ONE
    N = 2
    assert X.scan_pairs(lambda ma, mb: X.cover_Z(N, {(2,): ma, (1, 1): mb}, uniform=True), rng=range(2, 7)) == [(4, 4)]
    F = X.plog(X.so_Z(N, X.so_uniform(N)))
    Fc = X.plog(X.cover_Z(N, {}, uniform=True))
    from fractions import Fraction
    G = {k: F.get(k, X.ZERO) - X.conv(Fraction(1, 2)) * Fc.get(k, X.ZERO) for k in set(F) | set(Fc)}
    for E in (F, Fc, G):
        assert all(all(X.hatted_by_rho(E, k).values()) for k in (2, 3))
    F2c = X.plog2(X.cover_Z2(N, {}, uniform=True))
    assert all(all(X.hatted2_by_rho(F2c, *ij).values()) for ij in ((1, 1), (2, 1)))
    assert all(v.is_integer_laurent() for v in X.hatted(G, 3).values())
    w = s / r
    Tm = r * U / s**2
    assert is_laurent(W_SO(Tm**-6) - w * h_half_ref())
    assert sp.cancel(W_SO(Tm**-6) - w * h_half_ref() - w * g_ref()) == 0


def test_torus_family_uniform_normalization():
    """T(2,2m+1) の data は行列式で (q/t)^{δ n/2}（U(n) の det、δ = m）ずれる。一様な規格化 (t/q)^{δ|μ|/2} R_μ では
    固定 N = 2 の SO のベクトル表現が AS の Z_V·T_[1]^{−2(2m+1)}·unknot に一致し、大きな N で W_SO(T^{−2(2m+1)}) − (t/q)^{δ/2}½h は Laurent（m = 2）。"""
    import json
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "scripts"))
    import rk_k3_torus_family as Tf
    import rk_k5_trefoil_2box as T2
    from rcsknots.kauffman_knot import r as r_, s as s_

    def Rk(c, n):
        d = json.load(open(ROOT / f"data/torus/T2_5_{c}.json"))
        return sp.expand(sum(v * (s_**n)**(-int(k.split(",")[0])) * r_**(-int(k.split(",")[1])) * s_**(-int(k.split(",")[2])) for k, v in d["poly"].items()))
    m = 2
    assert sp.simplify(Rk("11", 2) - (r_ / s_)**(2 * m)) == 0
    Tm = r * U / s**2
    WSO = Z_V_T2(m) * Tm**(-2 * (2 * m + 1)) * Tf.unknot_SO
    lhs = (s_ / r_)**(2 * m) * Rk("1", 2)**2 * T2.unknot_K(2, (1,))
    assert sp.cancel(lhs - WSO.subs({r: r_, s: s_, U: s_**4})) == 0
    assert is_laurent(WSO - (s / r)**m * Tf.h_half(m))



def test_crosscap_criterion_uniform_normalization():
    """タスク B の判定式を一様な規格化で（scripts/rk_b_uniform.py）。T(2,2m+1) の [1]（m = 1, 2）：
    SO(2) の点で g^{true}_[1] = 1、𝐭 について偶、c = 2 = −(𝐚²−𝐪²)(𝐭²−1)𝐚^{−(4m+1)}𝐪^{−(2m−1)}𝐭^{−4m}Π_m（Π_m の係数は 1）、
    c = 1 の G 展開は ν^β ごとに同符号か (𝐭²−1)×同符号。固定 N = 2 の g_R（|R| ≤ 3）は一様な規格化のときだけ r → −r で不変。"""
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "scripts"))
    import rk_b_uniform as B
    for m in (1, 2):
        so2, even_t, c2, Pi = B.check_one(m)
        assert so2 == {1} and even_t and c2
        assert set(sp.Poly(Pi, B.Q, B.Tb, B.A).coeffs()) == {1}
    assert B.classify_c1(2) == {-10: "coherent+", -8: "tfactor+", -6: "coherent−", -4: "coherent+"}
    par = B.fixedN_parity(2)
    assert par["一様"] == (6, 6) and all(par["ĝ_R(−r) = (−1)^(|R|−1) ĝ_R(r)"].values())
    assert all(par[k][0] < par[k][1] for k in ("被覆 m=1・SO 一様", "SO K3'・被覆 m=1", "SO K3'・被覆一様"))


def test_T25_fixed_N_uniform_hatted():
    """T(2,5)（δ = 2）の固定 N = 2 の 3 箱を一様な規格化 data × (t/q)^{|μ|} で（scripts/rk_k5_torus_fixedN_hatted.py）。
    SO の f̂・被覆の f̂（対角）・ĝ が |R| ≤ 3 の全方向で Laurent、ĝ_R（|R| = 3）は整数、g_R は r → −r で不変、
    被覆の composite (2;1) の単項式は走査で予言どおり (t/q)^4。"""
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "scripts"))
    import rk_k5_trefoil_2box as T2
    import rk_k5_trefoil_3box as X
    from fractions import Fraction
    old_knot, old_delta = T2.KNOT, X.DELTA
    try:
        T2.KNOT, X.DELTA = "T2_5", 2
        X.RP.cache_clear()
        N, W = 2, X.W_
        rr = sp.Symbol("r", positive=True)
        F = X.plog(X.so_Z(N, X.so_uniform(N)))
        Fc = X.plog(X.cover_Z(N, {}, uniform=True))
        G = {k: F.get(k, X.ZERO) - X.conv(Fraction(1, 2)) * Fc.get(k, X.ZERO) for k in set(F) | set(Fc)}
        for E in (F, Fc, G):
            assert all(all(X.hatted_by_rho(E, k).values()) for k in (2, 3))
        assert all(v.is_integer_laurent() for v in X.hatted(G, 3).values())
        for k in (1, 2, 3):
            for v in X.to_schur(G, k).values():
                e = v.to_sympy()
                assert sp.cancel(e - e.subs(rr, -rr)) == 0
        # (2;1) の単項式を (t/q)^{β/2} に替えると、β = 8 のときだけ Z^cov が Laurent
        def Zc(m):
            return X.cover_Z(N, {(2,): m}, uniform=True)
        Z0, Z1 = Zc(X.ZERO), Zc(X.ONE)
        B = {k: Z1.get(k, X.ZERO) - Z0.get(k, X.ZERO) for k in set(Z0) | set(Z1)}
        hits = [b for b in range(5, 12) if all((Z0.get(k, X.ZERO) + W**b * B.get(k, X.ZERO)).is_laurent() for k in B)]
        assert hits == [8]
    finally:
        T2.KNOT, X.DELTA = old_knot, old_delta
        X.RP.cache_clear()


def test_double_twist_offset_and_41_fixed_N():
    """二重ツイスト結び目（scripts/rk_k5_41_data.py の data）：data の行列式表現のずれ (q/t)^c の c は [1] の superpolynomial の δ 次数
    𝐚 + 𝐪/2 − 𝐭（5_2 で 1、6_1 で 0、ともに種数 1）。T(2,7) では c = 3。
    4_1（c = 0）の固定 N = 2：SO 側の framing は走査で一意に 1、SO の f̂（|R| ≤ 3）は全方向で整数 Laurent。
    被覆 U(3) の 2 箱（随伴 [2,1] は KN の候補）：随伴の framing は一意に 1、ĝ（|R| = 2）は Laurent、g は (r,s) → (−r,−s) で不変。"""
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "scripts"))
    import rk_k5_trefoil_2box as T2
    import rk_k5_trefoil_3box as X
    import rk_k5_dt_fixedN as D
    from rk_k5_torus_fixedN_hatted import scan_rep
    from fractions import Fraction
    rr, ss = sp.symbols("r s", positive=True)
    old_knot, old_delta = T2.KNOT, X.DELTA
    try:
        for knot, m, n, c in (("5_2", 2, 1, 1), ("6_1", -2, 1, 0)):
            T2.KNOT = knot
            assert D.delta_grading(m, n) == [c]
            assert sp.simplify(T2.Rdata("11", 2) - (rr / ss)**(2 * c)) == 0
        T2.KNOT = "T2_7"
        assert sp.simplify(T2.Rdata("11", 2) - (rr / ss)**6) == 0
        T2.KNOT, X.DELTA = "4_1", 0
        X.RP.cache_clear()
        N = 2
        fr = X.so_uniform(N)
        assert scan_rep(lambda f: X.so_Z(N, f), fr, (3, 3), 6, width=3) == [0]
        F = X.plog(X.so_Z(N, fr))
        assert all(v.is_integer_laurent() for k in (2, 3) for v in X.hatted(F, k).values())
        assert scan_rep(lambda mm: X.cover_Z(N, mm, kmax=2, uniform=True), {(1,): X.ONE}, (1,), 0, width=3) == [0]
        Fc = X.plog(X.cover_Z(N, {}, kmax=2, uniform=True), kmax=2)
        G = {k: F.get(k, X.ZERO) - X.conv(Fraction(1, 2)) * Fc.get(k, X.ZERO) for k in set(F) | set(Fc) if sum(k) <= 2}
        assert all(X.hatted_by_rho(G, 2).values())
        for v in X.to_schur(G, 2).values():
            e = v.to_sympy()
            assert sp.cancel(e - e.subs({rr: -rr, ss: -ss}, simultaneous=True)) == 0
    finally:
        T2.KNOT, X.DELTA = old_knot, old_delta
        X.RP.cache_clear()

