"""refined Kauffman K2'(a)：双対 Cauchy 型の Hopf 絡み目（probe 各 1 個）の log Z~ の係数。"""
import sys
from fractions import Fraction as F
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rcsknots.askey_wilson import askey_wilson_monic  # noqa: E402
from rcsknots.kauffman_hopf import hopf_log_coeffs, mobius_f, qs, ts, Us  # noqa: E402
from rcsknots.kauffman_refined import koornwinder_eval  # noqa: E402
from rcsknots.kauffman_refined import q as _q, t as _t, U as _U  # noqa: E402

A = B = 2
L = hopf_log_coeffs(A, B)


def _finite_N_log(N, qv, sqv, tv):
    """有限 N で Z~ を Fraction で直接計算して log の係数を返す（閉じた形を使わない参照値）。"""
    aw = askey_wilson_monic(N, tv, F(1), F(-1), sqv, -sqv)
    tb = (1, -1, sp.Rational(sqv.numerator, sqv.denominator), -sp.Rational(sqv.numerator, sqv.denominator))
    sub = {_q: sp.Rational(qv.numerator, qv.denominator), _t: sp.Rational(tv.numerator, tv.denominator),
           _U: sp.Rational(tv.numerator, tv.denominator)**N}
    Z = {}
    for k in range(0, min(N, A) + 1):
        w = koornwinder_eval((1,) * k if k else (), tb).subs(sub)
        Wk = F(int(w.p), int(w.q))
        ky = {k + 2 * j: aw[N - k].get(N - k - 2 * j, F(0)) for j in range((A - k) // 2 + 1)}
        pz = {0: F(1)}
        for i in range(1, N + 1):
            xi = tv**(N - i) * (qv if i <= k else 1)
            X = xi + 1 / xi
            new = {}
            for b, v in pz.items():
                for db, f in ((0, F(1)), (1, -X), (2, F(1))):
                    if b + db <= B:
                        new[b + db] = new.get(b + db, F(0)) + v * f
            pz = new
        for a, va in ky.items():
            for b, vb in pz.items():
                Z[(a, b)] = Z.get((a, b), F(0)) + (-1)**k * Wk * va * vb
    X = {ab: v for ab, v in Z.items() if ab != (0, 0)}
    assert Z[(0, 0)] == 1

    def mul(P, Q):
        R = {}
        for (a1, b1), v1 in P.items():
            for (a2, b2), v2 in Q.items():
                if a1 + a2 <= A and b1 + b2 <= B:
                    R[(a1 + a2, b1 + b2)] = R.get((a1 + a2, b1 + b2), F(0)) + v1 * v2
        return R

    Lf = {}
    P = dict(X)
    for m in range(1, A + B + 1):
        for ab, v in P.items():
            Lf[ab] = Lf.get(ab, F(0)) + F((-1)**(m + 1), m) * v
        P = mul(P, X)
    return Lf


def test_symbolic_matches_finite_N():
    """U を形式変数として求めた L_{ab} が、U = t^N（N = 3, 4）での直接計算と一致する。"""
    qv, sqv, tv = F(9, 49), F(3, 7), F(5, 11)
    for N in (3, 4):
        Lf = _finite_N_log(N, qv, sqv, tv)
        sub = {qs: sp.Rational(9, 49), ts: sp.Rational(5, 11), Us: sp.Rational(5, 11)**N}
        for ab, v in Lf.items():
            assert L[ab].subs(sub) == sp.Rational(v.numerator, v.denominator), ab


def test_unknot_and_hopf_fundamental():
    """f_{a0} = 0 (a ≥ 2)、f_{10} = (1−U)(U+t)/(U(t−1))、
    f_{11} = −f_{10} z_q (ν~ − ν~^{−1})（z_q = q^{1/2} − q^{−1/2}、ν~ = U q^{1/2}/t）、f_{12} = f_{11} (U/t + t/U)。"""
    f10 = mobius_f(L, 1, 0)
    assert sp.cancel(f10 - (1 - Us) * (Us + ts) / (Us * (ts - 1))) == 0
    assert mobius_f(L, 2, 0) == 0
    zq = (qs - 1) / sp.sqrt(qs)
    nut = Us * sp.sqrt(qs) / ts
    f11 = mobius_f(L, 1, 1)
    assert sp.simplify(f11 + f10 * zq * (nut - 1 / nut)) == 0
    assert sp.cancel(mobius_f(L, 1, 2) - f11 * (Us / ts + ts / Us)) == 0


def test_multi_probe_consistency():
    """成分 1 に probe 無限個（冪和）の生成関数：Schur 基底の f_{(a),b} が 1 変数 probe の f_{ab} に一致し、
    unknot 部分 f_{R,0} は |R| ≥ 2 で 0。Macdonald 基底では f_{(1,1),1} の分母に (qt − 1) が現れる。"""
    from rcsknots.kauffman_hopf import hopf_Ztilde_multi, log_coeffs_multi, plethystic_f
    Lm = log_coeffs_multi(hopf_Ztilde_multi(2, 2), 2, 2)
    fs = plethystic_f(Lm, 2, 2, "schur")
    for a in (1, 2):
        for b in (0, 1, 2):
            assert sp.cancel(fs[((a,), b)].as_expr() - mobius_f(L, a, b)) == 0, (a, b)
    assert fs[((1, 1), 0)] == 0 and fs[((2,), 0)] == 0
    num, den = sp.fraction(sp.cancel(fs[((1, 1), 1)].as_expr()))
    assert sp.factor(den).has(ts - 1) and not sp.factor(den).has(qs * ts - 1)
    fm = plethystic_f(Lm, 2, 2, "macdonald")
    num, den = sp.fraction(sp.cancel(fm[((1, 1), 1)].as_expr()))
    assert sp.factor(den).has(qs * ts - 1)
