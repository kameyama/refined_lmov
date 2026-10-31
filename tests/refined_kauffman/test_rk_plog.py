"""n = 1 の母関数の plethystic log（rk_f_plog、theory §12.31、progress.md 続き 101）の機能テスト。
  (1) X の母関数の plethystic log は有限：(1 + 2t + qt) z − t² z²（V ⊕ t の Fock 空間に関係 t²z²）。
  (2) 上の系列の plethystic log は V z と、(q − t) t^{m−1−j′} c_{j′} z^{m+j′}（c₀ = 1、c₁ = V² − t(1 + qt)、c₂ から負の係数）。"""
import sys
import sympy as sp
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import rk_f_plog as PL

q, t, z = PL.q, PL.t, PL.z


def test_plog_X_finite():
    pl = PL.plog(PL.series_coeffs(PL.gen_X(), 5))
    assert pl[1] == sp.expand(1 + 2 * t + q * t) and pl[2] == -t**2 and all(pl[n] == 0 for n in range(3, 6))


def test_plog_upper_series_generators():
    m = 4
    pl = PL.plog(PL.series_coeffs(PL.gen_upper_N(m), m + 2))
    assert pl[1] == sp.expand(PL.V) and pl[2] == 0 and pl[3] == 0
    c = PL.generators_N(m, 2)
    assert c[0] == 1 and c[1] == sp.expand(PL.V**2 - t * (1 + q * t))
    assert min(sp.Poly(c[2], q, t).coeffs()) < 0
