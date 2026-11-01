"""n = 1 の N を [1] と随伴の a の行で書く母関数（rk_f_rows_gen、theory §12.33、続き 103）の機能テスト。
  (1) T(2,3)：母関数から作った n_k が cover_sectors から直接作った N の全ブロックと一致し、k > K で 0。
  (2) T(2,3)・T(2,5)・T(3,4)：G = (A − R²)/((q−1)(t−1)) は Laurent、g_0 = 0、g_1 = (2ρ_0ρ_1 − ρ_1²)/q、移り目 n_{2δ} = Sym^{2δ}V − (ρ_0 − ρ_1)²。
  (3) T(3,4) の S̃ = ρ_0 − ρ_1 = (1 + qt)(1 + t + t² + qt² + q²t²)（§12.32 の S̃）。
  (4) 一般の τ の語：gamma_subst が (m, k) の従来の語と一致し、T3_5 の第 1 列は (5, 3)。"""
import sys
import sympy as sp
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import rk_f_rows_gen as RG
import rk_torus_tag as TG
from rcsknots.refined_kauffman.daha_dn import DN, gamma_factors, gamma_subst, gamma_first_column
from rcsknots.refined_kauffman.flint_rf import R as Rr, ONE

q, t, y = RG.q, RG.t, RG.y


def _rho(R):
    P = sp.Poly(R, y)
    return {int(m[0]): c for m, c in zip(P.monoms(), P.coeffs())}


def test_trefoil_blocks_from_rows_match_direct():
    nk, G, A, R = RG.blocks_from_rows(1, 2)
    dn = RG.direct_blocks(1, 2)
    assert all(sp.cancel(nk.get(k, 0) - dn.get(k, 0)) == 0 for k in set(nk) | set(dn))
    assert all(sp.cancel(nk[k]) == 0 for k in (4, 5))


def test_g_relations_and_transition():
    for m, k in ((1, 2), (2, 2), (1, 3), (2, 3), (1, 4)):          # T(2,3)・T(2,5)・T(3,4)・T(3,7)・T(4,5)
        nk, G, A, R = RG.blocks_from_rows(m, k, kmax=2 * TG.delta(m, k))
        rho = _rho(R)
        assert sp.expand(G.coeff(y, 0)) == 0
        assert sp.cancel(G.coeff(y, 1) - (2 * rho[0] * rho[1] - rho[1]**2) / q) == 0
        de = TG.delta(m, k)
        assert sp.cancel(nk[2 * de] - RG.sym_V(2 * de) + (rho[0] - rho[1])**2) == 0


def test_deepest_block_is_h_s_minus_2():
    """n_K = −t^{K−s+2} h_{s−2}(q, t)：T(2,5)（s = 2）、T(3,4)（s = 3）、T(4,5)（s = 4）。"""
    for m, k, h in ((2, 2, 1), (1, 3, q + t), (1, 4, q**2 + q * t + t**2)):
        nk, G, A, R = RG.blocks_from_rows(m, k)
        K, s_ = TG.rs(m, k) - 3, TG.first_column(m, k)[1]
        assert sp.cancel(nk[K] + t**(K - s_ + 2) * h) == 0 and all(sp.cancel(nk[K + i]) == 0 for i in (1, 2))


def test_t34_square_factor_is_two_rows():
    rho = _rho(RG.R_poly(1, 3))
    St = (1 + q * t) * (1 + t + t**2 + q * t**2 + q**2 * t**2)
    assert sp.expand(rho[0] - rho[1] - St) == 0
    assert sp.expand(rho[2] - t**3) == 0                       # 3 行目（入らない）


def test_gamma_words():
    D = DN(3)
    e1 = D.e(0)
    for m, k in ((1, 2), (2, 2), (1, 3)):
        new = gamma_subst(gamma_factors(m, k), [("X", e1, 1)], e1)
        old = D.gamma_X_old(0, m, k)
        def norm(w):
            c = ONE; rest = []
            for op in w:
                if op[0] == "c":
                    c = c * op[1]
                elif op[0] == "X" and abs(op[2]) > 1:
                    rest += [("X", op[1], 1 if op[2] > 0 else -1)] * abs(op[2])
                else:
                    rest.append(op)
            return str(c.to_sympy()), rest
        assert norm(new) == norm(old)
        assert gamma_first_column(gamma_factors(m, k)) == (m * k + 1, k)
    assert gamma_first_column(gamma_factors("T3_5")) == (5, 3) and TG.rs("T3_5") == 15 and TG.delta("T3_5") == 4
