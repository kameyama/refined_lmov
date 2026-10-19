"""タスク C：D_N・GL_N の DAHA の多項式表現と DAHA-Jones 多項式（src/rcsknots/daha_dn.py）、一般の N の補間（scripts/rk_c_*.py）。"""
import sys
from pathlib import Path

import pytest
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

rr, ss, uu = sp.symbols("r s u", positive=True)
mir = lambda e: e.subs({rr: 1 / rr, ss: 1 / ss}, simultaneous=True)


def test_dn_relations_and_macdonald_eigen():
    """D_3：2 次・組紐・交差関係、Y の可換性、保存済みの D_N の Macdonald 多項式 K_(1,1)（U = t³）が Σ(Y + Y^{−1}) の固有関数。"""
    import random
    from rcsknots.daha_dn import DN, K_poly, macdonald_op, _add
    from rcsknots.flint_rf import R, S, conv, ONE
    random.seed(0)
    D = DN(3)
    f = {}
    for _ in range(4):
        _add(f, {tuple(random.randint(-2, 2) for _ in range(3)): conv(random.randint(1, 5))})
    eq = lambda a, b: not _add(dict(a), b, conv(-1))
    for i in range(3):
        assert not _add(_add(D.T(i, D.T(i, f)), D.T(i, f), -(S - 1 / S)), f, conv(-1))
    assert eq(D.T(0, D.T(1, D.T(0, f))), D.T(1, D.T(0, D.T(1, f))))
    assert eq(D.T(0, D.Xmul(D.e(0), D.T(0, f))), D.Xmul(D.reflect(0, D.e(0)), f))
    Y2 = lambda g: D.T(0, D.Y1(D.T(0, g, inv=True)), inv=True)
    assert eq(D.Y1(Y2(f)), Y2(D.Y1(f)))
    K = K_poly(D, (1, 1))
    ev = sum((S**(2 * D.rho[i]) * R**(2 * x) + S**(-2 * D.rho[i]) * R**(-2 * x) for i, x in enumerate((1, 1, 0))), conv(0))
    assert not _add(macdonald_op(D, K), K, -ev)


def test_dn_daha_jones_matches_data_and_as():
    """D_N の DAHA-Jones（T(3,2)）：N = 2 で Z = T_λ^{−6} × 鏡像（一様な規格化の data の W_λ）、N = 4 の [1] は鏡像の AS (6.15)。"""
    from rcsknots.daha_dn import daha_jones
    from rcsknots.kauffman_trefoil import Z615, r as rT, s as sT, U as UT
    import rk_k5_trefoil_2box as T
    Zv, Kv = daha_jones(2, (2,), 1)
    W = sp.cancel(T.Rdata("2", 2)**2 * (ss / rr)**4 * T.unknot_K(2, (2,)))
    assert sp.cancel(Zv.to_sympy() / mir(W) - 1 / (rr**24 * ss**24)) == 0
    Zv, Kv = daha_jones(4, (1,), 1)
    assert sp.cancel((Zv / Kv).to_sympy() - mir(Z615.subs({rT: rr, sT: ss, UT: ss**8}))) == 0


def test_gl_daha_jones_monomial_rule_and_adjoint():
    """GL_N の DAHA-Jones = 単項式 q^{−3(b,b)} t^{−6(b,ρ)+2|b|(N−1)} × 鏡像（一様な規格化の data）：GL_3 の [2,1]、U(3) の随伴 (1;1)。"""
    from rcsknots.daha_dn import GLN, gl_gamma_P, gl_P_poly, gl_adj_gamma
    from rk_c_daha_cover import mono
    import rk_k5_trefoil_2box as T
    D = GLN(3)
    Z = D.evaluate(gl_gamma_P(D, (2, 1), 1, 2)).to_sympy()
    K = D.evaluate(gl_P_poly(D, (2, 1))).to_sympy()
    W = (ss / rr)**3 * T.Rdata("21", 3)
    assert sp.cancel(Z / K / (mir(W) * mono(3, [2, 1, 0])) - 1) == 0
    Zadj = D.evaluate(gl_adj_gamma(D, 1, 2)).to_sympy()
    Wadj = (ss / rr)**3 * T.Rdata("21", 3) * T.P_A((2, 1), 3)
    assert sp.cancel(Zadj / (mir(Wadj) * mono(3, [1, 0, -1])) - 1) == 0


def test_interpolation_vector_reproduces_as():
    """D_N の DAHA-Jones の [1] を N = 2〜5 で計算して U で補間すると、N = 6 を予言し、AS (6.15) × T_[1]^{−6} に一致する。"""
    import rk_c_interp as I
    from rcsknots.flint_rf import from_sympy
    from rcsknots.kauffman_trefoil import Z615, r as rT, s as sT, U as UT
    import os, pickle
    from rcsknots.daha_dn import daha_jones
    U = I.U
    vals = {}
    for N in range(2, 7):
        Z, K = daha_jones(N, (1,), 1)
        lamp = [1] + [0] * (N - 1)
        Tl = rr * ss**(2 * (N - 1))
        vals[N] = sp.cancel(Tl**-6 * mir(Z.to_sympy()) / K.to_sympy())
    c = I.interpolate_generic({N: vals[N] for N in (2, 3, 4, 5)}, 2, -5, -2)
    assert (from_sympy(vals[6]) - I.predict_generic(c, 2, 6)).is_zero()
    fit = sp.expand(sum(v.to_sympy() * U**j for j, v in c.items()))
    Tm = rr * U / ss**2
    assert sp.expand(fit - sp.cancel(Z615.subs({rT: rr, sT: ss, UT: U}) * Tm**-6)) == 0


def test_large_N_two_box_crosscap():
    """大きな N の 2 箱のクロスキャップ部分（scripts/rk_c_2box_largeN.py、保存済みの補間を使う）：g_(2), g_(1,1) の分母は 2 × 単項式 × (1 + t)、
    g_(2) + g_(1,1) と (1 + t)(g_(2) − g_(1,1)) は整数 Laurent、U = t² で固定 N の g（一様な規格化）と一致。"""
    import pickle
    import rk_c_2box_largeN as L
    W, co = L.so_parts()
    f = L.so_f(W, co)
    fc = L.cover_f()
    g = {k: f[k] - L.half * fc[k] for k in f}
    for k in ((2,), (1, 1)):
        den = sp.factor(sp.fraction(sp.cancel(g[k].to_sympy()))[1])
        assert sp.cancel(den * 1 / (1 + ss**2)).is_polynomial(rr, ss, uu)
        assert not g[k].is_laurent() and (g[k] * (1 + L.S * L.S) * 2).is_integer_laurent()
    assert (g[(2,)] + g[(1, 1)]).is_integer_laurent()
    fx = pickle.load(open(ROOT / "data/refined_kauffman/T2_3_fixedN_3box_g.pkl", "rb"))
    for k in ((1,), (2,), (1, 1)):
        fixed = sp.sympify(fx[2]["g_schur"][str(k)], locals={"r": rr, "s": ss})
        assert sp.cancel(g[k].to_sympy().subs(uu, ss**4) - fixed) == 0


def test_large_N_two_box_marino_dictionary():
    """t = q で Mariño 0904.1088 の trefoil（標準 framing）の鏡像に一致（scripts/rk_b_marino_tq.py）：SO の W_λ と被覆の W_μ は比 1、
    g_(1) = −鏡像 g^Mar_[1]、奇数 d に直して unknot の値を引いた M(q^{−1})^{−1} g は −鏡像(ĝ^Mar_R)（(4.8)、誤植訂正済み）。"""
    import rk_b_marino_tq as M
    assert all(M.compare().values())


def test_large_N_two_box_kn_frames():
    """KN の二つの枠（q 枠：Macdonald の Cauchy 核と M(t)、t 枠：双対 Cauchy 核と M(q^{−1})）の hatted 不変量が、
    大きな N の 2 箱の SO・被覆・クロスキャップ部分で一致する（scripts/rk_b_frames.py）。"""
    import rk_b_frames as F
    import rk_c_2box_largeN as L
    from rcsknots.flint_rf3 import S
    W, co = L.so_parts()
    zt = 1 / S - S
    ft = {k: zt * v for k, v in L.so_f(W, co).items()}
    fct = {k: zt * v for k, v in L.cover_f().items()}
    fq, fcq = F.so_q(W, co), F.cover_q()
    gq = {k: fq[k] - L.half * fcq[k] for k in fq}
    gt = {k: ft[k] - L.half * fct[k] for k in ft}
    for A, B in ((fq, ft), (fcq, fct), (gq, gt)):
        hq, ht = F.hatted(A, ss), F.hatted(B, 1 / rr)
        assert all((hq[k] - ht[k]).is_zero() for k in hq)


@pytest.mark.results
def test_large_N_two_box_crosscap_structure():
    """保存済みの大きな N の 2 箱の g（data/refined_kauffman/daha/largeN_2box_m1.pkl）の hatted ĝ（scripts/rk_b_2box_structure.py）：
    c = 1 部分の V ↔ 1/V 対称は [1] だけで成立し 2 箱では不成立、refined だけの部分 ĝ − ĝ|_{t=q}(ν~) は全成分で (q − t)(ν² − t) で割り切れる。"""
    import rk_b_2box_structure as B
    G = B.ghat()
    assert B.c1_symmetric(G[(1,)]) and not B.c1_symmetric(G["sum"])
    assert all(B.divisible_universal(G[k]) for k in ((1,), (2,), (1, 1), "sum"))


@pytest.mark.results
def test_general_N_three_box_tq_and_T25_universal_factor():
    """保存済みの一般の N の trefoil の 3 箱（R_U_m1_3, R_U_m1_111）は t = q で NRZ の [3] と F_[3](−ν, t^{−1}) の鏡像に比 1。
    大きな N の T(2,5) の 2 箱の ĝ は、述べ直しの判定式の普遍因子 (q − t)(ν² − t) で全成分割り切れる。"""
    import rk_c_check_tq as C
    import rk_unrefined_check as rk
    import rk_b_2box_structure as B
    import rk_b_2box_structure_m as Bm
    nu, t = C.nu, C.t
    assert sp.cancel(C.R_tq((3,)) / C.mir(sp.cancel(rk.nrz_trefoil(3, nu, t, sp.Integer(-1))))) == 1
    assert sp.cancel(C.R_tq((1, 1, 1)) / C.mir(sp.cancel(rk.nrz_trefoil(3, -nu, 1 / t, sp.Integer(-1))))) == 1
    f, fc, g = Bm.load(2)
    G = Bm.hat(g)
    assert all(B.divisible_universal(G[k]) for k in ((1,), (2,), (1, 1), "sum"))
