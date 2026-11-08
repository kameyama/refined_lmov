"""T1（続き 107）：unrefined の随伴の大きな N の Rosso–Jones（rk_t1_composite_rj）が DAHA の R_adj|_{t=q} を再現する（T(2,3)、全行、単項式の倍を除く）。
ψ^r(s_[1;1]) の展開で中間の κ の項が消え、自明表現の係数が r − 1 になること。"""
import sys, pickle
import sympy as sp
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import rk_t1_composite_rj as T1

q, a = T1.q, T1.a
r_, s_ = sp.symbols("r s", positive=True)


def test_adjoint_rj_matches_daha_trefoil():
    Rad, ps = T1.adjoint_rj(2, 3)
    rows = T1.rows_a(Rad)
    radj = pickle.load(open("data/refined_kauffman/daha/Radj_A_m1.pkl", "rb"))
    daha = {int(j): sp.expand(sp.sympify(v, locals={"r": r_, "s": s_}).subs({r_: sp.sqrt(q), s_: sp.sqrt(q)})) for j, v in radj.items()}
    daha = {j: v for j, v in daha.items() if v != 0}
    assert sorted(rows) == sorted(daha)                       # 同じ ν の冪（ν^{-4}..ν^{-10}）
    ratios = {sp.cancel(rows[j] / daha[j]) for j in rows}
    assert len(ratios) == 1                                   # 全行で同じ単項式の倍
    rt = ratios.pop(); n_, d_ = sp.fraction(rt)
    assert sp.Poly(n_, q).is_monomial and sp.Poly(d_, q).is_monomial


def test_intermediate_kappa_terms_vanish():
    for r in (2, 3):
        ps = T1.psi_adjoint(r)
        assert all(sum(al) in (0, r) for (al, be) in ps)      # |α| = r（κ = ∅）か自明表現だけ
        assert ps[((), ())] == r - 1
