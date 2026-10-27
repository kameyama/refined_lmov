"""rk_f_twisted22_lam の機能テスト：冪和 → Schur の分解と、整数計画の条件ごとの重み。"""
import sys
from fractions import Fraction
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_f_twisted22_lam as L
import rk_b_twisted_multiplet as TM


def test_schur_split_of_products():
    q = Fraction(1, 4)
    chi = {"2": {(1, 1): 1, (2,): 1}, "11": {(1, 1): 1, (2,): -1}}
    for l in L.LAMS:
        for m in L.LAMS:
            # s_l(Y) s_m(Y′) = Σ_{ρ,σ} χ^l(ρ) χ^m(σ) p_ρ p_σ′ / (z_ρ z_σ)、z = 2
            a = {(rho, sig): q * chi[l][rho] * chi[m][sig] for rho in ((1, 1), (2,)) for sig in ((1, 1), (2,))}
            M = L.schur_split(a)
            assert all(M[k] == (1 if k == (l, m) else 0) for k in M)


def test_solve_weight_splits_variables():
    r, s, nu = TM.r, TM.s, TM.nu
    mult = {("a", "mono", 0, 0, 0): 1, ("b", "mono", 0, 0, 0): 1}
    basis = lambda key: TM.twisted(key[1:], "mono")
    # a だけで +1、和で 0 → b = −1
    conds = [(sp.Integer(1), {}, lambda k: sp.Integer(k[0] == "a")), (sp.Integer(0), {})]
    res = TM.solve(mult, basis, conds)
    assert res["feasible"] and not res["free"]
    assert res["e"] == {("a", "mono", 0, 0, 0): 1, ("b", "mono", 0, 0, 0): -1}
