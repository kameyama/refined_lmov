"""rk_f_twisted_tensor の機能テスト：条件の行の作り方（β 奇の行だけを使う指定）と、全状態が不変（e = ±c）の整数計画。"""
import sys
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_f_twisted_tensor as TT
import rk_b_twisted_multiplet as TM

r, s, nu = TT.r, TT.s, TT.nu


def test_rows_beta_odd_filter():
    mult = {("mono", 0, 0, 2): 1, ("mono", 0, 0, 4): 1}       # ν² は β = 1（奇）、ν⁴ は β = 2（偶）
    basis = lambda k: TM.twisted(k, "mono")
    conds = {"a": (nu**2 + nu**4, {}, lambda k: 1)}
    _, A, b = TT.rows(mult, basis, conds)
    _, A2, b2 = TT.rows(mult, basis, conds, beta_odd_only=("a",))
    assert len(b) == 2 and len(b2) == 1


def test_all_invariant_solve():
    mult = {("mono", 0, 0, 0): 1, ("mono", 2, 0, 0): -1}      # 1 − r²
    basis = lambda k: TM.twisted(k, "mono")
    res = TT.solve(mult, basis, [(1 + r**2, {})], all_inv=True)
    assert res["feasible"] and not res["free"]
    assert res["e"] == {("mono", 0, 0, 0): 1, ("mono", 2, 0, 0): 1}


def test_rule_milp_sign_rule_and_pairs():
    import numpy as np
    # 同じ群の 2 単項式、d = 0, 1：e = (σ, −σ)·w。条件 e_a = 1 → σ = 1、e_b = −1（一意）
    keys = ["a", "b"]
    mult = {"a": 1, "b": 1}
    res = TT.rule_milp(keys, mult, np.array([[1.0, 0.0]]), np.array([1.0]), lambda k: 0, lambda k: keys.index(k))
    assert res["e"] == {"a": 1, "b": -1} and res["free"] == [] and res["min_pairs"] == 0
    # 規則に反する条件 e_a + e_b = 2 は解なし
    assert TT.rule_milp(keys, mult, np.array([[1.0, 1.0]]), np.array([2.0]), lambda k: 0, lambda k: keys.index(k)) is None
    # c = 2 の単項式で e = 0：入れ替わる対が 1 個。対を許さなければ解なし
    res = TT.rule_milp(["c"], {"c": 2}, np.array([[1.0]]), np.array([0.0]), lambda k: 0, lambda k: 0)
    assert res["e"] == {"c": 0} and res["pairs"] == 1 and res["min_pairs"] == 1
    assert TT.rule_milp(["c"], {"c": 2}, np.array([[1.0]]), np.array([0.0]), lambda k: 0, lambda k: 0, pairs=False) is None
