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


def test_greedy_cells_and_mirror():
    # 重み q, t, qt の 3 方向の弦の積 [2]_q [2]_t [6]_qt（r² = q、s² = t）：貪欲法は 3 セルを選び、最も長い qt が最初
    box = {}
    for a in range(3):
        for b in range(3):
            for c in range(7):
                k = (2 * (a + c), 2 * (b + c), 0); box[k] = box.get(k, 0) + 1
    cand = [(x, y, 0, sg) for x in (0, 2, 4) for y in (0, 2, 4) if (x, y) != (0, 0) for sg in (1, -1)]
    cost = lambda P: sum(abs(v) for v in P.values())
    seq, _ = TT.greedy_cells(box, cand, cost)
    assert [x for x, _, _ in seq] == [(2, 2, 0, -1), (0, 2, 0, -1), (2, 0, 0, -1)]
    # 単体の形（次数の和 ≤ 4）の打ち切りでは qt は選ばれない
    simplex = {}
    for a in range(5):
        for b in range(5 - a):
            for c in range(5 - a - b):
                k = (2 * (a + c), 2 * (b + c), 0); simplex[k] = simplex.get(k, 0) + 1
    assert {x[:2] for x, _, _ in TT.greedy_cells(simplex, cand, cost)[0]} == {(0, 2), (2, 0)}
    # V = 1 + t + qt は鏡映 (q, t) → (1/t, 1/q) で自己双対、1 + t は双対でない
    assert TT.mirror_sign({(0, 0, 0): 1, (0, 2, 0): 1, (2, 2, 0): 1}) == {0: "+"}
    assert TT.mirror_sign({(0, 0, 0): 1, (0, 2, 0): 1}) == {0: "×"}


def test_choice_milp_selects_pattern():
    import numpy as np
    keys = ["a", "b"]
    mult = {"a": 1, "b": 1}
    A, b = np.eye(2), np.array([1.0, -1.0])
    same = [{"a": 1, "b": 1}, {"a": -1, "b": -1}]
    alt = [{"a": 1, "b": -1}, {"a": -1, "b": 1}]
    part = lambda k: "p"
    # 同符号の型だけでは e = (1, −1) を作れない（c が奇数なので w = 1）
    assert TT.choice_milp(keys, mult, A, b, part, {"p": same}) is None
    mp, chosen, e = TT.choice_milp(keys, mult, A, b, part, {"p": same + alt})
    assert mp == 0 and chosen["p"] == 2 and e == {"a": 1, "b": -1}
    # c = 2 で e = 0：入れ替わる対が 1 個
    mp, _, e = TT.choice_milp(["c"], {"c": 2}, np.array([[1.0]]), np.array([0.0]), lambda k: "p", {"p": [{"c": 1}]})
    assert mp == 1 and e == {"c": 0}


def test_choice_milp_ranges():
    import numpy as np
    # e_a + e_b = 0、c = (1, 1)：e = (1, −1) か (−1, 1)。対の数 0 のままどちらも可なので範囲は (−1, 1)
    keys = ["a", "b"]; mult = {"a": 1, "b": 1}
    pats = {k: [{k: 1}, {k: -1}] for k in keys}
    rg = TT.choice_milp(keys, mult, np.array([[1.0, 1.0]]), np.array([0.0]), lambda k: k, pats, max_pairs=0, ranges=keys)
    assert rg == {"a": (-1, 1), "b": (-1, 1)}
    # e_a = 1 を加えると一意
    rg = TT.choice_milp(keys, mult, np.array([[1.0, 1.0], [1.0, 0.0]]), np.array([0.0, 1.0]), lambda k: k, pats, max_pairs=0, ranges=keys)
    assert rg == {"a": (1, 1), "b": (-1, -1)}


def test_rule_milp_ranges():
    import numpy as np
    keys = ["a", "b"]; mult = {"a": 1, "b": 1}
    own = lambda k: k
    rg = TT.rule_milp(keys, mult, np.array([[1.0, 1.0]]), np.array([0.0]), own, lambda k: 0, max_pairs=0, ranges=keys)
    assert rg == {"a": (-1, 1), "b": (-1, 1)}
    rg = TT.rule_milp(keys, mult, np.array([[1.0, 1.0], [1.0, 0.0]]), np.array([0.0, 1.0]), own, lambda k: 0, max_pairs=0, ranges=keys)
    assert rg == {"a": (1, 1), "b": (-1, -1)}
    # c = 2、条件なし：対 0 なら e = ±2、対 1 まで許すと e = 0 も入るが範囲は (−2, 2) のまま
    assert TT.rule_milp(["c"], {"c": 2}, np.zeros((1, 1)), np.array([0.0]), own, lambda k: 0, max_pairs=0, ranges=["c"]) == {"c": (-2, 2)}


def test_rule_milp_samples():
    import numpy as np
    # e_a + e_b = 0、c = (1, 1)：解は (1, −1) と (−1, 1) の 2 個。標本はどちらかで、条件を満たす
    keys = ["a", "b"]; mult = {"a": 1, "b": 1}
    sm = TT.rule_milp(keys, mult, np.array([[1.0, 1.0]]), np.array([0.0]), lambda k: k, lambda k: 0, max_pairs=0, samples=(6, 0))
    assert len(sm) == 6 and all(x in ({"a": 1, "b": -1}, {"a": -1, "b": 1}) for x in sm)
    assert len({tuple(x.values()) for x in sm}) == 2


def test_rule_milp_flip():
    import numpy as np
    keys = ["a", "b"]; mult = {"a": 1, "b": 1}
    # e_a + e_b = 0：e0 = (1, −1) から両方の反転 (−1, 1) が可能
    val, e = TT.rule_milp(keys, mult, np.array([[1.0, 1.0]]), np.array([0.0]), lambda k: k, lambda k: 0, max_pairs=0, flip={"a": 1, "b": -1})
    assert val == -2 and e == {"a": -1, "b": 1}
    # e_a = 1 も課すと反転できず、最小値 = Σ|c| = 2（全部一意の証明）
    val, e = TT.rule_milp(keys, mult, np.array([[1.0, 1.0], [1.0, 0.0]]), np.array([0.0, 1.0]), lambda k: k, lambda k: 0, max_pairs=0, flip={"a": 1, "b": -1})
    assert val == 2 and e == {"a": 1, "b": -1}


def test_row_min():
    # c = (3, 1)、w = (1, −1)：S = 2 は対なしで (3, 1) の 1 通り、S = 0 は対 1 個で (±1, ±1)、S = 6 は届かない
    assert TT.row_min([3, 1], [1, -1], 2) == (0, [{3}, {1}])
    assert TT.row_min([3, -1], [1, -1], 0) == (1, [{-1, 1}, {-1, 1}])
    assert TT.row_min([3, 1], [1, -1], 6) == (None, None)


def test_row_min_matches_rule_milp():
    # 2 行（a, b, c と d, e）の条件で、行ごとの動的計画法と全体の整数計画の、対の最小と e の範囲が一致する
    import numpy as np
    rng = np.random.default_rng(1)
    keys = list("abcde"); rows_ = {"a": 0, "b": 0, "c": 0, "d": 1, "e": 1}
    for _ in range(3):
        mult = {k: int(rng.choice([-1, 1]) * rng.integers(1, 5)) for k in keys}
        w = {k: int(rng.choice([-1, 1])) for k in keys}
        e0 = {k: int(rng.choice(range(-abs(c), abs(c) + 1, 2))) for k, c in mult.items()}
        A = np.array([[float(w[k]) if rows_[k] == j else 0.0 for k in keys] for j in (0, 1)])
        b = A @ np.array([float(e0[k]) for k in keys])
        mp, rg = 0, {}
        for j in (0, 1):
            ks = [k for k in keys if rows_[k] == j]
            best, vals = TT.row_min([mult[k] for k in ks], [w[k] for k in ks], int(b[j]))
            mp += best; rg.update({k: (min(v), max(v)) for k, v in zip(ks, vals)})
        res = TT.rule_milp(keys, mult, A, b, lambda k: k, lambda k: 0, uniq=False)
        assert res["min_pairs"] == mp
        assert TT.rule_milp(keys, mult, A, b, lambda k: k, lambda k: 0, max_pairs=mp, ranges=keys) == rg


def test_hidden_milp():
    import numpy as np
    # c = (1, 0)：e_a = 1, e_b = 2 を求めると、b に見えない対が 1 個要る
    mult = {"a": 1}
    best, e, h = TT.hidden_milp(["a", "b"], mult, np.eye(2), np.array([1.0, 2.0]))
    assert (best, e, h) == (1, {"a": 1, "b": 2}, {"a": 0, "b": 1})
    # e_a + e_b = 3：Σh = 1 の解は (1, 2) と (3, 0)
    A, b = np.array([[1.0, 1.0]]), np.array([3.0])
    assert TT.hidden_milp(["a", "b"], mult, A, b)[0] == 1
    assert TT.hidden_milp(["a", "b"], mult, A, b, max_h=1, ranges=["a", "b"], hrange=["b"]) == {"a": (1, 3), "b": (0, 2), ("h", "b"): (0, 1)}


def test_pairs_onehot_matches_row_min():
    # 1 行の問題：0/1 変数の定式化と動的計画法の対の最小が一致する。届かない S は解なし
    import numpy as np
    keys = ["a", "b", "c"]; mult = {"a": 3, "b": -1, "c": 2}; w = [1, -1, 1]
    A = np.array([w], float)
    for S in (6, 2, 0, -1, 7):
        best, _ = TT.row_min([mult[k] for k in keys], w, S)
        st, val, _ = TT.pairs_onehot(keys, mult, A, np.array([float(S)]))
        assert val == best


def test_n2_quotient_rows_equivalent():
    # 割った形の条件の階数が元の条件と同じ（(i) 154、(ii) 142、(i)+(ii) 285）で、続き 79 の e が (i) を満たす
    import numpy as np, pickle, os
    Q = TT.n2_quotient()
    mult, _ = TT.n2_setup()
    keys = sorted(mult)
    ranks = [np.linalg.matrix_rank(TT.quotient_rows(keys, Q, use)[0]) for use in (("i22", "i4"), ("ii",), ("i22", "i4", "ii"))]
    assert ranks == [154, 142, 285]
    path = "data/refined_kauffman/daha/twisted_n2_pattern.pkl"
    if os.path.exists(path):
        e0 = pickle.load(open(path, "rb"))["e"]
        A, b = TT.quotient_rows(keys, Q, ("i22", "i4"))
        assert np.abs(A @ np.array([e0[k] for k in keys], float) - b).max() < 1e-9
