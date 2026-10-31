"""rk_f_twisted_tensor の機能テスト：条件の行の作り方（β 奇の行だけを使う指定）と、全状態が不変（e = ±c）の整数計画。"""
import sys
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import numpy as np
import pytest
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


def test_crosscap_split():
    # c = 1 は β ≡ n − 1 (mod 2)。和は元に戻り、n = 1 と n = 2 で偶奇の役が入れ替わる
    nu = TT.nu
    F = 3 * nu**-4 + 5 * nu**-3 - nu**2 + 7 * nu
    c1, c2 = TT.crosscap_split(F, 1)
    assert sp.expand(c1 - (3 * nu**-4 - nu**2)) == 0 and sp.expand(c2 - (5 * nu**-3 + 7 * nu)) == 0
    d1, d2 = TT.crosscap_split(F, 2)
    assert sp.expand(d1 - c2) == 0 and sp.expand(d2 - c1) == 0


def test_n1_crosscap_c2_vanishes_at_tq():
    # n = 1（trefoil）：(ii) の F = g_(1) + 1 の c = 2 部分は 0 でないが、t = q で消える（theory §12.18）
    import rk_b_knot_fixed_sector as K
    r, s = TT.r, TT.s
    c1, c2 = TT.crosscap_split(sp.cancel(K.load_g(1)[(1,)] + 1), 1)
    assert sp.cancel(c2) != 0 and sp.cancel(c2.subs(s, r)) == 0 and sp.cancel(c1.subs(s, r)) != 0


def test_n2_quotient_c_intersection():
    # (ii′) の割った形の目標は (i22) の目標と交点 q = t = −1 で全ブロック一致し、(ii) とは ν^{−24} で全体の符号が逆
    Q, Q0 = TT.n2_quotient("c"), TT.n2_quotient()
    R = sp.Symbol("R")
    poly = lambda d, c: sum((v * R**a for (a, b), v in d.items() if b == c), sp.Integer(0))
    for c in sorted({k[1] for k in Q["i22"]} | {k[1] for k in Q["ii"]}):
        assert sp.expand(poly(Q["i22"], c).subs(R, sp.I) * sp.I - poly(Q["ii"], c).subs(R, sp.I)) == 0
    assert sp.expand(poly(Q["ii"], -24) + poly(Q0["ii"], -24)) == 0


def test_n2_top_block_square_at_tq():
    # ν^{−12}：t = q で N_λ = −X_λ²（X は正係数）、(ii′) の目標 H′ = Σ_λ ψ²[X_λ]（X ⊗ X の入れ替えのトレース、theory §12.19）
    mult, _ = TT.n2_setup()
    Q = TT.n2_quotient("c")
    r = TT.r
    X = {"2": (r**2 + 1)**2 * (r**4 + 1)**2 / r**7, "11": (r**2 + 1)**2 * (r**4 + 1) / r**5}
    for lab in X:
        Nq = sum(int(v) * r**(k[2] + k[3]) for k, v in mult.items() if k[0] == lab and k[-1] == -12)
        assert sp.expand(Nq + X[lab]**2) == 0
    H = sum(v * r**a for (a, c), v in Q["ii"].items() if c == -12)
    assert sp.expand(H - X["2"].subs(r, r**2) - X["11"].subs(r, r**2)) == 0
    # 平方は ν^{−14} に広がらない：N_λ(q, q) は X_λ で割り切れない（theory §12.20）
    for lab in X:
        Nq = sum(int(v) * r**(k[2] + k[3]) for k, v in mult.items() if k[0] == lab and k[-1] == -14)
        assert not sp.Poly(sp.fraction(sp.cancel(Nq / X[lab]))[1], r).is_monomial


def test_top_square_from_g():
    # X_λ は g_λ の t = q での ν^{−3} の係数 G_λ と単項式の違いしかない：X_11 = q^{−2}G_2 = −q G_11、X_2 = (q + q^{−1})X_11（theory §12.20）
    import rk_b_knot_fixed_sector as K
    import rk_f_frozen_basis as FB
    r, s, nu = TT.r, TT.s, TT.nu
    g = K.load_g(1)
    G = {}
    for lab, key in (("2", (2,)), ("11", (1, 1))):
        n_, d_ = sp.fraction(sp.cancel(g[key].subs(s, r)))
        dn = sp.Poly(d_, nu).degree()
        G[lab] = sp.cancel(FB.nu_blocks(sp.expand(n_))[dn - 3] / (d_ / nu**dn))
    X11 = (r**2 + 1)**2 * (r**4 + 1) / r**5
    assert sp.cancel(X11 - G["2"] / r**4) == 0 and sp.cancel(X11 + r**2 * G["11"]) == 0


def test_pairs_onehot_extra():
    # 見えない対：c = 0 の単項式に e = ±2 を置くと費用 1。1 変数 c = 1 で e = 3 を求めれば費用 1
    keys = [("x",), ("y",)]
    mult = {("x",): 1}
    A, b = np.array([[1.0, 0.0], [0.0, 1.0]]), np.array([3.0, 2.0])
    assert TT.pairs_onehot(keys, mult, A, b)[1] is None
    assert TT.pairs_onehot(keys, mult, A, b, extra=1)[1] == 2


def test_n2_skew_parts():
    # (ii′) の目標は ψ²[[3]_q] で割り切れ、Y_(2) = (q + q^{−1})Φ/[3]_q、Y_(11) = Φ/[3]_q。積の部分 −Y_λ(ν)Y_λ(−ν) は
    # ν^{−12} で N_λ(q, q) に一致し、ν^{−14} では残りの係数の絶対値の和が小さい（theory §12.21）
    r = TT.r
    P = TT.n2_skew_parts()
    l1 = lambda e: sum(abs(x) for x in sp.Poly(sp.expand(e * r**80), r).coeffs()) if e != 0 else 0
    for lab in ("2", "11"):
        assert P["R"][lab][-12] == 0
        assert 5 * l1(P["R"][lab][-14]) < l1(P["N"][lab][-14])
        assert all(sp.Poly(sp.fraction(sp.cancel(y))[1], r).is_monomial for y in P["Y"][lab].values())
        # 偶奇：R_λ の係数は ν^{−14} まですべて偶数、ν^{−16} で初めて奇数が現れる（λ ごとの目標が e ≡ c と両立しなくなる）
        odd = lambda e: [x for x in sp.Poly(sp.expand(e * r**80), r).coeffs() if x % 2] if e != 0 else []
        assert not odd(P["R"][lab][-14]) and odd(P["R"][lab][-16])


def test_c1_c0_unrefined():
    # 合成表現の相関関数 ⟨P_α P_β^∨⟩ = W_{[α;β]} + c₁W_adj + c₀ の係数は、t = q で Koike の公式の値（c₁ = 1、c₀ = 1 または 0）になる。
    # 一般の t では A の有理式（分母に 1 − q^a A/t^b 型の因子）（theory §12.22）
    import rk_f_cover_fixed22 as CF
    u = sp.Symbol("u", positive=True)
    for key, want in ((((2,), (2,)), (1, 1)), (((2,), (1, 1)), (1, 0)), (((1, 1), (1, 1)), (1, 1))):
        c1, c0, _, _ = CF.decomposition_A(*key)
        f = lambda e: sp.cancel(e.to_sympy().subs(u, nu * s).subs(s, r))
        assert (f(c1), f(c0)) == want
        assert sp.fraction(sp.cancel(c1.to_sympy().subs(u, nu * s)))[1].has(nu)


@pytest.mark.results
def test_cover_parts_sum():
    # 保存済みの項ごとの寄与（rk_f_twisted22_lam.cover_parts）の和が N_λ に一致し、真部分集合の和は Laurent でない（theory §12.22）
    import itertools
    Np = TT.n2_parts_N()
    mult, _ = TT.n2_setup()
    names = ("comp", "c1", "c0", "low")
    for lab in ("2", "11"):
        N = sum(int(c) * r**k[2] * s**k[3] * nu**k[4] for k, c in mult.items() if k[0] == lab)
        assert sp.expand(sp.cancel(sum(Np[(n, lab)] for n in names)) - N) == 0
        for k in (1, 2, 3):
            for sub in itertools.combinations(names, k):
                e = sp.cancel(sum(Np[(n, lab)] for n in sub))
                assert not sp.Poly(sp.fraction(e)[1], r, s, nu).is_monomial


def test_rule_milp_fix():
    import numpy as np
    # e_a + e_b = 0、c = (1, 1)：e_a = 1 に固定すると e_b = −1 が一意。e_a = 1、e_b = 1 の固定は条件に反して解なし
    keys = ["a", "b"]; mult = {"a": 1, "b": 1}
    res = TT.rule_milp(keys, mult, np.array([[1.0, 1.0]]), np.array([0.0]), lambda k: k, lambda k: 0, fix={"a": 1})
    assert res["e"] == {"a": 1, "b": -1} and res["free"] == []
    assert TT.rule_milp(keys, mult, np.array([[1.0, 1.0]]), np.array([0.0]), lambda k: k, lambda k: 0, fix={"a": 1, "b": 1}) is None
    # c = 3 の単項式を e = 1 に固定：対 1 個（規則 e = σ sign(c) w、w = |c| − 2y）
    res = TT.rule_milp(["c"], {"c": 3}, np.zeros((1, 1)), np.array([0.0]), lambda k: 0, lambda k: 0, fix={"c": 1})
    assert res["e"] == {"c": 1} and res["pairs"] == 1


def test_principle_reproduces_unique_E():
    """状態の空間の同定（principle、続き 89）が T(2,3)・T(2,5) の一意な E（続き 75、(i)+(ii) から決めたもの）を
    条件を入力にせずに再現し、条件 (i)・(ii) を満たす。"""
    import io, contextlib, pickle
    for m in (1, 2):
        with contextlib.redirect_stdout(io.StringIO()):
            E, e, ok = TT.principle(m, verbose=False)
        assert ok == {"i": True, "ii": True}
        saved = pickle.load(open(f"data/refined_kauffman/daha/twisted_E_n1_m{m}.pkl", "rb"))
        Es = sp.sympify(saved["E"], locals={"r": TT.r, "s": TT.s, "nu": TT.nu})
        pre = TT.s**2 * (TT.t - TT.nu**2) * (TT.nu**2 + 1) * (TT.q + 1) * (TT.nu**2 * TT.q + 1)
        assert sp.expand(E * pre - Es) == 0


@pytest.mark.results
def test_principle_T29_consistent():
    """T(2,9)（続き 95、theory §12.29）：原理（R_{m+1} = t^{m−2}(V² + t²)、V 固定）で k ≠ 9 を固定し k = 9 の補正の対の符号だけを自由にすると
    (i)・(i)+(ii) に解があり、一意なのは 15 中 2 個。largeN_2box_m4.pkl が要る（約 10 分）。"""
    out = TT.principle_R(4)
    for use in (("i",), ("i", "ii")):
        rg = out[(6, use)]
        assert rg is not None
        assert sum(1 for v in rg.values() if v[0] == v[1]) == 2


def test_correction_generating_function():
    """補正 F_{j′} の母関数（続き 97、theory §12.30）：F₀ = 1、F₁ = V² + t²、F₂ = V²ψ²V + 2qt²V + t⁴ + q²t²（T(2,9)・T(2,11) の実測）、
    状態の数は C(2j′+3, 3)、係数は j′ ≤ 5 で非負。"""
    q, t = TT.q, TT.t
    V = 1 + t + q * t; pV = 1 + t**2 + q**2 * t**2
    assert TT.correction_F(0) == 1
    assert sp.expand(TT.correction_F(1) - V**2 - t**2) == 0
    assert sp.expand(TT.correction_F(2) - (V**2 * pV + 2 * q * t**2 * V + t**4 + q**2 * t**2)) == 0
    for jp in range(6):
        cs = sp.Poly(TT.correction_F(jp), TT.r, TT.s).coeffs()
        assert all(c > 0 for c in cs) and sum(cs) == sp.binomial(2 * jp + 3, 3)


def test_diag_X():
    """D_{j′} = ψ²X_{j′} の X（続き 98、theory §12.31）：X₁ = 1 + 2t + qt、X₂ = 1 + 2t + qt + 2t² + 2qt² + q²t²、状態の数 (j′+1)²、
    母関数 (1 + tz)/((1 − z)(1 − tz)(1 − qtz))。ψ²X₁ = ψ²V + t²、ψ²X₂ = ψ²V² + t⁴ − q²t²（T(2,9) の D₁, D₂）。残り (F − ψ²X)/2 は非負の整数。"""
    q, t, z = TT.q, TT.t, sp.Symbol("z")
    pV = 1 + t**2 + q**2 * t**2
    psi2 = lambda X: X.subs({TT.r: TT.r**2, TT.s: TT.s**2}, simultaneous=True)
    assert sp.expand(TT.diag_X(1) - (1 + 2 * t + q * t)) == 0 and sp.expand(TT.diag_X(2) - (1 + 2 * t + q * t + 2 * t**2 + 2 * q * t**2 + q**2 * t**2)) == 0
    assert sp.expand(psi2(TT.diag_X(1)) - pV - t**2) == 0 and sp.expand(psi2(TT.diag_X(2)) - pV**2 - t**4 + q**2 * t**2) == 0
    ser = sp.series((1 + t * z) / ((1 - z) * (1 - t * z) * (1 - q * t * z)), z, 0, 5).removeO()
    for jp in range(5):
        assert sp.expand(ser.coeff(z, jp) - TT.diag_X(jp)) == 0 and sum(sp.Poly(TT.diag_X(jp), TT.r, TT.s).coeffs()) == (jp + 1)**2
        Fd = sp.Poly(TT.correction_F(jp), TT.r, TT.s).as_dict(); Dd = sp.Poly(psi2(TT.diag_X(jp)), TT.r, TT.s).as_dict()
        assert all(Dd.get(m_, 0) <= Fd.get(m_, 0) and (Fd.get(m_, 0) - Dd.get(m_, 0)) % 2 == 0 for m_ in set(Fd) | set(Dd))


def test_correction_pairs_generating_function():
    """F = ψ²X ⊕ 2P の対の母関数（続き 99、theory §12.31 追記）：Σ P z^{j′} = t z (V′ + q t² z)/((1 − z)(1 − qtz)(1 − t²z)(1 − q²t²z))、V′ = 1 + q + qt。
    j′ ≤ 4 で F − ψ²X = 2P が成り立ち、P は非負。"""
    q, t, z = TT.q, TT.t, sp.Symbol("z")
    psi2 = lambda X: X.subs({TT.r: TT.r**2, TT.s: TT.s**2}, simultaneous=True)
    PGF = t * z * (1 + q + q * t + q * t**2 * z) / ((1 - z) * (1 - q * t * z) * (1 - t**2 * z) * (1 - q**2 * t**2 * z))
    ser = sp.series(PGF, z, 0, 5).removeO()
    for jp in range(5):
        P = sp.expand(ser.coeff(z, jp))
        assert sp.expand(TT.correction_F(jp) - psi2(TT.diag_X(jp)) - 2 * P) == 0
        assert all(c > 0 for c in sp.Poly(P, TT.r, TT.s).coeffs()) if jp else P == 0
