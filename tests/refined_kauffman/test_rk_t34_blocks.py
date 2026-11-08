"""T(3,4)（γ = τ₊τ₋³、k = 3）の n = 1 の N のブロック（続き 102、theory §12.32）の機能テスト。
  下の系列 k ≤ 5 は Sym^k V（V = 1 + t + qt）、移り目 k = 6 は Sym⁶V − S̃²、S̃ = (1 + qt)(1 + t + t² + qt² + q²t²) = Sym³V + (q − t)t²
  （T(3,4) の簡約 superpolynomial（11 状態）ではない）。上の端 k = 9 は −t⁸(q + t)。K = 9 = rs − 3。"""
import io, contextlib, sys
import sympy as sp
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import rk_f_twisted_tensor as TT
import rk_f_frozen_basis as FB

r, s, q, t = TT.r, TT.s, TT.q, TT.t


def test_t34_blocks():
    with contextlib.redirect_stdout(io.StringIO()):
        N, _ = TT._knot_N(1, 3)
    blocks = FB.nu_blocks(N)
    nk = {(-b - 6) // 2: sp.expand(blk * r**(4 + 2 * ((-b - 6) // 2))) for b, blk in blocks.items()}
    assert sorted(nk) == list(range(10))
    cells = lambda n: [(a_, b_, c_) for a_ in range(n + 1) for b_ in range(n + 1 - a_) for c_ in [n - a_ - b_]]
    sym = lambda n: sum(t**(b_ + c_) * q**c_ for a_, b_, c_ in cells(n))
    for k in range(6):
        assert sp.expand(nk[k] - sym(k)) == 0
    St = (1 + q * t) * (1 + t + t**2 + q * t**2 + q**2 * t**2)
    assert sp.expand(St - sym(3) - (q - t) * t**2) == 0
    assert sp.expand(nk[6] - (sym(6) - St**2)) == 0
    assert sp.expand(nk[9] + t**8 * (q + t)) == 0
    S = TT.superpoly_data(1, 3)
    assert sp.expand(S - St) != 0 and sum(sp.Poly(S, r, s).coeffs()) == 11 and sum(sp.Poly(sp.expand(St), r, s).coeffs()) == 10
