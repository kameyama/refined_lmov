r"""境界ゼロモードへの Ĩ の作用の規則の検証（theory_orientifold_index.md §11.3、progress.md 続き 44）。
  LMOV/KN：境界 h 個の曲線の境界ゼロモードは dθ_1..dθ_h（Σ dθ_i = 0）、t-brane 枠で 1 個あたり (1 − q^{−1})（規格化 q^{1/2}）。
  Ĩ は L の境界と I(L) の境界を入れ替えるので、対ごとに和 dθ+dθ' が偶、差 dθ−dθ' が奇。
  Ĩ を挿入したトレースでは奇のモードの因子が (1 + q^{−1}) に替わる。
  予言：Ĩ 不変な曲線の境界の対が n 組（境界 2n 個）なら、E の境界因子は z_q^{n−1} (z_q^+)^n、P（挿入なし）は z_q^{2n−1}。
  Hopf (2,2)・(2,4)：n = 2 → z_q (z_q^+)²。結び目の p₂ 方向：n = 1 → z_q^+（z_q 無し）。"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp

r, s, nu = sp.symbols("r s nu", positive=True)
zq, zqp = r - 1 / r, r + 1 / r


def multiplicity(e, f):
    """e を f で割り切れる回数。"""
    n = 0
    while True:
        q_, rem = sp.div(sp.Poly(sp.expand(e), r, s, nu), sp.Poly(sp.expand(f), r, s, nu))
        if not rem.is_zero:
            return n
        e = q_.as_expr(); n += 1


def hopf():
    d = pickle.load(open("data/refined_kauffman/hopf_E_fixed.pkl", "rb"))
    out = {}
    for ab in ((2, 2), (2, 4)):
        E = sp.sympify(d[ab], locals={"r": r, "s": s, "nu": nu})
        num = sp.fraction(sp.cancel(E))[0]
        out[ab] = (multiplicity(num, zq * r), multiplicity(num, zqp * r))   # r を掛けて多項式にする
    return out


if __name__ == "__main__":
    res = hopf()
    for ab, (mq, mp) in res.items():
        print(f"Hopf E_{ab}：z_q の冪 {mq}、z_q^+ の冪 {mp}（予言：1, 2）", flush=True)


def hopf_P():
    """挿入なしの固定セクターの指数 P = (1 − t) f^{cov}_{固定セクター}：予言は z_q^{2n−1}、z_q^+ 無し。"""
    import rk_k4_fixed_sector as H
    t = s**2
    out = {}
    for ab in ((2, 2), (2, 4)):
        sec = (ab[0] // 2, ab[0] // 2, ab[1] // 2, ab[1] // 2)
        P = sp.cancel((1 - t) * H.f4[sec])
        num = sp.fraction(P)[0]
        out[ab] = (multiplicity(num, zq * r), multiplicity(num, zqp * r))
    return out


def knot_P(m):
    """結び目の固定セクター（brane–image、n = 1）：P は z_q^1、E（整数計画の解、trefoil）は z_q^0 (z_q^+)^1。"""
    import io, contextlib
    import rk_b_cover_split as C
    u = sp.Symbol("u", positive=True)
    with contextlib.redirect_stdout(io.StringIO()):
        sec = C.cover_sectors(m)
    x = sp.cancel(sec["x"].to_sympy().subs(u, nu * s))
    num = sp.fraction(sp.cancel((1 - s**2) * x))[0]
    return multiplicity(num, zq * r), multiplicity(num, zqp * r)
