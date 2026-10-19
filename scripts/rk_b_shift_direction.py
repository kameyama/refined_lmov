r"""述べ直しの判定式（ĝ − ĝ|_{t=q}(ν → ν w^k) が (q − t)(ν² − t) で割り切れる、w = (q/t)^{1/2}）の、ずれの向き k を調べる（progress.md 続き 38）。
閉弦（AgSc §6、scripts/rk_b_closed_orientifold.py）では、凍結した RP² の電荷 x と向き付けられた状態の電荷 Q0 の比は x²/Q0 = w^{−1}。
開弦では、被覆（½h）の Kähler 変数は 𝐚_U² = ν² w^{−1}。凍結した状態の ν のずれ w^k と比べる。"""
import sys
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp
import rk_b_2box_structure as B

r, s, nu = B.r, B.s, B.nu


def divisible(e, k):
    X = sp.cancel(e - sp.cancel(e.subs(s, r)).subs(nu, nu * (r / s)**k))
    n_ = sp.Poly(sp.fraction(X)[0], nu, r)
    return {f"{f}": sp.rem(n_, sp.Poly(f, nu, r)).is_zero for f in (r**2 - s**2, nu**2 - s**2)}


def run(cases, ks=(-3, -2, -1, 0, 1, 2, 3)):
    for lab, e in cases:
        res = {k: all(divisible(e, k).values()) for k in ks}
        print(f"{lab}：割り切れる k = {[k for k, v in res.items() if v]}", flush=True)


if __name__ == "__main__":
    import rk_b_uniform as Bu
    import rk_b_2box_structure_m as Bm
    import rk_b_3box_structure as B3
    from rcsknots.kauffman_trefoil import r as r0, s as s0, U as U0
    cases = [(f"[1] T(2,{2*m+1})", sp.cancel(Bu.g_true(m).subs({r0: r, s0: s}).subs(U0, nu * s))) for m in (1, 2)]
    for m in (1, 2):
        h = Bm.hat(Bm.load(m)[2])
        cases += [(f"2 箱 T(2,{2*m+1}) ĝ_{k}", h[k]) for k in ((2,), (1, 1), "sum")]
    gh3 = B3.schur(B3.load("g"), 3, hat=True)
    cases += [(f"3 箱 trefoil ĝ_{R_}", gh3[R_]) for R_ in ((3,), (2, 1), (1, 1, 1))]
    run(cases)
