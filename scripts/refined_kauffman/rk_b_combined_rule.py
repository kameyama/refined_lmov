r"""(B) 組み合わせた規則の走査（theory §11.9）：固定セクター P = s z_q Pk（n = 1）、Pk = Σ c ν^β G^g W^j（G 展開）。
  E の候補 = s z_q^+ Σ ε(β, j) c ν^β X^g W^j。ε(β, j) = (−1)^{aβ/2 + b j/2}（a, b ∈ {0,1}；β, j は偶数）、X は取っ手の捻り。
  KN の仮定（S_R が BPS 状態に自明）のもとで Ĩ = ε_X e^{iπJ₂} が「t の指数と ν の指数の偶奇」で決まるなら、この族に入るはず。"""
import sys, io, contextlib, itertools
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_b_twisted_multiplet as M
import rk_b_sign_rule as SR
from rk_k4_crosscap import V, W, G

r, s, nu = M.r, M.s, M.nu


def run(m):
    import rk_b_knot_fixed_sector as K
    import rk_b_cover_split as C
    u = sp.Symbol("u", positive=True)
    g = {k: v.subs({K.r: r, K.s: s, K.nu: nu}) for k, v in K.load_g(m).items()}
    with contextlib.redirect_stdout(io.StringIO()):
        sec = C.cover_sectors(m)
    P = sp.cancel((1 - s**2) * sec["x"].to_sympy().subs(u, nu * s))
    Pk = sp.cancel(P / (s * (r - 1 / r)))
    mult = M.multiplets(Pk)
    zq, zt, zqp, ztp = r - 1 / r, s - 1 / s, r + 1 / r, s + 1 / s
    V2 = r * s
    handles = {"G": zq * zt, "Hq": zq, "Ht": zt, "Havg": (zq + zt) / 2, "Hqt": V2 - 1 / V2, "Hdiff": zqp - ztp, "A": zqp * zt, "B": zq * ztp, "Hplus": zqp * ztp}
    hits = []
    for (a, b), (xl, X) in itertools.product(itertools.product((0, 1), repeat=2), handles.items()):
        E = 0
        for (beta, gdeg, j), c in mult.items():
            eps = (-1)**((a * beta + b * j) // 2) if (a * beta + b * j) % 2 == 0 else None
            if eps is None:
                E = None; break
            E += eps * c * nu**beta * X**gdeg * (r / s)**(j // 2)
        if E is None:
            continue
        E = sp.expand(s * zqp * E)
        for sg in (1, -1):
            ci, cii = SR.check_conditions(sg * E, g[(2,)], g[(1,)])
            if ci or cii:
                hits.append((sg, (a, b), xl, ci, cii))
    print(f"T(2,{2*m+1})：多重項 {len(mult)}。(i) か (ii) を満たす (符号, (a,b), X, (i), (ii)) = {hits}", flush=True)


if __name__ == "__main__":
    for m in (1, 2):
        run(m)
