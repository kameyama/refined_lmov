r"""結び目の固定セクターの P = (1 − t)x と、必要な E の t = −1 での値 R = 2(1 + t) g_(2) を因数分解して比べる（progress.md 続き 34）。"""
import sys, io, contextlib; sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_b_knot_fixed_sector as K
import rk_b_cover_split as C
r, s, nu = K.r, K.s, K.nu
u = sp.Symbol("u", positive=True); t = s**2; I = sp.I
rr, nn = sp.symbols("rr nn")
m = 1
g = K.load_g(m)
with contextlib.redirect_stdout(io.StringIO()):
    sec = C.cover_sectors(m)
x = sp.cancel(sec["x"].to_sympy().subs(u, nu * s))
P = sp.cancel((1 - t) * x)
R = sp.cancel(2 * (1 + t) * g[(2,)])
ev = lambda e: sp.factor(sp.expand(sp.cancel(e.subs({r: rr, nu: nn})).subs(s, I)))
print("P(t=−1) =", ev(P))
print("R(t=−1) =", ev(R))
print("比 R/P (t=−1) =", sp.factor(sp.cancel(ev(R) / ev(P))))
