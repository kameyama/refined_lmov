r"""(a) の続き：固定セクター P = (1 − t)x = t^{1/2} z_q Pk、Pk = (t^{−1/2} − t^{1/2}) x / z_q（σ 対称、G 展開できる）を
Pk = Σ c_{β,g,j} ν^β G^g W^j と書き、ゼロモードの (−1)^F を和に替えた E = t^{1/2} z̃ Σ ε_β c_{β,g,j} ν^β (G^+)^g W^j
（G^+ = z_q^+ z_t^+、z̃ ∈ {z_q, z_q^+}、ε_β ∈ {1, (−1)^{β/2}, (−1)^{β}, …}）を t = −1 で R と比べる（progress.md 続き 36）。t = −1 では z_t^+ = 0 なので g = 0 の項だけが効く。"""
import sys, io, contextlib, itertools
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp
import rk_b_knot_fixed_sector as K
import rk_b_cover_split as C
from rk_k4_crosscap import V, W, G, in_G

r, s, nu = K.r, K.s, K.nu
u = sp.Symbol("u", positive=True)
t = s**2
I = sp.I
nuk = sp.Symbol("nu")


def main(m):
    g = K.load_g(m)
    with contextlib.redirect_stdout(io.StringIO()):
        sec = C.cover_sectors(m)
    x = sp.cancel(sec["x"].to_sympy().subs(u, nu * s))
    Pk = sp.cancel((1 / s - s) * x / (r - 1 / r))
    ev = sp.expand(sp.cancel(Pk.subs(nu, nuk)).subs({r: V * W, s: V / W}))
    P1 = sp.expand(in_G(ev))
    g0 = sp.expand(P1.subs(G, 0))                                   # g = 0 の項（ν, W の Laurent 多項式）
    R = sp.cancel(2 * (1 + t) * g[(2,)])
    Rm1 = sp.expand(sp.cancel(R.subs(s, I)))
    hits = []
    for zt_lab, zt in (("z_q", r - 1 / r), ("z_q^+", r + 1 / r)):
        for elab, eps in (("1", lambda b: 1), ("(−1)^{β/2}", lambda b: (-1)**(b // 2)), ("(−1)^β", lambda b: (-1)**b), ("i^β", lambda b: I**b)):
            terms = sp.Add.make_args(g0)
            Etw = sum(eps(int(sp.Poly(tm, nuk).degree()) if tm.has(nuk) else 0) * tm for tm in terms) if False else 0
            # ν の冪ごとに符号を付ける
            Etw = 0
            for tm in terms:
                b = sp.degree(tm, nuk) if tm.has(nuk) and not tm.as_coeff_exponent(nuk)[1] < 0 else tm.as_coeff_exponent(nuk)[1]
                Etw += eps(int(b)) * tm
            # W = (q/t)^{1/4}、V = (qt)^{1/4} に戻し、t = −1（s = i）で評価
            E_m1 = sp.expand(sp.cancel((s * zt * Etw).subs(nuk, nu).subs({W: sp.sqrt(r / s)})).subs(s, I))
            for c in (1, -1, I, -I):
                if sp.expand(c * E_m1 - Rm1) == 0:
                    hits.append((zt_lab, elab, c))
    print(f"T(2,{2*m+1})：t = −1 で R に一致する捻りの規則 (z̃, ε_β, 全体の係数)：{hits}", flush=True)
    print(f"   g = 0 の項 = {str(sp.factor(g0))[:300]}", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]))
