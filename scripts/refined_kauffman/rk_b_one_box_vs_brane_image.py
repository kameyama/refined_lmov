r"""(C) [1] の c = 1 部分（凍結した状態、σ 対称で G 展開できる）と、被覆の brane–image のセクター（Ĩ 不変な環状の弦、σ 対称）の G 展開を、
  ν^β のブロックごとに比べる（theory §11.9）。予想：[1] の凍結した状態（Möbius の帯）の被覆は brane–image の環状の弦なので、
  [1] の ν^β のブロックは brane–image の ν^{2β} のブロックと、G の次数・(𝐭² − 1) の因子・符号の型で対応する。"""
import sys, io, contextlib
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rk_k4_crosscap import V, W, G, in_G, to_VW, split_parity
import rk_b_uniform as Bu
from rcsknots.refined_kauffman.kauffman_trefoil import r as r0, s as s0, U as U0

nuk = sp.Symbol("nu", positive=True)


def blocks(P1):
    out = {}
    for b in sorted({tm.as_coeff_exponent(nuk)[1] for tm in sp.Add.make_args(P1)}):
        blk = sp.expand(P1.coeff(nuk, b))
        if blk == 0:
            continue
        cs = sp.Poly(sp.expand(blk * W**400), G, W).coeffs()
        sgn = "+" if all(c > 0 for c in cs) else "−" if all(c < 0 for c in cs) else "±"
        quo = sp.cancel(blk / (W**4 - 1))
        tf = sp.Poly(sp.fraction(quo)[1], G, W).is_monomial
        out[int(b)] = dict(deg=sp.Poly(blk, G).degree(), n=len(cs), sign=sgn, tfactor=tf)
    return out


def one_box(m):
    g = sp.cancel(Bu.g_true(m))
    E, _ = split_parity(to_VW(g))
    return blocks(sp.expand(in_G(E)))


def brane_image(m):
    import rk_b_cover_split as C
    with contextlib.redirect_stdout(io.StringIO()):
        sec = C.cover_sectors(m)
    from rcsknots.refined_kauffman.flint_rf3 import R, S
    zt_inv = (1 / S - S)
    e = ((zt_inv * sec["x"]) / (R - 1 / R)).to_sympy()
    u = sp.Symbol("u", positive=True)
    e = sp.cancel(e.subs(u, nuk * sp.Symbol("s", positive=True)))
    r_, s_ = sp.symbols("r s", positive=True)
    ev = sp.expand(e.subs({r_: V * W, s_: V / W}))
    P1 = in_G(ev)
    assert P1 is not None
    return blocks(sp.expand(P1))


if __name__ == "__main__":
    for m in (1, 2):
        a, b = one_box(m), brane_image(m)
        print(f"T(2,{2*m+1})：[1] の c=1 のブロック {a}", flush=True)
        print(f"   brane–image のブロック {b}", flush=True)
        print("   対応（[1] の β ↔ brane–image の 2β）：", {k: (a[k], b.get(2 * k)) for k in a}, flush=True)
