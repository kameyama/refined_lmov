r"""(a)(b) クロスキャップの多重項の変数の探索（progress.md 続き 31）。
ν^β ごとのブロックを (r, w)（r = q^{1/2}、w = (q/t)^{1/2}、t^{1/2} = r/w）で書き、
  w → 1/w（展開変数 w + 1/w）、w → −1/w（w − 1/w）、r → −1/r（z_q）、r → 1/r（q^{1/2} + q^{−1/2}）
のどの対称性をもつかを調べる。対象：大きな N の trefoil の 2 箱の ĝ_(2) + ĝ_(1,1) と、その refined だけの部分の商、[1] の g^true（T(2,3), T(2,5)）。"""
import sys
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp
import rk_b_2box_structure as B

r, s, nu = B.r, B.s, B.nu
w = sp.Symbol("w", positive=True)
MAPS = {"w→1/w": {w: 1 / w}, "w→−1/w": {w: -1 / w}, "r→−1/r": {r: -1 / r}, "r→1/r": {r: 1 / r}}


def blocks_rw(e):
    e = sp.cancel(e.subs(s, r / w))
    n, d = sp.fraction(e)
    out = {}
    for (k, *rest), c in sp.Poly(sp.expand(n), nu, r, w).terms():
        out[k] = out.get(k, 0) + c * r**rest[0] * w**rest[1]
    return {k: sp.cancel(v * nu**k / d) for k, v in out.items()}


def report(lab, e):
    bl = blocks_rw(e)
    res = {}
    for b, v in sorted(bl.items(), key=lambda kv: sp.degree(sp.fraction(kv[1])[0], nu) - sp.degree(sp.fraction(kv[1])[1], nu)):
        beta = sp.Poly(sp.fraction(sp.cancel(v))[0], nu).monoms()[0][0] - sp.Poly(sp.fraction(sp.cancel(v))[1], nu).monoms()[0][0]
        sy = []
        for name, mp in MAPS.items():
            img = sp.cancel(v.subs(mp, simultaneous=True))
            rat = sp.cancel(img / v)
            n_, d_ = sp.fraction(sp.factor(rat))
            if sp.Poly(n_, r, w).is_monomial and sp.Poly(d_, r, w).is_monomial:
                sy.append(name if rat == 1 else f"{name}×{rat}")
        res[beta] = sy or "なし"
    print(f"{lab}：ν^β ごとの対称性（単項式倍まで）{res}", flush=True)


if __name__ == "__main__":
    import rk_b_uniform as Bu
    from rcsknots.kauffman_trefoil import r as r0, s as s0, U as U0
    G = B.ghat()
    report("2 箱 ĝ_(2)+ĝ_(1,1)", G["sum"])
    _, Q = B.refined_only(G["sum"])
    report("2 箱 和の商", Q)
    for m in (1, 2):
        g = sp.cancel(Bu.g_true(m).subs({r0: r, s0: s}).subs(U0, nu * s))
        report(f"[1] T(2,{2*m+1}) g^true", g)
