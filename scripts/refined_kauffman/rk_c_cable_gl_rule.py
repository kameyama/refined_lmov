r"""ケーブル K1・K3 の GL_N の DAHA-Jones の規格化の規則を、普通の表現で確かめる（progress.md 続き 141）。
一様な規格化の W_b(N) = (t/q)^{δ|b|/2} × framed_knot の GN の値（鏡像、a^{1/2} = t^{N/2}）× P_b(t^ρ) に対し
  M_b := Z_b / 鏡像(W_b)、Z_b = {γ̂_comp(γ̂_pat(P_b)·1)}（daha_iter.iter_jones_gl）
が M_b = T_b^{−F} X(N)^{|b|}（T_b = q^{(b,b)/2} t^{(b,ρ)}、ρ = (N−1, …, 0)、F = p₂q₂）の形かを見る。
|b| = 0 の随伴では X が効かないので、この形が成り立てば随伴の規格化は T_adj^{−F} に決まる（rk_c_cable_lp.py gl）。"""
import sys
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rcsknots.refined_kauffman.daha_iter import iter_jones_gl
import rk_c_2box_largeN as L2
import rk_torus_tag as TG

r, s, u = sp.symbols("r s u", positive=True)
mir = lambda e: e.subs({r: 1 / r, s: 1 / s}, simultaneous=True)
LAMS = ((1,), (2,), (1, 1), (2, 1), (3,), (1, 1, 1))

if __name__ == "__main__":
    ok = True
    for cable in ("K1", "K3"):
        F, delta = TG.CABLE[cable]
        for N in (2, 3, 4):
            X = None
            for lam in LAMS:
                if len(lam) > N or (N == 4 and sum(lam) == 3):
                    continue
                Z, _ = iter_jones_gl(N, lam, cable)
                W = (L2.S / L2.R)**(delta * sum(lam)) * L2.Rdata_nu("".join(map(str, lam)), cable) * L2.Psym(lam)
                W = sp.cancel(W.to_sympy().subs(u, s**(N + 1)))                 # ν = U/s = s^N
                M = sp.factor(sp.cancel(Z.to_sympy() / mir(W)))
                b = list(lam) + [0] * (N - len(lam))
                bb = sum(x * x for x in b); br = sum(x * (N - 1 - i) for i, x in enumerate(b))
                rest = sp.factor(M / (r**(-F * bb) * s**(-2 * F * br)))
                n_, d_ = sp.fraction(rest)
                mono = sp.Poly(n_, r, s).is_monomial and sp.Poly(d_, r, s).is_monomial
                Xb = sp.simplify(rest**sp.Rational(1, sum(lam))) if mono else None
                if X is None:
                    X = Xb
                same = mono and sp.simplify(rest - X**sum(lam)) == 0
                ok &= same
                print(f"{cable} N={N} λ={lam}: M/T^(-F) = {rest}、X^|b| の形 {same}", flush=True)
            print(f"   {cable} N={N}：X(N) = {X}", flush=True)
    print("すべて T_b^{-F} X(N)^{|b|} の形：", ok)
