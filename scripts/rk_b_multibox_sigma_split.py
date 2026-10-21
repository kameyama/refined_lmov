r"""2 箱の refined c = 1（t = q で z 偶の部分 = hatted の ν 奇）を σ 偶部分 S と σ 奇部分 O に分け、
  O = (q − t)(qt − 1) K（§11.16 で確認した因子）として、S と K の整数性と G 展開のブロックの符号を見る（progress.md 続き 55）。
  σ：V → 1/V（r = VW, s = V/W）。G = V² + V⁻² − W² − W⁻²。"""
import sys
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp
from rk_b_multibox_sigma import load, sigma, c1, nterms, r_, s_, nu
from rk_k4_crosscap import V, W, G, in_G

c2 = lambda e: sp.cancel(e - c1(e))
nuk = sp.Symbol("nu", positive=True)


def to_vw(e):
    return sp.expand(sp.cancel(e).subs({r_: V * W, s_: V / W}, simultaneous=True))


def is_int_laurent(e):
    """整数係数の Laurent 多項式か（分母は係数 1 の単項式でなければならない。2 が残れば半整数）。"""
    n_, d_ = sp.fraction(sp.cancel(e))
    Pd = sp.Poly(d_, r_, s_, nu)
    if not Pd.is_monomial:
        return False
    cd = Pd.coeffs()[0]
    return all((c / cd).is_integer for c in sp.Poly(sp.expand(n_), r_, s_, nu).coeffs())


def blocks(P1):
    out = {}
    for b in sorted({tm.as_coeff_exponent(nuk)[1] for tm in sp.Add.make_args(P1)}):
        blk = sp.expand(P1.coeff(nuk, b))
        if blk == 0:
            continue
        cs = sp.Poly(sp.expand(blk * W**400), G, W).coeffs()
        sgn = "+" if all(c > 0 for c in cs) else "−" if all(c < 0 for c in cs) else "±"
        out[int(b)] = (sp.Poly(blk, G).degree(), len(cs), sgn)
    return out


if __name__ == "__main__":
    which = sys.argv[1:] or ["sum1", "diff1", "sum2", "diff2"]
    for m in (1, 2):
        g, f = load(m); g2, g11 = g["(2,)"], g["(1, 1)"]
        zq = r_ - 1 / r_
        print(f"== T(2,{2*m+1})", flush=True)
        cases = []
        if f"sum{m}" in which:
            cases.append(("ĝ2+ĝ11 の ν 奇（c=1）", c2((g2 + g11) / zq)))
        if f"diff{m}" in which:
            cases.append(("(1+t)(ĝ2−ĝ11) の ν 奇（c=1）", c2((1 + s_**2) * (g2 - g11) / (r_ + 1 / r_))))
        if f"even{m}" in which:
            cases.append(("ĝ2+ĝ11 の ν 偶（c=2）", c1((g2 + g11) / zq)))
        for lab, X in cases:
            S = sp.cancel((X + sigma(X)) / 2); O = sp.cancel((X - sigma(X)) / 2)
            K = sp.cancel(O / ((r_**2 - s_**2) * (r_**2 * s_**2 - 1)))
            print(f"  {lab}: S 整数 Laurent {is_int_laurent(S)}, O/((q−t)(qt−1)) 整数 Laurent {is_int_laurent(K)}, K の分母 {sp.factor(sp.fraction(sp.cancel(K))[1])}", flush=True)
            e = to_vw(S); P1 = in_G(e)
            print(f"     S の G 展開のブロック（β: (G 次数, 項数, 符号)）: {blocks(P1) if P1 is not None else '不可'}", flush=True)
            # K は σ 奇の因子を除いた後なので σ 対称：(qt − 1) = V⁴ − 1 → (V² − V⁻²) V²、(q − t) = V²(W² − W⁻²)
            Kp = sp.cancel(O / ((V**2 - V**-2) * (W**2 - W**-2)).subs({V: sp.sqrt(r_ * s_), W: sp.sqrt(r_ / s_)}))
            e2 = to_vw(Kp); P2 = in_G(e2)
            print(f"     O/((V²−V⁻²)(W²−W⁻²)) の G 展開のブロック: {blocks(P2) if P2 is not None else '不可（σ 非対称）'}", flush=True)
