r"""凍結した状態（クロスキャップ部分 C）の展開の仮定の検査（progress.md 続き 37、docs/refined_kauffman/theory_orientifold_index.md）。
  C = Σ_β ν^β [A_β(G, w) + X · B_β(G, w)]、A, B ∈ Z[G, w^{±}]（σ 対称）、X は c = 2 の余分の因子（t = q で z）。
  σ（q → t^{−1}, t → q^{−1}、ν 固定）で X ≠ σX なら B = (e − σe)/(X − σX)、A = e − X B が一意に決まる。
  A, B が Laurent で G 展開できるか、G 展開の係数の符号が ν^β ごとにそろうかを調べる。
  対象：[1] の g^true（T(2,3), T(2,5)、一様な規格化）、2 箱の ĝ_(2) + ĝ_(1,1)（trefoil, T(2,5)、E を含まない）、3 箱の ĝ_R（trefoil、E = 0）。"""
import sys
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp
import rk_b_2box_structure as B
from rk_k4_crosscap import V, W, G


def in_G(e, maxit=200):
    """V ↔ 1/V 対称で V の偶数冪だけの Laurent 多項式を G = V² + V^{−2} − W² − W^{−2} の多項式に。できなければ None。"""
    e = sp.expand(e)
    if e == 0:
        return sp.Integer(0)
    if sp.expand(e - e.subs(V, 1 / V)) != 0:
        return None
    if any(tm.as_coeff_exponent(V)[1] % 2 for tm in sp.Add.make_args(e)):
        return None
    out, rest = sp.Integer(0), e
    for _ in range(maxit):
        if rest == 0:
            return sp.expand(out)
        d = max(tm.as_coeff_exponent(V)[1] for tm in sp.Add.make_args(rest))
        if d == 0:
            return sp.expand(out + rest)
        c = rest.coeff(V, d)
        out += c * G**(d // 2)
        rest = sp.expand(rest - c * (V**2 + V**-2 - W**2 - W**-2)**(d // 2))
    return None

r, s, nu = B.r, B.s, B.nu
nuk = sp.Symbol("nu")
zq, zt = r - 1 / r, s - 1 / s
XS = {"z_q": zq, "z_t": zt, "(z_q−z_t)/2": (zq - zt) / 2, "(z_q+z_t)/2": (zq + zt) / 2}


def to_VW(e):
    return sp.expand(sp.cancel(e.subs(nu, nuk)).subs({r: V * W, s: V / W}))


def blocks_sign(P1):
    out = {}
    for b in sorted({tm.as_coeff_exponent(nuk)[1] for tm in sp.Add.make_args(P1)}):
        blk = sp.expand(P1.coeff(nuk, b))
        if blk == 0:
            continue
        cs = [c for _, c in sp.Poly(sp.expand(blk * W**400), G, W).terms()]
        out[b] = "+" if all(c > 0 for c in cs) else "−" if all(c < 0 for c in cs) else "±"
    return out


def test(lab, e):
    res = {}
    for xl, X in XS.items():
        sX = B.sigma(X)
        Bp = sp.cancel((e - B.sigma(e)) / (X - sX))
        Ap = sp.cancel(e - X * Bp)
        lau = all(sp.Poly(sp.fraction(v)[1], r, s, nu).is_monomial for v in (Ap, Bp))
        if not lau:
            res[xl] = "Laurent でない"; continue
        PA, PB = in_G(to_VW(Ap)), in_G(to_VW(Bp))
        if PA is None or PB is None:
            res[xl] = f"G 展開できない（A {PA is not None}, B {PB is not None}）"; continue
        sa, sb = blocks_sign(sp.expand(PA)), blocks_sign(sp.expand(PB))
        res[xl] = f"A の符号 {sum(v != '±' for v in sa.values())}/{len(sa)}、B の符号 {sum(v != '±' for v in sb.values())}/{len(sb)}"
    print(f"{lab}：{res}", flush=True)


if __name__ == "__main__":
    import rk_b_uniform as Bu
    import rk_b_2box_structure_m as Bm
    import rk_b_3box_structure as B3
    from rcsknots.kauffman_trefoil import r as r0, s as s0, U as U0
    for m in (1, 2):
        test(f"[1] T(2,{2*m+1}) g^true", sp.cancel(Bu.g_true(m).subs({r0: r, s0: s}).subs(U0, nu * s)))
    for m in (1, 2):
        f, fc, g = Bm.load(m)
        test(f"2 箱 T(2,{2*m+1}) ĝ_(2)+ĝ_(1,1)", Bm.hat(g)["sum"])
    g3 = B3.load("g")
    gh3 = B3.schur(g3, 3, hat=True)
    for R_ in ((3,), (2, 1), (1, 1, 1)):
        test(f"3 箱 trefoil ĝ_{R_}", gh3[R_])


def test_parity(lab, e):
    """e を V の偶奇（c = 1 と c = 2）に分け、e_even と e_odd/X が G 展開できるかを見る。"""
    ev = to_VW(e)
    even = sp.expand(sum(tm for tm in sp.Add.make_args(ev) if tm.as_coeff_exponent(V)[1] % 2 == 0))
    odd = sp.expand(ev - even)
    PE = in_G(even)
    res = {"c=1": "G 展開できる" + (f"、符号 {sum(v != '±' for v in blocks_sign(PE).values())}/{len(blocks_sign(PE))}" if PE not in (None, 0) else "") if PE is not None else "できない"}
    for xl, X in XS.items():
        Xv = to_VW(X)
        quo = sp.cancel(odd / Xv)
        if not sp.Poly(sp.fraction(quo)[1], V, W, nuk).is_monomial:
            res[f"c=2 ÷ {xl}"] = "割り切れない"; continue
        PO = in_G(sp.expand(quo))
        res[f"c=2 ÷ {xl}"] = "G 展開できない" if PO is None else f"G 展開できる、符号 {sum(v != '±' for v in blocks_sign(PO).values())}/{len(blocks_sign(PO))}"
    print(f"{lab}：{res}", flush=True)


def main_parity():
    import rk_b_uniform as Bu
    import rk_b_2box_structure_m as Bm
    import rk_b_3box_structure as B3
    from rcsknots.kauffman_trefoil import r as r0, s as s0, U as U0
    for m in (1, 2):
        test_parity(f"[1] T(2,{2*m+1}) g^true", sp.cancel(Bu.g_true(m).subs({r0: r, s0: s}).subs(U0, nu * s)))
    for m in (1, 2):
        f, fc, g = Bm.load(m)
        test_parity(f"2 箱 T(2,{2*m+1}) ĝ_(2)+ĝ_(1,1)", Bm.hat(g)["sum"])
    gh3 = B3.schur(B3.load("g"), 3, hat=True)
    for R_ in ((3,), (2, 1), (1, 1, 1)):
        test_parity(f"3 箱 trefoil ĝ_{R_}", gh3[R_])
