r"""K4'：refined の交差帽展開。g（Laurent 多項式）を
  g = P_1(G, w, ν) + X · P_2(G, w, ν)、G = (q^{1/2}−q^{−1/2})(t^{1/2}−t^{−1/2}) = v + 1/v − w − 1/w、v = (qt)^{1/2}、w = (q/t)^{1/2}
に分解する。v の整数冪の部分（v ↔ 1/v 対称なら G の多項式に書ける）を c=1、残りを X で割って c=2 とする。
対象：unknot の g = −1、Hopf の g_{11} = z_q(ν~ − 1/ν~)、trefoil [1] の g^{ref}（K3' の仮説の単項式）。
"""
import sys
sys.path.insert(0, "src")
import sympy as sp
from rcsknots.refined_kauffman.kauffman_trefoil import g_ref, r, s, U

V, W, nu, G = sp.symbols("V W nu G", positive=True)   # V = (qt)^{1/4}, W = (q/t)^{1/4}：r = V W, s = V/W


def to_VW(expr):
    """(r, s, U) → (V, W, ν)：r = V W、s = V/W、U = ν s。"""
    return sp.expand(expr.subs(U, nu * s).subs({r: V * W, s: V / W}))


def split_parity(e):
    """V の指数が偶数（v = V² の整数冪）の部分と奇数の部分に分ける。"""
    E = O = 0
    for term in sp.Add.make_args(e):
        if term.as_coeff_exponent(V)[1] % 2 == 0:
            E += term
        else:
            O += term
    return sp.expand(E), sp.expand(O)


def in_G(e):
    """V ↔ 1/V 対称な Laurent 多項式（V の偶数冪）を G = V² + V^{−2} − W² − W^{−2} の多項式にする。非対称なら None。"""
    if sp.expand(e - e.subs(V, 1 / V)) != 0:
        return None
    out = 0
    rest = sp.expand(e)
    while rest != 0:
        d = max(term.as_coeff_exponent(V)[1] for term in sp.Add.make_args(rest))
        if d == 0:
            out += rest
            break
        c = sp.expand(rest).coeff(V, d)
        out += c * G**(d // 2)
        rest = sp.expand(rest - c * (V**2 + V**-2 - W**2 - W**-2)**(d // 2))
    return sp.expand(out)


def crosscap_expand(g_rsU, name):
    e = to_VW(g_rsU)
    E, O = split_parity(e)
    print(f"== {name}")
    P1 = in_G(E)
    print("  c=1 部分（G, w, ν）:", P1 if P1 is not None else "v ↔ 1/v 非対称：G で書けない")
    if O == 0:
        print("  c=2 部分: 0"); return
    for xname, X in (("z_q", V * W - 1 / (V * W)), ("z_t", V / W - W / V), ("w−1/w", W**2 - W**-2)):
        quo = sp.cancel(O / X)
        num, den = sp.fraction(quo)
        if not sp.Poly(den, V, W, nu).is_monomial:
            print(f"  X = {xname}: 割り切れない"); continue
        quo = sp.expand(quo)
        P2 = in_G(quo)
        print(f"  X = {xname}: 商 = {quo}\n      G 展開: {P2 if P2 is not None else 'v ↔ 1/v 非対称'}")


if __name__ == "__main__":
    nu_t = U * r / s**2                       # ν~ = U q^{1/2}/t
    zq = r - 1 / r
    crosscap_expand(sp.Integer(-1), "unknot g_10 = −1")
    crosscap_expand(zq * (nu_t - 1 / nu_t), "Hopf g_11 = z_q(ν~ − 1/ν~)")
    crosscap_expand(g_ref(), "trefoil g^ref_[1]（仮説の単項式 T^{-6}(q/t)^{1/2}）")
    Q, T, A = sp.symbols("Q T A")
    print("== trefoil g^ref_[1] を DGR 変数で（𝐪=Q, 𝐭=T, 𝐚=A）:",
          sp.expand(sp.cancel(g_ref().subs({r: -T * Q, U: A * Q}).subs(s, Q))))
