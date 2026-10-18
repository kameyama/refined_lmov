r"""B：クロスキャップ部分の分類（report.md §8 の判定式）を一般の式（変数 r = q^{1/2}, s = t^{1/2}, ν）に適用する。
  c = 1：V = (qt)^{1/4} の偶数冪の部分（v の整数冪）。G = z_q z_t の多項式にして ν^β ごとに
         'coherent±'（係数が同符号）／'tfactor±'（(W⁴ − 1) = (𝐭² − 1) で割れて商が同符号）／'mixed' に分類。
  c = 2：V の奇数冪の部分。X ∈ {z_q, z_t, z_q + z_t, w − 1/w} で割り、商が v ↔ 1/v 対称（G で書ける）なら同じ分類をする。
r = V W、s = V/W（W = (q/t)^{1/4}）。"""
import sys
sys.path.insert(0, "scripts")
import sympy as sp
from rk_k4_crosscap import split_parity, in_G, V, W, nu, G

r, s = sp.symbols("r s", positive=True)


def to_VW_nu(e):
    return sp.expand(sp.cancel(e).subs({r: V * W, s: V / W}))


def classify_poly(P):
    """G, W, ν の Laurent 多項式を ν^β ごとに分類。"""
    out = {}
    P = sp.expand(P)
    for beta in sorted({t.as_coeff_exponent(nu)[1] for t in sp.Add.make_args(P)}):
        blk = sp.expand(P.coeff(nu, beta))
        cs = sp.Poly(sp.expand(blk * W**120), G, W).coeffs()
        if all(c > 0 for c in cs) or all(c < 0 for c in cs):
            out[beta] = "coherent" + ("+" if cs[0] > 0 else "−"); continue
        quo = sp.cancel(blk / (W**4 - 1))
        if sp.Poly(sp.fraction(quo)[1], G, W).is_monomial:
            qc = sp.Poly(sp.expand(quo * W**120), G, W).coeffs()
            if all(c > 0 for c in qc) or all(c < 0 for c in qc):
                out[beta] = "tfactor" + ("+" if qc[0] > 0 else "−"); continue
        out[beta] = "mixed"
    return out


def classify(e):
    """{'c1': 分類 or 'G で書けない', 'c2': {X: 分類} or ...}。"""
    E, O = split_parity(to_VW_nu(e))
    res = {}
    if E == 0:
        res["c1"] = "0"
    else:
        P1 = in_G(E)
        res["c1"] = classify_poly(P1) if P1 is not None else "G で書けない（v ↔ 1/v 非対称）"
    if O == 0:
        res["c2"] = "0"
    else:
        Xs = {"z_q": V * W - 1 / (V * W), "z_t": V / W - W / V, "w−1/w": W**2 - W**-2}
        c2 = {}
        for name, X in Xs.items():
            quo = sp.cancel(O / X)
            if not sp.Poly(sp.fraction(quo)[1], V, W, nu).is_monomial:
                continue
            P2 = in_G(sp.expand(quo))
            if P2 is not None:
                c2[name] = classify_poly(P2)
        res["c2"] = c2 or "どの X でも G で書けない"
    return res
