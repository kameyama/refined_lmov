r"""(c) 大きな N の被覆を brane（Y）と image（Y'）の 2 組の probe に分けた生成関数と M ⊗ M の hatted（progress.md 続き 31）。
  Z = 1 − W_1 (p_1 + p_1') + [W_(2) e_2(Y) + W_(1,1) P_(2)(Y; t, q) + (Y → Y')] + (W_adj + κ) p_1(Y) p_1(Y')（双対 Cauchy、W は一様な規格化、荷電共役で brane と image は同じ W）。
  Y = Y' で rk_c_2box_largeN.cover_f に戻る（main で確認）。M ⊗ M の hatted は冪和の成分を φ_ρ φ_σ で割る：φ_(1) = 1、φ_(1,1) = z_q、φ_(2) = q^{1/2} + q^{−1/2}。
  クロスキャップは、SO の p_1² 成分を z_q で、被覆の brane–image の部分を 1 で割る版 ĝ^split も作る。変数は (r, s, U)（rcsknots.flint_rf3）。"""
import sys
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import pickle
import sympy as sp
from fractions import Fraction
from rcsknots.flint_rf3 import RF, R, S, U, from_sympy, conv
import rk_c_2box_largeN as L

q, t = R * R, S * S
half = conv(Fraction(1, 2))


def cover_sectors():
    w = S / R
    Wc = {lam: w**sum(lam) * L.Rdata_nu("".join(map(str, lam))) * L.Psym(lam) for lam in ((1,), (1, 1), (2,))}
    radj = pickle.load(open("data/refined_kauffman/daha/Radj_A_m1.pkl", "rb"))
    r_, s_ = sp.symbols("r s", positive=True)
    nu = L.nu
    Radj = sum((from_sympy(sp.sympify(v, locals={"r": r_, "s": s_})) * nu**j for j, v in radj.items()), conv(0))
    A = nu * nu
    qn = (nu - 1 / nu) / (S - 1 / S)
    kappa = (1 - A) * (1 - q) / ((1 - t) * (1 - q * A / t))
    Wadj = Radj * (qn * qn - kappa)
    cp = (1 + t) * (1 - q) / (1 - t * q)                        # P_(2)(Y; t, q) = m_2 + c' m_11
    W1 = Wc[(1,)]
    # 冪和の係数（片側）：e_2 = (p1² − p2)/2、P_(2)(Y;t,q) = p2 + c'(p1² − p2)/2
    Lp11 = half * Wc[(2,)] + half * cp * Wc[(1, 1)] - half * W1 * W1
    Lp2 = -half * Wc[(2,)] + (1 - half * cp) * Wc[(1, 1)]
    Lx = Wadj + kappa - W1 * W1                                 # p1(Y) p1(Y') の係数（brane–image）
    f1 = -W1
    fp2 = Lp2 - half * f1.adams(2)
    return {"f1": f1, "p11": Lp11, "p2": fp2, "x": Lx}


def zq():
    return R - 1 / R


def main():
    sec = cover_sectors()
    fc = L.cover_f()
    # Y = Y' での検算：p1² 成分 = 2 Lp11 + Lx、p2 成分 = 2 fp2
    p11_single = (fc[(2,)] + fc[(1, 1)]) * half
    p2_single = (fc[(2,)] - fc[(1, 1)]) * half
    ok = ((p11_single - (sec["p11"] * 2 + sec["x"])).is_zero(), (p2_single - sec["p2"] * 2).is_zero(), (fc[(1,)] - sec["f1"] * 2).is_zero())
    print("Y = Y' で cover_f に一致（p1², p2, p1）:", ok, flush=True)
    hat = {"brane (2)": sec["p11"] / zq() + sec["p2"] / (R + 1 / R), "brane (1,1)": sec["p11"] / zq() - sec["p2"] / (R + 1 / R),
           "brane–image (1),(1)": sec["x"]}
    for k, v in hat.items():
        print(f"被覆の M⊗M hatted f̂_{k}：Laurent {v.is_laurent()}、整数 Laurent {v.is_integer_laurent()}", flush=True)
    pickle.dump({k: str(v.to_sympy()) for k, v in hat.items()}, open("data/refined_kauffman/daha/largeN_2box_m1_cover_split.pkl", "wb"))
    return hat, sec




def kn_check(hat):
    """KN の t 枠 (t^{-1/2} − t^{1/2}) f̂ にして、分母・σ 対称性・G 展開後の ν^β ごとの符号を見る。"""
    from rk_k4_crosscap import V, W, G, in_G
    r_, s_, u_ = sp.symbols("r s u", positive=True)
    nuk = sp.Symbol("nu")
    out = {}
    for k, v in hat.items():
        e = (1 / S - S) * v
        es = sp.cancel(e.to_sympy())
        den = sp.factor(sp.fraction(es)[1])
        ev = sp.expand(sp.cancel(es.subs(u_, nuk * s_)).subs({r_: V * W, s_: V / W}))
        P1 = in_G(ev) if e.is_laurent() else None
        if P1 is None:
            print(f"(1/s − s) f̂_{k}：分母 {den}、G 展開できない", flush=True)
            continue
        P1 = sp.expand(P1)
        res = {}
        for b in sorted({tm.as_coeff_exponent(nuk)[1] for tm in sp.Add.make_args(P1)}):
            cs = sp.Poly(sp.expand(sp.expand(P1.coeff(nuk, b)) * W**200), G, W).terms()
            res[b] = "そろう" if len({sp.sign(c) for _, c in cs}) == 1 else "混在"
        out[k] = res
        print(f"(1/s − s) f̂_{k}：分母 {den}、G 展開できる。ν^β ごとの符号：{res}", flush=True)
    return out


if __name__ == "__main__":
    hat, sec = main()
    kn_check(hat)
