r"""SO と Sp の 2 箱のクロスキャップ部分の関係の走査（progress.md 続き 52）。変数 (r, s, ν)。
  候補：M ∈ {id, σ (r→1/s, s→1/r), 鏡像 (r→1/r, s→1/s), q↔t (r↔s)} × ν ∈ {ν, −ν, 1/ν, −1/ν} × 転置 {なし, あり} × 符号 ±、単項式倍を許す。
  対象：g_R、hatted の p1² 成分 (g2+g11)/z_q、p2 成分 (g2−g11)/z_q^+、f_R（SO 側そのもの）。"""
import sys, pickle, itertools
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp

r_, s_, u_ = sp.symbols("r s u", positive=True)
nu = sp.Symbol("nu", positive=True)
loc = {"r": r_, "s": s_, "u": u_}
OUT = "data/refined_kauffman/daha"


def load():
    sp_d = pickle.load(open(f"{OUT}/largeN_2box_sp_m1.pkl", "rb"))
    so_d = pickle.load(open(f"{OUT}/largeN_2box_m1.pkl", "rb"))
    g_sp = {k: sp.cancel(sp.sympify(v, locals=loc).subs(u_, nu / s_)) for k, v in sp_d["g"].items()}
    g_so = {k: sp.cancel(sp.sympify(v, locals=loc).subs(u_, nu * s_)) for k, v in so_d["g"].items()}
    f_sp = {k: sp.cancel(sp.sympify(v, locals=loc).subs(u_, nu / s_)) for k, v in sp_d["f"].items()}
    f_so = {k: sp.cancel(sp.sympify(v, locals=loc).subs(u_, nu * s_)) for k, v in so_d["f"].items()}
    return g_sp, g_so, f_sp, f_so


MAPS = {"id": {}, "σ": {r_: 1 / s_, s_: 1 / r_}, "鏡像": {r_: 1 / r_, s_: 1 / s_}, "q↔t": {r_: s_, s_: r_}}
NUS = {"ν": {}, "−ν": {nu: -nu}, "1/ν": {nu: 1 / nu}, "−1/ν": {nu: -1 / nu}}


def is_mono(e):
    n_, d_ = sp.fraction(sp.cancel(e))
    try:
        return sp.Poly(n_, r_, s_, nu).is_monomial and sp.Poly(d_, r_, s_, nu).is_monomial
    except sp.PolynomialError:
        return False


def scan(lab, X_sp, X_so):
    """X_sp[k] と X_so[k'] の間で X_sp[k] = mono · X_so[k'](M, ν) となる組を探す。"""
    hits = []
    for mlab, msub in MAPS.items():
        for nlab, nsub in NUS.items():
            for k in X_sp:
                for kp in X_so:
                    tgt = X_so[kp].subs({**msub, **nsub}, simultaneous=True)
                    ratio = sp.cancel(X_sp[k] / tgt)
                    if is_mono(ratio):
                        hits.append((k, kp, mlab, nlab, sp.factor(ratio)))
    print(f"{lab}：{hits}", flush=True)
    return hits


if __name__ == "__main__":
    g_sp, g_so, f_sp, f_so = load()
    zq, zqp = r_ - 1 / r_, r_ + 1 / r_
    scan("g_R", g_sp, g_so)
    h_sp = {"p1²": sp.cancel((g_sp["(2,)"] + g_sp["(1, 1)"]) / zq), "p2": sp.cancel((g_sp["(2,)"] - g_sp["(1, 1)"]) / zqp)}
    h_so = {"p1²": sp.cancel((g_so["(2,)"] + g_so["(1, 1)"]) / zq), "p2": sp.cancel((g_so["(2,)"] - g_so["(1, 1)"]) / zqp)}
    scan("hatted", h_sp, h_so)
    scan("f_R（SO 側の不変量）", f_sp, f_so)
