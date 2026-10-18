r"""報告書（docs/refined_kauffman/tex/report.tex）に載せる計算結果の式を、保存済みのデータとコードから LaTeX で生成する。
出力：docs/refined_kauffman/tex/gen/*.tex（report.tex が \input する）。式を手で写さないための生成器で、生成時に既知の関係を assert で確かめる。"""
import sys, pickle
from pathlib import Path
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp

OUT = Path("docs/refined_kauffman/tex/gen")
BOLD = [(r"\ba", 1), (r"\bq", 1), (r"\bt", 1)]          # (𝐚, 𝐪, 𝐭)
QTU = [("q", 2), ("t", 2), ("U", 1)]                     # r = q^{1/2}, s = t^{1/2}, U


def laurent(expr, gens):
    """分母が単項式の有理式 → {指数: 係数}。"""
    n, d = sp.fraction(sp.cancel(sp.together(expr)))
    Pd = sp.Poly(d, *gens)
    assert Pd.is_monomial, Pd
    (dm, dc), = Pd.terms()
    return {tuple(a - b for a, b in zip(m, dm)): sp.Rational(c) / dc for m, c in sp.Poly(sp.expand(n), *gens).terms()}


def _exp(f):
    return str(f.p) if f.q == 1 else f"{f.p}/{f.q}"


def mono(ex, spec):
    out = []
    for e, (name, den) in zip(ex, spec):
        if e:
            f = sp.Rational(e, den)
            out.append(name if f == 1 else f"{name}^{{{_exp(f)}}}")
    return " ".join(out)


def coef(c, m, first):
    sgn = "-" if c < 0 else ("" if first else "+")
    a = abs(c)
    if a == 1 and m:
        body = m
    elif a.q == 1:
        body = f"{a.p}{(' ' + m) if m else ''}"
    else:
        body = rf"\tfrac{{{a.p}}}{{{a.q}}}{(' ' + m) if m else ''}"
    return f"{sgn}{body}"


def terms_lines(d, spec, per=6):
    keys = sorted(d, key=lambda k: tuple(-x for x in k))
    ts = [coef(d[k], mono(k, spec), i == 0) for i, k in enumerate(keys)]
    return [" ".join(ts[i:i + per]) for i in range(0, len(ts), per)] or ["0"]


def grouped(d, idx, gname, spec_rest, per=6, label=None):
    """d の idx 番目の変数の冪ごとに分けた align*（label は前に置く段落）。"""
    groups = {}
    for k, v in d.items():
        groups.setdefault(k[idx], {})[k[:idx] + k[idx + 1:]] = v
    rows = []
    for e in sorted(groups):
        ls = terms_lines(groups[e], spec_rest, per)
        g = gname if e == 1 else f"{gname}^{{{e}}}"
        rows.append(rf"[{g}]\colon\ & {ls[0]}")
        rows += [rf"& \quad {l}" for l in ls[1:]]
    return (label + "\n" if label else "") + align(" \\\\\n".join(rows))


def align(body):
    return "\\begin{align*}\n" + body + "\n\\end{align*}\n"


def write(name, text):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(text)
    print("wrote", OUT / name, flush=True)


# ---- 1. unrefined：Mariño (4.8) ----
def gen_unrefined():
    import rk_b_marino_tq as M
    rk, H, G, gM, ghM = M.marino_trefoil()
    body = []
    for R, lab in (((2,), r"\hat g_{[2]}"), ((1, 1), r"\hat g_{[1,1]}")):
        e = rk.z_expansion(ghM[R])                         # {(z の冪, ν の冪): 係数}
        d = {(k, Q): sp.Rational(c) for (k, Q), c in e.items()}
        body.append(grouped(d, 0, "z", [(r"\nu", 1)], per=6, label=rf"${lab}$（$z = t_M - t_M^{{-1}}$ の冪ごと）："))
    write("unrefined.tex", "\n".join(body))


# ---- 2. [1]：T(2,2m+1) の g^{true}_[1]、c=1 の G 展開、c=2 の Π_m ----
def gen_torus_one_box(ms=(1, 2, 3)):
    import rk_b_uniform as Bu
    from rk_k4_crosscap import to_VW, split_parity, in_G, V, W, nu, G
    Q, T, A = Bu.Q, Bu.Tb, Bu.A
    out = []
    for m in ms:
        g = Bu.g_true(m)
        d = laurent(Bu.bold(g), (A, Q, T))
        out.append(grouped(d, 0, r"\ba", BOLD[1:], per=6, label=rf"$m = {m}$（$T(2,{2*m+1})$、{len(d)} 項）：$g^{{\mathrm{{true}}}}_{{[1]}}$ の $\ba$ の冪ごとの係数"))
        E, _ = split_parity(to_VW(g))
        P1 = sp.expand(in_G(E))
        dG = laurent(P1, (nu, G, W))
        out.append(grouped(dG, 0, r"\nu", [("G", 1), ("W", 1)], per=6, label=rf"$m = {m}$：$c = 1$ 部分の $G$ 展開（$\nu^\beta$ ごと、$W = (q/t)^{{1/4}}$）"))
        so2, even_t, c2, Pi = Bu.check_one(m)
        assert so2 == {1} and even_t and c2
        dP = laurent(Pi, (A, Q, T))
        out.append(align(r"\Pi_{" + str(m) + r"} &= " + " \\\\\n&\\quad ".join(terms_lines(dP, BOLD, per=6))))
    write("torus_one_box.tex", "\n".join(out))


# ---- 3. DAHA：一般の N の R_λ(U) ----
def gen_daha():
    import rk_c_interp as I
    r, s = sp.symbols("r s", positive=True)
    c1 = I.interpolate(1, (1,), [2, 3, 4, 5], -5, -2)
    saved = pickle.load(open("data/refined_kauffman/daha/R_U_m1_2box.pkl", "rb"))
    Rl = {(1,): {j: v.to_sympy() for j, v in c1.items()}}
    for lam in ((2,), (1, 1)):
        Rl[lam] = {j: sp.sympify(v, locals={"r": r, "s": s}) for j, v in saved[lam].items()}
    out = []
    for lam, c in Rl.items():
        rows = []
        for j in sorted(c):
            d = laurent(c[j].subs({sp.Symbol("r", positive=True): r, sp.Symbol("s", positive=True): s}), (r, s))
            ls = terms_lines(d, QTU[:2], per=6)
            rows.append(rf"[U^{{{j}}}]\colon\ & {ls[0]}")
            rows += [rf"& \quad {l}" for l in ls[1:]]
        name = "(" + ",".join(map(str, lam)) + ")"
        out.append(rf"$R_{{{name}}}(U)$（$U^j$ の係数）：" + "\n" + align(" \\\\\n".join(rows)))
    radj = pickle.load(open("data/refined_kauffman/daha/Radj_A_m1.pkl", "rb"))
    rows = []
    for j in sorted(radj):
        d = laurent(sp.sympify(radj[j], locals={"r": r, "s": s}), (r, s))
        ls = terms_lines(d, QTU[:2], per=6)
        rows.append(rf"[A^{{{_exp(sp.Rational(j, 2))}}}]\colon\ & {ls[0]}")
        rows += [rf"& \quad {l}" for l in ls[1:]]
    out.append(r"被覆の随伴 $(1;1)$ の $R_{\mathrm{adj}}$（$A^{1/2} = \nu$ の冪の係数）：" + "\n" + align(" \\\\\n".join(rows)))
    write("daha.tex", "\n".join(out))


# ---- 4. 大きな N の 2 箱：ĝ、留数、商の符号 ----
def gen_largeN():
    import rk_b_2box_structure as B
    import rk_b_uniform as Bu
    r, s, nu = B.r, B.s, B.nu
    A, Q, T = sp.symbols("A Q T")
    bold = lambda e: sp.cancel(e.subs({r: -T * Q, s: Q, nu: A}, simultaneous=True))
    G = B.ghat()
    # ĝ_(1) = −g^{true}_[1]（m = 1）の確認
    gt = laurent(Bu.bold(Bu.g_true(1)), (Bu.A, Bu.Q, Bu.Tb))
    g1 = laurent(bold(G[(1,)]), (A, Q, T))
    assert g1 == {k: -v for k, v in gt.items()}
    out = []
    dsum = laurent(bold(G["sum"]), (A, Q, T))
    out.append(grouped(dsum, 0, r"\ba", BOLD[1:], per=6, label=rf"$\hat g_{{(2)}} + \hat g_{{(1,1)}}$（整数 Laurent、{len(dsum)} 項）："))
    D = sp.cancel((1 + s**2) * (G[(2,)] - G[(1, 1)]))
    dD = laurent(bold(D), (A, Q, T))
    out.append(grouped(dD, 0, r"\ba", BOLD[1:], per=6, label=rf"$(1 + t)(\hat g_{{(2)}} - \hat g_{{(1,1)}})$（整数 Laurent、{len(dD)} 項）："))
    _, Qs = B.refined_only(G["sum"])
    dQ = laurent(bold(Qs), (A, Q, T))
    out.append(grouped(dQ, 0, r"\ba", BOLD[1:], per=6, label=rf"和の商 $X/((q - t)(\nu^2 - t))$（{len(dQ)} 項）："))
    _, Q1 = B.refined_only(G[(1,)])
    dQ1 = laurent(bold(Q1), (A, Q, T))
    out.append(r"$\hat g_{(1)}$ の商：" + "\n" + align(r"X_{(1)}/((q - t)(\nu^2 - t)) &= " + " ".join(terms_lines(dQ1, BOLD, per=8))))
    write("largeN_2box.tex", "\n".join(out))
    # 符号の表
    rows = []
    for lab, k in ((r"\hat g_{(1)}", (1,)), (r"\hat g_{(2)} + \hat g_{(1,1)}", "sum")):
        sb = B.sign_blocks(G[k])
        rows.append(rf"${lab}$ & " + ", ".join(rf"$\ba^{{{e}}}$: {n}（{sg}）" for e, (n, sg) in sb.items()) + r" \\")
    write("sign_table.tex", "\\begin{center}\n\\begin{tabular}{lp{11.5cm}}\n\\toprule\n" + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n\\end{center}\n")
    # (1 + t) の留数
    rr, ss, uu = sp.symbols("r s u", positive=True)
    d = pickle.load(open("data/refined_kauffman/daha/largeN_2box_m1.pkl", "rb"))
    g = {eval(k): sp.sympify(v, locals={"r": rr, "s": ss, "u": uu}) for k, v in d["g"].items()}
    det = rr**2 - rr**-2
    Dr = sp.cancel((1 + ss**2) * ((rr * g[(2,)] + g[(1, 1)] / rr) - (g[(2,)] / rr + rr * g[(1, 1)])) / det)
    res = sp.cancel(Dr.subs(ss, sp.I) / (2 * sp.I))
    claim = -sp.I * (uu**2 - 1)**2 * (rr**2 * uu**2 - 1) * (rr**6 * uu**6 + rr**6 * uu**4 - rr**2 * uu**2 + 1) / (2 * rr**9 * uu**12)
    assert sp.cancel(res - claim) == 0
    print("留数の式を確認", flush=True)


if __name__ == "__main__":
    import rk_b_marino_tq as M
    assert all(M.compare().values())
    gen_unrefined()
    gen_torus_one_box()
    gen_daha()
    gen_largeN()
