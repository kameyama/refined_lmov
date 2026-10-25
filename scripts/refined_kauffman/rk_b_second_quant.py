r"""二次量子化の形での電荷ごとの符号（progress.md 続き 30）。
自由エネルギー log Z~ = Σ_d (1/d) Σ_R f_R(x^d) s_R(Y^d) を、AgSc 1202.4456 (3.13) 型の積 Z~ = Π (1 − (電荷の単項式)·y^α)^{−c} に書いたときの指数 c を調べる。
  - c は Y の単項対称式の基底の係数：c_μ = Σ_R K_{Rμ} f_R（K は Kostka 数。2 箱は c_(2) = f_(2)、c_(1,1) = f_(2) + f_(1,1)）。
  - f の分母 (1 − t^p) は t の冪級数（モードの塔）に展開する。展開した係数 c_{μ;β,k,a}（ν^β t^{k/2} q^{a/2} の係数）が電荷ごとの「状態数 × (−1)^F」。
対象：SO の f（orientifold の指数そのもの）と、比較のための被覆 f^cov。[1] は T(2,2m+1)（m = 1, 2, 3、一様な規格化）、2 箱は trefoil（大きな N）。
変数は (r, s, ν) = (q^{1/2}, t^{1/2}, U t^{−1/2})。"""
import sys, pickle, itertools
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp

r, s, nu = sp.symbols("r s nu", positive=True)


def laurent_terms(e, gens):
    n, d = sp.fraction(sp.cancel(e))
    Pd = sp.Poly(d, *gens)
    assert Pd.is_monomial, sp.factor(d)
    (dm, dc), = Pd.terms()
    return {tuple(a - b for a, b in zip(m, dm)): sp.Rational(c) / dc for m, c in sp.Poly(sp.expand(n), *gens).terms()}


def series_coeffs(e, extra=2):
    """e = L/(1 − s^{2p})（L は Laurent、p は 1 か 2）を s の冪級数に展開し、{(β, k, a): c}（ν^β s^k r^a）を返す。
    k は L の s の最高次 + extra·2p まで（そこから先は 2p 周期で一定）。"""
    for p in (1, 2):
        L = sp.cancel(e * (1 - s**(2 * p)))
        if sp.Poly(sp.fraction(L)[1], r, s, nu).is_monomial:
            break
    else:
        raise ValueError("分母が 1 − t または 1 − t² ではない")
    T = laurent_terms(L, (nu, s, r))
    ks = [k for (_, k, _) in T]
    kmax = max(ks) + extra * 2 * p
    out = {}
    for (b, j, a), c in T.items():
        for k in range(j, kmax + 1, 2 * p):
            out[(b, k, a)] = out.get((b, k, a), 0) + c
    return {key: c for key, c in out.items() if c != 0}, p, kmax


def analyze(lab, e):
    C, p, kmax = series_coeffs(e)
    assert all(c.q == 1 for c in C.values()), "整数でない指数"
    blocks = {}
    for (b, k, a), c in C.items():
        blocks.setdefault((b, k), []).append(c)
    mixed = sorted(key for key, cs in blocks.items() if not (all(c > 0 for c in cs) or all(c < 0 for c in cs)))
    # 単項式ごとの偶奇の符号則 sign = ε (−1)^{α a + γ k + δ β}
    rules = [(al, ga, de) for al, ga, de in itertools.product((0, 1), repeat=3)
             if len({sp.sign(c) * (-1)**(al * a + ga * k + de * b) for (b, k, a), c in C.items()}) == 1]
    # ν^β ごとに偶奇の符号則（(a, k) の一次式）がある β の数
    perbeta = {}
    for b in sorted({b for (b, _, _) in C}):
        sub = {(k, a): c for (bb, k, a), c in C.items() if bb == b}
        ok = [(al, ga) for al, ga in itertools.product((0, 1), repeat=2) if len({sp.sign(c) * (-1)**(al * a + ga * k) for (k, a), c in sub.items()}) == 1]
        perbeta[b] = ok[0] if ok else None
    print(f"{lab}：塔 1/(1 − t^{p})、s^k を k ≤ {kmax} まで。ブロック (β, k) は {len(blocks)} 個、そのうち符号が混在 {len(mixed)} 個 {mixed[:12]}{' …' if len(mixed) > 12 else ''}", flush=True)
    print(f"    単項式ごとの偶奇の符号則（全体）：{rules if rules else 'なし'}；ν^β ごとの符号則 (α, γ)：{perbeta}", flush=True)
    return C, blocks, mixed


def one_box(m):
    from rcsknots.refined_kauffman.kauffman_trefoil import Z_V_T2, unknot_SO, r as r0, s as s0, U as U0
    import rk_k3_torus_family as Tf
    mono = r0**(-(3 * m + 2)) * s0**(7 * m + 4) * U0**(-2 * (2 * m + 1)) * (s0 / r0)**m
    W = sp.cancel(Z_V_T2(m) * mono * unknot_SO)
    hh = sp.cancel((s0 / r0)**m * Tf.h_half(m))
    to = lambda e: sp.cancel(sp.sympify(e).subs({r0: r, s0: s}).subs(U0, nu * s).subs(sp.Symbol("nu"), nu))
    return to(-W), to(-2 * hh)                       # f_(1) = −W_SO、f^cov_(1) = −2·(t/q)^{δ/2}·½h


def two_box():
    u = sp.Symbol("u", positive=True)
    d = pickle.load(open("data/refined_kauffman/daha/largeN_2box_m1.pkl", "rb"))
    get = lambda key, k: sp.cancel(sp.sympify(d[key][str(k)], locals={"r": r, "s": s, "u": u}).subs(u, nu * s))
    f = {k: get("f", k) for k in ((1,), (2,), (1, 1))}
    fc = {k: get("fcov", k) for k in ((1,), (2,), (1, 1))}
    mono = lambda F: {"m_(2)": F[(2,)], "m_(1,1)": sp.cancel(F[(2,)] + F[(1, 1)])}
    return mono(f), mono(fc)


if __name__ == "__main__":
    for m in (1, 2, 3):
        f1, fc1 = one_box(m)
        analyze(f"[1] T(2,{2*m+1}) SO の f", f1)
        analyze(f"[1] T(2,{2*m+1}) 被覆の f^cov", fc1)
    fso, fcov = two_box()
    for k in fso:
        analyze(f"2 箱 trefoil SO の f の {k} の係数", fso[k])
        analyze(f"2 箱 trefoil 被覆の f^cov の {k} の係数", fcov[k])


def analyze_factored(lab, e):
    """ν^β ごとに、塔を外した L_β = (1 − t^p)·[e]_β から (q − t) の冪を最大限くくり出し、商を塔に展開して (β, k) ブロックの符号を見る。"""
    for p in (1, 2):
        L = sp.cancel(e * (1 - s**(2 * p)))
        if sp.Poly(sp.fraction(L)[1], r, s, nu).is_monomial:
            break
    T = laurent_terms(L, (nu, s, r))
    out, nmix, nblk = {}, 0, 0
    for b in sorted({b for (b, _, _) in T}):
        Lb = sum(c * s**k * r**a for (bb, k, a), c in T.items() if bb == b)
        kq = 0
        while True:
            quo = sp.cancel(Lb / (r**2 - s**2))
            if not sp.Poly(sp.fraction(quo)[1], r, s).is_monomial:
                break
            Lb, kq = quo, kq + 1
        C, _, _ = series_coeffs(sp.cancel(Lb * nu**b / (1 - s**(2 * p))))
        blocks = {}
        for (_, k, a), c in C.items():
            blocks.setdefault(k, []).append(c)
        mix = [k for k, cs in blocks.items() if not (all(c > 0 for c in cs) or all(c < 0 for c in cs))]
        sg = {"+" if all(c > 0 for cs in blocks.values() for c in cs) else "−" if all(c < 0 for cs in blocks.values() for c in cs) else "±"}
        out[b] = f"(q−t)^{kq}{'' if not mix else f' 混在 {len(mix)}/{len(blocks)}'}{sg.pop()}"
        nmix += len(mix); nblk += len(blocks)
    print(f"{lab}：(q − t) を ν^β ごとにくくり出した後の混在ブロック {nmix}/{nblk}；{out}", flush=True)


def two_box_hatted():
    """SO の f̂ = M(q^{-1})^{-1} f と被覆の f̂^cov（Schur 基底の hatted）。"""
    f, fc = two_box_raw()
    det = r**2 - r**-2
    hat = lambda F: {"(2)": sp.cancel((r * F[(2,)] + F[(1, 1)] / r) / det), "(1,1)": sp.cancel((F[(2,)] / r + r * F[(1, 1)]) / det)}
    return hat(f), hat(fc)


def two_box_raw():
    u = sp.Symbol("u", positive=True)
    d = pickle.load(open("data/refined_kauffman/daha/largeN_2box_m1.pkl", "rb"))
    get = lambda key, k: sp.cancel(sp.sympify(d[key][str(k)], locals={"r": r, "s": s, "u": u}).subs(u, nu * s))
    return ({k: get("f", k) for k in ((1,), (2,), (1, 1))}, {k: get("fcov", k) for k in ((1,), (2,), (1, 1))})


def main2():
    for m in (1, 2, 3):
        f1, fc1 = one_box(m)
        analyze_factored(f"[1] T(2,{2*m+1}) SO の f", f1)
    fh, fch = two_box_hatted()
    for k in fh:
        for lab, e in ((f"2 箱 SO の f̂_{k}", fh[k]), (f"2 箱 被覆の f̂^cov_{k}", fch[k])):
            try:
                analyze(lab, e)
                analyze_factored(lab, e)
            except (AssertionError, ValueError) as ex:
                print(f"{lab}：{type(ex).__name__} {ex}", flush=True)


def kn_G_check():
    """KN の t 枠 (t^{-1/2} − t^{1/2})·f̂ を (V, W, ν) に直し、σ 対称なら G 展開して ν^β ごと・(G^g, W^j) ごとの符号を見る。"""
    from rk_k4_crosscap import V, W, G, in_G
    nuk = sp.Symbol("nu")
    fh, fch = two_box_hatted()
    for lab, F in (("被覆", fch), ("SO", fh)):
        for k, e in F.items():
            e = sp.cancel(e * (1 / s - s))
            lap = sp.Poly(sp.fraction(e)[1], r, s, nu).is_monomial
            ev = sp.expand(e.subs(nu, nuk).subs({r: V * W, s: V / W}))
            P1 = in_G(ev)
            if P1 is None:
                print(f"{lab} の (1/s − s) f̂_{k}：Laurent {lap}、σ 対称でない（G 展開できない）", flush=True)
                continue
            P1 = sp.expand(P1)
            res = {}
            for b in sorted({tm.as_coeff_exponent(nuk)[1] for tm in sp.Add.make_args(P1)}):
                blk = sp.expand(P1.coeff(nuk, b))
                cs = sp.Poly(sp.expand(blk * W**200), G, W).terms()
                coh = len({sp.sign(c) for _, c in cs}) == 1
                # KN の (−1)^{2J_r}：W の指数の偶奇で符号が決まるか
                kn = len({sp.sign(c) * (-1)**(((m[1] - 200) // 2) % 2) for m, c in cs}) == 1 or len({sp.sign(c) * (-1)**(m[1] % 2) for m, c in cs}) == 1
                gw = len({(sp.sign(c), m[0]) for m, c in cs}) == len({m[0] for m, c in cs})   # g ごとに同符号
                res[b] = "そろう" if coh else ("g ごと" if gw else ("偶奇" if kn else "混在"))
            print(f"{lab} の (1/s − s) f̂_{k}：Laurent {lap}、G 展開できる。ν^β ごとの符号：{res}", flush=True)


