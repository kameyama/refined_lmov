"""核まで brane/image に分けた n = 2 の記号ブロックと、入れ替え Ĩ のトレースの候補（続き 117、theory §12.47）。
片側（brane b / image i）の記号：行 a2,a4,b4,…（rk_f_n2_rows_lab と同じ）、次元 D{λ}_{k}（Psym(λ) の ν^{−k} の係数）、
Macdonald→Schur の係数 M{λ}_{a}、測度の φ_ρ。composite の次元は D_α^b D_β^i − c₁(D₁^bD₁^i + (q−1)(t−1)C₁₁) − c₀（c₁, c₀ は定数）、
[α;1̄] は D_α^b D₁^i − c_α D₁^b（縮約は brane の 1 箱）。共有の核：κ（背景）、C₁₁（随伴の連結部分）、(1−t)/z_q、1/den。
トレース：G を含まない単項式で brane の内容 = image の内容のものだけ残し、片側の値に ψ²、共有の係数に規則 K0（そのまま）/K1（q → −q、ν^{−2} → −ν^{−2}）。
背景（記号なしの項）は別に扱う。実行（サンドボックス外）：PYTHONPATH=src:scripts/refined_kauffman python scripts/refined_kauffman/rk_f_n2_rows_ker.py [kmax]"""
import sys, pickle, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_f_n2_rows_sym as RS
from rk_f_n2_rows_sym import r, s, q, t, nu, Ser, rat_series, Psym, dimc, decomposition_flint, macdonald_P, P2, LAMS
import rk_f_n2_rows_lab as L

OUT = "data/refined_kauffman/daha"
LAMKEY = {(1,): "1", (2,): "2", (1, 1): "11"}


def build_ker(kmax=14, verbose=True):
    t0 = time.time()
    RS.KMAX = kmax
    qt1 = (q - 1) * (t - 1)
    # ---- 片側の記号と値 ----
    sym = {}          # 記号 → (側, 基底名, 値)
    def mk(name, side, val):
        x = sp.Symbol(f"{name}{side}")
        sym[x] = (side, name, val)
        return x
    rows = {"b": {}, "i": {}}
    val = RS.actual_rows()
    for base, d in (("b", L.BR), ("i", L.IM)):
        for k, x in d.items():
            sym[x] = (base, k, val[getattr(RS, k)])
            rows[base][k] = x
    Dser = {lam: rat_series(Psym(lam)) for lam in ((1,), (2,), (1, 1))}
    Dsym = {}
    for lam, ser in Dser.items():
        for side in ("b", "i"):
            for k, c in ser.c.items():
                Dsym[(lam, side, k)] = mk(f"D{LAMKEY[lam]}_{k}", side, c.as_expr())
    # Macdonald→Schur の係数と φ の記号も先に作る（SYMS を一度に決める）
    PY = {(): {(): sp.Integer(1)}}
    for n in (1, 2):
        for lam, v in macdonald_P(n, t, q, lambda c: sp.Rational(c.numerator, c.denominator)).items():
            PY[lam] = v
    Msym = {}
    for lam, d in PY.items():
        for a_, x in d.items():
            for side in ("b", "i"):
                Msym[(lam, a_, side)] = mk(f"M{''.join(map(str, lam)) or '0'}_{''.join(map(str, a_)) or '0'}", side, sp.sympify(x))
    rs = ((1, 1), (2,))
    zq = r - 1 / r
    phi = {(1, 1): zq, (2,): r + 1 / r}
    IPhi_ = {(rho, side): mk(f"IPHI{''.join(map(str, rho))}", side, 1 / phi[rho]) for rho in rs for side in ("b", "i")}
    RS.SYMS = list(sym.keys()) + [x for d in L.GS for x in d.values()]
    RS.ZERO = RS.P(0)
    P = RS.P
    D = {(lam, side): Ser({k: P(Dsym[(lam, side, k)]) for k in Dser[lam].c}) for lam in Dser for side in ("b", "i")}
    def Rser(side, lam):
        rr = rows[side]
        return {(1,): Ser({2: P(rr["a2"]), 4: P(rr["a4"])}), (2,): Ser({4: P(rr["b4"]), 6: P(rr["b6"]), 8: P(rr["b8"])}), (1, 1): Ser({4: P(rr["c4"]), 6: P(rr["c6"]), 8: P(rr["c8"])})}[lam]
    Gadj = Ser({j: P(v) for j, v in RS.gadj.items()})
    Gc = {(2,): Ser({j: P(v) for j, v in RS.g21.items()}), (1, 1): Ser({j: P(v) for j, v in RS.g111.items()})}
    Gb = {((2,), (2,)): Ser({j: P(v) for j, v in RS.g22.items()}), ((2,), (1, 1)): Ser({j: P(v) for j, v in RS.g211.items()}), ((1, 1), (1, 1)): Ser({j: P(v) for j, v in RS.g1111.items()})}
    # ---- 共有の核 ----
    A = nu**2
    qn = (nu - 1 / nu) / (s - 1 / s)
    kappa_e = (1 - A) * (1 - q) / ((1 - t) * (1 - q * A / t))
    kappa = rat_series(kappa_e)
    C11 = rat_series(sp.cancel((qn * qn - kappa_e - Psym((1,))**2) / qt1))
    Kadj = D[((1,), "b")] * D[((1,), "i")] + C11 * qt1
    W = {(lam, side): Rser(side, lam) * D[(lam, side)] for lam in ((1,), (2,), (1, 1)) for side in ("b", "i")}
    Wadj = (Rser("b", (1,)) * Rser("i", (1,)) + Gadj * qt1) * Kadj
    pair = {((), ()): Ser.const(1), ((1,), (1,)): Wadj + kappa}
    for lam in ((1,), (2,), (1, 1)):
        pair[(lam, ())] = W[(lam, "b")]
        pair[((), lam)] = W[(lam, "i")]
    # composite の次元 = D_α^b D_β^i + (q−1)(t−1) C_{αβ}（C は共有の連結部分。d_α・c₁・c₀ は ν に依存するので共有の級数にする）
    conn = {}
    for alpha in P2:
        dc = dimc(alpha, (1,))
        C = rat_series(sp.cancel((dc - Psym(alpha) * Psym((1,))) / qt1))
        conn[(alpha, (1,))] = C
        pair[(alpha, (1,))] = (Rser("b", alpha) * Rser("i", (1,)) + Gc[alpha] * qt1) * (D[(alpha, "b")] * D[((1,), "i")] + C * qt1) - C * Rser("b", (1,)) * qt1
        pair[((1,), alpha)] = (Rser("b", (1,)) * Rser("i", alpha) + Gc[alpha] * qt1) * (D[((1,), "b")] * D[(alpha, "i")] + C * qt1) - C * Rser("i", (1,)) * qt1
    for alpha in P2:
        for beta in P2:
            key = (alpha, beta) if (alpha, beta) != ((1, 1), (2,)) else ((2,), (1, 1))
            c1, c0, Kb = decomposition_flint(*key)
            C = rat_series(sp.cancel((Kb - Psym(alpha) * Psym(beta)) / qt1))
            conn[key] = C
            pair[(alpha, beta)] = (Rser("b", alpha) * Rser("i", beta) + Gb[key] * qt1) * (D[(alpha, "b")] * D[(beta, "i")] + C * qt1) + Wadj * rat_series(c1) + rat_series(c0)
    if verbose:
        print(f"  核と対の値：{time.time() - t0:.0f} 秒；連結部分 C の級数の最上の係数 = { {k: sp.factor(min(v.c.items())[1].as_expr()) for k, v in conn.items()} }", flush=True)
    # ---- Z と log（Macdonald→Schur の係数も片側の記号） ----
    T_ = lambda lam: tuple(sum(1 for x in lam if x > j) for j in range(lam[0])) if lam else ()
    Z = {}
    for (lam, mu), v in pair.items():
        sg = (-1)**(sum(lam) + sum(mu))
        for a_ in PY[T_(lam)]:
            for b_ in PY[T_(mu)]:
                Z[(a_, b_)] = Z.get((a_, b_), Ser()) + v * (sg * Msym[(T_(lam), a_, "b")] * Msym[(T_(mu), b_, "i")])
    def mul2(A_, B_):
        out = {}
        for (a1, a2_), x in A_.items():
            for (b1, b2), y in B_.items():
                if sum(a1) + sum(b1) <= 2 and sum(a2_) + sum(b2) <= 2:
                    k = (tuple(sorted(a1 + b1, reverse=True)), tuple(sorted(a2_ + b2, reverse=True)))
                    out[k] = out.get(k, Ser()) + x * y
        return out
    X = {k: v for k, v in Z.items() if k != ((), ())}
    Lg, Xn = {}, dict(X)
    for n in range(1, 5):
        for k, v in Xn.items():
            Lg[k] = Lg.get(k, Ser()) + v * sp.Rational((-1)**(n + 1), n)
        if n < 4:
            Xn = mul2(Xn, X)
    if verbose:
        print(f"  log Z：{time.time() - t0:.0f} 秒", flush=True)
    x11 = Lg[((1,), (1,))]
    CHI = {"2": {(1, 1): 1, (2,): 1}, "11": {(1, 1): 1, (2,): -1}}
    den = s * (t - nu**2) * (nu**2 - 1) * (nu**2 * q - 1)
    inv_den = rat_series(1 / den)
    N = {}
    for lab in LAMS:
        M = Ser()
        for rho in rs:
            for sig in rs:
                v = Lg.get((rho, sig), Ser())
                # 1/(φ_ρ^b φ_σ^i) は記号の逆数になるので、記号 IPHI（値 1/φ）を使う
                M = M + v * ((1 - t) * CHI[lab][rho] * CHI[lab][sig] / zq) * IPhi_[(rho, "b")] * IPhi_[(sig, "i")]
        N[lab] = M * inv_den
    adams_coef = {lab: sp.Rational(-1, 2) * (1 - t) * CHI[lab][(2,)]**2 / (zq * phi[(2,)]**2) for lab in LAMS}
    if verbose:
        print(f"  N_λ（Adams の項を除く）：{time.time() - t0:.0f} 秒", flush=True)
    return N, x11, adams_coef, inv_den, sym


def _level(name):
    """記号の ν の準位（ν^{−ℓ}）：行 a2 → 2、D1_k → k、M・IPHI → 0。"""
    if name[0] in "abc" and name[1:].isdigit():
        return int(name[1:])
    if name.startswith("D"):
        return int(name.split("_")[1])
    return 0


def trace(N, sym, kmax):
    """各 (λ, K) で、G を含まない単項式を brane/image の内容で分け、対角のトレース（片側の値に ψ²、共有の係数は K0：そのまま／
    K1：q → −q と (−1)^{(K−ℓ)/2}（ℓ は単項式の準位））と、非対角・G を含む・背景の部分の値を flint で返す。"""
    from rcsknots.refined_kauffman.flint_rf3 import RF, from_sympy
    gsyms = {x for d in L.GS for x in d.values()}
    vals = {x: from_sympy(sp.sympify(v)) for x, (side, name, v) in sym.items()}
    out = {}
    for lab in LAMS:
        for K, poly in N[lab].c.items():
            if K > kmax:
                continue
            tr0 = tr1 = off = withG = bg = RF.coerce(0)
            n_diag = n_off = n_G = 0
            for mon, c in poly.terms():
                ce = from_sympy(sp.sympify(c).as_expr() if hasattr(c, "as_expr") else sp.sympify(c))
                names = {"b": [], "i": []}; vb = vi = RF.coerce(1); hasG = False; lev = 0
                for x, e in zip(RS.SYMS, mon):
                    if not e:
                        continue
                    if x in gsyms:
                        hasG = True; continue
                    side, name, v = sym[x]
                    names[side] += [name] * e
                    lev += _level(name) * e
                    if side == "b":
                        vb = vb * vals[x]**e
                    else:
                        vi = vi * vals[x]**e
                full = ce * vb * vi
                if hasG:
                    withG = withG + full; n_G += 1; continue
                if not names["b"] and not names["i"]:
                    bg = bg + full; continue
                if sorted(names["b"]) == sorted(names["i"]):
                    n_diag += 1
                    tr0 = tr0 + ce * vb.adams(2)
                    # K1：共有の係数に準位の符号 (−1)^{(K−ℓ)/2} だけを掛ける（q → −q は係数が r の奇関数を含むので定義できない）
                    sgn = (-1)**(((K - lev) // 2) % 2) if (K - lev) % 2 == 0 else None
                    assert sgn is not None, (lab, K, lev)
                    tr1 = tr1 + ce * vb.adams(2) * sgn
                else:
                    off = off + full; n_off += 1
            cv = lambda x: x.to_sympy()
            out[(lab, K)] = {"tr_K0": cv(tr0), "tr_K1": cv(tr1), "off": cv(off), "withG": cv(withG), "bg": cv(bg), "n": (n_diag, n_off, n_G)}
            print(f"λ = {lab}、ν^-{K}：単項式（対角 {n_diag}、非対角 {n_off}、G {n_G}）、対角のトレース K0 = {sp.factor(out[(lab, K)]['tr_K0'])}"[:400], flush=True)
    return out


if __name__ == "__main__":
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    t0 = time.time()
    N, x11, adams_coef, inv_den, sym = build_ker(kmax)
    res = trace(N, sym, kmax)
    out = {"trace": {k: {kk: str(v) for kk, v in d.items()} for k, d in res.items()},
           "x11": {k: str(v.as_expr()) for k, v in x11.c.items()}, "adams_coef": {lab: str(v) for lab, v in adams_coef.items()},
           "inv_den": {k: str(v.as_expr()) for k, v in inv_den.c.items()}, "syms": {str(x): (side, name, str(v)) for x, (side, name, v) in sym.items()}}
    f = f"{OUT}/n2_rows_ker_m1_k{kmax}.pkl"
    pickle.dump(out, open(f, "wb"))
    print(f"保存：{f}（{time.time() - t0:.0f} 秒）", flush=True)
