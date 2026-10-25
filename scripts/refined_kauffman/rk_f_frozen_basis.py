r"""凍結した状態の多重項の基底の検査（theory §12、progress.md 続き 63）。
  凍結した状態 N = Σ̂/I（電荷 β）の被覆 Σ̂ は Ĩ 不変な曲線（電荷 2β）で、N のゼロモードは Σ̂ の Ĩ 不変なモードである。
    (a) I で入れ替わる取っ手の対 → G = z_q z_t（q モードと t モードが 1 個ずつ）：σ 対称、
    (b) I で自分自身に移る取っ手 → σ 対称でない部分の源、
    (c) 境界モード：入れ替わる境界の対ごとに和が残り、拘束で 1 個減る（§11.3 の n − 1 個）。hatted（M(q^{−1})^{−1}）はこれを除く。
  hatted の c の規則 β ≡ |R| − 1 (mod 2) は「境界モードを除いた z の冪 2g + c − 1 の偶奇」で、c = 1 ⟺ クロスキャップ 1 個（k 奇）。
  検査 1（ℓ 重みの非対称部分）：ブロック X_β = Σ_m c_m(w) q_ℓ^m（q_ℓ = (qt)^{1/2} = V²、w = (q/t)^{1/2} = W²、m は Lefschetz の重み）に対し
    d_m = c_m − c_{−m}（m > 0）が (1 − w²) × 同符号の多項式で、符号がブロック内の全 m で一致するか。
    X_β = A⁻_β + Σ_{m>0} d_m q_ℓ^m = A⁺_β − Σ_{m>0} d_m q_ℓ^{−m}（A^± は σ 対称で整数）。A^± の G 展開の符号も見る。
  検査 2（c = 2 の部分）：(𝐚² − 𝐪²)(𝐭² − 1) で割った商の符号。
  検査 3（任意、--ilp）：基底 {G^g w^j, z_q² G^g w^j, z_t² G^g w^j}（一次従属）で同符号の展開が存在するかの整数計画。
  対象：[1]（T(2,3), T(2,5)）、framed unknot の 2 箱の p₁² 方向（𝖿 = 1..4）、結び目の 2 箱の p₁² 方向（T(2,3), T(2,5)）、
        結び目の 3 箱の ĝ_R（T(2,3)）。"""
import sys, io, contextlib
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import numpy as np
import sympy as sp
from rk_k4_crosscap import V, W, G, in_G

r, s, nu = sp.symbols("r s nu", positive=True)
w = sp.Symbol("w")
t = s**2
zq, zt, zqp = r - 1 / r, s - 1 / s, r + 1 / r


def to_vw(e):
    return sp.expand(sp.cancel(e).subs({r: V * W, s: V / W}, simultaneous=True))


def sigma_vw(e):
    return sp.expand(e.subs(V, 1 / V))


def nu_blocks(X):
    """ν の冪ごとのブロック {β: 係数（r, s の Laurent 多項式）}。分母の ν の冪を基準にする。"""
    X = sp.cancel(X); n_, d_ = sp.fraction(X)
    Pd = sp.Poly(d_, nu); dn = Pd.monoms()[0][0] if Pd.is_monomial else 0
    P = sp.Poly(sp.expand(n_), nu)
    return {int(k[0] - dn): sp.cancel(c / (d_ / nu**dn)) for k, c in zip(P.monoms(), P.coeffs())}


def parity_split(X, even=True):
    b = nu_blocks(X)
    return sum((c * nu**k for k, c in b.items() if (k % 2 == 0) == even), sp.Integer(0))


def ell_table(blk):
    """ブロックを {(m, j): c}（q_ℓ^m w^j = V^{2m} W^{2j}）に。"""
    d = {}
    for tm in sp.Add.make_args(to_vw(blk)):
        c, mon = tm.as_coeff_Mul()
        ev = mon.as_coeff_exponent(V)[1] if mon != 1 else 0
        ew = mon.as_coeff_exponent(W)[1] if mon != 1 else 0
        assert ev % 2 == 0 and ew % 2 == 0 and c.is_integer, tm
        d[(int(ev) // 2, int(ew) // 2)] = int(c)
    return d


def sign_of(cs):
    cs = [c for c in cs if c != 0]
    if not cs:
        return "0"
    return "+" if all(c > 0 for c in cs) else "−" if all(c < 0 for c in cs) else "±"


def g_sign(e_vw):
    """σ 対称な Laurent 多項式（V, W）の G 展開の係数の符号。(𝐭² − 1) = W⁴ − 1 で割った版も。"""
    P1 = in_G(e_vw)
    if P1 is None:
        return "σ非対称", None
    cs = sp.Poly(sp.expand(P1 * W**400), G, W).coeffs()
    sg = sign_of(cs)
    quo = sp.cancel(e_vw / (W**4 - 1))
    sg2 = None
    if sp.Poly(sp.fraction(quo)[1], V, W).is_monomial:
        P2 = in_G(sp.expand(quo))
        sg2 = sign_of(sp.Poly(sp.expand(P2 * W**400), G, W).coeffs()) if P2 is not None else "σ非対称"
    return sg, sg2


def ell_asym(blk):
    """検査 1。戻り値：dict(sym, rows={m: (割り切れる, 符号, 商)}, uniform, Aminus_sign, Aplus_sign)。"""
    d = ell_table(blk)
    mmax = max((abs(m) for m, _ in d), default=0)
    rows = {}
    Am = sp.Integer(0); Ap = sp.Integer(0)
    for m in range(0, mmax + 1):
        cm = sum((c * w**j for (mm, j), c in d.items() if mm == m), sp.Integer(0))
        cmm = sum((c * w**j for (mm, j), c in d.items() if mm == -m), sp.Integer(0))
        if m == 0:
            Am += cm.subs(w, W**2); Ap += cm.subs(w, W**2)
            continue
        Am += sp.expand(cmm.subs(w, W**2) * (V**(2 * m) + V**(-2 * m)))
        Ap += sp.expand(cm.subs(w, W**2) * (V**(2 * m) + V**(-2 * m)))
        dm = sp.expand(cm - cmm)
        if dm == 0:
            continue
        q_, rem = sp.div(sp.Poly(sp.expand(dm * w**200), w), sp.Poly(1 - w**2, w))
        if rem.is_zero:
            rows[m] = (True, sign_of(q_.coeffs()), sp.expand(q_.as_expr() / w**200))
        else:
            rows[m] = (False, sign_of(sp.Poly(sp.expand(dm * w**200), w).coeffs()), dm)
    signs = {sg for ok, sg, _ in rows.values()}
    uniform = all(ok for ok, _, _ in rows.values()) and len(signs) <= 1 and "±" not in signs
    return dict(sym=not rows, rows=rows, uniform=uniform, asym_sign=(signs.pop() if len(signs) == 1 else "±"),
                Aminus=g_sign(sp.expand(Am)), Aplus=g_sign(sp.expand(Ap)))


def test_c1(label, X, verbose=False):
    print(f"== {label}", flush=True)
    out = {}
    for beta, blk in sorted(nu_blocks(X).items()):
        res = ell_asym(blk)
        out[beta] = res
        if res["sym"]:
            sg, sg2 = res["Aminus"]
            print(f"  ν^{beta}: σ 対称。G 展開の符号 {sg}" + (f"、/(𝐭²−1) で {sg2}" if sg2 else ""), flush=True)
        else:
            nrow = len(res["rows"])
            print(f"  ν^{beta}: σ 非対称。d_m（{nrow} 行）が (1−w²)×同符号で全行同じ符号: {res['uniform']}（符号 {res['asym_sign']}）；"
                  f" A⁻ の G 展開 {res['Aminus'][0]}" + (f"（/(𝐭²−1) で {res['Aminus'][1]}）" if res['Aminus'][1] else "")
                  + f"、A⁺ の G 展開 {res['Aplus'][0]}" + (f"（/(𝐭²−1) で {res['Aplus'][1]}）" if res['Aplus'][1] else ""), flush=True)
            if verbose:
                for m, (ok, sg, q_) in sorted(res["rows"].items(), reverse=True):
                    print(f"      m={m:2d}: d_m = " + (f"(1 − w²) × [{q_}]" if ok else f"{q_}（割り切れない）") + f"  符号 {sg}", flush=True)
    return out


def test_c2(label, X):
    """c = 2 の部分 X：(ν² − t)(q − t) で割った商の Laurent 性と、係数の符号（全体とブロックごと）。"""
    if sp.cancel(X) == 0:
        print(f"== {label}：0", flush=True); return
    Q = sp.cancel(X / ((nu**2 - t) * (r**2 - s**2)))
    lau = sp.Poly(sp.fraction(Q)[1], r, s, nu).is_monomial
    print(f"== {label}：(𝐚²−𝐪²)(𝐭²−1) で割り切れる {lau}", flush=True)
    if not lau:
        Q2 = sp.cancel(X / (r**2 - s**2)); lau2 = sp.Poly(sp.fraction(Q2)[1], r, s, nu).is_monomial
        print(f"   (𝐭²−1) だけで割り切れる {lau2}", flush=True)
        if not lau2:
            return
        Q = Q2
    cs = sp.Poly(sp.expand(sp.fraction(Q)[0]), r, s, nu).coeffs()
    print(f"   商の全係数の符号 {sign_of(cs)}（{len(cs)} 項）", flush=True)
    for beta, blk in sorted(nu_blocks(Q).items()):
        cs = sp.Poly(sp.expand(sp.fraction(sp.cancel(blk))[0]), r, s).coeffs()
        print(f"   ν^{beta}: {sign_of(cs)}（{len(cs)} 項）", flush=True)


# ---- 検査 3（整数計画、任意） -------------------------------------------------------------------------
E1 = {(1, 1): 1, (1, -1): 1, (0, 2): 1, (0, 0): -2, (0, -2): 1}       # z_q² + z_t² = G(w + w⁻¹) + (w − w⁻¹)²


def gw_poly(e):
    out = {}
    for tm in sp.Add.make_args(sp.expand(e)):
        c, mon = tm.as_coeff_Mul()
        dg = sp.Integer(mon.as_coeff_exponent(G)[1]) if mon != 1 else sp.Integer(0)
        jw = sp.Integer(mon.as_coeff_exponent(W)[1]) if mon != 1 else sp.Integer(0)
        assert jw % 2 == 0 and c.is_integer, (e, tm)
        out[(int(dg), int(jw // 2))] = out.get((int(dg), int(jw // 2)), 0) + int(c)
    return {k: v for k, v in out.items() if v}


def decompose(Xvw):
    """X = A_tot + z_t² B_tot（A_tot, B_tot ∈ Z[G, w^±]）。整数でなければ None。"""
    zt2 = sp.expand((V / W - W / V)**2); zq2 = sp.expand((V * W - 1 / (V * W))**2)
    B = sp.cancel((Xvw - sigma_vw(Xvw)) / (zt2 - zq2))
    if not sp.Poly(sp.fraction(B)[1], V, W).is_monomial:
        return None
    A = sp.expand(Xvw - zt2 * B)
    PA, PB = in_G(sp.expand(A)), in_G(sp.expand(B))
    if PA is None or PB is None:
        return None
    return gw_poly(PA), gw_poly(PB)


def coherent_expansion(A, B, sign):
    """X = (A_tot − e₁B_q) + z_q² B_q + z_t² (B_tot + B_q) の 3 つの係数がすべて sign 側になる B_q が存在するか。"""
    from scipy.optimize import milp, LinearConstraint, Bounds
    keys = set(A) | set(B)
    if not keys:
        return True
    gmax = max(g for g, _ in keys); jmin = min(j for _, j in keys) - 2; jmax = max(j for _, j in keys) + 2
    sup = [(g, j) for g in range(gmax + 1) for j in range(jmin, jmax + 1)]
    idx = {k: i for i, k in enumerate(sup)}
    allkeys = set(keys) | {(g + dg, j + dj) for (g, j) in sup for (dg, dj) in E1}
    rows, lo = [], []
    for k in allkeys:
        row = np.zeros(len(sup))
        for (dg, dj), c in E1.items():
            kk = (k[0] - dg, k[1] - dj)
            if kk in idx:
                row[idx[kk]] -= c
        rows.append(sign * row); lo.append(-sign * A.get(k, 0))
    for k in sup:
        row = np.zeros(len(sup)); row[idx[k]] = 1
        rows.append(sign * row); lo.append(-sign * B.get(k, 0))
        rows.append(sign * row.copy()); lo.append(0)
    res = milp(c=np.zeros(len(sup)), constraints=LinearConstraint(np.array(rows), np.array(lo, dtype=float), np.inf),
               integrality=np.ones(len(sup)), bounds=Bounds(-10**6, 10**6))
    return res.status == 0


def test_ilp(label, X):
    print(f"== {label}（整数計画）", flush=True)
    for beta, blk in sorted(nu_blocks(X).items()):
        d = decompose(to_vw(blk))
        if d is None:
            print(f"  ν^{beta}: 整数分解不可", flush=True); continue
        ok = [sl for sg, sl in ((1, "+"), (-1, "−")) if coherent_expansion(*d, sg)]
        print(f"  ν^{beta}: 同符号の展開 {ok or '無'}", flush=True)


# ---- 対象 ----------------------------------------------------------------------------------------
def knot_one_box(m):
    import rk_b_uniform as Bu
    from rcsknots.refined_kauffman.kauffman_trefoil import r as r0, s as s0, U as U0
    return sp.cancel(Bu.g_true(m).subs({r0: r, s0: s}).subs(U0, nu * s))


def knot_two_box(m):
    """結び目の 2 箱の p₁² 方向の hatted（E を含まない）：(g_(2) + g_(1,1))/z_q。"""
    from rk_b_multibox_sigma import load, r_, s_, nu as nu_
    g, f = load(m)
    sub = lambda e: sp.cancel(e.subs({r_: r, s_: s, nu_: nu}, simultaneous=True))
    return sp.cancel((sub(g["(2,)"]) + sub(g["(1, 1)"])) / zq)


def knot_three_box():
    """trefoil の 3 箱の hatted ĝ_R（E = 0 なので全体が凍結した状態）。"""
    import rk_b_3box_structure as B3
    with contextlib.redirect_stdout(io.StringIO()):
        gh = B3.schur(B3.load("g"), 3, hat=True)
    sub = lambda e: sp.cancel(e.subs({B3.r: r, B3.s: s, B3.nu: nu}, simultaneous=True))
    return {lam: sub(e) for lam, e in gh.items()}


def framed_unknot(x):
    """framed unknot（framing 𝖿 = x）の 2 箱の p₁² 方向の hatted（E を含まない）。"""
    import rk_b_framed_unknot as F
    with contextlib.redirect_stdout(io.StringIO()):
        g, xsec = F.run(x)
    sub = lambda e: sp.cancel(e.subs({F.r_: r, F.s_: s}).subs(F.u_, nu * s))
    return sp.cancel((sub(g[(2,)]) + sub(g[(1, 1)])) / zq)


if __name__ == "__main__":
    args = sys.argv[1:]
    ilp = "--ilp" in args; verbose = "-v" in args
    which = [a for a in args if not a.startswith("-")] or ["one", "unknot", "two", "three"]
    if "one" in which:
        for m in (1, 2):
            X = knot_one_box(m)
            test_c1(f"[1] T(2,{2*m+1}) の c=1（ν 偶）", parity_split(X, True), verbose)
            test_c2(f"[1] T(2,{2*m+1}) の c=2（ν 奇）", parity_split(X, False))
    if "unknot" in which:
        for x in (1, 2, 3, 4):
            X = framed_unknot(x)
            test_c1(f"framed unknot 𝖿={x} の 2 箱 p₁² hatted の c=1（ν 奇）", parity_split(X, False), verbose)
            test_c2(f"framed unknot 𝖿={x} の 2 箱 p₁² hatted の c=2（ν 偶）", parity_split(X, True))
    if "two" in which:
        for m in (1, 2):
            X = knot_two_box(m)
            test_c1(f"T(2,{2*m+1}) の 2 箱 p₁² hatted の c=1（ν 奇）", parity_split(X, False), verbose)
            test_c2(f"T(2,{2*m+1}) の 2 箱 p₁² hatted の c=2（ν 偶）", parity_split(X, True))
            if ilp:
                test_ilp(f"T(2,{2*m+1}) の 2 箱 p₁² hatted の c=1", parity_split(X, False))
    if "three" in which:
        for lam, X in knot_three_box().items():
            test_c1(f"T(2,3) の 3 箱 ĝ_{lam} の c=1（ν 偶）", parity_split(X, True), verbose)
            test_c2(f"T(2,3) の 3 箱 ĝ_{lam} の c=2（ν 奇）", parity_split(X, False))
            if ilp:
                test_ilp(f"T(2,3) の 3 箱 ĝ_{lam} の c=1", parity_split(X, True))
