r"""n = 2（trefoil）の N_λ を結び目のデータ R − 1 の項ごとに分ける（theory §12.26、progress.md 続き 92）。
  結び目のデータ R を 4 つの群に分ける：W（W_λ の R_1, R_2, R_11）、adj（R_adj）、c（composite [α;1] の R_c）、b（composite [α;β] の R_b）。
  群の部分集合 S ごとに、S の群だけ結び目の値、他は基準値 b（b = 1：unknot、b = 0：Wilson loop なしの背景）にした N_λ(S) を計算し
  （枠 w = (t/q)^{1/2} は共通）、包除で T(S) = Σ_{S′ ⊆ S} (−1)^{|S| − |S′|} N_λ(S′) にする。N_λ = Σ_S T(S)。
  b = 1 の T(∅) は枠付き unknot の参照 N^{univ}_λ（rk_f_n2_univ）、b = 0 の T(∅) は κ・c₀ だけの背景（n = 1 の κ の項に当たる）。
  別に、すべての R を b + ε(R − b) にして ε の多項式（次数 ≤ 4）に分け、R − b の次数ごとの寄与 N_d も出す。
  terms [b]  ：T(S) と N_d を計算して保存（data/refined_kauffman/daha/n2_N_rterms_b{b}_m1.pkl）。
  analyze [b]：各 T(S)・N(S)・N_d が Laurent か、分母の因子、ν^{−k} の形式級数のブロック（分母と項数）、ブロックごとの打ち切り。
  background・fock：基準 0 の背景の閉じた形（6 セルの Fock 空間）と結び目の部分 K = N − T(∅) のブロック。
  hyptop・hyptop2：上端 2 ブロック（λ = (2) の ν^{−8}, ν^{−10}）での Ĩ の仮説の族の照合（(i)、(i)+(ii)）。
  hyp12・hyp12ii・hyp12b・hyp12bii・hyp12c・hyp12ci：ν^{−12}（両 λ）への延長（theory §12.27、続き 93）。
使い方：python rk_f_n2_rterms.py terms|analyze|background|fock|hyptop|hyptop2|hyp12|hyp12ii|hyp12b|hyp12bii|hyp12c|hyp12ci [1|0]"""
import sys, os, time, pickle, itertools
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from fractions import Fraction
import rk_f_cover_fixed22 as CF
import rk_c_2box_largeN as L
from rcsknots.refined_kauffman.flint_rf3 import RF, CTX, conv
from rk_f_n2_univ import hatted_N

OUT = "data/refined_kauffman/daha"
GROUPS = ("W", "adj", "c", "b")
LAMS = ("2", "11")
_ORIG = (CF.load_poly_A, L.Rdata_nu)
r_, s_, u_ = sp.symbols("r s u", positive=True)


def _group(path):
    name = os.path.basename(path)
    return "adj" if name.startswith("Radj") else "b" if name.startswith("Rcomp22") else "c"


def set_groups(on, eps=1, base=1):
    """on に入る群は base + eps (R − base)、他は base。"""
    load0, rdata0 = _ORIG
    e, b = conv(eps), conv(base)
    CF.load_poly_A = lambda path, unit: b + e * (load0(path, unit) - b) if _group(path) in on else b
    L.Rdata_nu = (lambda c, knot="T2_3": b + e * (rdata0(c, knot) - b)) if "W" in on else (lambda c, knot="T2_3": b)


def subsets():
    return [tuple(S) for k in range(5) for S in itertools.combinations(GROUPS, k)]


def name(S):
    return "+".join(S) or "∅"


def to_sym(x):
    cv = lambda p: sp.Poly.from_dict({k: int(v) for k, v in p.to_dict().items()}, r_, s_, u_).as_expr() if not p.is_zero() else sp.Integer(0)
    return cv(x.n) / cv(x.d)


def rf_dump(x):
    return ({k: int(v) for k, v in x.n.to_dict().items()}, {k: int(v) for k, v in x.d.to_dict().items()})


def rf_load(nd):
    return RF(CTX.from_dict(nd[0]), CTX.from_dict(nd[1]))


def terms(base=1):
    t0 = time.time()
    NS = {}
    for S in subsets():
        set_groups(S, base=base)
        Ns, _ = hatted_N()
        NS[S] = Ns
        print(f"N({name(S)})：Laurent (2) {Ns['2'].is_laurent()}、(11) {Ns['11'].is_laurent()}、0 {Ns['2'].is_zero()}（{time.time() - t0:.0f} 秒）", flush=True)
    T = {}
    for S in subsets():
        T[S] = {}
        for lab in LAMS:
            acc = CF.ZERO
            for k in range(len(S) + 1):
                for Sp in itertools.combinations(S, k):
                    acc = acc + conv((-1)**(len(S) - k)) * NS[Sp][lab]
            T[S][lab] = acc
    for lab in LAMS:
        tot = sum((T[S][lab] for S in subsets()), CF.ZERO)
        print(f"λ = {lab}：Σ_S T(S) = N {(tot - NS[GROUPS][lab]).is_zero()}", flush=True)
    vals = {}
    for e in range(5):
        set_groups(GROUPS, eps=e, base=base)
        vals[e], _ = hatted_N()
    Vm = sp.Matrix([[sp.Integer(e)**d for d in range(5)] for e in range(5)]).inv()
    Nd = {d: {lab: sum((conv(Fraction(int(Vm[d, e].p), int(Vm[d, e].q))) * vals[e][lab] for e in range(5)), CF.ZERO) for lab in LAMS} for d in range(5)}
    set_groups(GROUPS, eps=5, base=base)
    chk, _ = hatted_N()
    for lab in LAMS:
        tot = sum((conv(5**d) * Nd[d][lab] for d in range(5)), CF.ZERO)
        print(f"λ = {lab}：ε = 5 で次数 ≤ 4 の多項式が再現 {(tot - chk[lab]).is_zero()}、N_d は Laurent "
              + " ".join(f"d={d}:{Nd[d][lab].is_laurent()}" for d in range(5)), flush=True)
    out = {"T": {(S, lab): rf_dump(T[S][lab]) for S in subsets() for lab in LAMS},
           "N": {(S, lab): rf_dump(NS[S][lab]) for S in subsets() for lab in LAMS},
           "Nd": {(d, lab): rf_dump(Nd[d][lab]) for d in range(5) for lab in LAMS}}
    f = f"{OUT}/n2_N_rterms_b{base}_m1.pkl"
    pickle.dump(out, open(f, "wb"))
    print(f"保存：{f}（{time.time() - t0:.0f} 秒）", flush=True)


def load(base=1):
    d = pickle.load(open(f"{OUT}/n2_N_rterms_b{base}_m1.pkl", "rb"))
    return {k: {kk: rf_load(v) for kk, v in d[k].items()} for k in d}


def u_coeffs(p):
    """fmpz_mpoly（r, s, u）→ {u の冪: RF（r, s）}。"""
    out = {}
    for (a, b, c), v in p.to_dict().items():
        out[c] = out.get(c, {})
        out[c][(a, b, 0)] = int(v)
    return {c: RF(CTX.from_dict(dd)) for c, dd in out.items()}


def series_u(x, kmax):
    """x（RF）の u = ∞ での形式級数 Σ_k c_k u^{−k}（k ≤ kmax）。{k: RF（r, s）}。"""
    nc, dc = u_coeffs(x.n), u_coeffs(x.d)
    dn, dd = max(nc), max(dc)
    nc = {dn - k: v for k, v in nc.items()}; dc = {dd - k: v for k, v in dc.items()}   # u^{-j} の係数
    c = {}
    for k in range(kmax + 1 + (dd - dn)):
        v = nc.get(k, CF.ZERO) - sum((dc[j] * c[k - j] for j in dc if j >= 1 and k - j in c), CF.ZERO)
        c[k] = v / dc[0]
    return {k + dd - dn: v for k, v in c.items() if k + dd - dn <= kmax}


def nu_block(c_u, k):
    """u^{−k} の係数 → ν^{−k} の係数（ν = u/s）。"""
    return c_u * CF.S**(-k)


def den_factor(x):
    return str(sp.factor(to_sym(RF(x.d)))).replace("**", "^")


def analyze(base=1, kmax=26):
    D = load(base)
    T, N, Nd = D["T"], D["N"], D["Nd"]
    print(f"基準 R = {base}", flush=True)
    print("(1) 各 T(S)・N(S)・N_d：Laurent か、分子の項数、分母", flush=True)
    for lab in LAMS:
        print(f"  λ = {lab}", flush=True)
        for S in subsets():
            x = T[(S, lab)]
            print(f"    T({name(S)})：" + ("0" if x.is_zero() else f"Laurent {x.is_laurent()}、分子 {len(x.n.coeffs())} 項、分母 {den_factor(x)}"), flush=True)
        print("    N(S) が Laurent：" + "、".join(f"{name(S)}:{N[(S, lab)].is_laurent()}" for S in subsets()), flush=True)
        for d in range(5):
            x = Nd[(d, lab)]
            print(f"    N_{d}（R − {base} の次数 {d}）：" + ("0" if x.is_zero() else f"Laurent {x.is_laurent()}、分子 {len(x.n.coeffs())} 項、分母 {den_factor(x)}"), flush=True)
    print(f"(2) ν^(-k) の形式級数のブロック（k ≤ {kmax}）：項ごとに [分母の因子 / 分子の項数]、和の項数", flush=True)
    ser = {}
    for lab in LAMS:
        for S in subsets():
            if not T[(S, lab)].is_zero():
                ser[(S, lab)] = {k: nu_block(v, k) for k, v in series_u(T[(S, lab)], kmax).items()}
        for d in range(5):
            if not Nd[(d, lab)].is_zero():
                ser[(d, lab)] = {k: nu_block(v, k) for k, v in series_u(Nd[(d, lab)], kmax).items()}
    pickle.dump({k: {kk: rf_dump(v) for kk, v in d.items()} for k, d in ser.items()}, open(f"{OUT}/n2_N_rterms_series_b{base}_m1.pkl", "wb"))
    for lab in LAMS:
        for kind, keys in (("T(S)", [(S, lab) for S in subsets()]), ("N_d", [(d, lab) for d in range(5)])):
            print(f"  λ = {lab}、{kind}", flush=True)
            for k in range(0, kmax + 1, 2):
                row, tot = [], CF.ZERO
                for key in keys:
                    c = ser.get(key, {}).get(k, CF.ZERO)
                    if c.is_zero():
                        continue
                    tot = tot + c
                    lab_ = name(key[0]) if isinstance(key[0], tuple) else f"d{key[0]}"
                    row.append(f"{lab_}[{den_factor(c) if not c.is_laurent() else 'L'}/{len(c.n.coeffs())}]")
                print(f"    ν^-{k}：和 {len(tot.n.coeffs())} 項（Laurent {tot.is_laurent()}）  " + " ".join(row), flush=True)


# ---- 基準 R = 0 の背景の項と結び目の部分 ----------------------------------------------------------------
def rs_poly(x, nu_sym=None):
    """RF（r, s；u を含まない）→ sympy の多項式（Laurent なら単項式で割る）。"""
    e = to_sym(x)
    return sp.expand(sp.cancel(e))


def monos(x):
    """Laurent な RF → {(a, b): c}（r^a s^b、分母の単項式を引く）。"""
    assert x.is_laurent()
    (dm, dc), = x.d.to_dict().items()
    return {(a - dm[0], b - dm[1]): Fraction(int(c), int(dc)) for (a, b, _), c in x.n.to_dict().items()}


def background(kmax=26):
    """基準 R = 0 の背景 T(∅)（κ・c₀ だけ）の閉じた形とブロック、結び目の部分 N − T(∅) のブロック、結び目の項の部分和の Laurent 性。"""
    D = load(0)
    T, N = D["T"], D["N"]
    nu = sp.Symbol("nu", positive=True)
    for lab in LAMS:
        x = T[((), lab)]
        e = sp.factor(to_sym(x).subs(u_, nu * s_))
        print(f"λ = {lab}：背景 T(∅) = {e}", flush=True)
    print("(1) 背景のブロック（ν^-k の係数、r^a s^b の単項式で、係数の絶対値の和と項数、(1 + qt) で割り切れるか）", flush=True)
    bg = {}
    for lab in LAMS:
        ser = {k: nu_block(v, k) for k, v in series_u(T[((), lab)], kmax).items()}
        bg[lab] = ser
        for k in sorted(ser):
            c = ser[k]
            if c.is_zero():
                continue
            m = monos(c)
            P = sp.Poly.from_dict({kk: sp.Rational(v.numerator, v.denominator) for kk, v in m.items()}, r_, s_) if m else None
            div = sp.rem(P, sp.Poly(1 + r_**2 * s_**2, r_, s_)) == 0 if P is not None else None
            amin = min(a for a, b in m); bmin = min(b for a, b in m)
            shifted = {(a - amin, b - bmin): v for (a, b), v in m.items()}
            print(f"  λ = {lab}、ν^-{k}：{len(m)} 項、Σ|c| {sum(abs(v) for v in m.values())}、係数の集合 {sorted(set(m.values()))}、(1+qt) | {div}、"
                  f"r^{amin} s^{bmin} × {sp.factor(sp.Poly.from_dict({kk: sp.Rational(v.numerator, v.denominator) for kk, v in shifted.items()}, r_, s_).as_expr()) if k <= 14 else '…'}", flush=True)
    print("(2) 結び目の部分 K = N − T(∅) のブロック", flush=True)
    kn = {}
    for lab in LAMS:
        K = N[(GROUPS, lab)] - T[((), lab)]
        ser = {k: nu_block(v, k) for k, v in series_u(K, kmax).items()}
        kn[lab] = ser
        for k in sorted(ser):
            c = ser[k]
            if c.is_zero():
                continue
            m = monos(c)
            amin = min(a for a, b in m); bmin = min(b for a, b in m)
            shifted = {(a - amin, b - bmin): v for (a, b), v in m.items()}
            P = sp.Poly.from_dict({kk: sp.Rational(v.numerator, v.denominator) for kk, v in shifted.items()}, r_, s_)
            print(f"  λ = {lab}、ν^-{k}：{len(m)} 項、Σ|c| {sum(abs(v) for v in m.values())}、r^{amin} s^{bmin} × {sp.factor(P.as_expr()) if k <= 14 else str(P.as_expr())[:150]}", flush=True)
    pickle.dump({"bg": {(lab, k): rf_dump(v) for lab, d in bg.items() for k, v in d.items()},
                 "knot": {(lab, k): rf_dump(v) for lab, d in kn.items() for k, v in d.items()}}, open(f"{OUT}/n2_N_bg_knot_m1.pkl", "wb"))
    print("(3) 結び目の項（W, adj, b, W+adj, W+c）の部分和がブロックごとに Laurent になる組", flush=True)
    keys = [S for S in subsets() if S and not T[(S, "2")].is_zero()]
    ser = {}
    for lab in LAMS:
        for S in keys:
            ser[(S, lab)] = {k: nu_block(v, k) for k, v in series_u(T[(S, lab)], kmax).items()}
    for lab in LAMS:
        for n_ in range(1, len(keys)):
            for sub in itertools.combinations(keys, n_):
                ok = all(sum((ser[(S, lab)].get(k, CF.ZERO) for S in sub), CF.ZERO).is_laurent() for k in range(8, kmax + 1, 2))
                if ok:
                    print(f"  λ = {lab}：{' + '.join(name(S) for S in sub)} はブロックごとに Laurent", flush=True)
    print("  （全 5 項の和 = K は Laurent）", flush=True)


# ---- 背景の閉じた形と、上端のブロックでの Ĩ の仮説の照合 ----------------------------------------------------
def fock():
    """背景 T(∅)（基準 R = 0）の閉じた形：−r^7 s^3 (q − t) ν^6 / [(A − t)(Aq − 1)(Aq − t)(Aq − t²)(Aq² − t)(A²q² − t²)]（λ = (2)）、
    λ = (1,1) は分子が r^3 s^7 (q − t) ν^2。A = ν²。セル（ν^{−2} あたり）：t, q^{−1}, tq^{−1}, t²q^{−1}, tq^{−2}、ν^{−4} あたり t²q^{−2}。
    各ブロックが (1 + qt) で割り切れるかも確かめる。"""
    D = load(0)
    q, t, nu, A = CF.q, CF.t, CF.nu, CF.A
    den = (A - t) * (A * q - 1) * (A * q - t) * (A * q - t**2) * (A * q**2 - t) * (A**2 * q**2 - t**2)
    for lab, numer in (("2", -CF.R**7 * CF.S**3 * (q - t) * nu**6), ("11", -CF.R**3 * CF.S**7 * (q - t) * nu**2)):
        x = D["T"][((), lab)]
        print(f"λ = {lab}：T(∅) = −{'r^7 s^3' if lab == '2' else 'r^3 s^7'} (q − t) ν^{6 if lab == '2' else 2} / [(A − t)(Aq − 1)(Aq − t)(Aq − t²)(Aq² − t)(A²q² − t²)]：{(x - numer / den).is_zero()}", flush=True)
        ser = series_u(x, 30)
        one_qt = sp.Poly(1 + r_**2 * s_**2, r_, s_)
        divs = []
        for k in sorted(ser):
            c = nu_block(ser[k], k)
            if c.is_zero():
                continue
            m = monos(c)
            amin, bmin = min(a for a, b in m), min(b for a, b in m)
            P = sp.Poly.from_dict({(a - amin, b - bmin): sp.Rational(v.numerator, v.denominator) for (a, b), v in m.items()}, r_, s_)
            divs.append((k, sp.rem(P, one_qt).is_zero))
        print(f"   ブロックが (1 + qt) で割り切れる：{divs}", flush=True)
    # 2 つの λ の比：T(∅)_(2) / T(∅)_(11) = (r/s)^4 ν^4 = q² A² / t²
    print(f"T(∅)_(2) / T(∅)_(11) = q² A² / t²：{(D['T'][((), '2')] - D['T'][((), '11')] * q**2 * A**2 / t**2).is_zero()}", flush=True)


def hyp_top(with_ii=False):
    """λ = (2) の上端 2 ブロックで Ĩ の仮説を (i)（t = −1 の p_2²・p_4）の一次方程式と照合する。
    ν^{−8}：N = −r^{−7}s³(q − t)（背景だけ、2 状態）。ν^{−10}：N = −r^{−11}s³(q − t)[B₁₀ + K₁₀]、B₁₀ = q²t + qt² + qt + q + t（背景の 5 セル）、
    K₁₀ = (1 + qt)²（結び目の部分、(1 + qt) ⊗ (1 + qt)）。
    仮説：分子の 2 状態に ε_q, ε_t、背景のセルに ε₁..ε₅（対角）、K₁₀ には入れ替え（トレース σ(1 + q²t²)）または対角（σ(1 + ε₆ qt)²）。
    各仮説の E を、E = pre Σ e_m m（pre = s(t − ν²)(ν² + 1)(ν²q + 1)、rk_f_twisted_tensor._n2_model と同じ）の e に直し、
    ν^{−8}, ν^{−10} の単項式だけを含む行で A e = b を検査する。hyptop2 は (ii)（t = q の和）の行も入れる。"""
    import numpy as np
    import rk_f_twisted_tensor as TT
    import rk_b_twisted_multiplet as TM
    mult, conds = TT.n2_setup()
    pre = TT.s * (TT.t - TT.nu**2) * (TT.nu**2 + 1) * (TT.nu**2 * TT.q + 1)
    use = ("i22", "i4", "ii") if with_ii else ("i22", "i4")
    keys, A, b = TT.rows(mult, lambda key: pre * TM.twisted(key[1:], "mono"), {k: conds[k] for k in use})
    print(f"条件 {use}", flush=True)
    top = {key: mult[key] for key in keys if key[0] == "2" and key[4] in (-8, -10)}
    assert not any(key[0] == "11" and key[4] in (-8, -10) for key in keys)
    idx = {key: i for i, key in enumerate(keys)}
    cols = [idx[key] for key in top]
    rows_ = [i for i in range(A.shape[0]) if all(j in cols or A[i, j] == 0 for j in range(A.shape[1]))]
    print(f"ν^-8, ν^-10（λ = (2)）の単項式 {len(top)} 個、それだけを含む (i) の行 {len(rows_)} 本（全 {A.shape[0]} 本）", flush=True)
    q, t = r_**2, s_**2
    B10 = q**2 * t + q * t**2 + q * t + q + t
    N8 = sp.expand(-r_**-7 * s_**3 * (q - t)); N10 = sp.expand(-r_**-11 * s_**3 * (q - t) * (B10 + (1 + q * t)**2))
    chk = {(-8, a, b_): c for (a, b_), c in sp.Poly(N8 * r_**20 * s_**20, r_, s_).as_dict().items()}
    chk = {(k, a - 20, b_ - 20): c for (k, a, b_), c in chk.items()}
    chk.update({(-10, a - 20, b_ - 20): c for (a, b_), c in sp.Poly(N10 * r_**20 * s_**20, r_, s_).as_dict().items()})
    ok = all(int(mult[key]) == chk.get((key[4], key[2], key[3]), 0) for key in top) and len(chk) == len(top)
    print(f"N の ν^-8, ν^-10 が −r^-7 s^3 (q − t)、−r^-11 s^3 (q − t)(B₁₀ + (1 + qt)²) に一致：{ok}", flush=True)
    st = pickle.load(open(f"{OUT}/twisted_n2_pattern.pkl", "rb"))
    rg = st["rg"]
    good = []
    for eq, et in itertools.product((1, -1), repeat=2):
        vac = eq * q - et * t
        for eps in itertools.product((1, -1), repeat=5):
            Bt = eps[0] * q**2 * t + eps[1] * q * t**2 + eps[2] * q * t + eps[3] * q + eps[4] * t
            Kopts = [("swap", sg, None, sg * (1 + q**2 * t**2)) for sg in (1, -1)] + [("diag", sg, e6, sg * (1 + e6 * q * t)**2) for sg in (1, -1) for e6 in (1, -1)]
            for kind, sg, e6, Kt in Kopts:
                E8 = sp.expand(-r_**-7 * s_**3 * vac); E10 = sp.expand(-r_**-11 * s_**3 * vac * (Bt + Kt))
                e = {}
                for k_, E in ((-8, E8), (-10, E10)):
                    for (a, b_), c in sp.Poly(sp.expand(E * r_**20 * s_**20), r_, s_).as_dict().items():
                        e[("2", "mono", a - 20, b_ - 20, k_)] = int(c)
                if any(key not in top for key in e) or any(abs(e.get(key, 0)) > abs(int(mult[key])) or (e.get(key, 0) - int(mult[key])) % 2 for key in top):
                    continue                                     # |e| ≤ |c| と偶奇（N に無い単項式は c = 0）
                ev = np.zeros(len(keys))
                for key, v in e.items():
                    ev[idx[key]] = v
                res = np.abs(A[rows_] @ ev - b[rows_]).max()
                if res < 1e-6:
                    inrg = sum(1 for key in top if key in rg and rg[key][0] <= e.get(key, 0) <= rg[key][1])
                    good.append((eq, et, eps, kind, sg, e6, inrg))
    print(f"(i) の行を満たす仮説 {len(good)} 個（ε_q, ε_t, セルの符号 ε₁..ε₅ (q²t, qt², qt, q, t), K の型, σ, ε₆, pattern の範囲に入る単項式の数/{sum(1 for key in top if key in rg)}）", flush=True)
    for g in good:
        print("   ", g, flush=True)
    uni = {key: rg[key] for key in top if key in rg and rg[key][0] == rg[key][1]}
    print(f"pattern（t = −1、対が最小）で一意な単項式 {len(uni)}/{len(top)}：", {(k[2], k[3], k[4]): (int(mult[k]), v[0]) for k, v in sorted(uni.items())}, flush=True)



# ---- ν^{−12}（両 λ が結合する）への延長 ----------------------------------------------------------------------------
def _n2_rows(with_ii, kmin):
    """n2_setup の行のうち、ν^{−8}..ν^{kmin} の単項式だけを含む行。(keys, A, b, idx, sel, rows_)。"""
    import rk_f_twisted_tensor as TT
    import rk_b_twisted_multiplet as TM
    mult, conds = TT.n2_setup()
    pre = TT.s * (TT.t - TT.nu**2) * (TT.nu**2 + 1) * (TT.nu**2 * TT.q + 1)
    use = ("i22", "i4", "ii") if with_ii else ("i22", "i4")
    keys, A, b = TT.rows(mult, lambda key: pre * TM.twisted(key[1:], "mono"), {k: conds[k] for k in use})
    idx = {key: i for i, key in enumerate(keys)}
    sel = [key for key in keys if key[4] >= kmin]
    cols = {idx[key] for key in sel}
    rows_ = [i for i in range(A.shape[0]) if all(j in cols or A[i, j] == 0 for j in range(A.shape[1]))]
    return mult, keys, A, b, idx, sel, rows_, use


def _e_of(expr, lab, k):
    """sympy の Laurent 多項式（r, s）→ {key: e}。"""
    e = {}
    for (a, b_), c in sp.Poly(sp.expand(expr * r_**40 * s_**40), r_, s_).as_dict().items():
        e[(lab, "mono", a - 40, b_ - 40, k)] = int(c)
    return e


def _top_E():
    """上端の仮説（hyp_top の 2 通りが与える共通の E）：E₈, E₁₀（λ = (2)）。"""
    q, t = r_**2, s_**2
    return {**_e_of(r_**-7 * s_**3 * (q - t), "2", -8),
            **_e_of(r_**-11 * s_**3 * (q - t) * (q**2 * t + q * t**2 - q * t - q - t - 1 - q**2 * t**2), "2", -10)}


def _ranges(mult, keys, A, b, rows_, idx, unknown, fixed):
    """A e = b（rows_ の行）、|e| ≤ |c|、e ≡ c (mod 2) の整数解の集合で、unknown の各 e の範囲。解が無ければ None。
    変数 y_m ∈ {0..|c_m|}、e_m = c_m − 2 sign(c_m) y_m。"""
    import numpy as np
    from scipy.optimize import milp, LinearConstraint, Bounds
    un = list(unknown)
    c = np.array([float(mult[k]) for k in un]); sg = np.sign(c); ca = np.abs(c)
    Ar = A[rows_][:, [idx[k] for k in un]]
    rhs = b[rows_].copy()
    for k, v in fixed.items():
        rhs -= A[rows_][:, idx[k]] * v
    # Ar (c − 2 sg y) = rhs → Ar diag(−2 sg) y = rhs − Ar c
    M = Ar * (-2 * sg); beq = rhs - Ar @ c
    cons = LinearConstraint(M, beq - 1e-9, beq + 1e-9)
    bnd = Bounds(np.zeros(len(un)), ca)
    out = {}
    for i, k in enumerate(un):
        cc = np.zeros(len(un)); cc[i] = -2 * sg[i]
        lo_ = milp(cc, constraints=cons, bounds=bnd, integrality=np.ones(len(un)))
        if lo_.status != 0:
            return None
        hi_ = milp(-cc, constraints=cons, bounds=bnd, integrality=np.ones(len(un)))
        out[k] = (round(c[i] + lo_.fun), round(c[i] - hi_.fun))
    return out


def _poly_of(e, lab, k):
    """{key: e} の λ = lab、ν^k の部分 → sympy（r, s）。"""
    return sp.expand(sum(v * r_**key[2] * s_**key[3] for key, v in e.items() if key[0] == lab and key[4] == k))


def hyp12(with_ii=False):
    """ν^{−12}（λ = (2) の 3 段目、λ = (1,1) の最上段）への延長。ν^{−8}, ν^{−10} の e は上端の仮説（hyp_top）に固定し、
    ν^{−8}..ν^{−12} の単項式だけを含む行で、ν^{−12} の両 λ の e の範囲（整数計画）を求める。
    次に λ = (1,1) の構造的な仮説の族（K₁₂ = −Y⊗Y − (q − t)·t·(1 + qt)⊗(1 + qt)、Y = (1 + q)(1 + t)(1 + qt)、背景 = 真空）と
    λ = (2) の背景（6 セルの Sym² に誘導される符号）を当て、残る λ = (2) の結び目の部分の e を読む。"""
    mult, keys, A, b, idx, sel, rows_, use = _n2_rows(with_ii, -12)
    q, t = r_**2, s_**2
    fixed = _top_E()
    assert all(key in idx for key in fixed)
    unknown = [key for key in sel if key[4] == -12]
    n2 = sum(1 for k in unknown if k[0] == "2"); n11 = len(unknown) - n2
    print(f"条件 {use}：ν^-8..-12 の単項式 {len(sel)} 個（ν^-12：λ = (2) {n2}、λ = (1,1) {n11}）、それだけを含む行 {len(rows_)} 本（全 {A.shape[0]} 本）", flush=True)
    rg = _ranges(mult, keys, A, b, rows_, idx, unknown, fixed)
    if rg is None:
        print("  上端の E を固定すると ν^-12 の解が無い", flush=True); return
    uni = {k: v[0] for k, v in rg.items() if v[0] == v[1]}
    print(f"  解あり。ν^-12 で一意な単項式 {len(uni)}/{len(unknown)}（λ = (2) {sum(1 for k in uni if k[0] == '2')}/{n2}、λ = (1,1) {sum(1 for k in uni if k[0] == '11')}/{n11}）", flush=True)
    for lab in LAMS:
        print(f"  λ = {lab} の範囲 (r^a s^b: c → [min, max])：", {(k[2], k[3]): (int(mult[k]), v) for k, v in sorted(rg.items()) if k[0] == lab and v[0] != v[1]}, flush=True)
    # ---- λ = (1,1) の仮説の族 ----
    # N₁₂(11) = r^-13 s^3 [−Y² − (q − t) t (1 + qt)² − q t² (q − t)]（最後が背景の真空）
    Y_f = ((1 + q, 1 + q**2), (1 + t, 1 + t**2), (1 + q * t, 1 + q**2 * t**2))   # (因子, ψ² の像)
    def tr(fs, choice):                                                               # choice ∈ {swap, +, −} ごとの F⊗F のトレース
        out = sp.Integer(1)
        for (F, F2), ch in zip(fs, choice):
            out *= F2 if ch == "swap" else _flip(F, ch) ** 2
        return out
    cand = []
    chs = ("swap", "+", "-")
    for ch_Y in itertools.product(chs, repeat=3):
        for ch_R in chs:
            for sY, sR, eq, et in itertools.product((1, -1), repeat=4):
                E = r_**-13 * s_**3 * (-sY * tr(Y_f, ch_Y) - sR * (q - t) * t * tr((Y_f[2],), (ch_R,)) - q * t**2 * (eq * q - et * t))
                e = _e_of(E, "11", -12)
                if any(k not in rg for k in e) or any(abs(e.get(k, 0)) > abs(int(mult[k])) or (e.get(k, 0) - int(mult[k])) % 2 for k in rg if k[0] == "11"):
                    continue
                if all(rg[k][0] <= e.get(k, 0) <= rg[k][1] for k in rg if k[0] == "11"):
                    cand.append(((ch_Y, ch_R, sY, sR, eq, et), e))
    print(f"  λ = (1,1) の仮説の族（Y の 3 因子と補正の (1 + qt) に swap/+/−、符号 σ_Y, σ_R、真空 ε_q, ε_t）で範囲に入るもの {len(cand)} 個", flush=True)
    seenE = {}
    for par, e in cand:
        key = tuple(sorted(e.items()))
        seenE.setdefault(key, []).append(par)
    print(f"    異なる E の数 {len(seenE)}", flush=True)
    # 各 λ=(1,1) の E を固定して λ=(2) の範囲
    for key, pars in seenE.items():
        e11 = dict(key)
        rg2 = _ranges(mult, keys, A, b, rows_, idx, [k for k in unknown if k[0] == "2"], {**fixed, **{k: e11.get(k, 0) for k in unknown if k[0] == "11"}})
        if rg2 is None:
            print(f"    {pars[0]}（ほか {len(pars) - 1}）：λ = (2) の解なし", flush=True); continue
        uni2 = {k: v[0] for k, v in rg2.items() if v[0] == v[1]}
        print(f"    {pars[0]}（ほか {len(pars) - 1}）：λ = (2) の一意な単項式 {len(uni2)}/{n2}；E₁₂(11) = {sp.factor(_poly_of(e11, '11', -12))}", flush=True)
        if len(uni2) == n2:
            E2 = _poly_of(uni2, "2", -12)
            print(f"      E₁₂(2) = {sp.factor(E2)}", flush=True)
        else:
            print(f"      λ = (2) の範囲：", {(k[2], k[3]): (int(mult[k]), v) for k, v in sorted(rg2.items()) if v[0] != v[1]}, flush=True)


def _flip(F, ch):
    """因子 F（1 + m）の対角な符号：+ なら 1 + m、− なら 1 − m。"""
    m = sp.expand(F - 1)
    return 1 + m if ch == "+" else 1 - m



def _E11_12():
    """hyp12 で λ = (2) に解を残す唯一の λ = (1,1) の ν^{−12} の E：全因子に入れ替え、σ_Y = σ_R = −1、真空 ε_q = ε_t = −1。"""
    q, t = r_**2, s_**2
    return _e_of(r_**-13 * s_**3 * ((1 + q**2) * (1 + t**2) * (1 + q**2 * t**2) + (q - t) * t * (1 + q**2 * t**2) + q * t**2 * (q - t)), "11", -12)


def hyp12b(with_ii=False):
    """λ = (2) の ν^{−12}：背景 B₁₂ = −r^{−15}s³(q − t)[Sym²(5 セル) + q²t²] には上端の仮説の符号（真空 (−,−)、セルの対角な符号の 2 通り、
    ν^{−4} のセル q²t² に ε₆）を誘導し、結び目の部分 K₁₂ = −r^{−17}s³(1 + qt)P₂₃ の e_K は |e_K| ≤ |c_K|・偶奇だけで自由にして、
    ν^{−8}..ν^{−12} の行（λ = (1,1) は _E11_12 に固定）で e_K の範囲を求める。"""
    import numpy as np
    from scipy.optimize import milp, LinearConstraint, Bounds
    mult, keys, A, b, idx, sel, rows_, use = _n2_rows(with_ii, -12)
    q, t = r_**2, s_**2
    fixed = {**_top_E(), **_E11_12()}
    D = load(0)
    Bq = sp.expand(to_sym(nu_block(series_u(D["T"][((), "2")], 12)[12], 12)))
    Kq = sp.expand(to_sym(nu_block(series_u(D["N"][(GROUPS, "2")] - D["T"][((), "2")], 12)[12], 12)))
    cells = (q**2 * t, q * t**2, q * t, q, t)
    assert sp.expand(Bq + r_**-15 * s_**3 * (q - t) * (sum(ci * cj for i, ci in enumerate(cells) for cj in cells[i:]) + q**2 * t**2)) == 0
    cK = _e_of(Kq, "2", -12)
    print(f"条件 {use}：行 {len(rows_)} 本。λ = (2) ν^-12：背景 13 単項式、結び目の部分 {len(cK)} 単項式（N は {sum(1 for k in sel if k[0] == '2' and k[4] == -12)}）", flush=True)
    unK = sorted(cK)
    for vname, eps in (("diag (+,+,+,−,−)", (1, 1, 1, -1, -1)), ("swap (+,+,−,−,−)", (1, 1, -1, -1, -1))):
        for e6 in (1, -1):
            Eb = sp.expand(r_**-15 * s_**3 * (q - t) * (sum(eps[i] * eps[j] * ci * cj for i, ci in enumerate(cells) for j, cj in enumerate(cells) if j >= i) + e6 * q**2 * t**2))
            eb = _e_of(Eb, "2", -12)
            # 行：A (e_fixed + e_b + e_K) = b。e_K = c_K − 2 sg y
            c = np.array([float(cK[k]) for k in unK]); sg = np.sign(c); ca = np.abs(c)
            Ar = A[rows_][:, [idx[k] for k in unK]]
            rhs = b[rows_].copy()
            for k, v in {**fixed, **eb}.items():
                rhs = rhs - A[rows_][:, idx[k]] * v
            M = Ar * (-2 * sg); beq = rhs - Ar @ c
            cons = LinearConstraint(M, beq - 1e-9, beq + 1e-9); bnd = Bounds(np.zeros(len(unK)), ca)
            rg = {}; ok = True
            for i, k in enumerate(unK):
                cc = np.zeros(len(unK)); cc[i] = -2 * sg[i]
                lo_ = milp(cc, constraints=cons, bounds=bnd, integrality=np.ones(len(unK)))
                if lo_.status != 0:
                    ok = False; break
                hi_ = milp(-cc, constraints=cons, bounds=bnd, integrality=np.ones(len(unK)))
                rg[k] = (round(c[i] + lo_.fun), round(c[i] - hi_.fun))
            if not ok:
                print(f"  背景の符号 {vname}、ε₆ = {e6:+d}：解なし", flush=True); continue
            uni = {k: v[0] for k, v in rg.items() if v[0] == v[1]}
            print(f"  背景の符号 {vname}、ε₆ = {e6:+d}：解あり、一意 {len(uni)}/{len(unK)}、対の数の範囲（Σ(|c|−|e|)/2）", flush=True)
            print("    範囲 (r^a s^b: c_K → [min, max])：", {(k[2], k[3]): (cK[k], v) for k, v in sorted(rg.items())}, flush=True)
            if len(uni) == len(unK):
                print(f"    E_K = {sp.factor(_poly_of(uni, '2', -12))}", flush=True)



QQ, TT = sp.symbols("Q T")


def _qt_poly(expr):
    """r, s の多項式（指数が偶数）→ Q = r², T = s² の多項式。"""
    d = sp.Poly(sp.expand(expr * r_**40 * s_**40), r_, s_).as_dict()
    return sp.Poly.from_dict({((a - 40) // 2, (b - 40) // 2): v for (a, b), v in d.items()}, QQ, TT)


def hyp12c(with_ii=True):
    """λ = (2) の ν^{−12} の結び目の部分に構造を入れる：K₁₂ = −r^{−17}s³[Y⊗Y + (q − t)R]、Y(q, q) = (1 + q)²(1 + q²)² を満たす
    因子 (1 + q), (1 + t), (1 + qt), (1 + q²), (1 + t²) の積 Y ごとに R = (K + Y²)/(q − t) を取り、
    E_K = σ_Y tr(Y⊗Y)（因子ごとに swap/+/−）+ (q − t) Σ_m e_m m（|e_m| ≤ |R_m|、e_m ≡ R_m）とする。
    上端・λ = (1,1) を固定し、背景（swap, ε₆ = +1）の E_b を加えた E₁₂(2) の全単項式で |e| ≤ |c|（N に無い単項式では e = 0）を課し、
    ν^{−8}..ν^{−12} の行で e_m の範囲を求める。"""
    import numpy as np
    from scipy.optimize import milp, LinearConstraint, Bounds
    mult, keys, A, b, idx, sel, rows_, use = _n2_rows(with_ii, -12)
    q, t = r_**2, s_**2
    fixed = {**_top_E(), **_E11_12()}
    D = load(0)
    Kq = sp.expand(to_sym(nu_block(series_u(D["N"][(GROUPS, "2")] - D["T"][((), "2")], 12)[12], 12)))
    Kp = sp.expand(Kq * r_**17 / s_**3)                                   # −(1 + qt)P₂₃
    cells = (q**2 * t, q * t**2, q * t, q, t); eps = (1, 1, -1, -1, -1)
    eb = _e_of(r_**-15 * s_**3 * (q - t) * (sum(eps[i] * eps[j] * ci * cj for i, ci in enumerate(cells) for j, cj in enumerate(cells) if j >= i) + q**2 * t**2), "2", -12)
    rhs0 = b[rows_].copy()
    for k, v in fixed.items():
        rhs0 = rhs0 - A[rows_][:, idx[k]] * v
    facs = ((1 + q, 1 + q**2), (1 + t, 1 + t**2), (1 + q * t, 1 + q**2 * t**2), (1 + q**2, 1 + q**4), (1 + t**2, 1 + t**4))
    target = sp.expand((1 + q)**2 * (1 + q**2)**2)
    print(f"条件 {use}：行 {len(rows_)} 本。K₁₂(2) = −(1 + qt)P₂₃ = −[Y⊗Y + (q − t)R] の形", flush=True)
    ys = [c for n_ in range(1, 5) for c in itertools.combinations_with_replacement(range(5), n_)
          if sp.expand(sp.prod([facs[i][0] for i in c]).subs(t, q) - target) == 0]
    print(f"  Y の候補 {len(ys)} 個", flush=True)
    key_of = lambda a, b_: ("2", "mono", 2 * a - 17, 2 * b_ + 3, -12)       # q^a t^b（r^-17 s^3 を基準）
    for c in ys:
        Y = sp.prod([facs[i][0] for i in c])
        R = sp.div(_qt_poly(-Kp - Y**2), sp.Poly(QQ - TT, QQ, TT))
        assert R[1].is_zero
        Rm = {m: int(v) for m, v in R[0].as_dict().items()}
        R = R[0].as_expr(); mons = sorted(Rm)
        cR = np.array([float(Rm[m]) for m in mons]); sg = np.sign(cR); ca = np.abs(cR)
        feas = []
        for ch in itertools.product(("swap", "+", "-"), repeat=len(c)):
            trY = sp.Integer(1)
            for i, ch_ in zip(c, ch):
                trY *= facs[i][1] if ch_ == "swap" else _flip(facs[i][0], ch_) ** 2
            for sY in (1, -1):
                eY = _e_of(r_**-17 * s_**3 * sY * trY, "2", -12)
                # 各単項式 key の e = k_key + w_key·y（y は対の変数、e_m = R_m − 2 sg_m y_m）
                kc = {}; wv = {}
                for k, v in list(eb.items()) + list(eY.items()):
                    kc[k] = kc.get(k, 0) + v
                for i, (a, b_) in enumerate(mons):
                    for k, sgn in ((key_of(a + 1, b_), 1), (key_of(a, b_ + 1), -1)):
                        kc[k] = kc.get(k, 0) + sgn * Rm[(a, b_)]
                        wv.setdefault(k, np.zeros(len(mons)))[i] += sgn * (-2 * sg[i])
                allk = sorted(kc)
                Wm = np.array([wv.get(k, np.zeros(len(mons))) for k in allk]); kv = np.array([float(kc[k]) for k in allk])
                cab = np.array([abs(float(mult[k])) if k in idx else 0.0 for k in allk])
                Mrow = sum((np.outer(A[rows_][:, idx[k]], Wm[j]) for j, k in enumerate(allk) if k in idx), np.zeros((len(rows_), len(mons))))
                beq = rhs0 - sum((A[rows_][:, idx[k]] * kv[j] for j, k in enumerate(allk) if k in idx), np.zeros(len(rows_)))
                cons = [LinearConstraint(Mrow, beq - 1e-9, beq + 1e-9), LinearConstraint(Wm, -cab - kv - 1e-9, cab - kv + 1e-9)]
                bnd = Bounds(np.zeros(len(mons)), ca)
                res = milp(np.zeros(len(mons)), constraints=cons, bounds=bnd, integrality=np.ones(len(mons)))
                if res.status != 0:
                    continue
                rg = {}
                for i, m in enumerate(mons):
                    cc = np.zeros(len(mons)); cc[i] = -2 * sg[i]
                    lo_ = milp(cc, constraints=cons, bounds=bnd, integrality=np.ones(len(mons)))
                    hi_ = milp(-cc, constraints=cons, bounds=bnd, integrality=np.ones(len(mons)))
                    rg[m] = (round(cR[i] + lo_.fun), round(cR[i] - hi_.fun))
                nun = sum(1 for v in rg.values() if v[0] == v[1])
                # 補正の入れ替わる対 Σ y の最小と、最小での範囲
                pm = milp(np.ones(len(mons)), constraints=cons, bounds=bnd, integrality=np.ones(len(mons)))
                pmin = round(pm.fun)
                cons2 = cons + [LinearConstraint(np.ones((1, len(mons))), -np.inf, pmin + 1e-9)]
                rg2 = {}
                for i, m in enumerate(mons):
                    cc = np.zeros(len(mons)); cc[i] = -2 * sg[i]
                    lo_ = milp(cc, constraints=cons2, bounds=bnd, integrality=np.ones(len(mons)))
                    hi_ = milp(-cc, constraints=cons2, bounds=bnd, integrality=np.ones(len(mons)))
                    rg2[m] = (round(cR[i] + lo_.fun), round(cR[i] - hi_.fun))
                feas.append((ch, sY, nun, rg, pmin, rg2))
        print(f"  Y = {sp.factor(Y)}：R = {sp.factor(R)}（{len(mons)} 状態、Σ|R| {sum(abs(v) for v in Rm.values())}）、可能な (swap/±, σ_Y) {len(feas)}/{2 * 3**len(c)}", flush=True)
        for ch, sY, nun, rg, pmin, rg2 in feas:
            print(f"     {ch} σ_Y = {sY:+d}：R の e が一意 {nun}/{len(mons)}、範囲 {{(a,b): (R, [min,max])}} = {{{', '.join(f'{m}: ({Rm[m]}, {v})' for m, v in rg.items() if v[0] != v[1])}}}", flush=True)
            nun2 = sum(1 for v in rg2.values() if v[0] == v[1])
            print(f"       補正の対の最小 {pmin}（Σ|R| {sum(abs(v) for v in Rm.values())}）、最小で一意 {nun2}/{len(mons)}、範囲 {{{', '.join(f'{m}: ({Rm[m]}, {v})' for m, v in rg2.items() if v[0] != v[1])}}}", flush=True)
            if nun2 == len(mons):
                ER = sum(v[0] * QQ**m[0] * TT**m[1] for m, v in rg2.items())
                print(f"       対が最小の E_R = {sp.factor(ER)}", flush=True)

if __name__ == "__main__":
    b = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    {"terms": terms, "analyze": analyze, "background": lambda b: background(), "fock": lambda b: fock(), "hyptop": lambda b: hyp_top(), "hyptop2": lambda b: hyp_top(True), "hyp12": lambda b: hyp12(), "hyp12ii": lambda b: hyp12(True), "hyp12b": lambda b: hyp12b(), "hyp12bii": lambda b: hyp12b(True), "hyp12c": lambda b: hyp12c(), "hyp12ci": lambda b: hyp12c(False)}[sys.argv[1]](b)
