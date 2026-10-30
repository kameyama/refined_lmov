r"""n = 2（trefoil）の N_λ を結び目のデータ R − 1 の項ごとに分ける（theory §12.26、progress.md 続き 92）。
  結び目のデータ R を 4 つの群に分ける：W（W_λ の R_1, R_2, R_11）、adj（R_adj）、c（composite [α;1] の R_c）、b（composite [α;β] の R_b）。
  群の部分集合 S ごとに、S の群だけ結び目の値、他は基準値 b（b = 1：unknot、b = 0：Wilson loop なしの背景）にした N_λ(S) を計算し
  （枠 w = (t/q)^{1/2} は共通）、包除で T(S) = Σ_{S′ ⊆ S} (−1)^{|S| − |S′|} N_λ(S′) にする。N_λ = Σ_S T(S)。
  b = 1 の T(∅) は枠付き unknot の参照 N^{univ}_λ（rk_f_n2_univ）、b = 0 の T(∅) は κ・c₀ だけの背景（n = 1 の κ の項に当たる）。
  別に、すべての R を b + ε(R − b) にして ε の多項式（次数 ≤ 4）に分け、R − b の次数ごとの寄与 N_d も出す。
  terms [b]  ：T(S) と N_d を計算して保存（data/refined_kauffman/daha/n2_N_rterms_b{b}_m1.pkl）。
  analyze [b]：各 T(S)・N(S)・N_d が Laurent か、分母の因子、ν^{−k} の形式級数のブロック（分母と項数）、ブロックごとの打ち切り。
使い方：python rk_f_n2_rterms.py terms|analyze [1|0]"""
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


if __name__ == "__main__":
    b = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    {"terms": terms, "analyze": analyze, "background": lambda b: background(), "fock": lambda b: fock(), "hyptop": lambda b: hyp_top(), "hyptop2": lambda b: hyp_top(True)}[sys.argv[1]](b)
