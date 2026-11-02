"""n = 2（trefoil）の G の行の解析（続き 111、theory §12.41）。§12.40 の記号版のブロック（n2_rows_sym_m1.pkl）を使う。
  A：ν^{−12} の結び目の部分 K₁₂ = −(s³/r)[X_λ² + 2X_λZ + cZ²] の平方の補完。c の候補ごとに Z を解き、判別式が完全平方になる c を探す
     （λ = (1,1) では §12.27 の Y = (1+q)(1+t)(1+qt) に当たる Z = Y − X が c = q/t で出るはず。λ = (2) の Y の選択（18 通り）を行の構造で決める）。
  B：ν^{−14} の行の 2 次形式を −(s³/r)[2X₀X₁ + β₁X₀²]（X₁ = 2 行目の同じ型の線形結合）に当てはめる（n = 1 の −x^{2δ}R(qx)²/((1−x)(1−tx)(1−qtx)) の n = 2 版）。
  C：G の線形項の核の ν^{−12} → ν^{−14} の移り（q の規則の n = 2 版）。
  D：Z と行の関係（Z / a₂a₄ など）。
実行：PYTHONPATH=src:scripts/refined_kauffman python scripts/refined_kauffman/rk_f_n2_grows.py [A|B|C|D|all]"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_f_n2_rows_sym as RS
from rk_f_n2_rows_sym import r, s, q, t, a2, a4, b4, b6, b8, c4, c6, c8, gadj, g21, g111, g22, g211, g1111, SYMS, LAMS

PKL = "data/refined_kauffman/daha/n2_rows_sym_m1.pkl"
LOC = {str(x): x for x in SYMS}; LOC.update({"r": r, "s": s})
_D = None


def blocks():
    """記号版のブロック {(λ, k): expr}、x11 {k: expr}、adams_coef {λ: expr}、inv_den {k: expr}。"""
    global _D
    if _D is None:
        d = pickle.load(open(PKL, "rb"))
        cv = lambda e: sp.sympify(e, locals=LOC)
        _D = {"N": {k: cv(v) for k, v in d["N"].items()}, "x11": {k: cv(v) for k, v in d["x11"].items()},
              "adams_coef": {k: cv(v) for k, v in d["adams_coef"].items()}, "inv_den": {k: cv(v) for k, v in d["inv_den"].items()}}
    return _D


def X_of(lab):
    """§12.40 の X_λ（行の記号）。"""
    if lab == "2":
        return a2**2 / ((q - 1) * (t - 1)) - b4 / ((q - 1) * (q * t - 1)) - t * c4 / ((t - 1) * (q * t - 1))
    return a2**2 / ((q - 1) * (t - 1)) - q * b4 / ((q - 1) * (q * t - 1)) - c4 / ((t - 1) * (q * t - 1))


def adams_block(lab, k, val):
    """Adams の項 −½ψ²(x₁₁)·(1−t)/(z_q φ_(2)²)·χχ/den の ν^{−k} のブロック（実際の値で）。"""
    D = blocks()
    tot = 0
    for j, v in D["x11"].items():
        m = k - 2 * j
        if m in D["inv_den"]:
            vj = sp.expand(v.subs(val))
            tot += vj.subs({r: r**2, s: s**2}, simultaneous=True) * D["inv_den"][m]
    return sp.expand(sp.cancel(tot * D["adams_coef"][lab]))


def background(lab, k):
    import rk_f_n2_rterms as NR
    D = NR.load(0)
    ser = NR.series_u(D["T"][((), lab)], k)
    return sp.expand(NR.to_sym(NR.nu_block(ser[k], k))) if k in ser else sp.Integer(0)


def knot_block(lab, k, val):
    """結び目の部分 K_k = N_k（実際の値、Adams を含む）− 背景。"""
    D = blocks()
    Nk = sp.expand(sp.cancel(D["N"][(lab, k)].subs(val)))
    return sp.expand(sp.cancel(Nk + adams_block(lab, k, val) - background(lab, k)))


def is_square(e):
    """e（r, s の有理式）が完全平方なら平方根を返す。"""
    f = sp.factor(e)
    rt = sp.sqrt(f)
    rt = sp.powdenest(rt, force=True)
    return rt if sp.cancel(rt**2 - e) == 0 and not rt.has(sp.Pow) or (sp.cancel(rt**2 - e) == 0 and all(p.exp.is_integer for p in rt.atoms(sp.Pow))) else None


def part_A(val):
    print("== A：ν^{−12} の平方の補完 K₁₂ = −(s³/r)[X² + 2XZ + cZ²]", flush=True)
    out = {}
    for lab in LAMS:
        X = sp.cancel(X_of(lab).subs(val))
        K = knot_block(lab, 12, val)
        Dp = sp.cancel(-(K + s**3 / r * X**2) * r / s**3)          # Δ′ = 2XZ + cZ²
        print(f"λ = {lab}：X = {sp.factor(X)}", flush=True)
        print(f"   Δ′ = K₁₂·(−r/s³) − X² = {sp.factor(Dp)}", flush=True)
        for cname, c in (("q/t", q / t), ("t/q", t / q), ("1", 1), ("q", q), ("t", t), ("1/q", 1 / q), ("1/t", 1 / t), ("qt", q * t), ("1/(qt)", 1 / (q * t)), ("q²/t", q**2 / t), ("t²/q", t**2 / q)):
            disc = sp.cancel(X**2 + c * Dp)
            rt = is_square(disc)
            if rt is not None:
                for sg in (1, -1):
                    Z = sp.cancel((-X + sg * rt) / c)
                    Y = sp.cancel(X + Z)
                    print(f"   c = {cname}：完全平方。Z = {sp.factor(Z)}、Y = X + Z = {sp.factor(Y)}", flush=True)
                    out.setdefault(lab, []).append((cname, c, Z, Y))
            else:
                print(f"   c = {cname}：判別式は完全平方でない（{str(sp.factor(disc))[:80]}）", flush=True)
    return out


def quad_monos():
    tops = [a2**2, b4, c4]
    nexts = [a2 * a4, b6, c6]
    return tops, nexts


def part_B():
    print("== B：ν^{−14} の 2 次形式を −(s³/r)[2X₀X₁ + β₁X₀²] に当てはめる（記号、Adams を除く）", flush=True)
    D = blocks()
    tops, nexts = quad_monos()
    monos = sorted({tops[i] * tops[j] for i in range(3) for j in range(i, 3)} | {x * y for x in tops for y in nexts}, key=str)
    for lab in LAMS:
        N14 = sp.Poly(sp.expand(D["N"][(lab, 14)]), *SYMS, domain=RS.DOM)
        X0 = X_of(lab)
        be, al, bb, ga, wa, wb, wc = sp.symbols("beta1 alpha beta gamma w_a w_b w_c")
        unk = [be, al, bb, ga, wa, wb, wc]
        X1 = al * a2**2 + wa * a2 * a4 + bb * b4 + wb * b6 + ga * c4 + wc * c6
        Qf = sp.Poly(sp.expand(-(s**3 / r) * (2 * X0 * X1 + be * X0**2)), *SYMS, domain=sp.QQ.frac_field(r, s, *unk))
        eqs = []
        for e in monos:
            cN = sp.sympify(N14.coeff_monomial(e)); cN = cN.as_expr() if hasattr(cN, "as_expr") else cN
            cQ = sp.sympify(Qf.coeff_monomial(e)); cQ = cQ.as_expr() if hasattr(cQ, "as_expr") else cQ
            eqs.append((e, sp.cancel(cQ - cN)))
        print(f"λ = {lab}：2 次形式の単項式 {len(monos)} 個（N₁₄ に無いもの {sum(1 for e, _ in eqs if N14.coeff_monomial(e) == 0)}）、未知数 7", flush=True)
        sol = sp.solve([eq for _, eq in eqs], unk, dict=True)
        if not sol:
            print("   解なし。w を X₀ の重みに固定して単項式ごとの残差（= 0 なら整合）：", flush=True)
            fix = {wa: 2 / ((q - 1) * (t - 1)), wb: -(q if lab == "11" else 1) / ((q - 1) * (q * t - 1)), wc: -(1 if lab == "11" else t) / ((t - 1) * (q * t - 1))}
            for e, eq in eqs:
                print(f"     [{e}] {sp.factor(eq.subs(fix))}", flush=True)
            # 次の行との交差項だけで w を解く
            cross = [eq for e, eq in eqs if any(e.has(y) for y in nexts)]
            print("   交差項（tops × nexts）だけの解：", sp.solve(cross, [wa, wb, wc, be, al, bb, ga], dict=True), flush=True)
        else:
            for k_, v in sol[0].items():
                print(f"   {k_} = {sp.factor(v)}", flush=True)


def part_C():
    print("== C：G の線形項の核：ν^{−12} と ν^{−14} の係数の比", flush=True)
    D = blocks()
    for lab in LAMS:
        N12 = sp.Poly(sp.expand(D["N"][(lab, 12)]), *SYMS, domain=RS.DOM)
        N14 = sp.Poly(sp.expand(D["N"][(lab, 14)]), *SYMS, domain=RS.DOM)
        print(f"λ = {lab}：", flush=True)
        for fam, name in ((gadj, "gadj"), (g21, "g21"), (g111, "g111"), (g22, "g22"), (g211, "g211"), (g1111, "g1111")):
            js = sorted(fam)
            c12 = {j: sp.sympify(N12.coeff_monomial(fam[j])).as_expr() for j in js}
            c14 = {j: sp.sympify(N14.coeff_monomial(fam[j])).as_expr() for j in js}
            for j in js:
                if c12.get(j, 0) != 0 or c14.get(j, 0) != 0:
                    nxt = j + 2
                    r1 = sp.factor(sp.cancel(c14[j] / c12[j])) if c12[j] != 0 else None
                    r2 = sp.factor(sp.cancel(c14[nxt] / c12[j])) if (nxt in c14 and c12[j] != 0) else None
                    print(f"   {name}_{j}：ν^-12 の核 {sp.factor(c12[j])}；ν^-14 で同じ行の核/ν^-12 の核 = {r1}；次の行 {name}_{nxt} の核/ν^-12 の核 = {r2}", flush=True)
        # a2·g の項
        for e in (a2 * g21[8], a2 * g111[8], a2**2 * gadj[6], a2 * g21[10], a2 * g111[10], a2**2 * gadj[8], a4 * g21[8], a4 * g111[8], a2 * a4 * gadj[6]):
            c1 = sp.sympify(N12.coeff_monomial(e)).as_expr(); c2 = sp.sympify(N14.coeff_monomial(e)).as_expr()
            print(f"   [{e}] ν^-12：{sp.factor(c1)}；ν^-14：{sp.factor(c2)}", flush=True)


def part_D(val, outA):
    print("== D：Z と行の関係", flush=True)
    for lab, lst in outA.items():
        for cname, c, Z, Y in lst:
            print(f"λ = {lab}、c = {cname}：Z/(a₂a₄) = {sp.factor(sp.cancel(Z / (a2 * a4).subs(val)))}、Z/a₂² = {sp.factor(sp.cancel(Z / (a2**2).subs(val)))}、"
                  f"Z/b₆ = {sp.factor(sp.cancel(Z / val[b6]))}、Z/c₆ = {sp.factor(sp.cancel(Z / val[c6]))}、Y/a₂ = {sp.factor(sp.cancel(Y / val[a2]))}", flush=True)


# ---- flint による高速版（sympy の cancel が ν^{−14} で止まるため） ----
from rcsknots.refined_kauffman.flint_rf3 import RF, from_sympy as _fs, ONE as _ONE


def _rf(e):
    return _fs(sp.sympify(e))


def rows_rf(val):
    return {k: _rf(v) for k, v in val.items()}


def block_rf(lab, k, valrf, poly=None):
    """記号ブロックに行の値（RF）を代入（Adams・背景を含まない）。"""
    P = poly if poly is not None else sp.Poly(sp.expand(blocks()["N"][(lab, k)]), *SYMS, domain=RS.DOM)
    tot = RF.coerce(0)
    for mon, c in P.terms():
        v = _rf(sp.sympify(c).as_expr() if hasattr(c, "as_expr") else c)
        for sym, e in zip(SYMS, mon):
            if e:
                v = v * valrf[sym]**e
        tot = tot + v
    return tot


def adams_rf(lab, k, valrf):
    D = blocks()
    tot = RF.coerce(0)
    for j, v in D["x11"].items():
        m = k - 2 * j
        if m in D["inv_den"]:
            Pj = sp.Poly(sp.expand(v), *SYMS, domain=RS.DOM)
            vj = RF.coerce(0)
            for mon, c in Pj.terms():
                w = _rf(sp.sympify(c).as_expr() if hasattr(c, "as_expr") else c)
                for sym, e in zip(SYMS, mon):
                    if e:
                        w = w * valrf[sym]**e
                vj = vj + w
            tot = tot + vj.adams(2) * _rf(D["inv_den"][m])
    return tot * _rf(D["adams_coef"][lab])


def background_rf(lab, k):
    import rk_f_n2_rterms as NR
    D = NR.load(0)
    ser = NR.series_u(D["T"][((), lab)], k)
    return NR.nu_block(ser[k], k) if k in ser else RF.coerce(0)


def knot_rf(lab, k, valrf):
    return block_rf(lab, k, valrf) + adams_rf(lab, k, valrf) - background_rf(lab, k)


def to_sym(x):
    return x.to_sympy().subs({sp.Symbol("r", positive=True): r, sp.Symbol("s", positive=True): s})


def part_E(val):
    """ν^{−12}・ν^{−14} の結び目の部分を X² の項・(q−t) の行の項・G の項（＋結び目の Adams）に分け、G の項が平方を補完するか（flint 版）。"""
    print("== E：結び目の部分の分解（flint）：K_k = −(s³/r)[X の 2 次形式] + (q−t) の行の項 + (G の項 + Adams_knot)", flush=True)
    valrf = rows_rf(val)
    zero = {k: RF.coerce(0) for k in valrf}
    gsyms = set(gadj.values()) | set(g21.values()) | set(g111.values()) | set(g22.values()) | set(g211.values()) | set(g1111.values())
    tops, nexts = quad_monos()
    out = {}
    for lab in LAMS:
        X0 = sp.cancel(X_of(lab).subs(val)); X0rf = _rf(X0)
        X1 = sp.cancel(X_of(lab).subs({a2: sp.sqrt(2 * a2 * a4), b4: b6, c4: c6}, simultaneous=True).subs(val)); X1rf = _rf(X1)
        for k in (12, 14):
            P = sp.Poly(sp.expand(blocks()["N"][(lab, k)]), *SYMS, domain=RS.DOM)
            cls = {"quad": RF.coerce(0), "qt": RF.coerce(0), "G": RF.coerce(0), "const": RF.coerce(0)}
            for mon, c in P.terms():
                e = sp.Poly.from_dict({mon: 1}, *SYMS).as_expr()
                v = _rf(sp.sympify(c).as_expr() if hasattr(c, "as_expr") else c)
                for sym, ex in zip(SYMS, mon):
                    if ex:
                        v = v * valrf[sym]**ex
                fs = e.free_symbols
                if fs & gsyms:
                    cls["G"] = cls["G"] + v
                elif e == 1:
                    cls["const"] = cls["const"] + v
                elif any(e == x * y for x in tops for y in tops) or any(e == x * y for x in tops for y in nexts) or any(e == x * y for x in nexts for y in nexts):
                    cls["quad"] = cls["quad"] + v
                else:
                    cls["qt"] = cls["qt"] + v
            ad = adams_rf(lab, k, valrf); ad0 = adams_rf(lab, k, zero); bg = background_rf(lab, k)
            K = cls["quad"] + cls["qt"] + cls["G"] + cls["const"] + ad - bg
            Gp = cls["G"] + (ad - ad0)
            print(f"λ = {lab}、ν^-{k}：背景の照合（定数 + 背景の Adams − 背景 = 0）：{(cls['const'] + ad0 - bg).is_zero()}；K は Laurent：{K.is_laurent()}", flush=True)
            print(f"   (q−t) の行の項 = {sp.factor(to_sym(cls['qt']))}", flush=True)
            print(f"   G の項 + Adams_knot：Laurent {Gp.is_laurent()}", flush=True)
            if k == 12:
                Y2 = X0rf * X0rf - Gp * (R_ / (S_ * S_ * S_))
                e = to_sym(Y2); rt = is_square(e)
                print(f"   Y² := X² − (r/s³)(G の項 + Adams_knot) = {sp.factor(e)}；完全平方 {rt is not None}；Y = {sp.factor(rt) if rt is not None else None}；t=q：{sp.factor(e.subs(s, r))}", flush=True)
                C = -cls["qt"] * (R_ / (S_ * S_ * S_))
                print(f"   補正 C（K = −(s³/r)[Y² + C]）= {sp.factor(to_sym(C))}", flush=True)
                out[(lab, "Y2")] = e; out[(lab, "C")] = to_sym(C)
            else:
                # ν^{−14}：2 次形式 − (−(s³/r)·2X₀X₁) の残り、G の項の値、(q−t) の項
                quad_rest = cls["quad"] + X0rf * X1rf * 2 * (S_ * S_ * S_ / R_)
                print(f"   2 次形式 + (s³/r)·2X₀X₁ = {sp.factor(to_sym(quad_rest))}", flush=True)
                print(f"   G の項 + Adams_knot = {sp.factor(to_sym(Gp)) if Gp.is_laurent() else '非 Laurent'}", flush=True)
                tot14 = to_sym(K)
                print(f"   K₁₄ = {sp.factor(tot14)}；t=q：{sp.factor(tot14.subs(s, r))}", flush=True)
                out[(lab, "K14")] = tot14; out[(lab, "quad_rest")] = to_sym(quad_rest); out[(lab, "G14")] = Gp; out[(lab, "qt14")] = to_sym(cls["qt"])
    return out


R_ = RF.coerce(0) + _rf(r)
S_ = RF.coerce(0) + _rf(s)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    val = RS.actual_rows()
    outA = None
    if what in ("A", "all", "D"):
        outA = part_A(val)
    if what in ("B", "all"):
        part_B()
    if what in ("C", "all"):
        part_C()
    if what in ("D", "all"):
        part_D(val, outA)
    if what in ("E", "all"):
        part_E(val)


