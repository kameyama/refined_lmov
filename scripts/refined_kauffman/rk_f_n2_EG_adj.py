"""G_adj の行の Ĩ（続き 114、theory §12.44）。
  part A：n = 1 の trefoil の E（twisted_E_n1_m1.pkl）を §12.35 の形 E(x) = F^ε(1 − x²ψ²R(q²x²)) + x·E_G(x) に分け、Ĝ_j = [E_G(1+tx)]_j（j = 1, 2, 3）を
          q の規則 g_j(−q,t)(−q)^j と比べる（対の補正 D_j）。
  part B：n = 2 の λ = (2) の ν^{−12} の補正の (D) 部分 (q−t)q⁵(V′−t)gadj6 に、E_D = (q−t)·κ^ε·Ĝ（κ^ε：核 V′−t の 4 セルの符号、Ĝ：g₁ の Ĩ の候補）
          を課し、(A)(B)(C) の可能な組（§12.43）と合わせて (i)・(i)+(ii) の可否を見る。
実行：PYTHONPATH=src:scripts/refined_kauffman python scripts/refined_kauffman/rk_f_n2_EG_adj.py [A|B|all]"""
import sys, itertools, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_f_rows_gen as RG
import rk_f_frozen_basis as FB

q, t, y = RG.q, RG.t, RG.y
r_, s_ = RG.r_, RG.s_
OUT = "data/refined_kauffman/daha"


def _rs_dict(e):
    d = sp.Poly(sp.expand(e * r_**60 * s_**60), r_, s_).as_dict()
    return {(a - 60, b - 60): int(v) for (a, b), v in d.items()}


def n1_blocks_E(m=1):
    """n = 1 の N と E のブロック（q, t の多項式、rows_gen の規格化）。E の規格化は |e| ≤ |c|・偶奇を満たすものを選ぶ。"""
    import rk_b_cover_split as C
    import io, contextlib
    u = sp.Symbol("u", positive=True)
    nu = FB.nu
    with contextlib.redirect_stdout(io.StringIO()):
        sec = C.cover_sectors(m, 2)
    tt, qq = s_**2, r_**2
    P = sp.cancel((1 - tt) * sp.cancel(sec["x"].to_sympy().subs(u, nu * s_)))
    Nrs = sp.cancel(P / (s_**2 * (tt - nu**2) * (nu**2 - 1) * (qq - 1) * (nu**2 * qq - 1)))
    blN = FB.nu_blocks(Nrs)
    saved = pickle.load(open(f"{OUT}/twisted_E_n1_m{m}.pkl", "rb"))
    # e の辞書：鍵 "('mono', a, b, c)" = r^a s^b ν^c（N と同じ規格化）、値 e
    Ee = sp.Integer(0)
    for blk_, inner in saved["e"].items():
        for key, v in inner.items():
            _, a, b, c = eval(key)
            Ee += v * r_**a * s_**b * nu**c
    # 古い規格化：e の単項式は N·pre₁（pre₁ = r s (t − ν²)(ν² + 1)(ν²q + 1)）のもの。E = E_old/pre₁ が §12.35 の規格化（N と同じブロック）。
    pre1 = r_ * s_ * (s_**2 - nu**2) * (nu**2 + 1) * (nu**2 * r_**2 + 1)
    Ee = sp.cancel(Ee / pre1)
    assert len(sp.Poly(sp.fraction(Ee)[1], r_, s_, nu).terms()) == 1, sp.fraction(Ee)[1]      # 分母は単項式（Laurent）
    blE = FB.nu_blocks(sp.expand(Ee))
    assert set(blE) <= set(blN), (sorted(blE), sorted(blN))
    for b in blE:
        cN, cE = _rs_dict(blN[b]), _rs_dict(blE[b])
        for mon, v in cE.items():
            c = cN.get(mon, 0)
            assert abs(v) <= abs(c) and (v - c) % 2 == 0, (b, mon, v, c)
    cv = lambda blk, b: sp.expand((blk * r_**(4 + 2 * ((-b - 6) // 2))).subs({r_: sp.sqrt(q), s_: sp.sqrt(t)}))
    N = {(-b - 6) // 2: cv(blk, b) for b, blk in blN.items()}
    E = {(-b - 6) // 2: cv(blk, b) for b, blk in blE.items()}
    return N, E, "e の辞書"


def part_A(m=1):
    N, E, name = n1_blocks_E(m)
    de = 1
    print(f"n = 1 の trefoil：E の規格化 = {name}、ブロック k = {sorted(E)}", flush=True)
    G, A, R = RG.G_poly(m)
    rho = RG._coeffs(R, y); g = RG._coeffs(G, y)
    print(f"  R = {R}、G = {sp.factor(G)}", flush=True)
    # F^ε のセル 1, t, qt に ε = (−1)^{a+b}
    def Fe(kk):
        return sp.expand(sum((-1)**(a_ + b_) * t**(b_ + c_) * q**c_ for a_ in range(kk + 1) for b_ in range(kk + 1 - a_) for c_ in [kk - a_ - b_])) if kk >= 0 else sp.Integer(0)
    psi = {2 * i: sp.expand(c.subs({q: q**2, t: t**2}, simultaneous=True)) for i, c in rho.items()}
    K = max(N)
    base = {kk: sp.expand(Fe(kk) - sum(psi[a] * q**a * Fe(kk - 2 * de - a) for a in psi if kk - 2 * de - a >= 0)) for kk in range(K + 1)}
    EG = {kk: sp.expand(E.get(kk, 0) - base[kk]) for kk in range(K + 1)}
    print(f"  E − base（E_G x^{2 * de - 1} の部分）：{ {kk: sp.factor(v) for kk, v in EG.items()} }", flush=True)
    Gh = {}
    for j in range(1, K - 2 * de + 3):
        kk = 2 * de - 1 + j
        if kk > K:
            break
        Gh[j] = sp.expand(EG[kk] + t * EG.get(kk - 1, 0))
        rule = sp.expand(g.get(j, 0).subs(q, -q) * (-q)**j) if j in g else sp.Integer(0)
        D = sp.expand((Gh[j] - rule) / 2)
        print(f"  j = {j}（k = {kk}）：g_j = {sp.factor(g.get(j, 0))}、Ĝ_j = {sp.factor(Gh[j])}、q の規則 = {sp.factor(rule)}、D_j = (Ĝ_j − 規則)/2 = {sp.factor(D)}", flush=True)
    return Gh, g


def part_B():
    import rk_f_n2_E_Y2 as EY
    from rk_f_n2_E_Y2 import solve, psi2, sym2, _rows12, _qt_poly
    Gh, g = part_A(1)
    qq, tt = EY.q, EY.t
    K2, EB, Yh, R, parts = EY.setup()
    g1 = -tt * (2 + tt + 2 * qq * tt) / qq**3          # gadj6 = g₁
    assert sp.expand(g.get(1, 0).subs({q: qq, t: tt}) / qq**2 - g1) == 0        # rows_gen の g₁ = q²·gadj6
    cells = [sp.Integer(1), qq, qq * tt, -tt]          # V′ − t
    G1hat_n1 = sp.expand(Gh[1].subs({q: qq, t: tt}) / qq**3)   # n = 1 の Ĩ(q·g₁^{RG})/q³ = Ĩ(gadj6)（trefoil では −t²/q³ = −ρ₁²/q：入れ替え）
    Gcands = {"g₁（恒等）": g1, "−g₁": -g1, "g₁(−q,t)（q の規則）": g1.subs(qq, -qq), "−g₁(−q,t)": -g1.subs(qq, -qq),
              "n=1 の Ĝ₁（入れ替え、−ρ₁²/q）": G1hat_n1, "−n=1 の Ĝ₁": -G1hat_n1, "0": sp.Integer(0)}
    a2q = 1 + qq * tt; PB = sp.expand((1 + qq + qq * tt) * (1 + qq * tt))
    rules_A = {"ψ²": sp.expand(qq**2 * tt * psi2(a2q)), "恒等": parts["A"], "−恒等": -parts["A"]}
    rule_B = sp.expand(tt * psi2(PB)); rule_C = sp.Integer(0)
    for with_ii in (False, True):
        S = _rows12(with_ii)
        print(f"== part B：条件 {S[-1]}。(B) ψ²、(C) 0、(A) ∈ {{ψ², ±恒等}}、(D) = (q−t)·q⁵·κ^ε·Ĝ", flush=True)
        n_ok = 0
        for (nG, Gc), eps in itertools.product(Gcands.items(), itertools.product((1, -1), repeat=4)):
            kap = sum(e * c for e, c in zip(eps, cells))
            ED = sp.expand(qq**5 * kap * Gc)
            for nA, eA in rules_A.items():
                fixed = [psi2(Yh), EB, sp.expand((qq - tt) * (eA + rule_B + rule_C + ED))]
                out = solve(S, fixed, sp.Integer(0), with_min=False)
                if out is not None:
                    n_ok += 1
                    print(f"   可：Ĝ = {nG}、κ の符号 {eps}（セル 1, q, qt, −t）、(A) {nA}", flush=True)
        print(f"   可能な組 {n_ok}/{len(Gcands) * 16 * len(rules_A)}", flush=True)
        # (A)(B)(C) を自由（R の (A)(B)(C) 部分の e を自由）にして (D) だけ規則
        print("   (A)(B)(C) の e を自由にして (D) だけ規則：", flush=True)
        free_ABC = sp.expand(parts["A"] + parts["B"] + parts["C"])
        n_ok = 0
        for (nG, Gc), eps in itertools.product(Gcands.items(), itertools.product((1, -1), repeat=4)):
            kap = sum(e * c for e, c in zip(eps, cells))
            ED = sp.expand(qq**5 * kap * Gc)
            out = solve(S, [psi2(Yh), EB, sp.expand((qq - tt) * ED)], free_ABC, with_min=True)
            if out is not None:
                n_ok += 1
                nun, n, rg, pmin, rg2, pb = out
                print(f"   可：Ĝ = {nG}、κ の符号 {eps}：(A)(B)(C) の e が一意 {nun}/{n}、対の最小 {pmin}", flush=True)
        print(f"   可能な組 {n_ok}/{len(Gcands) * 16}", flush=True)


def part_C():
    """全部入れ替えの解を明示する：E₁₂(2) = r^{−17}s³[ψ²Ŷ + E_B + (q−t)(E_A + E_B′ + E_C + E_D)]、
    E_A = q²t·ψ²(1+qt)、E_B′ = t·ψ²((1+q+qt)(1+qt))、E_C = 0、E_D = −q²t²(1 + q − qt − t)（Ĝ(g₁) = −ρ₁²/q、核 V′−t の符号は F^ε 型 (+,+,−)）。"""
    import rk_f_n2_E_Y2 as EY
    from rk_f_n2_E_Y2 import solve, psi2, _rows12, _qt_poly
    qq, tt = EY.q, EY.t
    K2, EB, Yh, R, parts = EY.setup()
    EA = sp.expand(qq**2 * tt * psi2(1 + qq * tt)); EBp = sp.expand(tt * psi2((1 + qq + qq * tt) * (1 + qq * tt)))
    ED = sp.expand(-qq**2 * tt**2 * (1 + qq - qq * tt - tt))
    corr = sp.expand(EA + EBp + ED)
    print(f"全部入れ替えの E の補正の部分 Σ e_m m（R に対する e）= {corr}", flush=True)
    Rm = {m: int(v) for m, v in _qt_poly(R).as_dict().items()}
    em = {m: int(v) for m, v in _qt_poly(corr).as_dict().items()}
    print(f"   R の単項式ごとの (R_m, e_m)：{ {m: (Rm.get(m, 0), em.get(m, 0)) for m in sorted(set(Rm) | set(em))} }", flush=True)
    print(f"   |e| ≤ |R| と偶奇：{all(abs(em.get(m, 0)) <= abs(Rm.get(m, 0)) and (em.get(m, 0) - Rm.get(m, 0)) % 2 == 0 for m in set(Rm) | set(em))}、入れ替わる対の数 Σ(|R| − |e|)/2 = {sum(abs(Rm.get(m, 0)) - abs(em.get(m, 0)) for m in set(Rm) | set(em)) // 2}", flush=True)
    for with_ii in (False, True):
        S = _rows12(with_ii)
        out = solve(S, [psi2(Yh), EB, sp.expand((qq - tt) * corr)], sp.Integer(0), with_min=False)
        print(f"   条件 {S[-1]}：{'満たす' if out is not None else '満たさない'}", flush=True)
    E12 = sp.expand(psi2(Yh) + EB + (qq - tt) * corr)
    print(f"   E₁₂(2)·r¹⁷/s³ = {E12}", flush=True)
    print(f"   E₁₂(2) at t = q：{sp.factor(E12.subs(tt, qq))}；N₁₂(2) at t=q：{sp.factor((-(K2) + 0).subs(tt, qq))}（K₂ は結び目の部分）", flush=True)
    return E12


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what == "A":
        part_A(1)
    elif what == "C":
        part_C()
    else:
        part_B()
