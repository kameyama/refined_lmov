"""E への適用（続き 113、theory §12.43）：行から決めた Y_(2)（§12.42）で λ = (2) の ν^{−12} の Ĩ を再検査する。
K₁₂(2) = −r^{−17}s³[Ŷ² + (q − t)R]、Ŷ = (1+q)(1+t)(1+qt)² + qt(q − t)、R = q⁸C′（15 状態）。
  part 1：平方の部分に入れ替え（E_Y = σ_Y ψ²Ŷ）、補正 R の e を自由にして、(i)・(i)+(ii) の行で e の範囲（hyp12c と同じ設定）。
  part 2：補正を行の部分 (A) q²t·a₂²q²（(1,1) 型）、(B) t(V′a₂q)²、(C) −2q³ a₂S̃q²、(D) q⁵(V′−t)gadj6 に分け、
          (A)(B) に Sym²/ψ²/恒等、(C) に 0（入れ替え）/±恒等 の規則を課し、(D) を自由にして可否と一意性を見る。
実行：PYTHONPATH=src:scripts/refined_kauffman python scripts/refined_kauffman/rk_f_n2_E_Y2.py [1|2|all]"""
import sys, itertools
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import numpy as np
import sympy as sp
from scipy.optimize import milp, LinearConstraint, Bounds
import rk_f_n2_rterms as NR
from rk_f_n2_rterms import r_, s_, QQ, TT, _e_of, _qt_poly, _KEY12, _rows12, _k12_setup

q, t = r_**2, s_**2
Q, T = QQ, TT


def psi2(P):
    return sp.expand(P.subs({r_: r_**2, s_: s_**2}, simultaneous=True))


def sym2(P):
    return sp.expand((P**2 + psi2(P)) / 2)


def setup():
    K2, B2, EB, _ = _k12_setup()                      # r^{−17}s³ で規格化（q, t の多項式）
    Yh = sp.expand((1 + q) * (1 + t) * (1 + q * t)**2 + q * t * (q - t))
    Rp = sp.div(_qt_poly(-K2 - Yh**2), sp.Poly(Q - T, Q, T))
    assert Rp[1].is_zero
    R = sp.expand(Rp[0].as_expr().subs({Q: q, T: t}))
    # 行の部分（q⁸ C′）
    a2 = (1 + q * t) / q; S = (1 + t + q * t) / q; V = 1 + q + q * t
    gadj6 = -t * (2 + t + 2 * q * t) / q**3
    parts = {"A": sp.expand(q**8 * t * a2**2 / q**4), "B": sp.expand(q**8 * t * V**2 * a2**2 / q**6),
             "C": sp.expand(q**8 * (-2 * a2 * S) / q**3), "D": sp.expand(q**8 * (V - t) * gadj6 / q**3)}
    assert sp.expand(sum(parts.values()) - R) == 0
    return K2, EB, Yh, R, parts


def e_dict(poly):
    """q, t の多項式（整数係数）→ {key: e}（r^{−17}s³ 基準の ν^{−12}、λ = (2)）。"""
    return _e_of(r_**-17 * s_**3 * poly, "2", -12)


def solve(S, fixed_polys, free_poly, with_min=True):
    """fixed_polys：e が決まっている部分（q, t の多項式の和、(q−t) を掛ける前のものは呼び出し側で掛けておく）。
    free_poly：(q − t) × Σ_m R_m m の R（e_m = R_m − 2 sg y_m を自由にする）。
    戻り：None（不可）または (一意な数, 全数, 範囲, 対の最小, 最小での範囲, 偶奇の破れ)。"""
    mult, keys, A, b, idx, rows_, Ar, rhs0, use = S
    kc = {}
    for P in fixed_polys:
        for k, v in e_dict(P).items():
            kc[k] = kc.get(k, 0) + v
    Rm = {m: int(v) for m, v in _qt_poly(free_poly).as_dict().items()} if free_poly != 0 else {}
    mons = sorted(Rm)
    cR = np.array([float(Rm[m]) for m in mons]); sg = np.sign(cR); ca = np.abs(cR)
    wv = {}
    for i, (a, b_) in enumerate(mons):
        for k, sgn in ((_KEY12(a + 1, b_), 1), (_KEY12(a, b_ + 1), -1)):
            kc[k] = kc.get(k, 0) + sgn * Rm[(a, b_)]
            wv.setdefault(k, np.zeros(len(mons)))[i] += sgn * (-2 * sg[i])
    allk = sorted(kc)
    # 偶奇：固定部分 + R の寄与（y = 0 のとき）が c と同じ偶奇か
    parity_bad = [k for k in allk if (k in idx and (kc[k] - mult[k]) % 2 != 0) or (k not in idx and kc[k] % 2 != 0)]
    Wm = np.array([wv.get(k, np.zeros(len(mons))) for k in allk]) if mons else np.zeros((len(allk), 0))
    kv = np.array([float(kc[k]) for k in allk])
    cab = np.array([abs(float(mult[k])) if k in idx else 0.0 for k in allk])
    Mrow = sum((np.outer(Ar[:, idx[k]], Wm[j]) for j, k in enumerate(allk) if k in idx), np.zeros((len(rows_), len(mons))))
    beq = rhs0 - sum((Ar[:, idx[k]] * kv[j] for j, k in enumerate(allk) if k in idx), np.zeros(len(rows_)))
    if not mons:
        ok = np.all(np.abs(beq) < 1e-9) and np.all(np.abs(kv) <= cab + 1e-9) and not parity_bad
        return (0, 0, {}, 0, {}, parity_bad) if ok else None
    cons = [LinearConstraint(Mrow, beq - 1e-9, beq + 1e-9), LinearConstraint(Wm, -cab - kv - 1e-9, cab - kv + 1e-9)]
    bnd = Bounds(np.zeros(len(mons)), ca)
    res = milp(np.zeros(len(mons)), constraints=cons, bounds=bnd, integrality=np.ones(len(mons)))
    if res.status != 0 or parity_bad:
        return None
    def ranges(cons_):
        rg = {}
        for i, m in enumerate(mons):
            cc = np.zeros(len(mons)); cc[i] = -2 * sg[i]
            lo_ = milp(cc, constraints=cons_, bounds=bnd, integrality=np.ones(len(mons)))
            hi_ = milp(-cc, constraints=cons_, bounds=bnd, integrality=np.ones(len(mons)))
            rg[m] = (round(cR[i] + lo_.fun), round(cR[i] - hi_.fun))
        return rg
    rg = ranges(cons)
    nun = sum(1 for v in rg.values() if v[0] == v[1])
    if not with_min:
        return nun, len(mons), rg, None, None, parity_bad
    pm = milp(np.ones(len(mons)), constraints=cons, bounds=bnd, integrality=np.ones(len(mons)))
    pmin = round(pm.fun)
    rg2 = ranges(cons + [LinearConstraint(np.ones((1, len(mons))), -np.inf, pmin + 1e-9)])
    return nun, len(mons), rg, pmin, rg2, parity_bad


def fmt(rg, Rm):
    return "{" + ", ".join(f"{m}: ({Rm[m]}, {v})" for m, v in rg.items() if v[0] != v[1]) + "}"


def part1():
    K2, EB, Yh, R, parts = setup()
    Rm = {m: int(v) for m, v in _qt_poly(R).as_dict().items()}
    print(f"Ŷ = {sp.factor(Yh)}、R = {sp.factor(R)}（{len(Rm)} 状態、Σ|R| = {sum(abs(v) for v in Rm.values())}）", flush=True)
    for with_ii in (False, True):
        S = _rows12(with_ii)
        print(f"== part 1：条件 {S[-1]}、行 {len(S[5])} 本。平方の部分 E_Y = σ_Y ψ²Ŷ、補正 R は自由", flush=True)
        for name, EY in (("ψ²Ŷ", psi2(Yh)), ("Ŷ²（恒等）", sp.expand(Yh**2))):
            for sY in (1, -1):
                out = solve(S, [sY * EY, EB], R)
                if out is None:
                    print(f"   E_Y = {sY:+d}·{name}：解なし", flush=True)
                    continue
                nun, n, rg, pmin, rg2, pb = out
                print(f"   E_Y = {sY:+d}·{name}：R の e が一意 {nun}/{n}、範囲 {fmt(rg, Rm)}", flush=True)
                print(f"      対の最小 {pmin}、最小で一意 {sum(1 for v in rg2.values() if v[0] == v[1])}/{n}、範囲 {fmt(rg2, Rm)}", flush=True)


def part2():
    K2, EB, Yh, R, parts = setup()
    a2q = 1 + q * t; PB = sp.expand((1 + q + q * t) * (1 + q * t)); S_ = 1 + t + q * t
    rules_A = {"Sym²": sp.expand(q**2 * t * sym2(a2q)), "ψ²": sp.expand(q**2 * t * psi2(a2q)), "恒等": parts["A"], "−恒等": -parts["A"]}
    rules_B = {"Sym²": sp.expand(t * sym2(PB)), "ψ²": sp.expand(t * psi2(PB)), "恒等": parts["B"], "−恒等": -parts["B"]}
    rules_C = {"0（入れ替え）": sp.Integer(0), "恒等": parts["C"], "−恒等": -parts["C"]}
    assert sp.expand(rules_A["恒等"] - q**2 * t * a2q**2) == 0 and sp.expand(rules_B["恒等"] - t * PB**2) == 0
    for with_ii in (False, True):
        S = _rows12(with_ii)
        print(f"== part 2：条件 {S[-1]}。E_Y = +ψ²Ŷ、(A)(B)(C) に規則、(D) = q⁵(V′−t)gadj6 の e は自由（{len(_qt_poly(parts['D']).as_dict())} 状態）", flush=True)
        feas = []
        for (na, ea), (nb, eb_), (nc, ec) in itertools.product(rules_A.items(), rules_B.items(), rules_C.items()):
            fixed = [psi2(Yh), EB, sp.expand((q - t) * (ea + eb_ + ec))]
            out = solve(S, fixed, parts["D"])
            if out is None:
                continue
            nun, n, rg, pmin, rg2, pb = out
            feas.append((na, nb, nc, nun, n, rg, pmin, rg2))
        print(f"   可能な規則の組 {len(feas)}/{len(rules_A) * len(rules_B) * len(rules_C)}", flush=True)
        Dm = {m: int(v) for m, v in _qt_poly(parts["D"]).as_dict().items()}
        for na, nb, nc, nun, n, rg, pmin, rg2 in feas:
            print(f"   (A) {na}, (B) {nb}, (C) {nc}：(D) の e が一意 {nun}/{n}、範囲 {fmt(rg, Dm)}；対の最小 {pmin}、最小で一意 {sum(1 for v in rg2.values() if v[0] == v[1])}/{n} {fmt(rg2, Dm)}", flush=True)
        # (D) も規則で：gadj6 = g₁ = (2ρ₀ρ₁ − ρ₁²)/q → 候補 −(V′−t)·[2ρ₀ρ₁（入れ替え → 0 か ±）− ρ₁²（恒等）]
        print("   (D) の規則の候補（(A)(B)(C) は可能な組をすべて試す）：", flush=True)
        rho0, rho1 = (1 + q * t) / q, -t / q
        for nd, ed in (("恒等", parts["D"]), ("−恒等", -parts["D"]), ("ρ₀ρ₁ を 0、ρ₁² 恒等", sp.expand(q**8 * (1 + q + q * t - t) * (-rho1**2 / q) / q**3)),
                       ("ρ₀ρ₁ を 0、ρ₁² を −", sp.expand(q**8 * (1 + q + q * t - t) * (rho1**2 / q) / q**3)), ("0", sp.Integer(0))):
            for na, nb, nc, *_ in feas:
                fixed = [psi2(Yh), EB, sp.expand((q - t) * (rules_A[na] + rules_B[nb] + rules_C[nc] + ed))]
                out = solve(S, fixed, sp.Integer(0), with_min=False)
                if out is not None:
                    print(f"      (D) {nd} と (A) {na}, (B) {nb}, (C) {nc}：全体が (i){'+(ii)' if with_ii else ''} を満たす", flush=True)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("1", "all"):
        part1()
    if what in ("2", "all"):
        part2()
