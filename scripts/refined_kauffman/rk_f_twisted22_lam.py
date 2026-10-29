r"""n = 2（brane 2・image 2）の E を brane の表現 λ ∈ {(2), (1,1)} で分けて整数計画で解き直す（theory §12.9、progress.md 続き 73）。
  被覆の分解：f^{cov} の (2,2) の成分（Adams の ½ψ²[x₁₁] p_2 p_2′ を除く）の冪和の係数 a_{ρσ}（p_ρ(Y) p_σ(Y′)）。
    Schur のまま λ で分けると、P_λ = (1 − t) M_{λλ} は z_q でしか割り切れない（z_q³ は和にだけ現れる）。
    冪和の方向ごとに境界の因子が違う（KN の hatted：φ_(1,1) = z_q、φ_(2) = z_q^+ = q^{1/2} + q^{−1/2}）ためで、
    hatted の P̂_λμ = Σ_{ρ,σ} χ_λ(ρ) χ_μ(σ) (1 − t) a_{ρσ} / (z_q φ_ρ φ_σ) を使う（z_q は brane と image を結ぶ境界のモード、
    brane–image (1),(1) の P ∝ z_q と同じ）。p_1² p_1′² では z_q · z_q · z_q = z_q³（§12.6 の B_unt）を再現する。
    不変な状態は P̂_λλ にだけ入り、(2),(1,1) の組は Ĩ で (1,1),(2) と入れ替わる対。
  捻ったトレース【仮説】：不変な状態（brane λ、hatted の電荷 m、符号 ε）は g に
    ½ε m s z_q^+ ψ²[Σ_ρ φ_ρ χ_λ(ρ) p_ρ / z_ρ] / (1 + t) を与える（ψ²φ_(1,1) = z_q z_q^+、ψ²φ_(2) = q + q^{−1}）。
    E_λ = Σ ε m とおくと
      4(1 + t) g_{p_2²} = s z_q (z_q^+)² (E_(2) + E_(1,1))（§12.7 の B_tw を再現）、
      4(1 + t) g_{p_4}  = s z_q^+ (q + q^{−1}) (E_(2) − E_(1,1))。
  条件：(i) t = −1 で上の二つが 4(1 + t) f^{SO}_{(2,2)}、4(1 + t) f^{SO}_{(4)} に等しい（f^{cov} の分母は (1 − t) だけ）、
        (ii) t = q で s z_q (z_q^+)² (E_(2) + E_(1,1)) = −(1 + q) ψ²[g_(2) + g_(1,1)]（§12.7、和だけ）。
  E_λ = Σ_m e^λ_m m（m は P̂_λλ/s の単項式、|e^λ_m| ≤ |c^λ_m|、e ≡ c mod 2）。
使い方：python rk_f_twisted22_lam.py cover（被覆の分解を保存）、python rk_f_twisted22_lam.py solve [ν のブロック ...]。"""
import sys, os, pickle, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp

OUT = "data/refined_kauffman/daha"
COVER = f"{OUT}/largeN_cover_fixed22_lam_m1.pkl"
LAMS = ("2", "11")
CHI = {"2": {(1, 1): 1, (2,): 1}, "11": {(1, 1): 1, (2,): -1}}


def schur_split(a):
    """冪和の係数 a[(ρ, σ)]（ρ, σ ∈ {(1,1), (2,)}）→ Schur の M[(λ, μ)]。χ^(2) = (1, 1)、χ^(1,1) = (1, −1)（ρ = (1,1), (2)）。"""
    return {(l, m): sum(CHI[l][rho] * CHI[m][sig] * v for (rho, sig), v in a.items()) for l in LAMS for m in LAMS}


def cover():
    import rk_f_cover_fixed22 as CF
    t0 = time.time()
    P, x11, Lg = CF.fixed_sector()
    print(f"被覆の log Z：{time.time() - t0:.0f} 秒", flush=True)
    t, zq = CF.t, CF.R - 1 / CF.R
    phi = {(1, 1): zq, (2,): CF.R + 1 / CF.R}
    a = {(rho, sig): Lg.get((rho, sig), CF.ZERO) for rho in ((1, 1), (2,)) for sig in ((1, 1), (2,))}
    a[((2,), (2,))] = a[((2,), (2,))] - CF.half * x11.adams(2)          # Adams の項を除く
    M = schur_split(a)
    print("Y ↔ Y′ の対称性 M_{(2),(1,1)} = M_{(1,1),(2)}:", (M[("2", "11")] - M[("11", "2")]).is_zero(), flush=True)
    print("和の検算 ΣM = 4 a_{(1,1),(1,1)}:", (sum(M.values(), CF.ZERO) - a[((1, 1), (1, 1))] * 4).is_zero(), flush=True)
    ah = {k: (1 - t) * v / (zq * phi[k[0]] * phi[k[1]]) for k, v in a.items()}
    for k, v in ah.items():
        print(f"(1 − t) a_{k} / (z_q φ φ) は Laurent: {v.is_laurent()}、整数係数: {v.is_integer_laurent()}", flush=True)
    Mh = schur_split(ah)
    out = {}
    for lab, key in (("2", ("2", "2")), ("11", ("11", "11")), ("off", ("2", "11"))):
        out["hat" + lab] = str(CF.to_rs_nu(Mh[key])[0])
        out[lab] = str(CF.to_rs_nu((1 - t) * M[key])[0])
    pickle.dump(out, open(COVER, "wb"))
    print(f"保存：{COVER}（{time.time() - t0:.0f} 秒）", flush=True)


PARTS = f"{OUT}/largeN_cover_fixed22_lam_parts_m1.pkl"


def cover_parts():
    """被覆の (2,2) 成分を項ごとに分ける（theory §12.22）。|α| = |β| = 2 の ⟨P_α P_β^∨⟩ = R_b K_b + c₁ W_adj + c₀ は log Z の (2,2) 成分に線形に入る。
    comp：R_b K_b（合成表現、Koike の σ = ∅）、c1：c₁ W_adj（σ = (1)）、c0：c₀（σ = (2), (1,1)）、low：log が引く低い次数の積と Adams の項。
    各項の hatted な成分 P̂_λλ を cover() と同じ規格化で作り、和が cover() の値に一致することを確かめて保存する。"""
    import rk_f_cover_fixed22 as CF
    t0 = time.time()
    pair = CF.pair_values()
    t, zq = CF.t, CF.R - 1 / CF.R
    phi = {(1, 1): zq, (2,): CF.R + 1 / CF.R}
    kappa = (1 - CF.A) * (1 - CF.q) / ((1 - t) * (1 - CF.q * CF.A / t))
    Wadj = pair[((1,), (1,))] - kappa
    big = [k for k in pair if sum(k[0]) == 2 and sum(k[1]) == 2]
    parts = {"comp": {}, "c1": {}, "c0": {}}
    for k in big:
        key = k if k != ((1, 1), (2,)) else ((2,), (1, 1))
        tag = f"{''.join(map(str, key[0]))}_{''.join(map(str, key[1]))}"
        Rb = CF.load_poly_A(os.path.join(OUT, f"Rcomp22_A_m1_{tag}.pkl"), 2)
        c1, c0, Kb, _ = CF.decomposition_A(*key)
        parts["comp"][k], parts["c1"][k], parts["c0"][k] = Rb * Kb, c1 * Wadj, c0
        assert (parts["comp"][k] + parts["c1"][k] + parts["c0"][k] - pair[k]).is_zero()
    rs = ((1, 1), (2,))
    a = {}
    for name, pp in parts.items():
        Z = CF.cover_Z(pp)
        a[name] = {(rho, sig): Z.get((rho, sig), CF.ZERO) for rho in rs for sig in rs}
    low = {k: v for k, v in pair.items() if k not in big}
    Lg = CF.log2(CF.cover_Z(low))
    x11 = Lg[((1,), (1,))]
    a["low"] = {(rho, sig): Lg.get((rho, sig), CF.ZERO) for rho in rs for sig in rs}
    a["low"][((2,), (2,))] = a["low"][((2,), (2,))] - CF.half * x11.adams(2)
    print(f"項ごとの (2,2) 成分：{time.time() - t0:.0f} 秒", flush=True)
    _, _, Lfull = CF.fixed_sector()
    tot = {k: sum((a[n][k] for n in a), CF.ZERO) for k in a["low"]}
    tot[((2,), (2,))] = tot[((2,), (2,))] + CF.half * x11.adams(2)
    print("項の和が log Z の (2,2) 成分に一致:", all((tot[k] - Lfull.get(k, CF.ZERO)).is_zero() for k in tot), flush=True)
    out = {}
    for name in a:
        ah = {k: (1 - t) * v / (zq * phi[k[0]] * phi[k[1]]) for k, v in a[name].items()}
        Mh = schur_split(ah)
        for lab in LAMS:
            out[(name, lab)] = str(CF.to_rs_nu(Mh[(lab, lab)])[0])
    pickle.dump(out, open(PARTS, "wb"))
    print(f"保存：{PARTS}（{time.time() - t0:.0f} 秒）", flush=True)


# ---- 整数計画 ---------------------------------------------------------------------------------------
_r = sp.Symbol("r", positive=True)                     # rk_f_twisted22.r と同じ記号
W22 = (_r - 1 / _r) * (_r + 1 / _r)**2                    # p_2² 方向の捻った境界の因子 z_q (z_q^+)²（§12.7 の B_tw）
W4 = (_r + 1 / _r) * (_r**2 + _r**-2)                     # p_4 方向 z_q^+ (q + q^{−1})


def load():
    import rk_f_twisted22 as T
    r, s, nu = T.r, T.s, T.nu
    d = pickle.load(open(COVER, "rb"))
    P = {k: sp.sympify(d["hat" + k], locals={"r": r, "s": s, "nu": nu}) for k in ("2", "11", "off")}
    u = sp.Symbol("u", positive=True)
    f = pickle.load(open(f"{OUT}/largeN_4box_so_m1.pkl", "rb"))["f"]
    t = s**2
    R = {k: sp.cancel(4 * (1 + t) * sp.sympify(f[k], locals={"r": r, "s": s, "u": u}).subs(u, nu * s))
         for k in ("(2, 2)", "(4,)")}
    return P, R["(2, 2)"], R["(4,)"]


CACHE = "data/cache/twisted22_lam_blocks.pkl"


def blocks():
    """ν のブロック（係数、ν^k を除く）：P̂_λ/s、(i) の目標 R22・R4、(ii) の目標 Eq。SO の f の約分が重いので data/cache に保存する。"""
    import rk_f_twisted22 as T
    import rk_f_frozen_basis as FB
    r, s, nu, I = T.r, T.s, T.nu, sp.I
    loc = {"r": r, "s": s, "nu": nu}
    if os.path.exists(CACHE):
        d = pickle.load(open(CACHE, "rb"))
        return {n: {k: sp.sympify(v, locals=loc) for k, v in b.items()} for n, b in d.items()}
    P, R22, R4 = load()
    _, F = T.sectors()[2]
    Eq = sp.cancel(-(1 + r**2) * T.psi2(F))
    bl = {}
    for lab in ("2", "11", "off"):
        Pk = sp.cancel(P[lab] / s)
        lau = sp.Poly(sp.fraction(Pk)[1], r, s, nu).is_monomial
        bl["P" + lab] = FB.nu_blocks(Pk)
        print(f"P̂_{lab} / s は Laurent: {lau}、ν のブロック {sorted(bl['P' + lab])}", flush=True)
    for name, X, w in (("R22", R22, W22), ("R4", R4, W4)):
        bl[name] = FB.nu_blocks(X)
        q_ = sp.cancel(sp.cancel(X.subs(s, I)) / (I * w))
        print(f"{name}|_(t=−1) / (s × 捻った境界の因子) は Laurent: {sp.Poly(sp.fraction(q_)[1], r, nu).is_monomial}", flush=True)
    bl["Eq"] = FB.nu_blocks(Eq)
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    pickle.dump({n: {k: str(v) for k, v in b.items()} for n, b in bl.items()}, open(CACHE, "wb"))
    return bl


def setup(k, bl, use=("i22", "i4", "ii"), tr=False):
    """ν^k のブロックの変数と条件。use で条件を選ぶ（i22：t = −1 の p_2²、i4：t = −1 の p_4、ii：t = q の和）。
    tr = True：Ĩ が brane の λ を image の λᵀ に移すとする仮説。不変な状態は (λ, λᵀ) の成分（P̂_off）に入り、
    (λ, λ) の成分は (λᵀ, λᵀ) との対になる。"""
    import rk_b_twisted_multiplet as TM
    import rk_f_twisted22 as T
    r, s, nu, I = T.r, T.s, T.nu, sp.I
    get = lambda name: bl[name].get(k, sp.Integer(0)) * nu**k
    mult = {}
    for lab in LAMS:
        for key, c in TM.monomials(get("Poff" if tr else "P" + lab)).items():
            mult[(lab,) + key] = c
    basis = lambda key: s * TM.twisted(key[1:], "mono")
    conds = {"i22": (get("R22"), {s: I}, lambda key: W22),
             "i4": (get("R4"), {s: I}, lambda key: CHI[key[0]][(2,)] * W4),
             "ii": (get("Eq"), {s: r}, lambda key: W22)}
    return mult, basis, [conds[c] for c in use]


def run(ks=None, tr=False):
    import rk_b_twisted_multiplet as TM
    import rk_f_twisted22 as T
    s, I = T.s, sp.I
    bl = blocks()
    allk = sorted(set(bl["Poff"]) if tr else set(bl["P2"]) | set(bl["P11"]))
    print("仮説：" + ("Ĩ は λ → λᵀ（不変な状態は P̂_off）" if tr else "Ĩ は λ → λ（不変な状態は P̂_λλ）"), flush=True)
    extra = sorted((set(bl["R22"]) | set(bl["R4"])) - set(allk))
    vanish = all(sp.expand(sp.cancel(bl[n].get(k, 0)).subs(s, I)) == 0 for n in ("R22", "R4") for k in extra)
    print(f"(i) の目標だけにあるブロック {extra}（t = −1 で 0：{vanish}）", flush=True)
    out = {}
    for k in ks or allk:
        t0 = time.time()
        mult, basis, conds = setup(k, bl, tr=tr)
        res = TM.solve(mult, basis, conds, uniq=True)
        nl = {lab: sum(1 for key in mult if key[0] == lab) for lab in LAMS}
        if res["feasible"]:
            fl = {lab: sum(1 for key in res["free"] if key[0] == lab) for lab in LAMS}
            msg = f"解あり、一意でない変数 (2): {fl['2']}/{nl['2']}、(1,1): {fl['11']}/{nl['11']}"
        else:
            msg = f"解なし（変数 (2): {nl['2']}、(1,1): {nl['11']}）"
        print(f"  ν^{k}：{msg}（{time.time() - t0:.0f} 秒）", flush=True)
        out[k] = res
    return out


def linear_ok(mult, basis, conds):
    """整数性と範囲を外した一次方程式 A e = b が解をもつか（最小二乗の残差で判定）。"""
    import numpy as np
    import rk_b_twisted_multiplet as TM
    keys = sorted(mult)
    bs = [basis(k) for k in keys]
    A, b = [], []
    for target, sub, w in conds:
        a_, b_ = TM.rows_for([x * w(k) for k, x in zip(keys, bs)], target, sub)
        A += a_; b += b_
    A, b = np.array(A), np.array(b)
    x = np.linalg.lstsq(A, b, rcond=None)[0]
    return float(np.abs(A @ x - b).max()) < 1e-6


def diagnose(ks):
    """解の無いブロックで、条件の組ごとに整数計画と一次方程式の解の有無を調べる。"""
    import rk_b_twisted_multiplet as TM
    bl = blocks()
    subsets = (("i22",), ("i4",), ("ii",), ("i22", "i4"), ("i22", "ii"), ("i4", "ii"), ("i22", "i4", "ii"))
    for k in ks:
        row = []
        for use in subsets:
            mult, basis, conds = setup(k, bl, use)
            res = TM.solve(mult, basis, conds, uniq=False)
            lin = linear_ok(mult, basis, conds)
            row.append(f"{'+'.join(use)}: 整数計画 {'○' if res['feasible'] else '×'} 一次 {'○' if lin else '×'}")
        print(f"  ν^{k}：" + "、".join(row), flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "cover":
        cover()
    elif sys.argv[1] == "parts":
        cover_parts()
    elif sys.argv[1] == "diagnose":
        diagnose([int(a) for a in sys.argv[2:]])
    else:
        tr = "--tr" in sys.argv
        out = run([int(a) for a in sys.argv[2:] if a != "--tr"] or None, tr)
        pickle.dump(out, open(f"{OUT}/twisted_E_lam{'_tr' if tr else ''}_n2_m1.pkl", "wb"))
