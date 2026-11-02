"""n = 2（brane 2・image 2、trefoil）の x セクターを composite の因子化で分ける（続き 109、theory §12.39）。
§12.22 の cover_parts（comp：R_b K_b、c1、c0、low）の comp を、§12.36 の R_{[α;β]} = R_α R_β + (q − 1)(t − 1) G_{αβ} で
  RR：R_α R_β K_{[α;β]}（brane⊗image の積 × composite の次元）、G：(q − 1)(t − 1) G_{αβ} K_{[α;β]}（連結部分）
に分け、各項の hatted な λ 成分（λ = (2), (1,1)）の ν のブロックを出す。n = 1 の N = Fock(1 − x^{2δ}R²) + G 部分 の n = 2 版を探す。
実行：PYTHONPATH=src:scripts/refined_kauffman python scripts/refined_kauffman/rk_f_n2_factor_parts.py"""
import sys, os, pickle, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_f_cover_fixed22 as CF
import rk_f_twisted22_lam as TL
import rk_c_2box_largeN as L
import rk_f_frozen_basis as FB

OUT = TL.OUT
LAMS = TL.LAMS
PARTS2 = f"{OUT}/n2_factor_parts_m1.pkl"


def parts_split():
    t0 = time.time()
    pair = CF.pair_values()
    t, zq = CF.t, CF.R - 1 / CF.R
    q = CF.q
    phi = {(1, 1): zq, (2,): CF.R + 1 / CF.R}
    kappa = (1 - CF.A) * (1 - q) / ((1 - t) * (1 - q * CF.A / t))
    Wadj = pair[((1,), (1,))] - kappa
    big = [k for k in pair if sum(k[0]) == 2 and sum(k[1]) == 2]
    w = CF.S / CF.R                                     # (t/q)^{δ/2}、δ = 1
    Rdat = {lam: w**sum(lam) * L.Rdata_nu("".join(map(str, lam)), "T2_3") for lam in ((1, 1), (2,))}
    parts = {"RR": {}, "G": {}, "c1": {}, "c0": {}}
    for k in big:
        key = k if k != ((1, 1), (2,)) else ((2,), (1, 1))
        tag = f"{''.join(map(str, key[0]))}_{''.join(map(str, key[1]))}"
        Rb = CF.load_poly_A(os.path.join(OUT, f"Rcomp22_A_m1_{tag}.pkl"), 2)
        c1, c0, Kb, _ = CF.decomposition_A(*key)
        RR = Rdat[key[0]] * Rdat[key[1]]
        G = (Rb - RR) / ((q - 1) * (t - 1))
        assert G.is_laurent(), f"G_{tag} が Laurent でない"
        parts["RR"][k], parts["G"][k] = RR * Kb, (q - 1) * (t - 1) * G * Kb
        parts["c1"][k], parts["c0"][k] = c1 * Wadj, c0
        assert (parts["RR"][k] + parts["G"][k] + parts["c1"][k] + parts["c0"][k] - pair[k]).is_zero()
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
        Mh = TL.schur_split(ah)
        for lab in LAMS:
            out[(name, lab)] = Mh[(lab, lab)]                       # RF（r, s, u）のまま
    import rk_f_n2_rterms as NR
    pickle.dump({k: NR.rf_dump(v) for k, v in out.items()}, open(PARTS2, "wb"))
    print(f"保存：{PARTS2}（{time.time() - t0:.0f} 秒）", flush=True)
    return out


def analyze(out=None, kmax=26):
    """各項の λ 成分を u = ∞ の形式級数で ν^{−k} のブロックに分け（flint）、状態数と Laurent 性を表示する。"""
    import rk_f_n2_rterms as NR
    if out is None:
        out = {k: NR.rf_load(v) for k, v in pickle.load(open(PARTS2, "rb")).items()}
    names = ["RR", "G", "c1", "c0", "low"]
    r_, s_ = NR.r_, NR.s_
    states = lambda e: sum(abs(v) for v in sp.expand(e).as_coefficients_dict().values())
    for lab in LAMS:
        print(f"== λ = {lab}", flush=True)
        ser = {name: NR.series_u(out[(name, lab)], kmax) for name in names}
        allk = sorted({k for name in names for k in ser[name]})
        tot = {k: sum((ser[n].get(k, CF.ZERO) for n in names), CF.ZERO) for k in allk}
        for k in allk:
            cells = {}
            for name in names:
                c = ser[name].get(k)
                if c is None or c.is_zero():
                    cells[name] = 0; continue
                cells[name] = states(NR.to_sym(NR.nu_block(c, k))) if c.is_laurent() else "×"
            tk = tot[k]
            cells["全体"] = states(NR.to_sym(NR.nu_block(tk, k))) if tk.is_laurent() else "×"
            print(f"   ν^-{k}：{cells}", flush=True)
        for name in ("RR", "G"):
            ks = [k for k in sorted(ser[name]) if not ser[name][k].is_zero()][:3]
            for k in ks:
                c = ser[name][k]
                print(f"   {name} の ν^-{k}：" + (str(sp.factor(NR.to_sym(NR.nu_block(c, k)))) if c.is_laurent() else "Laurent でない：分母 " + NR.den_factor(c))[:300], flush=True)


if __name__ == "__main__":
    out = parts_split() if "--analyze" not in sys.argv else None
    analyze(out)
