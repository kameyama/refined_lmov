"""n = 1 の E を、行の母関数（theory §12.33）の読みで一般のトーラス結び目に当てる（続き 104）。
下の系列 k < 2δ：Sym^k V のセル 1^a t^b (qt)^c に ε = (−1)^{a+b}（principle と同じ）。
移り目 k = 2δ：Sym^{2δ}V（同じ ε）− ψ²S̃、S̃ = ρ_0 − ρ_1（[1] の最上 2 行。2 本撚りでは簡約 superpolynomial）。
上の系列 k > 2δ：e を自由（|e| ≤ |c|、偶奇）にして、条件 (i)・(i)+(ii) の解の有無、e の一意性、対の最小を整数計画で求める。
実行：PYTHONPATH=src:scripts/refined_kauffman python scripts/refined_kauffman/rk_f_rows_E.py [m [--k=K] | T3_5] ..."""
import sys, io, contextlib
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import numpy as np
import sympy as sp
from scipy.optimize import milp, LinearConstraint, Bounds
import rk_f_twisted_tensor as TT
import rk_f_frozen_basis as FB
import rk_f_rows_gen as RG
import rk_torus_tag as TG
import rk_b_twisted_multiplet as TM

r, s, nu = TT.r, TT.s, TT.nu
q, t = TT.q, TT.t


def g_rule_rows(m, k, j0, de, K):
    """G の q の規則（theory §12.35）：Ĝ_j = g_j(−q, t)(−q)^j（j ≥ j0）を、E のブロックの線形制約
    E_k − base_k + t(E_{k−1} − base_{k−1}) = Ĝ_j（k = 2δ − 1 + j）として返す：{(k, a, b): (右辺の係数)}。
    base = F^ε(x)(1 − x^{2δ}ψ²R(q²x²))、F^ε = 1/((1 + x)(1 + tx)(1 − qtx))。"""
    G, A, R = RG.G_poly(m, k)
    rho = RG._coeffs(R, RG.y); g = RG._coeffs(G, RG.y) if G != 0 else {}
    T_ = max(g) if g else 0
    qq, tt = RG.q, RG.t
    def Fe(kk):
        return sp.expand(sum((-1)**(a_ + b_) * tt**(b_ + c_) * qq**c_ for a_ in range(kk + 1) for b_ in range(kk + 1 - a_) for c_ in [kk - a_ - b_])) if kk >= 0 else sp.Integer(0)
    psi = {2 * i: sp.expand(c.subs({qq: qq**2, tt: tt**2}, simultaneous=True)) for i, c in rho.items()}
    base = {}
    for kk in range(K + 3):
        base[kk] = sp.expand(Fe(kk) - sum(psi[a] * qq**a * Fe(kk - 2 * de - a) for a in psi if kk - 2 * de - a >= 0))
    eqs = {}
    for j in range(j0, max(T_, K - 2 * de + 1) + 3):
        kk = 2 * de - 1 + j
        if kk > K + 2:
            break
        Gh = sp.expand(g.get(j, 0).subs(qq, -qq) * (-qq)**j) if j in g else sp.Integer(0)
        rhs = sp.expand(Gh + base[kk] + tt * base[kk - 1])
        eqs[(kk, j)] = rhs
    return eqs


def principle_rows(m, k=2, fix_upper=None, verbose=True, g_rule_from=None):
    """fix_upper = {k: E_k（q, t の多項式）} で上のブロックを固定できる（他は自由）。
    g_rule_from = j0：G の q の規則（Ĝ_j = g_j(−q,t)(−q)^j、j ≥ j0）を線形制約として加える（§12.35）。
    戻り値：{use: (解の有無, 一意な数, 自由な数, 対の最小, 範囲)}。"""
    N, allc = TT._knot_N(m, k)
    mult = TM.monomials(N)
    blocks_ = FB.nu_blocks(N)
    de, K = TG.delta(m, k), TG.rs(m, k) - 3
    Rr = RG.R_poly(m, k)
    rho = {int(mm[0]): c for mm, c in zip(sp.Poly(Rr, RG.y).monoms(), sp.Poly(Rr, RG.y).coeffs())}
    St = sp.expand(rho[0] - rho[1]).subs({RG.q: q, RG.t: t})
    psi2 = lambda X: X.subs({r: r**2, s: s**2}, simultaneous=True)
    cells = lambda kk: [(a_, b_, c_) for a_ in range(kk + 1) for b_ in range(kk + 1 - a_) for c_ in [kk - a_ - b_]]
    sym = lambda kk: sum(t**(b_ + c_) * q**c_ for a_, b_, c_ in cells(kk))
    E = sp.Integer(0); free_nu = set()
    for b_, blk in blocks_.items():
        kk = (-b_ - 6) // 2
        nk = sp.expand(blk * r**(4 + 2 * kk))
        if kk < 2 * de:
            assert sp.expand(nk - sym(kk)) == 0, f"k={kk} ≠ Sym^k V"
            Ek = sum((-1)**(a_ + b_) * t**(b_ + c_) * q**c_ for a_, b_, c_ in cells(kk))
        elif kk == 2 * de:
            assert sp.expand(nk - (sym(kk) - St**2)) == 0, "移り目 ≠ Sym^{2δ}V − S̃²"
            Ek = sum((-1)**(a_ + b_) * t**(b_ + c_) * q**c_ for a_, b_, c_ in cells(kk)) - psi2(St)
        elif fix_upper is not None and kk in fix_upper:
            Ek = fix_upper[kk]
        else:
            free_nu.add(b_); continue
        E += sp.expand(Ek) * nu**b_ * r**(-(4 + 2 * kk))
    E = sp.expand(E)
    e = TM.monomials(E)
    bad = [k_ for k_ in e if k_ not in mult]
    viol = [k_ for k_ in e if k_ in mult and (abs(e[k_]) > abs(mult[k_]) or (e[k_] - mult[k_]) % 2)]
    name = TG.knot_name(m, k)
    print(f"== {name}（δ = {de}、K = {K}）：下と移り目を規則で固定（台の外 {len(bad)}、|e| ≤ |c|・偶奇の違反 {len(viol)}）、上の {len(free_nu)} ブロックは自由", flush=True)
    pre = s**2 * (t - nu**2) * (nu**2 + 1) * (q + 1) * (nu**2 * q + 1)
    out = {}
    for use in (("i",), ("i", "ii")):
        keys, A, b = TT.rows(mult, lambda k_: pre * TM.twisted(k_, "mono"), {nm: allc[nm] for nm in use})
        un = [k_ for k_ in keys if k_[-1] in free_nu]
        fx = {k_: int(e.get(k_, 0)) for k_ in keys if k_[-1] not in free_nu}
        idx = {k_: i for i, k_ in enumerate(keys)}
        if not un:
            ev = np.array([float(e.get(k_, 0)) for k_ in keys]); res = A @ ev - b
            ok = bool(np.max(np.abs(res)) < 1e-9)
            print(f"   条件 {use}：全ブロック固定、{'満たす' if ok else '満たさない'}", flush=True); out[use] = (ok, None, None, None, None); continue
        c = np.array([float(mult[k_]) for k_ in un]); sg = np.sign(c); ca = np.abs(c)
        if g_rule_from is not None:
            # 追加の等式：各 (k, q^a t^b) について e[k; a, b] + e[k−1; a, b−1] = rhs の係数（右辺に無い単項式は 0）
            eqs = g_rule_rows(m, k, g_rule_from, de, K)
            kof = lambda k_: (-k_[-1] - 6) // 2
            qe = lambda k_: (k_[1] + 4 + 2 * kof(k_)) // 2
            te = lambda k_: k_[2] // 2
            rows_add, rhs_add = [], []
            for (kk, j), rhs_poly in eqs.items():
                cd = sp.expand(rhs_poly).as_coefficients_dict()
                monos = set()
                for k_ in keys:
                    if kof(k_) == kk: monos.add((qe(k_), te(k_)))
                    if kof(k_) == kk - 1: monos.add((qe(k_), te(k_) + 1))
                for mono_, v in cd.items():
                    pm = sp.Poly(mono_, RG.q, RG.t).monoms()[0] if mono_ != 1 else (0, 0)
                    monos.add((int(pm[0]), int(pm[1])))
                for (a_, b_) in monos:
                    row = np.zeros(len(keys))
                    for k_ in keys:
                        if kof(k_) == kk and qe(k_) == a_ and te(k_) == b_: row[idx[k_]] += 1
                        if kof(k_) == kk - 1 and qe(k_) == a_ and te(k_) == b_ - 1: row[idx[k_]] += 1
                    val = float(cd.get(RG.q**a_ * RG.t**b_, 0))
                    rows_add.append(row); rhs_add.append(val)
            if rows_add:
                A = np.vstack([A, np.array(rows_add)]); b = np.concatenate([b, np.array(rhs_add)])
        Ar = A[:, [idx[k_] for k_ in un]]
        rhs = b - sum(A[:, idx[k_]] * v for k_, v in fx.items())
        M = Ar * (-2 * sg); beq = rhs - Ar @ c
        cons = LinearConstraint(M, beq - 1e-9, beq + 1e-9); bnd = Bounds(np.zeros(len(un)), ca); integ = np.ones(len(un))
        if milp(np.zeros(len(un)), constraints=cons, bounds=bnd, integrality=integ).status != 0:
            print(f"   条件 {use}：解なし", flush=True); out[use] = (False, None, None, None, None); continue
        rg = {}
        for i, k_ in enumerate(un):
            cc = np.zeros(len(un)); cc[i] = -2 * sg[i]
            lo_ = milp(cc, constraints=cons, bounds=bnd, integrality=integ); hi_ = milp(-cc, constraints=cons, bounds=bnd, integrality=integ)
            rg[k_] = (round(c[i] + lo_.fun), round(c[i] - hi_.fun))
        pm = milp(np.ones(len(un)), constraints=cons, bounds=bnd, integrality=integ)
        uni = {k_: v[0] for k_, v in rg.items() if v[0] == v[1]}
        print(f"   条件 {use}：解あり。上の e は一意 {len(uni)}/{len(un)}、対の最小 {round(pm.fun)}", flush=True)
        if verbose:
            for b_ in sorted(free_nu, reverse=True):
                kk = (-b_ - 6) // 2
                tosym = lambda d: sum(v * q**((k_[1] + 4 + 2 * kk) // 2) * t**((k_[2]) // 2) for k_, v in d.items() if k_[-1] == b_)
                nfree = sum(1 for k_, v in rg.items() if v[0] != v[1] and k_[-1] == b_)
                print(f"     k={kk}：一意な e の和 = {sp.factor(tosym(uni))}、一意でない単項式 {nfree}"[:400], flush=True)
        out[use] = (True, len(uni), len(un) - len(uni), round(pm.fun), rg)
    return out


if __name__ == "__main__":
    args = sys.argv[1:] or ["3", "1 --k=3"]
    kk = next((int(a[4:]) for a in args if a.startswith("--k=")), 2)
    j0 = next((int(a[5:]) for a in args if a.startswith("--j0=")), None)
    for a in args:
        if a.startswith("--"):
            continue
        principle_rows(int(a) if a.isdigit() else a, kk, g_rule_from=j0, verbose="--quiet" not in args)
