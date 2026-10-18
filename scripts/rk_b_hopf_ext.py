r"""拡張版の S3：両成分とも多 probe の Hopf（|R| ≤ AY、|S| ≤ AZ）の g_{R,S} = f_{R,S} − ½(Σ N N f^{cov})_{R,S} と、
Adams の予言 E_{R,S} = Σ_{ρ,σ} c_{2;ρ}^R c_{2;σ}^S ε_{ρσ}、hatted の不変量（M^{−1}(q^{−1}) を各成分に）。
使い方：python scripts/rk_b_hopf_ext.py AY AZ（AZ ≤ 2 なら S1 の K_λ（|λ| = 4）を待たずに走る）。"""
import sys, pickle, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp
from rcsknots.kauffman_hopf2 import so_f_general
from rcsknots.kauffman_hopf import qs, ts, Us
from rcsknots.kauffman_cover import qs_, ss_, nus_
from rcsknots.partitions import partitions, z as z_lam
from rcsknots.characters import chi
from rcsknots.bps import M_matrix
from rk_b_classify import classify, r, s
from rk_k4_crosscap import nu

t = s**2
nut = sp.Symbol("nut", positive=True)
loc = {"r": r, "s": s, "nu": nu}


def adams2(rho):
    n = sum(rho)
    return {R: c for R in partitions(2 * n)
            if (c := sum(sp.Rational(chi(rho, m), z_lam(m)) * chi(R, tuple(2 * x for x in m)) for m in partitions(n)))}


def mono(e):
    return sp.Poly(sp.fraction(sp.cancel(e))[1], r, s, nu).is_monomial


if __name__ == "__main__":
    AY, AZ = int(sys.argv[1]), int(sys.argv[2])
    t0 = time.time()
    fSO = {k: v.subs({qs: r**2, ts: s**2, Us: nu * s}) for k, v in so_f_general(AY, AZ).items()}
    print(f"SO 側 so_f_general({AY},{AZ}): {time.time()-t0:.0f}s", flush=True)
    q_, s_, n_ = sp.symbols("q s nu")
    fcLR = {k: sp.sympify(v, locals={"q": q_, "s": s_, "nu": n_}).subs({q_: r**2, s_: s, n_: nu})
            for k, v in pickle.load(open("data/refined_kauffman/hopf_fcovLR_A4.pkl", "rb")).items()}
    g = {k: sp.cancel(v - fcLR.get(k, 0) / 2) for k, v in fSO.items()}
    pickle.dump({k: str(v) for k, v in g.items()}, open(f"data/refined_kauffman/hopf_g_two_multi_{AY}{AZ}.pkl", "wb"))
    odd_bad = [k for k, v in g.items() if (sum(k[0]) % 2 or sum(k[1]) % 2) and not mono(v)]
    print("奇数のセクターで Laurent でないもの:", odd_bad, flush=True)
    # ε_{ρσ}（|ρ|, |σ| ≤ 1 の組み合わせのうち既知のもの）
    E1 = {k: sp.sympify(v, locals=loc) for k, v in pickle.load(open("data/refined_kauffman/hopf_E_fixed.pkl", "rb")).items()}
    e11_1 = sp.sympify(pickle.load(open("data/refined_kauffman/hopf_eps_multi.pkl", "rb"))["eps_11_1"], locals=loc)
    eps = {((1,), (1,)): E1[(2, 2)], ((2,), (1,)): E1[(2, 4)], ((1, 1), (1,)): e11_1,
           ((1,), (2,)): E1[(2, 4)], ((1,), (1, 1)): e11_1}
    C = {}
    for (R, S_), v in g.items():
        if sum(R) % 2 or sum(S_) % 2 or not sum(R) or not sum(S_):
            continue
        rhos = [rho for rho in partitions(sum(R) // 2)]
        sigs = [sg for sg in partitions(sum(S_) // 2)]
        if any((rho, sg) not in eps for rho in rhos for sg in sigs):
            continue
        E = sum(adams2(rho).get(R, 0) * adams2(sg).get(S_, 0) * eps[(rho, sg)] for rho in rhos for sg in sigs)
        C[(R, S_)] = sp.cancel(v - E / (2 * (1 + t)))
        print(f"  C_{{{R},{S_}}} = g − ½E/(1+t)（Adams）: Laurent = {mono(C[(R, S_)])}", flush=True)
    # hatted Ĉ（|R| = 4, |S| = 2 のブロック）
    for nR, nS in ((4, 2), (2, 4)):
        PR, PS = partitions(nR), partitions(nS)
        if not all((R, S_) in C for R in PR for S_ in PS):
            continue
        MR, MS = M_matrix(nR, 1 / r).inv(), M_matrix(nS, 1 / r).inv()
        Ch = MR * sp.Matrix([[C[(R, S_)] for S_ in PS] for R in PR]) * MS.T
        for i, R in enumerate(PR):
            for j, S_ in enumerate(PS):
                v = sp.cancel(Ch[i, j])
                w = sp.cancel(v.subs(nu, nut * s / r))
                print(f"  Ĉ_{{{R},{S_}}}: Laurent = {mono(v)}、ν~ 固定で t 非依存 = {not w.has(s)}", flush=True)
                print(f"      = {sp.factor(w)}", flush=True)
