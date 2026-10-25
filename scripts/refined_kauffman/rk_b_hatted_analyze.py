r"""hatted の Ĉ_{R,S} を「t = q の部分（ν~ で表す）＋ refined だけの部分」に分け、判定式（report.md §8 の述べ直し）を検査する：
  (i) t = q の部分が ν~ の冪ごとに符号がそろうか、
  (ii) refined だけの部分が普遍的な因子 (ν~² − q)(q − t)（太字で 𝐭²𝐪²(𝐚² − 𝐪²)(𝐭² − 1)）か (q − t) で割り切れ、商が ν~ の冪ごとに符号がそろうか。
入力：data/refined_kauffman/hopf_g_two_multi_{AY}{AZ}.pkl（scripts/rk_b_hopf_ext.py）。"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rcsknots.partitions import partitions, z as z_lam
from rcsknots.characters import chi
from rcsknots.bps import M_matrix
from rk_b_classify import r, s
from rk_k4_crosscap import nu
from rk_b_hopf_ext import adams2

t = s**2
nut = sp.Symbol("nut", positive=True)
loc = {"r": r, "s": s, "nu": nu}


def coherent_by_power(e, var):
    """var の冪ごとに（他の変数の）係数の符号がそろうか。"""
    e = sp.expand(e)
    n, d = sp.fraction(sp.cancel(e))
    P = sp.Poly(sp.expand(n), var, r, s)
    blocks = {}
    for (k, *_), c in P.terms():
        blocks.setdefault(k, []).append(c)
    return all(all(c > 0 for c in cs) or all(c < 0 for c in cs) for cs in blocks.values())


def strip_zq(e):
    """z_q = r − 1/r で割れる限り割る（z_q 自身が符号の混在した因子なので、判定の前に除く）。"""
    zq = r - 1 / r
    k = 0
    while e != 0:
        quo = sp.cancel(e / zq)
        if not sp.Poly(sp.fraction(quo)[1], r, s, nut).is_monomial:
            break
        e, k = quo, k + 1
    return e, k


def analyze(X, name):
    w = sp.cancel(X.subs(nu, nut * s / r))
    w0 = sp.cancel(w.subs(s, r))
    diff = sp.cancel(w - w0)
    w0s, k0 = strip_zq(w0)
    out = {"t=q": f"z_q^{k0} × （ν~ の冪ごとに符号がそろう = {coherent_by_power(w0s, nut)}）"}
    if diff == 0:
        out["refined"] = "0"
    else:
        for fname, F in (("(ν~²−q)(q−t)", (nut**2 - r**2) * (r**2 - s**2)), ("(q−t)", r**2 - s**2)):
            quo = sp.cancel(diff / F)
            if sp.Poly(sp.fraction(quo)[1], r, s, nut).is_monomial:
                qs_, kq = strip_zq(quo)
                out["refined"] = f"{fname} × z_q^{kq} × （ν~ の冪ごとに符号がそろう = {coherent_by_power(qs_, nut)}）"
                break
        else:
            out["refined"] = "普遍的な因子で割れない"
    print(f"  {name}: {out}", flush=True)
    return out


if __name__ == "__main__":
    AY, AZ = int(sys.argv[1]), int(sys.argv[2])
    g = {k: sp.sympify(v, locals=loc) for k, v in pickle.load(open(f"data/refined_kauffman/hopf_g_two_multi_{AY}{AZ}.pkl", "rb")).items()}
    E1 = {k: sp.sympify(v, locals=loc) for k, v in pickle.load(open("data/refined_kauffman/hopf_E_fixed.pkl", "rb")).items()}
    e11_1 = sp.sympify(pickle.load(open("data/refined_kauffman/hopf_eps_multi.pkl", "rb"))["eps_11_1"], locals=loc)
    eps = {((1,), (1,)): E1[(2, 2)], ((2,), (1,)): E1[(2, 4)], ((1, 1), (1,)): e11_1, ((1,), (2,)): E1[(2, 4)], ((1,), (1, 1)): e11_1}
    for nR, nS in ((2, 2), (4, 2)):
        PR, PS = partitions(nR), partitions(nS)
        C, E = {}, {}
        for R in PR:
            for S_ in PS:
                Ev = sum(adams2(rho).get(R, 0) * adams2(sg).get(S_, 0) * eps[(rho, sg)] for rho in partitions(nR // 2) for sg in partitions(nS // 2))
                E[(R, S_)] = Ev
                C[(R, S_)] = sp.cancel(g[(R, S_)] - Ev / (2 * (1 + t)))
        MR, MS = M_matrix(nR, 1 / r).inv(), M_matrix(nS, 1 / r).inv()
        Ch = MR * sp.Matrix([[C[(R, S_)] for S_ in PS] for R in PR]) * MS.T
        Eh = MR * sp.Matrix([[E[(R, S_)] for S_ in PS] for R in PR]) * MS.T
        print(f"== |R| = {nR}, |S| = {nS}", flush=True)
        for i, R in enumerate(PR):
            for j, S_ in enumerate(PS):
                analyze(sp.cancel(Ch[i, j]), f"Ĉ_{{{R},{S_}}}")
        for i, R in enumerate(PR):
            for j, S_ in enumerate(PS):
                print(f"  Ê_{{{R},{S_}}} = {sp.factor(sp.cancel(sp.cancel(Eh[i, j]).subs(nu, nut * s / r)))}", flush=True)
