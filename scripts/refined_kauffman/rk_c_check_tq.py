r"""① の一般の N の R_λ(U) の t = q での検算（progress.md 続き 32）：trefoil の R_λ(U)|_{t=q, U = ν q^{1/2}} を、NRZ 1310.2240 の [r] 彩色の Kauffman 多項式
（scripts/rk_unrefined_check.nrz_trefoil、正規化済み）の鏡像（ν → 1/ν, t_M → 1/t_M）と比べる。[1^r] は F_{[1^r]}(ν, t) = F_{[r]}(±ν, t^{−1}) で作る（r = 2 は Mariño (4.8) で確認済み）。"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_unrefined_check as rk

r, s = sp.symbols("r s", positive=True)
nu, t = rk.nu, rk.t
U = sp.Symbol("U")


def R_tq(lam):
    c = pickle.load(open(f"data/refined_kauffman/daha/R_U_m1_{''.join(map(str, lam))}.pkl", "rb"))
    e = sum(sp.sympify(v, locals={"r": r, "s": s}) * U**j for j, v in c.items())
    return sp.cancel(e.subs(s, r).subs(U, nu * r).subs(r, t))


mir = lambda e: sp.cancel(e.subs({nu: 1 / nu, t: 1 / t}, simultaneous=True))

if __name__ == "__main__":
    for lam in ((1,), (2,), (1, 1), (3,), (1, 1, 1)):
        n = sum(lam)
        cands = [rk.nrz_trefoil(n, nu, t, sp.Integer(-1))] if len(lam) == 1 else [rk.nrz_trefoil(n, sg * nu, 1 / t, sp.Integer(-1)) for sg in (1, -1)]
        res = []
        for F in cands:
            rat = sp.factor(sp.cancel(R_tq(lam) / mir(sp.cancel(F))))
            res.append(rat)
        print(f"λ={lam}：R_λ|_{{t=q}} / 鏡像(F_λ) = {res}", flush=True)
