r"""タスク C (b)：E を含まない 2 箱の方向（Schur の p_1² 方向 X = g_(2)+g_(1,1) と、P_λ(Y; t, q) 基底の P_(1,1) の係数 x_(1,1)）について、冪和の定数倍型の基底（÷ φ）で Laurent 性・σ 対称性・𝐚 ブロックの符号・述べ直しの普遍因子を調べる（progress.md 続き 28）。"""
import sys, pickle; sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_b_2box_structure as B
r, s, nu = B.r, B.s, B.nu
q, t = r**2, s**2
u = sp.Symbol("u", positive=True)
d = pickle.load(open("data/refined_kauffman/daha/largeN_2box_m1.pkl", "rb"))
g = {eval(k): sp.cancel(sp.sympify(v, locals={"r": r, "s": s, "u": u}).subs(u, nu * s)) for k, v in d["g"].items()}
X = sp.cancel(g[(2,)] + g[(1, 1)])
x11 = sp.cancel(g[(1, 1)] + (q - t) / (1 - q * t) * g[(2,)])          # P_λ(Y; t, q) 基底の P_(1,1) の係数
Y11 = sp.cancel((1 - q * t) * x11)
zq, zt = r - 1 / r, s - 1 / s
phis = {"1": 1, "z_q（Mariño の M(q⁻¹)）": zq, "z_t（M(t)）": zt, "z_q z_t": zq * zt, "q^{1/2}+q^{-1/2}": r + 1 / r, "t^{1/2}+t^{-1/2}": s + 1 / s,
        "q^{1/2}−t^{1/2}": r - s, "(qt)^{1/2}−1": r * s - 1, "1−q": 1 - q, "1−t": 1 - t, "1−qt": 1 - q * t, "q−t": q - t}
def mono(e):
    n_, d_ = sp.fraction(sp.factor(e))
    return sp.Poly(n_, r, s, nu).is_monomial and sp.Poly(d_, r, s, nu).is_monomial
for name, F in (("X = g_(2)+g_(1,1)（p₁² 方向）", X), ("(1−qt)·x_(1,1)（P(Y;t,q) 基底）", Y11)):
    print("==", name, "：Laurent", sp.Poly(sp.fraction(F)[1], r, s, nu).is_monomial, flush=True)
    for lab, ph in phis.items():
        e = sp.cancel(F / ph)
        if not sp.Poly(sp.fraction(e)[1], r, s, nu).is_monomial:
            continue
        c1 = sum((v for k, v in B.blocks(e).items() if k[0] == 0), sp.Integer(0))
        symm = c1 != 0 and mono(sp.cancel(B.sigma(c1) / c1))
        sb = B.sign_blocks(e) if False else None
        # 𝐚 ブロックごとの符号（e 自身）
        A, Q, T = sp.symbols("A Q T")
        n_, d_ = sp.fraction(sp.cancel(e.subs({r: -T * Q, s: Q, nu: A}, simultaneous=True)))
        sh = sp.Poly(d_, A).degree()
        by = {}
        for (ea, _, _), cf in sp.Poly(sp.expand(n_), A, Q, T).terms():
            by.setdefault(ea - sh, []).append(cf)
        coh = sum(1 for v in by.values() if all(c > 0 for c in v) or all(c < 0 for c in v))
        div = B.divisible_universal(e)
        print(f"  ÷ {lab}：Laurent、c=1 が σ で（単項式倍を除いて）対称 {symm}、𝐚 ブロックの符号がそろう {coh}/{len(by)}、述べ直しの普遍因子 {div}", flush=True)
print("== 普遍因子 (q−t)(ν²−t) で割った商の 𝐚 ブロックごとの符号")
for name, e in (("ĝ_sum = X/z_q", X / zq), ("(1−qt)x_(1,1)", Y11), ("(1−qt)x_(1,1)/z_q", Y11 / zq), ("(1−qt)x_(1,1)/(q^{1/2}+q^{-1/2})", Y11 / (r + 1 / r))):
    sb = B.sign_blocks(sp.cancel(e))
    coh = sum(1 for n_, sg in sb.values() if sg != "混在")
    print(f"  {name}：そろう {coh}/{len(sb)}", {k: f"{n_}{'' if sg == '混在' else sg}" for k, (n_, sg) in sb.items()}, flush=True)
