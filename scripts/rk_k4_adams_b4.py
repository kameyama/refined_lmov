r"""K4'(a)：成分 1 に probe 無限個（|R| ≤ 4）、成分 2 は 1 変数で b ≤ 4。Adams の予言を b = 4 で検査し、
ε_{(1,1),(2)}（E_[1^4],4 の固定セクター ((1,1),(1,1)) × (2,2)）を決め、ε_{(ρ),(2)}/ε_{(ρ),(1)} の比から E_44 = ε_{(2),(2)} を予言して検査する。
計算（約 15 分）は data/refined_kauffman/{hopf_fcov_multi_A4B4, hopf_f_multi_A4B4}.pkl に保存し、あれば読み込む。"""
import sys, os, time, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp
from rcsknots.characters import chi
from rcsknots.partitions import partitions, z as z_lam

r, s, nu = sp.symbols("r s nu", positive=True)
nut = sp.Symbol("nut")
t = s**2
loc = {"r": r, "s": s, "nu": nu}
A_, B = 4, 4
FC, FS = "data/refined_kauffman/hopf_fcov_multi_A4B4.pkl", "data/refined_kauffman/hopf_f_multi_A4B4.pkl"
if not os.path.exists(FC):
    from rcsknots.kauffman_cover import cover_f_multi, qs_, ss_, nus_
    t0 = time.time()
    fc = {k: v.subs({qs_: r**2, ss_: s, nus_: nu}) for k, v in cover_f_multi(A_, B).items()}
    pickle.dump({k: str(v) for k, v in fc.items()}, open(FC, "wb"))
    print(f"cover_f_multi({A_},{B}): {time.time()-t0:.0f}s", flush=True)
fc_raw = pickle.load(open(FC, "rb"))
if not os.path.exists(FS):
    from rcsknots.kauffman_hopf import hopf_Ztilde_multi, log_coeffs_multi, plethystic_f, qs, ts, Us
    t0 = time.time()
    f = plethystic_f(log_coeffs_multi(hopf_Ztilde_multi(A_, B), A_, B), A_, B, "schur")
    pickle.dump({k: str(v.as_expr().subs({qs: r**2, ts: s**2, Us: nu * s})) for k, v in f.items()}, open(FS, "wb"))
    print(f"SO multi ({A_},{B}): {time.time()-t0:.0f}s", flush=True)
fSO_raw = pickle.load(open(FS, "rb"))
# 後処理は有理関数体 Q(r,s,ν) の上で、必要な鍵だけ行う（sympy の cancel を式に使うと b = 4 で終わらない：2026-10-18 に 10 時間停止）
from sympy import QQ
from sympy.polys.fields import field
Kr, R_, S_, N_ = field("r,s,nu", QQ)
plain = {"r": sp.Symbol("r"), "s": sp.Symbol("s"), "nu": sp.Symbol("nu")}
to_pos = {plain["r"]: r, plain["s"]: s, plain["nu"]: nu}
conv = lambda st: Kr.from_expr(sp.sympify(st, locals=plain))
expr = lambda v: v.as_expr().subs(to_pos)
T = S_**2


def lr(r1, r2):
    out = {}
    for R in partitions(sum(r1) + sum(r2)):
        c = sum(sp.Rational(chi(r1, m1), z_lam(m1)) * sp.Rational(chi(r2, m2), z_lam(m2)) * chi(R, tuple(sorted(m1 + m2, reverse=True)))
                for m1 in partitions(sum(r1)) for m2 in partitions(sum(r2)))
        if c:
            out[R] = c
    return out


def adams2(rho):
    """p_2∘s_ρ = Σ_R c_R s_R。"""
    n = sum(rho)
    return {R: c for R in partitions(2 * n)
            if (c := sum(sp.Rational(chi(rho, m), z_lam(m)) * chi(R, tuple(2 * x for x in m)) for m in partitions(n)))}


need = [(R, 4) for n in (2, 4) for R in partitions(n)] + [(R, 2) for n in (1, 2) for R in partitions(n)] + [((1,), 1)]
t0 = time.time()
g = {}
for (R, b) in need:
    h = Kr.zero
    for (r1, r2, k, l), st in fc_raw.items():
        if k + l != b or sum(r1) + sum(r2) != sum(R):
            continue
        c = lr(r1, r2).get(R, 0) if sum(r1) + sum(r2) else 1
        if c:
            h += Kr(c) * conv(st) / 2
    g[(R, b)] = conv(fSO_raw[(R, b)]) - h
print(f"g（必要な鍵 {len(need)} 個）: {time.time()-t0:.0f}s", flush=True)


def ad(v):
    return Kr.from_expr(v.as_expr().subs({plain["r"]: plain["r"]**2, plain["s"]: plain["s"]**2, plain["nu"]: plain["nu"]**2}, simultaneous=True))


godd = {}
for key in sorted(need, key=lambda k: (sum(k[0]), k[1])):
    R, b = key
    v = g[key]
    if sum(R) % 2 == 0 and b % 2 == 0:
        for rho in partitions(sum(R) // 2):
            c = adams2(rho).get(R, 0)
            if c and (rho, b // 2) in godd:
                v = v + Kr(sp.Rational(c, 2)) * ad(godd[(rho, b // 2)])
    godd[key] = v
mono = lambda v: len(v.denom.terms()) == 1
G = lambda R, b: g[(tuple(R), b)]
print("Adams の予言（b = 4）:")
for name, combo in (("g_[4]+g_[3,1]", G((4,), 4) + G((3, 1), 4)), ("g_[2,1,1]+g_[1^4]", G((2, 1, 1), 4) + G((1, 1, 1, 1), 4)),
                    ("g_[2,2]−g_[4]−g_[1^4]", G((2, 2), 4) - G((4,), 4) - G((1, 1, 1, 1), 4)), ("g_[2]+g_[1,1]", G((2,), 4) + G((1, 1), 4))):
    print(f"  {name}: Laurent = {mono(combo)}", flush=True)
for key in ((((4,), 4)), ((1, 1, 1, 1), 4), ((2,), 4)):
    print(f"  {key}: 単独では Laurent = {mono(g[key])}（偶×偶なので False が期待値）", flush=True)
from rk_k4_fixed_sector import solve_generic
Efix = {k: sp.sympify(v, locals=loc) for k, v in pickle.load(open("data/refined_kauffman/hopf_E_fixed.pkl", "rb")).items()}
print("  E_[4],4 = E_44（1 変数）との整合：g_[4],4 − ½E_44/(1+t) は 1 変数の E_44 が一意でないので省略", flush=True)
eps = {}
for rho, R in (((2,), (4,)), ((1, 1), (1, 1, 1, 1))):
    b = 4
    sec = (rho, rho, 2, 2)
    sign = adams2(rho)[R]
    Pexpr = expr((1 - T) * conv(fc_raw[sec]))
    Rexpr = expr(2 * (1 + T) * g[(R, b)] * Kr(sign))
    Eq = sp.cancel(expr(2 * (1 + T) * (g[(R, b)] - godd[(R, b)]) * Kr(sign)).subs(s, r))
    t0 = time.time()
    res = solve_generic(Pexpr, Rexpr, Eq)
    tag = f"ε_{rho},(2)"
    print(f"{tag}: 解あり {res['feasible']}、状態 {res['n']}、一意でない変数 {len(res.get('free', []))}（{time.time()-t0:.0f}s）", flush=True)
    if res["feasible"]:
        eps[(rho, b)] = res["E"]
        print("   =", sp.factor(sp.cancel(res["E"].subs(nu, nut * s / r))) if len(str(res["E"])) < 3000 else "(長い)", flush=True)
# ε の解は b = 4 で一意でないので保存しない
e11_1 = sp.sympify(pickle.load(open("data/refined_kauffman/hopf_eps_multi.pkl", "rb"))["eps_11_1"], locals=loc)
if ((1, 1), 4) in eps:
    Y = sp.cancel(eps[((1, 1), 4)] / e11_1)
    print("ε_(1,1),(2)/ε_(1,1),(1) =", sp.factor(sp.cancel(Y.subs(nu, nut * s / r))), flush=True)
    E44g = sp.cancel(Efix[(2, 4)] * Y)       # 予言：ε_(2),(2) = ε_(2),(1) × 同じ比
    Rg = expr(2 * (1 + T) * g[((4,), 4)])
    Eq4 = sp.cancel(expr(2 * (1 + T) * (g[((4,), 4)] - godd[((4,), 4)])).subs(s, r))
    print("予言 E_44 = ε_(2),(1)·[ε_(1,1),(2)/ε_(1,1),(1)]：Laurent =", sp.Poly(sp.fraction(E44g)[1], r, s, nu).is_monomial,
          " t=−1 =", sp.cancel(sp.cancel(Rg - E44g).subs(s, sp.I)) == 0, " t=q =", sp.cancel(E44g.subs(s, r) - Eq4) == 0, flush=True)
print("完了", flush=True)
