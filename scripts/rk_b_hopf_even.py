r"""最小版：両成分とも多 probe（|R|, |S| ≤ 2）の Hopf のクロスキャップ部分と、その hatted 不変量。
  SO 側 f_{R,S}（kauffman_hopf2.so_f2）、被覆 f^{cov}（cover_f2）、g_{R,S} = f_{R,S} − ½Σ N^R_{ρρ'}N^S_{σσ'} f^{cov}_{ρρ'σσ'}。
  偶数×偶数：E_{R,S} = c_{2;(1)}^R c_{2;(1)}^S ε（Adams、p_2 = s_2 − s_11）、ε = E_22（1 変数で一意、hopf_E_fixed.pkl）、C = g − ½E/(1+t)。
  hatted：ĝ = M^{−1}(q^{−1}) ⊗ M^{−1}(q^{−1}) g（Ĉ, Ê も同様）。分類は scripts/rk_b_classify.py。
出力：data/refined_kauffman/hopf_g_two_multi.pkl。"""
import sys, pickle, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp
from rcsknots.kauffman_hopf2 import so_f2, cover_f2
from rcsknots.kauffman_hopf import qs, ts, Us
from rcsknots.kauffman_cover import qs_, ss_, nus_
from rcsknots.characters import chi
from rcsknots.partitions import partitions, z as z_lam
from rcsknots.bps import M_matrix
from rk_b_classify import classify, r, s
from rk_k4_crosscap import nu

t = s**2
nut = sp.Symbol("nut", positive=True)


def lr(r1, r2):
    out = {}
    for R in partitions(sum(r1) + sum(r2)):
        cval = sum(sp.Rational(chi(r1, m1), z_lam(m1)) * sp.Rational(chi(r2, m2), z_lam(m2)) * chi(R, tuple(sorted(m1 + m2, reverse=True)))
                   for m1 in partitions(sum(r1)) for m2 in partitions(sum(r2)))
        if cval:
            out[R] = cval
    return out


def LR(a, b):
    return lr(a, b) if sum(a) + sum(b) else {(): 1}


t0 = time.time()
fSO = {k: v.subs({qs: r**2, ts: s**2, Us: nu * s}) for k, v in so_f2(2).items()}
fc = {k: v.subs({qs_: r**2, ss_: s, nus_: nu}) for k, v in cover_f2(2).items()}
print(f"計算 {time.time()-t0:.0f}s", flush=True)
half = {}
for (r1, r2, s1, s2), v in fc.items():
    for R, a in LR(r1, r2).items():
        for S_, b in LR(s1, s2).items():
            half[(R, S_)] = half.get((R, S_), 0) + a * b * v / 2
g = {k: sp.cancel(v - half.get(k, 0)) for k, v in fSO.items()}
pickle.dump({k: str(v) for k, v in g.items()}, open("data/refined_kauffman/hopf_g_two_multi.pkl", "wb"))
mono = lambda e: sp.Poly(sp.fraction(sp.cancel(e))[1], r, s, nu).is_monomial
print("対称性 g_{R,S} = g_{S,R}:", all(sp.cancel(g[(R, S_)] - g[(S_, R)]) == 0 for (R, S_) in g if (S_, R) in g))
for (R, S_), v in sorted(g.items(), key=lambda kv: (sum(kv[0][0]) + sum(kv[0][1]), kv[0])):
    even = sum(R) % 2 == 0 and sum(S_) % 2 == 0 and sum(R) and sum(S_)
    if v != 0:
        print(f"  g_{{{R},{S_}}}: Laurent = {mono(v)}" + ("（偶×偶）" if even else ""))
E = {k: sp.sympify(v, locals={"r": r, "s": s, "nu": nu}) for k, v in pickle.load(open("data/refined_kauffman/hopf_E_fixed.pkl", "rb")).items()}
eps = E[(2, 2)]
sgn = {(2,): 1, (1, 1): -1}
C = {}
for R in ((2,), (1, 1)):
    for S_ in ((2,), (1, 1)):
        C[(R, S_)] = sp.cancel(g[(R, S_)] - sgn[R] * sgn[S_] * eps / (2 * (1 + t)))
        print(f"  C_{{{R},{S_}}} = g − ½(±ε)/(1+t): Laurent = {mono(C[(R, S_)])}", flush=True)
# hatted
Minv = M_matrix(2, 1 / r).inv()
parts = partitions(2)


def hat2(X):
    Mx = sp.Matrix([[X[(R, S_)] for S_ in parts] for R in parts])
    H = Minv * Mx * Minv.T
    return {(R, S_): sp.cancel(H[i, j]) for i, R in enumerate(parts) for j, S_ in enumerate(parts)}


Ch = hat2(C)
Eh = hat2({(R, S_): sgn[R] * sgn[S_] * eps for R in parts for S_ in parts})
for k in Ch:
    w = sp.cancel(Ch[k].subs(nu, nut * s / r))
    print(f"  Ĉ_{k}: Laurent = {mono(Ch[k])}、ν~ 固定で t 非依存 = {not w.has(s)}、分類 = {classify(Ch[k]) if mono(Ch[k]) else '—'}", flush=True)
    print(f"      = {sp.factor(w)}", flush=True)
for k in Eh:
    print(f"  Ê_{k} = {sp.factor(sp.cancel(Eh[k].subs(nu, nut * s / r)))}", flush=True)
# 奇数のセクター（|R| = 2, |S| = 1 と |R| = 1, |S| = 2）の hatted
for S_ in ((1,),):
    vec = Minv * sp.Matrix([g[(R, S_)] for R in parts])
    for R, v in zip(parts, vec):
        print(f"  ĝ_{{{R},{S_}}} = {sp.factor(sp.cancel(sp.cancel(v).subs(nu, nut * s / r)))}（既存の ĝ_{{R,[1]}} と比較）")

print("== 対角の Ĉ の refined だけの部分（ν~ を固定し、t = q との差）")
for k in (((2,), (2,)), ((1, 1), (1, 1))):
    w = sp.cancel(Ch[k].subs(nu, nut * s / r))
    w0 = sp.cancel(w.subs(s, r))
    diff = sp.factor(sp.cancel(w - w0))
    print(f"  Ĉ_{k}：t=q の値 = {sp.factor(w0)}")
    print(f"      差（refined だけ）= {diff}")
    quo = sp.cancel(diff / ((r**2 - s**2) / r**2))
    print(f"      差 / (1 − t/q) = {sp.factor(quo)}  Laurent = {sp.Poly(sp.fraction(quo)[1], r, s, nut).is_monomial}")
