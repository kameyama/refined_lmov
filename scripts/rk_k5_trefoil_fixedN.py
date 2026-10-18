r"""(b)：N = 2, 3 の refined SO(2N) trefoil（1^k、k ≤ N）の自由エネルギー f_a（a ≤ 2N）を厳密に求め、整数係数の Laurent 多項式かを見る。
W_{1^k} = m_k R K（scripts/rk_k5_framing_scan.py の setup）。framing の単項式は、Z~ の極が打ち消す条件（極の上でしか決まらない）に、t = q で 1 の条件を課して決める（conventions.md「固定 N の framing」）：
  N = 2：m_2 = t/q（q t = 1 の上で s^4）。N = 3：m_2 = (t/q)^2（q t³ = 1 の上で s^16）、m_3 = m_1、m_3' = (t/q)^2 m_1（q t = 1 の上で s^8 m_1）。
  m_1 は固定 N では自由：K3' の値 (t/q)^{1/2}（仮説 T_[1]^{−6}(q/t)^{1/2} を data の鏡像に対して書いたもの；N = 3 で確認）と
  hatted の値 t/q（rk_k5_trefoil_2box_cover.py）の両方で計算する。保存するのは K3' の値のもの。
  2026-10-18 の初版は t = q の条件を課さず、m_2 = s^4（N=2）・s^16（N=3）、m_3' = s^8 m_1 を使った。"""
import sys, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp
from rcsknots.kauffman_knot import free_energy, unknot_W, r, s
from rcsknots.kauffman_trefoil import Z615, unknot_SO, r as rT, s as sT, U as UT
import rk_k5_framing_scan as F

w = s / r


def framings(N, m1):
    return {2: [m1, w**2], 3: [m1, w**4, m1, w**4 * m1]}[N]


out = {}
for N, m1 in ((2, w), (3, w), (2, w**2), (3, w**2)):
    parts = F.setup(N)
    Ms, i, W = framings(N, m1), 0, [sp.Integer(1)]
    for k in range(1, N + 1):
        e = 0
        for part in parts[k]:
            e += Ms[i] * part; i += 1
        W.append(sp.cancel(e))
    if N == 3:
        k3 = sp.cancel((Z615 * unknot_SO * rT**-5 * sT**11 * UT**-6).subs({rT: r, sT: s, UT: s**6}))
        print("N=3：W_1 == K3' の仮説の W_SO（(6.15)·T^{−6}(q/t)^{1/2}·unknot）:", sp.cancel(W[1] - k3) == 0)
    t0 = time.time()
    f = free_energy(N, W, 2 * N)
    print(f"N = {N}、m_1 = {m1}（{time.time()-t0:.0f}s）:")
    for a, v in f.items():
        num, den = sp.fraction(sp.cancel(v))
        P = sp.Poly(sp.expand(num), r, s)
        D = sp.Poly(den, r, s)
        print(f"  f_{a}: Laurent = {D.is_monomial}、整数係数 = {D.is_monomial and all((c / D.LC()).is_integer for c in P.coeffs())}、項数 {len(P.terms())}")
    if m1 == w:
        out[N] = {a: str(v) for a, v in f.items()}
import pickle
pickle.dump(out, open("data/refined_kauffman/trefoil_fixedN_f.pkl", "wb"))
