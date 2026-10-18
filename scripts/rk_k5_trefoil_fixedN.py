r"""(b)：N = 2, 3 の refined SO(2N) trefoil（1^k、k ≤ N）の自由エネルギー f_a（a ≤ 2N）を厳密に求め、整数係数の Laurent 多項式かを見る。
W_{1^k} = m_k R K（scripts/rk_k5_framing_scan.py の setup）。framing の単項式は、Z~ の極が打ち消す条件（極の上でしか決まらない）を満たす最も簡単なもの：
  m_1 = r^{−1} s（K3' の仮説 T_[1]^{−6}(q/t)^{1/2} を data の鏡像に対して書いたもの；N = 3 で確認）、
  N = 2：m_2 = s^4（q t = 1 の上で s^4）。N = 3：m_2 = s^16（q t³ = 1 の上で s^16）、m_3 = m_1、m_3' = s^8 m_1（q t = 1 の上で）。"""
import sys, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp
from rcsknots.kauffman_knot import free_energy, unknot_W, r, s
from rcsknots.kauffman_trefoil import Z615, unknot_SO, r as rT, s as sT, U as UT
import rk_k5_framing_scan as F

M = {2: [r**-1 * s, s**4], 3: [r**-1 * s, s**16, r**-1 * s, s**8 * r**-1 * s]}
out = {}
for N in (2, 3):
    parts = F.setup(N)
    Ms, i, W = M[N], 0, [sp.Integer(1)]
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
    print(f"N = {N}（{time.time()-t0:.0f}s）:")
    for a, v in f.items():
        num, den = sp.fraction(sp.cancel(v))
        P = sp.Poly(sp.expand(num), r, s)
        print(f"  f_{a}: Laurent = {sp.Poly(den, r, s).is_monomial}、整数係数 = {all(c.is_integer for c in P.coeffs())}、項数 {len(P.terms())}")
    out[N] = {a: str(v) for a, v in f.items()}
import pickle
pickle.dump(out, open("data/refined_kauffman/trefoil_fixedN_f.pkl", "wb"))
