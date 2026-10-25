r"""K3'：trefoil [1]。AS (6.17) と (6.15) の一致、t = q での Mariño との関係、refined の g^{ref}_[1] の形。"""
import sys
sys.path.insert(0, "src")
import sympy as sp
from rcsknots.refined_kauffman.kauffman_trefoil import Z_V_T2, Z615, r, s, U, bq, bt, ba, marino_F31, g_ref, W_SO, h_half_ref, is_laurent

print("(6.17) m=1 − (6.15) =", sp.simplify(Z_V_T2(1) - Z615))
print("t=q: AS/F_M(ν→1/ν, z→−z) =", sp.factor(sp.cancel(Z615.subs(r, s) / marino_F31(s / U, 1 / s - s))))
for mono, name in ((r**-6 * s**12 * U**-6, "T_[1]^{-6}"), (r**-5 * s**11 * U**-6, "T_[1]^{-6} (q/t)^{1/2}")):
    g = g_ref(mono)
    print(f"単項式 {name}: g^ref Laurent = {is_laurent(g)}")
g = g_ref()
Q, T, A = sp.symbols("Q T A")   # 太字 𝐪, 𝐭, 𝐚（𝐚 = U/s）
gb = sp.expand(sp.cancel(g.subs({r: -T * Q, U: A * Q}).subs(s, Q)))
print("g^ref_[1](3_1) (r,s,U) =", sp.factor(g))
print("g^ref_[1](3_1) 太字 =", gb)
print("t = q (𝐭 = −1):", sp.expand(gb.subs(T, -1)))
print("W_SO 太字 =", sp.factor(sp.cancel(W_SO().subs({r: -T * Q, U: A * Q}).subs(s, Q))))
print("½h^ref 太字 =", sp.factor(sp.cancel(h_half_ref().subs({r: -T * Q, U: A * Q}).subs(s, Q))))
