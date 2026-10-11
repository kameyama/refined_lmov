"""（試作）3 箱の公式を、ノートブック new derivation of link rf_qt ver2-4.nb と同じく unreduced の C の短い式として取り出す。

\\hat f_ρ = Σ_i k_{ρ,i}(Q, T)·（C の単項式）。k は f^q の各単項式の係数を factor し、M(t)⁻¹ を掛けて factor したもの（全体の通分はしない）。
使い方：PYTHONPATH=src python scripts/dgr_formula_short.py
"""
import sys, time; sys.path.insert(0, "scripts")
import sympy as sp
t0 = time.time()
from rcsknots.bps import C, C1, M_matrix, Q, T, f_formulas
from rcsknots.partitions import partitions
F = f_formulas(3); print("f_formulas(3)", f"{time.time()-t0:.1f}s", flush=True)
parts3 = partitions(3)
Minv = M_matrix(3).inv().applyfunc(sp.factor); print("M⁻¹", f"{time.time()-t0:.1f}s", flush=True)
syms = [C((3,)), C((2, 1)), C((1, 1, 1)), C((2,)), C((1, 1)), C((1,)), C1(2), C1(3)]
mons = [C((3,)), C((2, 1)), C((1, 1, 1)), C((2,)) * C((1,)), C((1, 1)) * C((1,)), C((1,)) ** 3, C((1,)) * C1(2), C1(3)]
fq = [sp.expand(F[mu][0]) for mu in parts3]
# 各 f^q_μ の単項式の係数（Q, T の有理式）
coef = []
for e in fq:
    P = sp.Poly(e, *syms)
    d = {}
    for mon, c in P.terms():
        d[mon] = sp.factor(c)
    coef.append(d)
print("f^q の係数", f"{time.time()-t0:.1f}s", {i: len(d) for i, d in enumerate(coef)}, flush=True)
for i, rho in enumerate(parts3):
    tot = {}
    for j in range(3):
        for mon, c in coef[j].items():
            tot[mon] = tot.get(mon, 0) + Minv[i, j] * c
    tot = {mon: sp.factor(sp.together(c)) for mon, c in tot.items()}
    tot = {mon: c for mon, c in tot.items() if c != 0}
    print(rho, f"{time.time()-t0:.1f}s")
    for mon, c in tot.items():
        print("   ", mon, ":", c)
