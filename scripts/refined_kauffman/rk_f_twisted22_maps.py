r"""rk_f_twisted22 の refined の ψ² 型の候補を、t = q で ψ² に戻る写像 r → r^a s^b、s → r^c s^d（a + b = c + d = 2）と (s/r)^j のずれで走査する（theory §12.7）。"""
import sys, itertools
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_f_twisted22 as T
import rk_b_twisted_multiplet as TM
import rk_f_frozen_basis as FB
r, s, nu = T.r, T.s, T.nu
exps = [(2, 0), (1, 1), (0, 2)]
for n in (1, 2):
    P, F = T.sectors()[n]
    Bu, Bt = T.zq**(2 * n - 1), T.zq**(n - 1) * T.zqp**n
    Pk = sp.cancel(P / (s * Bu))
    blocks = sorted(FB.nu_blocks(Pk))
    res = {}
    for (a, b), (c, d) in itertools.product(exps, exps):
        Fm = F.subs({r: r**a * s**b, s: r**c * s**d, nu: nu**2}, simultaneous=True)
        Et = sp.cancel(-(1 + r**2) * Fm)
        good = []
        for k in blocks:
            mult = TM.monomials(T.block(Pk, k))
            cand = sp.cancel(T.block(Et, k) / (s * Bt))
            if not sp.Poly(sp.fraction(cand)[1], r, s, nu).is_monomial:
                continue
            for j in range(-30, 31):
                ce = TM.monomials(sp.expand(cand * (s / r)**j))
                if all(m in mult and abs(ce[m]) <= abs(mult[m]) and (ce[m] - mult[m]) % 2 == 0 for m in ce) and \
                   all(mult[m] % 2 == 0 for m in mult if m not in ce):
                    good.append((k, j)); break
        res[((a, b), (c, d))] = good
        print(f"n={n} 写像 r→r^{a}s^{b}, s→r^{c}s^{d}：条件を満たすブロック {len(good)}/{len(blocks)} {good}", flush=True)
