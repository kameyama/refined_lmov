r"""タスク C の (b)：大きな N の 2 箱（trefoil、一様な規格化）の t = q での Mariño 0904.1088 との辞書。
  - SO の W_λ と被覆の W_μ（普通の表現）は、t = q で Mariño の G_λ, H_μ（標準 framing）の鏡像（ν → 1/ν, t → 1/t）に比 1 で一致する。
  - クロスキャップ：g^{odd}_R = g_R + (−1)^{ℓ(R)−1}·½ g_(1)(x²)（全 d から奇数 d へ）から unknot の値（g^{odd}(unknot) = ∓½）を引き、
    M(q^{−1})^{−1} で hatted にすると、Mariño の ĝ_R の鏡像の −1 倍（R はそのまま）。g_(1) = −鏡像(g^{Mar}_[1])。
  Mariño の変数：t_M = q^{1/2} = r、ν_M = ν = U/t^{1/2}（t = q で U = ν r）。入力は scripts/rk_unrefined_check.py（tests/test_rk_unrefined.py で (4.8) と照合済み）。"""
import sys
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp

r_, s_, u_ = sp.symbols("r s u", positive=True)


def marino_trefoil():
    import rk_unrefined_check as rk
    from rcsknots import homfly
    nu, t, z = rk.nu, rk.t, rk.z
    sub = {homfly.a: nu**2, homfly.Q: t}
    H = {lam: rk.simp(homfly.torus_reduced(lam, 2, 3, -1).subs(sub) * rk.H_unk[lam]) for lam in ((1,), (2,), (1, 1))}
    Hadj = rk.simp(rk.H_adj_unk * (4*nu**4 - 4*nu**6 + nu**8 + z**2*(4*nu**4 - 7*nu**6 + 2*nu**8 + nu**10) + z**4*(nu**4 - 2*nu**6 + nu**8)))
    F31 = 2*nu**2 - nu**4 + z*(-nu**3 + nu**5) + z**2*(nu**2 - nu**4)
    G = {(1,): rk.simp(F31 * rk.G_unk[(1,)]),
         (2,): rk.simp(rk.nrz_trefoil(2, nu, t, sp.Integer(-1)) * rk.G_unk[(2,)]),
         (1, 1): rk.simp(rk.nrz_trefoil(2, -nu, 1/t, sp.Integer(-1)) * rk.G_unk[(1, 1)])}
    g1, g2, g11, _, _ = rk.reformulated(H, Hadj, G)
    gh2, gh11 = rk.hat2(g2, g11)
    return rk, H, G, {(1,): g1, (2,): g2, (1, 1): g11}, {(2,): gh2, (1, 1): gh11}


def compare():
    import rk_c_2box_largeN as L
    rk, H, G, gM, ghM = marino_trefoil()
    nu, t = rk.nu, rk.t
    tq = lambda e: sp.cancel(e.subs(s_, r_).subs(u_, nu * r_).subs(r_, t))
    mir = lambda e: sp.cancel(e.subs({nu: 1/nu, t: 1/t}, simultaneous=True))
    out = {}
    W, co = L.so_parts()
    out["SO W_λ = 鏡像 G_λ"] = all(sp.cancel(tq(W[k].to_sympy()) - mir(G[k])) == 0 for k in W)
    w = L.S / L.R
    out["被覆 W_μ = 鏡像 H_μ"] = all(sp.cancel(tq((w**sum(k) * L.Rdata_nu("".join(map(str, k))) * L.Psym(k)).to_sympy()) - mir(H[k])) == 0
                                for k in ((1,), (1, 1), (2,)))
    f, fc = L.so_f(W, co), L.cover_f()
    g = {k: tq((f[k] - L.half * fc[k]).to_sympy()) for k in f}
    out["g_(1) = −鏡像 g^Mar_[1]"] = sp.cancel(g[(1,)] + mir(gM[(1,)])) == 0
    half = sp.Rational(1, 2)
    g1x2 = g[(1,)].subs({t: t**2, nu: nu**2}, simultaneous=True)
    god = {(2,): g[(2,)] + half * g1x2 + half, (1, 1): g[(1, 1)] - half * g1x2 - half}
    det = t**2 - t**-2                                          # M(q^{-1}) = [[t, −1/t], [−1/t, t]]（bps.M_matrix(2, 1/t)）
    gh = {(2,): (t * god[(2,)] + god[(1, 1)] / t) / det, (1, 1): (god[(2,)] / t + t * god[(1, 1)]) / det}
    out["ĝ^odd − unknot = −鏡像 ĝ^Mar"] = all(sp.cancel(gh[k] + mir(ghM[k])) == 0 for k in gh)
    return out


if __name__ == "__main__":
    for k, v in compare().items():
        print(k, v, flush=True)
