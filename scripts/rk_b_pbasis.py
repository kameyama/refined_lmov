r"""タスク C (b)：自由エネルギーの右辺を Schur から Macdonald P_λ(Y; t, q)（双対 Cauchy の probe 側の基底）に替えたときの検査（progress.md 続き 29）。
  K：大きな N の trefoil の 2 箱。hatted（M(q^{-1})^{-1}）にしてから P 基底に替える版と、P 基底の係数を z_q で割る版（続き 28）。
  H：Hopf の多 probe の Ĉ_{R,S}（|R| ≤ 4, |S| ≤ 2、E を引いた C を hatted にしてから両成分を P 基底に替える）に述べ直しの判定式を当てる。
  F：固定 N = 2, 3 の trefoil の 3 箱（E = 0 のセクター）の g, ĝ を P 基底に替えたときの分母と整数性。"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp
from fractions import Fraction
from rcsknots.macdonald import macdonald_P
from rcsknots.characters import chi
from rcsknots.partitions import partitions, z as zlam

r, s = sp.symbols("r s", positive=True)
q, t = r**2, s**2
conv = lambda f: sp.Rational(f.numerator, f.denominator)


def schur_to_P(n):
    """C[R][λ]：s_R(Y) = Σ_λ C_{Rλ} P_λ(Y; t, q)。"""
    lams = partitions(n)
    P = macdonald_P(n, t, q, conv)
    Pm = sp.Matrix([[P[l].get(nu, 0) for nu in lams] for l in lams])          # 行 λ、列 ν（冪和）
    S = sp.Matrix([[sp.Rational(chi(R, nu), zlam(nu)) for nu in lams] for R in lams])
    C = (S * Pm.inv()).applyfunc(sp.cancel)
    return lams, C


def p_to_P(n):
    """D[ν][λ]：p_ν = Σ_λ D_{νλ} P_λ(Y; t, q)。"""
    lams = partitions(n)
    P = macdonald_P(n, t, q, conv)
    Pm = sp.Matrix([[P[l].get(nu, 0) for nu in lams] for l in lams])
    return lams, Pm.inv().applyfunc(sp.cancel)


def nonmono_den(e, gens):
    """分母の単項式でない部分。"""
    d = sp.fraction(sp.cancel(e))[1]
    c, fs = sp.factor_list(d)
    return sp.prod([f**k for f, k in fs if not sp.Poly(f, *gens).is_monomial])


def part_K():
    import rk_b_2box_structure as B
    nu = B.nu
    G = B.ghat()
    lams, C = schur_to_P(2)
    x = {l: sp.cancel(sum(C[i, j] * G[R] for i, R in enumerate(lams))) for j, l in enumerate(lams)}
    for l in lams:
        D = nonmono_den(x[l], (r, s, nu))
        print(f"K：hatted → P 基底の P_{l} の係数：単項式以外の分母 {sp.factor(D)}", flush=True)
    e = sp.cancel(x[(1, 1)] * nonmono_den(x[(1, 1)], (r, s, nu)))
    sb = B.sign_blocks(e)
    coh = sum(1 for n_, sg in sb.values() if sg != "混在")
    print(f"K：hatted → P 基底の P_(1,1) の係数（分母を払う）：普遍因子 {B.divisible_universal(e)}、商の 𝐚 ブロックの符号 {coh}/{len(sb)}", flush=True)


def part_H():
    import rk_b_hatted_analyze as H
    from rcsknots.bps import M_matrix
    from rk_b_hopf_ext import adams2
    rr, ss, nu, nut = H.r, H.s, H.nu, H.nut
    loc = {"r": rr, "s": ss, "nu": nu}
    gm = {k: sp.sympify(v, locals=loc) for k, v in pickle.load(open("data/refined_kauffman/hopf_g_two_multi_42.pkl", "rb")).items()}
    Eraw = {k: sp.sympify(v, locals=loc) for k, v in pickle.load(open("data/refined_kauffman/hopf_E_fixed.pkl", "rb")).items()}
    e11 = sp.sympify(pickle.load(open("data/refined_kauffman/hopf_eps_multi.pkl", "rb"))["eps_11_1"], locals=loc)
    eps = {((1,), (1,)): Eraw[(2, 2)], ((2,), (1,)): Eraw[(2, 4)], ((1, 1), (1,)): e11, ((1,), (2,)): Eraw[(2, 4)], ((1,), (1, 1)): e11}
    tt = ss**2
    sub = {r: rr, s: ss}
    rows = []
    for nR, nS in ((2, 2), (4, 2)):
        PR, PS = partitions(nR), partitions(nS)
        Cd = {}
        for R in PR:
            for S_ in PS:
                Ev = sum(adams2(rho).get(R, 0) * adams2(sg).get(S_, 0) * eps[(rho, sg)] for rho in partitions(nR // 2) for sg in partitions(nS // 2))
                Cd[(R, S_)] = sp.cancel(gm[(R, S_)] - Ev / (2 * (1 + tt)))
        MR, MS = M_matrix(nR, 1 / rr).inv(), M_matrix(nS, 1 / rr).inv()
        Ch = MR * sp.Matrix([[Cd[(R, S_)] for S_ in PS] for R in PR]) * MS.T
        _, CR = schur_to_P(nR)
        _, CS = schur_to_P(nS)
        CR, CS = CR.subs(sub), CS.subs(sub)
        Y = CR.T * Ch * CS                                     # Σ Ĉ_{RS} s_R s_S = Σ Y_{λμ} P_λ P_μ
        for i, l in enumerate(PR):
            for j, m in enumerate(PS):
                e = sp.cancel(Y[i, j])
                D = nonmono_den(e, (rr, ss, nu))
                e = sp.cancel(e * D)
                out = H.analyze(e, f"P 基底の Ĉ_{{{l},{m}}}（× {sp.factor(D)}）")
                rows.append((l, m, out))
    return rows


def p_to_P_rf(n):
    """flint の有理関数で D[ν][λ]：p_ν = Σ_λ D_{νλ} P_λ(Y; t, q)（Cramer の公式）。"""
    from rcsknots.flint_rf import R, S, conv as fc
    lams = partitions(n)
    P = macdonald_P(n, S * S, R * R, fc)
    Pm = [[P[l].get(nu, fc(0)) for nu in lams] for l in lams]
    def det(M):
        if len(M) == 1:
            return M[0][0]
        return sum(((-1)**j * M[0][j] * det([row[:j] + row[j + 1:] for row in M[1:]]) for j in range(len(M))), fc(0))
    dP = det(Pm)
    k = len(lams)
    inv = [[None] * k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            minor = [row[:i] + row[i + 1:] for a, row in enumerate(Pm) if a != j]
            inv[i][j] = (-1)**(i + j) * det(minor) / dP                    # (Pm^{-1})_{ij}
    return lams, inv


def part_F():
    import rk_k5_trefoil_3box as X
    lams, D = p_to_P_rf(3)
    for N in (2, 3):
        F = X.plog(X.so_Z(N, X.so_uniform(N)))
        Fc = X.plog(X.cover_Z(N, {}, uniform=True))
        Gd = {k: F.get(k, X.ZERO) - X.conv(Fraction(1, 2)) * Fc.get(k, X.ZERO) for k in set(F) | set(Fc)}
        for lab, vec in (("g", {rho: Gd.get(rho, X.ZERO) for rho in lams}), ("ĝ", {rho: Gd.get(rho, X.ZERO) / X.phi(rho) for rho in lams})):
            for j, l in enumerate(lams):
                y = sum((vec[nu_] * D[i][j] for i, nu_ in enumerate(lams)), X.ZERO)
                den = sp.factor(sp.fraction(sp.cancel(y.to_sympy()))[1])
                print(f"F：N = {N}、{lab} の P_{l}(Y;t,q) の係数：分母 {den}、整数 Laurent {y.is_integer_laurent()}", flush=True)
            if lab == "g":
                sch = X.to_schur(Gd, 3)
                print(f"F：N = {N}、比較：Schur 基底の g_R は整数 Laurent {all(v.is_integer_laurent() for v in sch.values())}", flush=True)


if __name__ == "__main__":
    part_K()
    part_H_K()
    part_F()


def part_H_K():
    """Hopf の 2 箱 × 2 箱：E を引いた C（hatted にしない）を両成分とも P 基底に替え、各成分を z_q で割る（K と同じ順序）。"""
    import rk_b_hatted_analyze as H
    from rk_b_hopf_ext import adams2
    rr, ss, nu = H.r, H.s, H.nu
    loc = {"r": rr, "s": ss, "nu": nu}
    gm = {k: sp.sympify(v, locals=loc) for k, v in pickle.load(open("data/refined_kauffman/hopf_g_two_multi_42.pkl", "rb")).items()}
    Eraw = {k: sp.sympify(v, locals=loc) for k, v in pickle.load(open("data/refined_kauffman/hopf_E_fixed.pkl", "rb")).items()}
    e11 = sp.sympify(pickle.load(open("data/refined_kauffman/hopf_eps_multi.pkl", "rb"))["eps_11_1"], locals=loc)
    eps = {((1,), (1,)): Eraw[(2, 2)], ((2,), (1,)): Eraw[(2, 4)], ((1, 1), (1,)): e11, ((1,), (2,)): Eraw[(2, 4)], ((1,), (1, 1)): e11}
    tt = ss**2
    P2 = partitions(2)
    Cd = sp.Matrix([[sp.cancel(gm[(R, S_)] - sum(adams2(a).get(R, 0) * adams2(b).get(S_, 0) * eps[(a, b)] for a in partitions(1) for b in partitions(1)) / (2 * (1 + tt)))
                     for S_ in P2] for R in P2])
    _, C2 = schur_to_P(2)
    C2 = C2.subs({r: rr, s: ss})
    Y = C2.T * Cd * C2
    zq = rr - 1 / rr
    for i, l in enumerate(P2):
        for j, m in enumerate(P2):
            e = sp.cancel(Y[i, j] / zq**2)
            D = nonmono_den(e, (rr, ss, nu))
            H.analyze(sp.cancel(e * D), f"K 順序の P 基底の C_{{{l},{m}}}/z_q²（× {sp.factor(D)}）")
