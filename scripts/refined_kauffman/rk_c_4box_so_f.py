r"""大きな N の trefoil の SO の自由エネルギー f を Y の冪和の次数 4 まで組み立てる（theory §12.7 の条件 (i)、progress.md 続き 72）。
  rk_c_3box_largeN の SO 側を次数 4 に広げたもの：Z~(Y) = Σ_ν (−1)^{|ν|} s_ν(Y) ⟨s_{ν'}(x, x^{−1})⟩、⟨S_ρ⟩ = Σ_λ Sinv[ρ][λ] W_λ、
  W_λ = R_λ(U) K_λ(t^ρ)（R_λ は R_U_m1_{λ}.pkl、4 箱は rk_c_4box_largeN で決めたもの）、K_λ の S_μ 展開は data/refined_kauffman/koornwinder_DN。
  f は全 d の Möbius 反転 f = Σ_d μ(d)/d ψ_d(log Z~)（ψ_d：r, s, u → d 乗、p_k → p_{dk}）。
  検算：次数 ≤ 3 の成分が largeN_3box_m1.pkl の f と全係数一致すること。
保存：data/refined_kauffman/daha/largeN_4box_so_m1.pkl（{"f": {冪和の分割: 式の文字列}}、変数 (r, s, u)）。"""
import sys, os, pickle, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from fractions import Fraction
from rcsknots.refined_kauffman.flint_rf3 import from_sympy, conv
from rcsknots.characters import chi
from rcsknots.partitions import partitions, z as zlam
import rk_c_3box_largeN as T3

OUT = "data/refined_kauffman/daha"
r_, s_, u_ = T3.r_, T3.s_, T3.u_
ZERO, ONE = T3.ZERO, T3.ONE
KMAX = 4
PARTS = [()] + [lam for n in range(1, KMAX + 1) for lam in partitions(n)]


def K_matrix():
    """S_μ = Σ_λ Sinv[μ][λ] K_λ（|λ|, |μ| ≤ 4）。K_λ の S_μ 展開は偶奇で分かれるので、次数の偶奇ごとに逆にする。"""
    qs, ts, Us = sp.symbols("q t U")
    C = {(): {(): ONE}}
    for lam in PARTS[1:]:
        raw = pickle.load(open("data/refined_kauffman/koornwinder_DN/K_" + "_".join(map(str, lam)) + ".pkl", "rb"))
        C[lam] = {mu: from_sympy(sp.sympify(v, locals={"q": qs, "t": ts, "U": Us}).subs({qs: r_**2, ts: s_**2, Us: u_})) for mu, v in raw.items()}
    Sinv = {}
    for par in (0, 1):
        ps = [lam for lam in PARTS if sum(lam) % 2 == par]
        assert all(sum(mu) % 2 == par for lam in ps for mu in C[lam])
        n = len(ps)
        M = [[C[lam].get(mu, ZERO) for mu in ps] + [ONE if i == j else ZERO for j in range(n)] for i, lam in enumerate(ps)]
        for col in range(n):
            piv = next(i for i in range(col, n) if not M[i][col].is_zero())
            M[col], M[piv] = M[piv], M[col]
            inv = 1 / M[col][col]
            M[col] = [x * inv for x in M[col]]
            for i in range(n):
                if i != col and not M[i][col].is_zero():
                    f = M[i][col]
                    M[i] = [x - f * y for x, y in zip(M[i], M[col])]
        for i, mu in enumerate(ps):
            Sinv[mu] = {lam: M[i][n + j] for j, lam in enumerate(ps) if not M[i][n + j].is_zero()}
    return Sinv


def so_W():
    from rcsknots.refined_kauffman.kauffman_hopf2 import K_rho
    from rcsknots.refined_kauffman.kauffman_hopf import qs, ts, Us
    sub = {qs: r_**2, ts: s_**2, Us: u_}
    W = {(): ONE}
    for lam in PARTS[1:]:
        c = pickle.load(open(os.path.join(OUT, f"R_U_m1_{''.join(map(str, lam))}.pkl"), "rb"))
        Rf = sum((from_sympy(sp.sympify(v, locals={"r": r_, "s": s_})) * T3.U**j for j, v in c.items()), ZERO)
        kr = K_rho(lam)
        W[lam] = Rf * from_sympy(sp.sympify(kr.as_expr() if hasattr(kr, "as_expr") else kr).subs(sub))
    return W


def so_Z(W, Sinv):
    Sexp = {rho: sum((c * W[lam] for lam, c in Sinv[rho].items()), ZERO) for rho in PARTS}
    Zp = {(): ONE}
    for n in range(1, KMAX + 1):
        for lam in partitions(n):
            lamT = tuple(sum(1 for x in lam if x > j) for j in range(lam[0]))
            val = sum((conv(Fraction(chi(lamT, rho), zlam(rho))) * Sexp[rho] for rho in partitions(n)), ZERO)
            Zp = T3.padd(Zp, T3.schur_p(lam), conv((-1)**n) * val)
    return Zp


def plog(Zp):
    """f = Σ_d μ(d)/d ψ_d(log Z~)（次数 ≤ 4）。"""
    X = {k: v for k, v in Zp.items() if k}
    Lg, Xn = {}, dict(X)
    for n in range(1, KMAX + 1):
        Lg = T3.padd(Lg, Xn, conv(Fraction((-1)**(n + 1), n)))
        Xn = T3.pmul(Xn, X, KMAX)
    mob = {1: 1, 2: -1, 3: -1, 4: 0}
    f = {}
    for d, mu in mob.items():
        if mu == 0:
            continue
        for k, v in Lg.items():
            if d * sum(k) <= KMAX:
                kk = tuple(sorted((d * x for x in k), reverse=True))
                f[kk] = f.get(kk, ZERO) + conv(Fraction(mu, d)) * v.adams(d)
    return f


def main():
    t0 = time.time()
    Sinv = K_matrix()
    print(f"K の逆 {time.time()-t0:.0f}s", flush=True)
    W = so_W()
    print(f"W_λ {time.time()-t0:.0f}s", flush=True)
    F = plog(so_Z(W, Sinv))
    print(f"f {time.time()-t0:.0f}s", flush=True)
    old = pickle.load(open(os.path.join(OUT, "largeN_3box_m1.pkl"), "rb"))["f"]
    ok = {k: (F[eval(k)] - from_sympy(sp.sympify(v, locals={"r": r_, "s": s_, "u": u_}))).is_zero() for k, v in old.items()}
    print(f"次数 ≤ 3 の f が largeN_3box_m1 と一致：{ok}", flush=True)
    for k in partitions(4):
        den = sp.factor(sp.fraction(sp.cancel(F[k].to_sympy()))[1])
        print(f"f_{k}：分母 {den}", flush=True)
    pickle.dump({"f": {str(k): str(v.to_sympy()) for k, v in F.items()}}, open(os.path.join(OUT, "largeN_4box_so_m1.pkl"), "wb"))
    print(f"保存 {time.time()-t0:.0f}s", flush=True)
    return F


if __name__ == "__main__":
    main()
