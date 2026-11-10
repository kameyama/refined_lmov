r"""大きな N の trefoil の 4 箱のクロスキャップ部分 g = f − ½f^{cov}（Y の冪和の次数 4 まで、progress.md 続き 139）。変数 (r, s, u)、flint_rf3。
SO：rk_c_4box_so_f（D_N の K_λ、|λ| ≤ 4、R_U_m1_{λ}）。
被覆 U(N_A)：Z~ = Σ_{λ,μ} (−1)^{|λ|+|μ|} ⟨P_λ(x) P_μ(x^{−1})⟩ P_{λ'}(Y; t, q) P_{μ'}(Y; t, q)（|λ| + |μ| ≤ 4、rk_c_3box_largeN と同じ）。
  ⟨P_λ P_μ^∨⟩：(|λ|, |μ|) ≤ (2, 2) は rk_f_cover_fixed22.pair_values（普通の表現、随伴 + κ、(α;1̄) |α| = 2、(α;β̄) |α| = |β| = 2）。
  普通の表現 |λ| = 3, 4：W_λ = (t/q)^{|λ|/2} R_λ(KN の data) P_λ(t^ρ)。
  (α;1̄)、|α| = 3：W_{(α;1)} + Σ_μ c_{α,μ} W_μ、W_{(α;1)} = R_{(α;1)}(A^{1/2}) P_{(α;1)}(t^ρ)
    （R は rk_c_fast_composite3 の値を rk_c_bounded_save.run_fixed で補間した Rcomp_A_m1_{α}.pkl）、
    c_{α,μ} は縦帯の Pieri の ψ'（rk_c_fast_composite3.psi_A）。荷電共役で ⟨P_1 P_α^∨⟩ = ⟨P_α P_1^∨⟩。
自由エネルギーは全 d の Möbius 反転（rk_c_4box_so_f.plog）。
検算：次数 ≤ 3 の g が largeN_3box_m1.pkl と全係数一致。被覆の f の固定 N との照合は rk_c_4box_cover_fixedN.py。
保存：data/refined_kauffman/daha/largeN_4box_m1.pkl（{"f"|"fcov"|"g": {冪和の分割: (分子の辞書, 分母の辞書)}}、flint_rf3 の (r, s, u) の指数 → 整数係数）。"""
import sys, os, pickle, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from fractions import Fraction
from rcsknots.refined_kauffman.flint_rf3 import RF, R, S, U, CTX, from_sympy, conv
from rcsknots.refined_kauffman.kauffman_refined import dim_composite
from rcsknots.macdonald import macdonald_P
from rcsknots.partitions import partitions
import rk_c_2box_largeN as L
import rk_c_3box_largeN as T3
import rk_c_4box_so_f as S4
import rk_f_cover_fixed22 as F22
import rk_c_fast_composite3 as C3

OUT = "data/refined_kauffman/daha"
q, t = R * R, S * S
nu = U / S
A = nu * nu
ZERO, ONE = conv(0), conv(1)
half = conv(Fraction(1, 2))
KMAX = 4
PARTS = [()] + [lam for n in range(1, KMAX + 1) for lam in partitions(n)]


def ser(x):
    return ({tuple(e): int(c) for e, c in x.n.to_dict().items()}, {tuple(e): int(c) for e, c in x.d.to_dict().items()})


def deser(p):
    return RF(CTX.from_dict(p[0]), CTX.from_dict(p[1]), normalize=False)


def cover_pairs():
    pair = F22.pair_values()
    w = S / R
    for n in (3, 4):
        for lam in partitions(n):
            v = w**n * L.Rdata_nu("".join(map(str, lam))) * L.Psym(lam)
            pair[(lam, ())] = pair[((), lam)] = v
    for alpha in partitions(3):
        # R_{(α;1)}(A^{1/2})：rk_c_bounded_save.run_fixed（r の次数ごとに冪の範囲を絞る補間、使わない N_A で全係数一致）で保存したもの
        Rc = F22.load_poly_A(os.path.join(OUT, f"Rcomp_A_m1_{''.join(map(str, alpha))}.pkl"), 1)
        dimc = F22.sym_rf(dim_composite(alpha, (1,)))
        cs = {mu: C3.psi_A(alpha, i, q, t, A) for i, mu in C3.lower(alpha)}
        # 主特殊化の恒等式（ev は乗法的）：P_α(t^ρ) P_1(t^ρ) = K_{(α;1)} + Σ_μ c_μ P_μ(t^ρ)
        assert (L.Psym(alpha) * L.Psym((1,)) - dimc - sum((c * L.Psym(mu) for mu, c in cs.items()), ZERO)).is_zero(), alpha
        v = Rc * dimc + sum((c * pair[(mu, ())] for mu, c in cs.items()), ZERO)
        pair[(alpha, (1,))] = pair[((1,), alpha)] = v
    return pair


def cover_Z(pair):
    PY = {(): {(): ONE}}
    for n in range(1, KMAX + 1):
        PY.update(macdonald_P(n, t, q, lambda f: conv(f)))      # P_λ(Y; t, q) の冪和表示
    T_ = lambda lam: tuple(sum(1 for x in lam if x > j) for j in range(lam[0])) if lam else ()
    Zp = {(): ONE}
    for lam in PARTS:
        for mu in PARTS:
            n = sum(lam) + sum(mu)
            if n == 0 or n > KMAX:
                continue
            term = T3.pmul(PY[T_(lam)], PY[T_(mu)], KMAX)
            Zp = T3.padd(Zp, term, conv((-1)**n) * pair[(lam, mu)])
    return Zp


def main():
    t0 = time.time()
    Sinv = S4.K_matrix()
    F = S4.plog(S4.so_Z(S4.so_W(), Sinv))
    print(f"SO の f {time.time()-t0:.0f}s", flush=True)
    pair = cover_pairs()
    print(f"被覆の ⟨P_λ P_μ^∨⟩ {time.time()-t0:.0f}s", flush=True)
    Fc = S4.plog(cover_Z(pair))
    print(f"被覆の f {time.time()-t0:.0f}s", flush=True)
    G = {k: F.get(k, ZERO) - half * Fc.get(k, ZERO) for k in set(F) | set(Fc)}
    old = pickle.load(open(os.path.join(OUT, "largeN_3box_m1.pkl"), "rb"))
    for key, X in (("f", F), ("fcov", Fc), ("g", G)):
        ok = all((X[eval(k)] - from_sympy(sp.sympify(v, locals={"r": T3.r_, "s": T3.s_, "u": T3.u_}))).is_zero() for k, v in old[key].items())
        print(f"次数 ≤ 3 の {key} が largeN_3box_m1 と一致：{ok}", flush=True)
    for k in partitions(4):
        for lab, X in (("f", F), ("fcov", Fc), ("g", G)):
            den = sp.factor(sp.fraction(sp.cancel(X[k].to_sympy()))[1]) if not X[k].is_laurent() else "単項式"
            print(f"{lab}_{k}：分母 {den}", flush=True)
    pickle.dump({key: {str(k): ser(v) for k, v in X.items()} for key, X in (("f", F), ("fcov", Fc), ("g", G))},
                open(os.path.join(OUT, "largeN_4box_m1.pkl"), "wb"))
    print(f"保存 {time.time()-t0:.0f}s", flush=True)
    return F, Fc, G


if __name__ == "__main__":
    main()
