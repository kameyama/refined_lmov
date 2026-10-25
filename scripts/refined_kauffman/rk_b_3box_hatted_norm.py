r"""trefoil の 3 箱の hatted の規格化（theory §11.17、progress.md 続き 56）。
  注意：largeN_3box_m1.pkl の "g"・"f"・"fcov" は冪和基底（キー ρ は p_ρ(Y) の係数）。Schur 成分は g_λ = Σ_ρ χ_λ(ρ) G[ρ]。
  (1) t = q で、冪和成分と Schur 成分の ν 偶奇 × q↔1/q 偶奇。
  (2) hatted ĝ_λ = M(q^{-1})^{-1} g_λ（Schur）の分母：refined、t = q（U 一般と U = q^N）。
  (3) refined の ν 偶・奇の σ 対称性（Schur と hatted）。"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rk_b_multibox_sigma import sigma, c1, nterms, r_, s_, nu, loc, u_
from rcsknots.bps import M_matrix, partitions, chi

c2 = lambda e: sp.cancel(e - c1(e))
OUT = "data/refined_kauffman/daha"


def load3():
    d = pickle.load(open(f"{OUT}/largeN_3box_m1.pkl", "rb"))
    P = {tuple(eval(k)) if k != "(1,)" else (1,): sp.cancel(sp.sympify(v, locals=loc)) for k, v in d["g"].items()}
    parts = partitions(3)
    S = {lam: sp.cancel(sum((chi(lam, rho) * P[rho] for rho in parts), sp.Integer(0))) for lam in parts}
    return P, S


def parity_table(e, lab):
    x = sp.cancel(e.subs(s_, r_))
    out = []
    for plab, part in (("ν偶", c1(x)), ("ν奇", c2(x))):
        Sy = sp.cancel((part + part.subs(r_, 1 / r_)) / 2); A = sp.cancel((part - part.subs(r_, 1 / r_)) / 2)
        out.append(f"{plab}: 対称 {0 if Sy == 0 else nterms(Sy)} 反対称 {0 if A == 0 else nterms(A)}")
    print(f"  {lab}: " + "；".join(out), flush=True)


if __name__ == "__main__":
    P, S = load3()
    parts = partitions(3)
    print("== (1) t = q の偶奇（u = U_D = νs を ν に直したもの）", flush=True)
    sub = lambda e: sp.cancel(e.subs(u_, nu * s_))
    for rho in parts:
        parity_table(sub(P[rho]), f"冪和 p_{rho}")
    for lam in parts:
        parity_table(sub(S[lam]), f"Schur g_{lam}")
    print("== (2) hatted ĝ_λ = M(q^{-1})^{-1} g_λ の分母", flush=True)
    Mi = M_matrix(3, 1 / r_).inv()
    vec = sp.Matrix([S[lam] for lam in parts])
    hat = Mi * vec
    for lam, h in zip(parts, hat):
        h = sp.cancel(h); htq = sp.cancel(h.subs(s_, r_))
        print(f"  ĝ_{lam}: refined の分母 {sp.factor(sp.fraction(h)[1])}；t=q（U 一般）{sp.factor(sp.fraction(htq)[1])}；"
              + "；".join(f"U=q^{N}: {sp.factor(sp.fraction(sp.cancel(htq.subs(u_, r_**(2*N))))[1])}" for N in (2, 4)), flush=True)
        parity_table(sub(h), f"   hatted ĝ_{lam}")
    print("== (3) refined の σ 対称性（ν 固定）", flush=True)
    for lam in parts:
        for plab, part in (("ν偶", c1(sub(S[lam]))), ("ν奇", c2(sub(S[lam])))):
            a = sp.cancel(part - sigma(part)); print(f"  Schur g_{lam} の {plab}: σ 対称 {a == 0}（非対称 {0 if a == 0 else nterms(a)} 項）", flush=True)
    for lam, h in zip(parts, hat):
        for plab, part in (("ν偶", c1(sub(h))), ("ν奇", c2(sub(h)))):
            a = sp.cancel(part - sigma(part)); print(f"  hatted ĝ_{lam} の {plab}: σ 対称 {a == 0}（非対称 {0 if a == 0 else nterms(a)} 項）", flush=True)
