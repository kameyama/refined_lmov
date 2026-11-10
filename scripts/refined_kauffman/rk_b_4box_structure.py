r"""大きな N の trefoil の 4 箱のクロスキャップ部分（data/refined_kauffman/daha/largeN_4box_m1.pkl、rk_c_4box_largeN_g.py）の検査（progress.md 続き 139）。
偶数のセクターなので g = ½E/(1+t) + C で、(1 + t) の極は E の方向（冪和に偶数の部分をもつ方向）だけにあるはず（全 d の構造の予言）。
検査（rk_b_3box_structure と同じ項目）：
  冪和の成分 g_ρ・hatted ĝ_ρ = g_ρ/φ_ρ の分母（(1 + t) の極の位置）。
  Schur の成分 ĝ_R（|R| = 4）と、(1 + t) を払った (1 + t)ĝ_R：整数 Laurent か、c = 1 部分の σ 対称性、
  述べ直しの判定式の普遍因子 (q − t)(ν² − t)、商の 𝐚 ブロックの符号。"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rcsknots.characters import chi
from rcsknots.partitions import partitions
import rk_b_2box_structure as B
import rk_b_3box_structure as B3
import rk_c_4box_largeN_g as G4

r, s, nu = B.r, B.s, B.nu
u = sp.Symbol("u", positive=True)


def load(key="g"):
    d = pickle.load(open("data/refined_kauffman/daha/largeN_4box_m1.pkl", "rb"))[key]
    return {eval(k): sp.cancel(G4.deser(v).to_sympy().subs(u, nu * s)) for k, v in d.items() if sum(eval(k)) == 4}


def report(lab, e):
    num, den = sp.fraction(sp.cancel(e))
    lau = sp.Poly(den, r, s, nu).is_monomial
    line = f"{lab}：分母 {sp.factor(den)}"
    if lau:
        line += f"、整数係数 {all(c.is_integer for c in sp.Poly(num, r, s, nu).coeffs())}、c=1 が σ 対称 {B.c1_symmetric(e)}、普遍因子で割り切れる {B.divisible_universal(e)}"
    print(line, flush=True)
    if lau:
        sb = B.sign_blocks(e)
        coh = sum(1 for n_, sg in sb.values() if sg != "混在")
        print(f"   商の 𝐚 ブロックの符号：{coh}/{len(sb)} {sb}", flush=True)


if __name__ == "__main__":
    g = load("g")
    for rho in partitions(4):
        print(f"g_{rho}：分母 {sp.factor(sp.fraction(g[rho])[1])}", flush=True)
    gs, gh = B3.schur(g, 4), B3.schur(g, 4, hat=True)
    for R_ in partitions(4):
        report(f"g_{R_}", gs[R_])
        report(f"ĝ_{R_}", gh[R_])
        if not sp.Poly(sp.fraction(gh[R_])[1], r, s, nu).is_monomial:
            report(f"(1+t)ĝ_{R_}", sp.cancel((1 + s**2) * gh[R_]))
