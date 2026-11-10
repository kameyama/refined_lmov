r"""大きな N の trefoil の 3 箱のクロスキャップ部分（data/refined_kauffman/daha/largeN_3box_m1.pkl、scripts/rk_c_3box_largeN.py）の検査（progress.md 続き 33）。
奇数のセクターなので E = 0 で、g_R（|R| = 3）は (1 + t) の極をもたないはず。hatted ĝ_R = Σ_ρ χ_R(ρ) g_ρ / φ_ρ（φ_ρ は M(q^{−1}) の冪和基底での固有値）。
検査：整数 Laurent、c = 1 部分の σ 対称性、述べ直しの判定式（ĝ − ĝ|_{t=q}(ν → ν~) が (q − t)(ν² − t) で割り切れるか）、商の 𝐚 ブロックの符号、二次量子化の指数。
T(2,2m+1)：python rk_b_3box_structure.py m（largeN_3box_m{m}.pkl、続き 139）。"""
import sys, os, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rcsknots.characters import chi
from rcsknots.partitions import partitions
from rcsknots.bps import B as Bfun
import rk_b_2box_structure as B
import rk_b_second_quant as Q

r, s, nu = B.r, B.s, B.nu
u = sp.Symbol("u", positive=True)


def phi(rho):
    return sp.cancel(sum(chi(sg, rho) * Bfun(sg, 1 / r) for sg in partitions(sum(rho))))


M = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 1


_ASM = {}


def load(key="g"):
    if M >= 4:                                       # 保存した文字列は sympy が読めない大きさなので、組み立てをやり直す（続き 139）
        if not _ASM:
            import rk_c_3box_check_fixedN as CF
            _ASM["g"], _ASM["f"] = CF.assemble(M)
        return {k: sp.cancel(v.to_sympy().subs(u, nu * s)) for k, v in _ASM[key].items() if sum(k) == 3}
    path = f"data/refined_kauffman/daha/largeN_3box_m{M}.pkl"
    if not os.path.exists(path):                     # g だけを保存したもの（rk_c_3box_g_only.py）
        path = f"data/refined_kauffman/daha/largeN_3box_m{M}_g.pkl"
    d = pickle.load(open(path, "rb"))
    if key not in d:
        return None
    return {eval(k): sp.cancel(sp.sympify(v, locals={"r": r, "s": s, "u": u}).subs(u, nu * s)) for k, v in d[key].items()}


def schur(F, n, hat=False):
    return {R_: sp.cancel(sum(chi(R_, rho) * F.get(rho, 0) / (phi(rho) if hat else 1) for rho in partitions(n))) for R_ in partitions(n)}


if __name__ == "__main__":
    g = load("g")
    gs, gh = schur(g, 3), schur(g, 3, hat=True)
    for R_ in partitions(3):
        for lab, e in (("g", gs[R_]), ("ĝ", gh[R_])):
            d_ = sp.factor(sp.fraction(e)[1])
            lau = sp.Poly(sp.fraction(e)[1], r, s, nu).is_monomial
            print(f"{lab}_{R_}：分母 {d_}" + (f"、c=1 が σ 対称 {B.c1_symmetric(e)}、普遍因子で割り切れる {B.divisible_universal(e)}" if lau else ""), flush=True)
        if sp.Poly(sp.fraction(gh[R_])[1], r, s, nu).is_monomial:
            sb = B.sign_blocks(gh[R_])
            coh = sum(1 for n_, sg in sb.values() if sg != "混在")
            print(f"   ĝ_{R_} の商の 𝐚 ブロックの符号：{coh}/{len(sb)} {sb}", flush=True)
    f = load("f")
    if f is None:
        print("f が保存されていないので二次量子化の検査は飛ばす", flush=True)
        sys.exit(0)
    fh = schur(f, 3, hat=True)
    for R_ in partitions(3):
        try:
            Q.analyze(f"二次量子化：SO の f̂_{R_}", fh[R_])
        except (AssertionError, ValueError) as ex:
            print(f"二次量子化：SO の f̂_{R_}：{type(ex).__name__} {str(ex)[:100]}", flush=True)
