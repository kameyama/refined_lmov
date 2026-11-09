r"""② T(2,2m+1) の大きな N の 2 箱（data/refined_kauffman/daha/largeN_2box_m{m}.pkl、scripts/rk_c_2box_largeN_m.py）に、trefoil で行った検査をかける（progress.md 続き 32）。
  極と整数性、c = 1 部分の σ 対称性、述べ直しの判定式（ĝ − ĝ|_{t=q}(ν → ν~) が (q − t)(ν² − t) で割り切れるか）、商の 𝐚 ブロックの符号、
  二次量子化の指数（hatted の SO の f̂ を t の冪級数に展開）、被覆の brane–image の分割（KN の t 枠で G 展開）。
使い方：python scripts/rk_b_2box_structure_m.py m"""
import sys, pickle, io, contextlib
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_b_2box_structure as B
import rk_b_second_quant as Q

r, s, nu = B.r, B.s, B.nu
u = sp.Symbol("u", positive=True)


def load(m):
    import rk_torus_tag as TG
    d = pickle.load(open(f"data/refined_kauffman/daha/largeN_2box_{TG.tag(m)}.pkl", "rb"))     # m が文字列ならケーブルなど（例 K1、続き 141）
    get = lambda key: {eval(k): sp.cancel(sp.sympify(v, locals={"r": r, "s": s, "u": u}).subs(u, nu * s)) for k, v in d[key].items()}
    return get("f"), get("fcov"), get("g")


def hat(F):
    det = r**2 - r**-2
    h2 = sp.cancel((r * F[(2,)] + F[(1, 1)] / r) / det)
    h11 = sp.cancel((F[(2,)] / r + r * F[(1, 1)]) / det)
    return {(1,): F[(1,)], (2,): h2, (1, 1): h11, "sum": sp.cancel(h2 + h11)}


def main(m):
    f, fc, g = load(m)
    G = hat(g)
    print(f"== {m if isinstance(m, str) else f'T(2,{2*m+1})'} の 2 箱、大きな N", flush=True)
    for k in ((1,), (2,), (1, 1), "sum"):
        e = G[k]
        print(f"ĝ_{k}：分母 {sp.factor(sp.fraction(e)[1])}、c=1 が σ 対称 {B.c1_symmetric(e)}、refined だけの部分が (q−t)(ν²−t) で割り切れる {B.divisible_universal(e)}", flush=True)
    for k in ((1,), "sum"):
        sb = B.sign_blocks(G[k])
        coh = sum(1 for n_, sg in sb.values() if sg != "混在")
        print(f"ĝ_{k} の商の 𝐚 ブロックの符号：{coh}/{len(sb)} {sb}", flush=True)
    print("-- 二次量子化の指数（hatted、t の冪級数）", flush=True)
    for lab, F in (("SO の f̂", hat(f)), ("被覆の f̂^cov", hat(fc))):
        for k in ((2,), (1, 1)):
            Q.analyze(f"{lab}_{k}", F[k])
    f1 = f[(1,)]
    Q.analyze(f"[1] SO の f（大きな N の DAHA から）", f1)


if __name__ == "__main__":
    main(int(sys.argv[1]) if sys.argv[1].isdigit() else sys.argv[1])
