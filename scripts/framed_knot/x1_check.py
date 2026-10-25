"""別研究 [2,1] の X1：KNTZ の refined 𝐁 の実装（rcsknots.double_twist_refined）の照合。

1. KNTZ のパッケージの値（scripts/framed_knot/x1_kntz_package.wls → data/framed_knot/reference/kntz_package/values.txt）と全係数一致（𝐅、𝒫、非 refined の H）。
2. 対称表現 [r]（r = 1, 2）：GNSSS 1512.07883 の四重次数の式（scripts/framed_knot/t3_pretzel.py の quadruple、第 2 期の入力）と一致。
   反対称表現 [1^r]：GGS 1304.3481 のミラー対称性（higher_v17.tex :1294-1299、t_r ↔ t_c）を介して GNSSS と一致。
3. t = q（σ = 1, a = A）で X0 の非 refined の式（KNTZ の長方形、Fang–Zhou の定理）に戻る。
4. (Y6) t_c 次数の三重次数 P(𝐚, 𝐪, 𝐭) = 𝒫(𝐚, 𝐐 = 𝐪, 𝐭_r = 1/𝐪, 𝐭_c = 𝐭𝐪)（第 2 期 T3 の規則）を KN の小文字の rCS に直し、
   三葉結び目 K_{1,1} で EHA の値（data/framed_knot/torus/T2_3_*.json、[2,2] はこのスクリプトで GN の構成から計算）と比べる。
使い方：python scripts/framed_knot/x1_check.py
"""
import json
import os
import sys
import time
from pathlib import Path

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rcsknots import double_twist as U  # noqa: E402
from rcsknots.double_twist_refined import A, Fref, Qb, a, poincare, poincare_pre, q, sig, t, tc, tr  # noqa: E402

PKG = Path("data/framed_knot/reference/kntz_package/values.txt")


def mma(s):
    """Mathematica の InputForm（a, Q, tr, tc, A, q, t）→ sympy。"""
    s = s.replace("Sqrt[", "sqrt(").replace("[", "(").replace("]", ")").replace("^", "**")
    # 非 refined の H には QBinomial[n, k, q²] が評価されずに残る（Wolfram Engine 15.0）。ガウスの二項係数として評価する
    return sp.sympify(s, locals={"a": a, "Q": Qb, "tr": tr, "tc": tc, "A": A, "q": q, "t": t, "QBinomial": U.gauss})


def load_package(keep=lambda kind, args: True):
    """パッケージの値 [(種類, 引数, 値)]。keep(種類, 引数) が真の行だけ値を読む。"""
    out = []
    for line in open(PKG):
        kind, args, val = line.rstrip("\n").split("\t")
        args = sp.sympify(args.replace("{", "[").replace("}", "]"))
        if keep(kind, args):
            out.append((kind, args, mma(val)))
    return out


def check_package(log=print):
    bad = 0
    for kind, args, val in load_package():
        t0 = time.time()
        if kind == "F":
            lam, m, s = tuple(args[0]), args[1], args[2]
            ours = Fref(lam, m, s)
        elif kind == "P":
            ours = poincare(*args)
        else:
            r, s, m, n = args
            ours = U.H_rect(r, s, m, n).subs({U.A: A, U.q: q})
        ok = sp.expand(sp.cancel(sp.together(ours - val))) == 0
        bad += not ok
        log(f"{kind} {list(args)}：{'一致' if ok else '不一致'}（{time.time() - t0:.1f}s）")
    return bad


def tq_limit(r, s, m, n):
    """置換の前の式で t = q, σ = 1, a = A とし、非 refined の H_rect と比べる。"""
    pre = poincare_pre(r, s, m, n).subs({t: q, sig: 1, a: A}, simultaneous=True)
    ref = U.H_rect(r, s, m, n).subs({U.A: A, U.q: q})
    return sp.expand(sp.cancel(sp.together(pre - ref))) == 0


def triply(r, s, m, n):
    """t_c 次数の三重次数（DGR の太字 𝐚, 𝐪, 𝐭 を a, Q, t_c の記号で返す：𝐪 = Q, 𝐭 = t_c）。

    GGS 1304.3481（higher_v17.tex :1244-1248）：R 行の長方形では 𝐪 = R𝐐 − 𝐭_r + 𝐭_c。よって 𝐐 → 𝐪^R, 𝐭_r → 1/𝐪, 𝐭_c → 𝐭𝐪（R = s）。
    第 2 期 T3 の規則（𝐐 = 𝐪）は R = 1 の場合。"""
    P = poincare(r, s, m, n)
    return sp.expand(P.subs({Qb: Qb ** s, tr: 1 / Qb, tc: tc * Qb}, simultaneous=True))


def bold_dict(expr):
    e = sp.expand(expr * (a * Qb * tc) ** 200)
    return {(i - 200, j - 200, k - 200): int(c) for (i, j, k), c in sp.Poly(e, a, Qb, tc).terms()}


def main():
    print("1. パッケージとの照合", flush=True)
    bad = check_package(lambda s: print("  ", s, flush=True))
    print("   不一致の数", bad, flush=True)
    print("3. t = q の極限", flush=True)
    for rs in ((2, 1), (1, 2), (2, 2), (3, 1)):
        for mn in ((1, 1), (-1, 1), (2, 1), (-2, 1), (2, 2)):
            print("  ", rs, mn, "一致" if tq_limit(*rs, *mn) else "不一致", flush=True)


if __name__ == "__main__":
    main()


def eha_trefoil(lam):
    """T(2,3) の rCS_λ（KN の (A, Q, T) の指数、reduced、framing 0）。保存値があれば読み、無ければ GN の構成で計算して (N2') で正規化する。"""
    name = "".join(map(str, lam))
    path = Path(f"data/framed_knot/torus/T2_3_{name}.json")
    if path.exists():
        d = json.loads(path.read_text())
        return {tuple(map(int, k.split(","))): v for k, v in d["poly"].items()}
    from fractions import Fraction
    from rcsknots.exact import reconstruct
    from rcsknots.partitions import n_stat, conjugate as conj
    raw = reconstruct(lambda G: G.torus_reduced_apoly(lam, 2, -3), log=lambda s: None)
    # (N2')：M_λ = M_□^{|λ|} T_λ^{6}、M_□ = −a q^{3/2} t^{−3/2}（data/framed_knot/torus/T2_3_1.json）、T_λ = q^{n(λ')} t^{−n(λ)}
    k = sum(lam)
    sign = (-1) ** k
    ea, eq, et = Fraction(k), Fraction(3 * k, 2) + 6 * n_stat(conj(lam)), Fraction(-3 * k, 2) - 6 * n_stat(lam)
    out = {(int(2 * (x + ea)), int(2 * (y + eq)), int(2 * (z + et))): sign * c for (x, y, z), c in raw.items()}
    meta = {"knot": "T2_3", "color": list(lam), "gn_slope": [2, -3], "normalization": "(N2') M_lambda = M_box^|lambda| T_lambda^6 (x1_check.py)",
            "variables": "A=a^(1/2), Q=q^(1/2), T=t^(1/2) of 1703.05408; exponents of (A,Q,T)",
            "poly": {",".join(map(str, kk)): v for kk, v in sorted(out.items())}}
    path.write_text(json.dumps(meta, indent=1))
    return out


def kntz_trefoil_small(lam):
    """KNTZ の 𝒫(K_{1,1}) を t_c 次数の三重次数（第 2 期 T3 の規則）にし、KN の小文字の指数に直す（scripts/framed_knot/t4_colored.bold_to_small）。"""
    from t4_colored import bold_to_small
    r, s = lam[0], len(lam)
    return bold_to_small(bold_dict(triply(r, s, 1, 1)))


def y6_trefoil(log=print):
    bad = 0
    for lam in ((1,), (2,), (1, 1), (2, 2)):
        ours, ref = kntz_trefoil_small(lam), eha_trefoil(lam)
        ok = ours == ref
        if not ok:
            ratio = {tuple(x - y for x, y in zip(k1, k2)) for k1 in ours for k2 in ref}
            log(f"   T(2,3) {lam}：不一致（項数 {len(ours)} / {len(ref)}）")
        else:
            log(f"   T(2,3) {lam}：EHA と一致（{len(ref)} 項）")
        bad += not ok
    return bad


def y6_twist(ps, log=print):
    """ツイスト結び目 K_{p,1} の [1], [2], [1,1]（t_c 次数の三重次数、太字）を、第 2 期 T3 の入力（NRZS の式＋EGP で決めた δ'、scripts/framed_knot/t3_twist.py）と比べる。"""
    from t3_twist import inputs
    from t4_colored import tc_colored
    from t4_egp import egp_solve
    bad = 0
    for p in ps:
        P1, gens = inputs(p)
        d = next(x / 2 for x in range(-24, 25) if egp_solve(gens, P1, [x / 2]) is not None)
        P2, P11 = tc_colored(gens, lambda x: d)
        ok = [bold_dict(triply(*rs, p, 1)) == ref for rs, ref in (((1, 1), P1), ((2, 1), P2), ((1, 2), P11))]
        bad += not all(ok)
        log(f"   K_{{{p},1}}：[1], [2], [1,1] {['一致' if x else '不一致' for x in ok]}")
    return bad
