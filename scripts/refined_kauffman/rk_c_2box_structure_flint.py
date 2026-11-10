r"""rk_b_2box_structure_m.py と同じ検査を flint の有理関数（flint_rf3.RF）で行う版（続き 141。ケーブル K1 の g は sympy では遅い）。
  g, f, f^cov は rk_c_2box_largeN_m の関数で組み直す。変数は (r, s, u) → (r, s, ν)（u = ν s、第 3 変数を ν と読み替える）。
  検査：hatted ĝ の分母、c = 1 部分の σ 対称性（σ：r → 1/s, s → 1/r）、X = ĝ − ĝ|_{t=q}(ν → ν (q/t)^{1/2}) が (q − t)(ν² − t) で割り切れるか、
  商の 𝐚 ブロックの符号（r = −𝐭𝐪, s = 𝐪, ν = 𝐚）、二次量子化の指数（rk_b_second_quant.analyze と同じ集計）。
使い方：python scripts/rk_c_2box_structure_flint.py K1（照合用に m = 2 なども可）"""
import sys, itertools
from fractions import Fraction
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import flint
from rcsknots.refined_kauffman.flint_rf3 import RF, CTX, R, S

half = RF.coerce(Fraction(1, 2))


def emap(e, M, sign=None):
    """単項式の指数を線形写像 M（(a,b,c) → 3 つ組、負も可）で写す。sign(a,b,c) は係数に掛ける ±1。"""
    def one(p):
        d = {}
        for ex, c in p.to_dict().items():
            k = M(*ex)
            d[k] = d.get(k, 0) + int(c) * (sign(*ex) if sign else 1)
        lo = [min(k[i] for k in d) for i in range(3)]
        return CTX.from_dict({tuple(k[i] - lo[i] for i in range(3)): c for k, c in d.items() if c}), lo
    n, ln = one(e.n)
    d, ld = one(e.d)
    sh = [ln[i] - ld[i] for i in range(3)]
    mono = lambda v: CTX.from_dict({tuple(max(x, 0) for x in v): 1})
    return RF(n * mono(sh), d * mono([-x for x in sh]))


to_nu = lambda e: emap(e, lambda a, b, c: (a, b + c, c))                  # u = ν s
sigma = lambda e: emap(e, lambda a, b, c: (-b, -a, c))
at_q = lambda e: emap(e, lambda a, b, c: (a + b, 0, c))                   # s = r
nu_shift = lambda e: emap(e, lambda a, b, c: (a + c, b - c, c))           # ν → ν r/s
NU = RF(CTX.gens()[2])


def hat(F):
    det = R * R - 1 / (R * R)
    h2 = (R * F[(2,)] + F[(1, 1)] / R) / det
    h11 = (F[(2,)] / R + R * F[(1, 1)]) / det
    return {(1,): F[(1,)], (2,): h2, (1, 1): h11, "sum": h2 + h11}


def c1_symmetric(e):
    (a0, b0, c0), = [list(e.d.to_dict())[0]]
    vp = (a0 + b0) % 2
    part = CTX.from_dict({ex: c for ex, c in e.n.to_dict().items() if (ex[0] + ex[1] - vp) % 2 == 0})
    out = []
    for par in (0, 1):                                                     # ν の偶奇ごと
        p = CTX.from_dict({ex: c for ex, c in part.to_dict().items() if (ex[2] - c0) % 2 == par})
        if not p.is_zero():
            v = RF(p, e.d)
            out.append((v - sigma(v)).is_zero())
    return all(out)


def refined_only(e):
    X = e - nu_shift(at_q(e))
    return X, X / ((R * R - S * S) * (NU * NU - S * S))


def divisible_universal(e):
    X, _ = refined_only(e)
    return all(X.num_divisible(f) for f in (R * R - S * S, NU * NU - S * S))


def laurent_terms(e, order):
    """Laurent な e の項 {指数（order の並び）: Fraction}。"""
    assert e.is_laurent(), "Laurent でない"
    (dm, dc), = e.d.to_dict().items()
    return {tuple(int(ex[i] - dm[i]) for i in order): Fraction(int(c), int(dc)) for ex, c in e.n.to_dict().items()}


def sign_blocks(e):
    _, Qt = refined_only(e)
    B = emap(Qt, lambda a, b, c: (a, a + b, c), sign=lambda a, b, c: (-1)**a)   # (T, Q, A) = (𝐭, 𝐪, 𝐚)
    by = {}
    for (ea, _, _), c in laurent_terms(B, (2, 1, 0)).items():
        by.setdefault(ea, []).append(c)
    return {k: (len(v), "+" if all(c > 0 for c in v) else "−" if all(c < 0 for c in v) else "混在") for k, v in sorted(by.items())}


def series_coeffs(e, extra=2):
    for p in (1, 2):
        L = e * (1 - S**(2 * p))
        if L.is_laurent():
            break
    else:
        raise ValueError("分母が 1 − t または 1 − t² ではない")
    T = laurent_terms(L, (2, 1, 0))                                        # (β, k, a) = (ν, s, r)
    kmax = max(k for (_, k, _) in T) + extra * 2 * p
    out = {}
    for (b, j, a), c in T.items():
        for k in range(j, kmax + 1, 2 * p):
            out[(b, k, a)] = out.get((b, k, a), 0) + c
    return {key: c for key, c in out.items() if c != 0}, p, kmax


def analyze(lab, e):
    C, p, kmax = series_coeffs(e)
    assert all(c.denominator == 1 for c in C.values()), "整数でない指数"
    blocks = {}
    for (b, k, a), c in C.items():
        blocks.setdefault((b, k), []).append(c)
    mixed = sorted(key for key, cs in blocks.items() if not (all(c > 0 for c in cs) or all(c < 0 for c in cs)))
    sg = lambda c: 1 if c > 0 else -1
    rules = [(al, ga, de) for al, ga, de in itertools.product((0, 1), repeat=3)
             if len({sg(c) * (-1)**(al * a + ga * k + de * b) for (b, k, a), c in C.items()}) == 1]
    perbeta = {}
    for b in sorted({b for (b, _, _) in C}):
        sub = {(k, a): c for (bb, k, a), c in C.items() if bb == b}
        ok = [(al, ga) for al, ga in itertools.product((0, 1), repeat=2) if len({sg(c) * (-1)**(al * a + ga * k) for (k, a), c in sub.items()}) == 1]
        perbeta[b] = ok[0] if ok else None
    print(f"{lab}：塔 1/(1 − t^{p})、s^k を k ≤ {kmax} まで。ブロック (β, k) は {len(blocks)} 個、そのうち符号が混在 {len(mixed)} 個 {mixed[:12]}{' …' if len(mixed) > 12 else ''}", flush=True)
    print(f"    単項式ごとの偶奇の符号則（全体）：{rules if rules else 'なし'}；ν^β ごとの符号則 (α, γ)：{perbeta}", flush=True)


def denom(e):
    fac = flint.fmpz_mpoly_ctx.get(("r", "s", "nu"), "lex")
    c, fs = fac.from_dict(e.d.to_dict()).factor()
    return " ".join([str(c)] + [f"({f})^{k}" if k > 1 else f"({f})" for f, k in fs])


def main(m, k=2):
    import time
    import rk_c_2box_largeN as L
    import rk_c_2box_largeN_m as M
    import rk_torus_tag as TG
    t0 = time.time()
    W, co = M.so_parts(m, k)
    f = L.so_f(W, co)
    fc = M.cover_f(m, TG.knot_name(m, k), k)
    g = {key: f[key] - half * fc[key] for key in f}
    f, fc, g = ({key: to_nu(v) for key, v in F.items()} for F in (f, fc, g))
    print(f"== {m if isinstance(m, str) else f'T(2,{2*m+1})'} の 2 箱、大きな N（組み立て {time.time() - t0:.0f} 秒）", flush=True)
    G = hat(g)
    for key in ((1,), (2,), (1, 1), "sum"):
        e = G[key]
        print(f"ĝ_{key}：分母 {denom(e)}、c=1 が σ 対称 {c1_symmetric(e)}、refined だけの部分が (q−t)(ν²−t) で割り切れる {divisible_universal(e)}", flush=True)
    print("ĝ_(2) + ĝ_(1,1) は整数 Laurent", G["sum"].is_integer_laurent(), flush=True)
    for key in ((1,), "sum"):
        sb = sign_blocks(G[key])
        coh = sum(1 for n_, sg in sb.values() if sg != "混在")
        print(f"ĝ_{key} の商の 𝐚 ブロックの符号：{coh}/{len(sb)} {sb}", flush=True)
    print("-- 二次量子化の指数（hatted、t の冪級数）", flush=True)
    for lab, F in (("SO の f̂", hat(f)), ("被覆の f̂^cov", hat(fc))):
        for key in ((2,), (1, 1)):
            analyze(f"{lab}_{key}", F[key])
    analyze("[1] SO の f（大きな N の DAHA から）", f[(1,)])
    print(f"計 {time.time() - t0:.0f} 秒", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]) if sys.argv[1].isdigit() else sys.argv[1])
