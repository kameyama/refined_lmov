"""仮説 H_u（D12）：N_λ := rCS^{red}_λ × Π_{□≠(1,1)}（unknot の分子の太字の因子）の太字の係数は非負（単項式倍と全体の符号を除く）。

根拠（report.md D12）：GN/EHA の unreduced の値は f̃_λ[1 − a] で、s_μ[1 − a] は hook (d−s, 1^s) で (1 − a)(−a)^s、他は 0。
よって unreduced/(1 − a) の a の行は hook の Schur 係数 c_s（符号 (−1)^s は太字で (𝐚²𝐭)^s に吸収される）。
reduced は unreduced ÷ unknot_λ で、unknot_λ の分子 Π_□(t^{l'} − a q^{a'}) は太字で Π_□(𝐪^{2l'} + 𝐚²𝐪^{2a'}𝐭^{2a'+1})（全て正）、
分母 Π_□(1 − q^{a(□)} t^{l(□)+1}) が負の符号を作りうる。そこで N_λ = rCS^{red}_λ × Π_{□≠(0,0)}(分子の因子) を調べる。
[1] は N = rCS_{[1]} で、Mellit の定理（P̃_{n,m}·1 の Schur 正値性）から正。[2], [1,1] も reduced が正なので自明。[2,1] が本命。

使い方：PYTHONPATH=src:scripts python scripts/dgr_hu.py check [結び目 ...]      既知の reduced から N_λ（|λ| ≤ 3）の符号を数える
        PYTHONPATH=src:scripts python scripts/dgr_hu.py convention T2_3 T2_5     GN の hook の Schur 係数と N の a の行の照合（有理数の評価点）
        PYTHONPATH=src:scripts python scripts/dgr_hu.py highbox m n N            T(m,n) の N 箱の rCS_λ を GN から補間して N_λ の符号を数える
        PYTHONPATH=src:scripts python scripts/dgr_hu.py series [結び目 ...]     仮説 H_u″：N_λ の a の行 ÷ Π_□(1 − q^{a(□)} t^{l(□)+1}) の冪級数（P 正規化の hook の Schur 係数）の符号
        PYTHONPATH=src:scripts python scripts/dgr_hu.py series-highbox m n N   同じ検査を N 箱で
記録：logs/dgr_hu.log、data/dgr_lmov/hu_check.json。
"""
import json
import sys
import time
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from framing_check import load_reduced  # noqa: E402
from t4_colored import bold_to_small  # noqa: E402
from rcsknots.bps import to_bold  # noqa: E402
from rcsknots.partitions import arm, boxes, leg, partitions  # noqa: E402

LOG = Path("logs/dgr_hu.log")
OUT = Path("data/dgr_lmov/hu_check.json")
TORUS = {"T2_3": (2, 3), "T2_5": (2, 5), "T2_7": (2, 7), "T2_9": (2, 9), "T2_11": (2, 11), "T3_4": (3, 4), "T3_5": (3, 5), "T3_7": (3, 7), "T4_5": (4, 5)}
ALGEBRAIC = ["4_1", "K1", "K3", "K6a", "K7", "K8", "K9"]
NONALG = ["K2", "K4", "K5a", "K5b", "K5c", "K6b"]


def log(s):
    print(s, flush=True)
    with LOG.open("a") as f:
        f.write(s + "\n")


def bold(P):
    return {k: v for k, v in to_bold(P, 1 if any(k[0] % 2 for k in P) else 0).items() if v}


def mul(A, B):
    out = defaultdict(int)
    for k1, v1 in A.items():
        for k2, v2 in B.items():
            out[(k1[0] + k2[0], k1[1] + k2[1], k1[2] + k2[2])] += v1 * v2
    return {k: v for k, v in out.items() if v}


def unknot_numerator_factors(lam):
    """unknot_λ の分子の因子（太字）。箱 (i, j)（0 始まり、i = l'、j = a'）の因子 t^i − a q^j = 𝐪^{2i} + 𝐚²𝐪^{2j}𝐭^{2j+1}。箱 (0,0) は除く。"""
    out = []
    for i, j in boxes(lam):
        if (i, j) == (0, 0):
            continue
        out.append({(0, 2 * i, 0): 1, (2, 2 * j, 2 * j + 1): 1})
    return out


def N_lambda(red_small, lam):
    N = bold(red_small)
    for f in unknot_numerator_factors(lam):
        N = mul(N, f)
    return N


def sign_summary(N):
    """全体の符号（多数派）を掛けた後の負の項の数と、a の行ごとの符号の一様性。"""
    pos = sum(1 for v in N.values() if v > 0)
    neg = sum(1 for v in N.values() if v < 0)
    s = 1 if pos >= neg else -1
    rows = defaultdict(lambda: [0, 0])
    for (a, _, _), v in N.items():
        rows[a][0 if s * v > 0 else 1] += 1
    uniform = {a: (p == 0 or n == 0) for a, (p, n) in rows.items()}
    return {"terms": len(N), "neg": min(pos, neg), "rows": {str(a): tuple(pn) for a, pn in sorted(rows.items())},
            "rows_uniform": all(uniform.values()), "bad_rows": sorted(a for a, u in uniform.items() if not u)}


def check(knots):
    res = {}
    for K in knots:
        try:
            R = load_reduced(K, 3)
        except AssertionError as e:
            log(f"{K}: {e}")
            continue
        res[K] = {}
        for lam in sorted(R, key=lambda l: (sum(l), l)):
            N = N_lambda(R[lam], lam)
            sm = sign_summary(N)
            red = sign_summary(bold(R[lam]))
            res[K][str(lam)] = {"N": sm, "red_neg": red["neg"], "red_terms": red["terms"]}
            log(f"{K:6s} λ={str(lam):10s} reduced 負 {red['neg']:4d}/{red['terms']:5d} → N 負 {sm['neg']:4d}/{sm['terms']:5d} 行の符号一様 {sm['rows_uniform']} 非一様な行 {sm['bad_rows']}")
    return res


# ---- 規約の照合：GN の hook の Schur 係数 ----------------------------------------------------------

def schur_hook_coeffs(G, lam, n, m):
    """f̃_λ = P_λ[p_k ↦ P̃_{kn,km}]·1 の hook (d−s, 1^s) の Schur 係数 c_s（s = 0..d−1）。c_μ = Σ_ρ f_ρ χ^μ(ρ)。"""
    from rcsknots.characters import chi
    vec = G.apply_symmetric(G.P(sum(lam))[lam], n, m, {(): G.one})
    f = G.to_p(vec)
    d = sum(lam) * n
    out = []
    for s in range(d):
        mu = (d - s,) + (1,) * s
        c = G.zero
        for rho, fr in f.items():
            x = chi(mu, rho)
            if x:
                c = c + fr * x
        out.append(c)
    return out


def convention(K, points=((Fraction(3, 7), Fraction(5, 11)), (Fraction(2, 9), Fraction(7, 13)), (Fraction(5, 3), Fraction(4, 7)))):
    """T(n, m) の [2,1]：N の a の各行（小文字 a の冪で並べる）と (−1)^s c_s Π_□(1 − q^{a(□)} t^{l(□)+1}) の比が、全ての行で同じ単項式になるか。"""
    from rcsknots.gn import GN
    n, m = TORUS[K]
    R = load_reduced(K, 3)
    lam = (2, 1)
    N = N_lambda(R[lam], lam)
    Ns = bold_to_small(N)                                     # (e_A, e_Q, e_T)
    As = sorted({k[0] for k in Ns})
    assert all((a - As[0]) % 2 == 0 for a in As), "A の指数の偶奇がそろわない"
    results = []
    for Q, T in points:
        q, t = Q * Q, T * T
        G = GN(q, t, Fraction(1), lambda c: Fraction(c))
        cs = schur_hook_coeffs(G, lam, n, -m)                 # 向き：T(2,2p+1) は GN の傾き (2, −(2p+1))（conventions §6）
        norm = Fraction(1)
        for i, j in boxes(lam):
            norm *= 1 - q ** arm(lam, i, j) * t ** (leg(lam, i, j) + 1)
        rows = {}
        for (ea, eq, et), c in Ns.items():
            rows[ea] = rows.get(ea, Fraction(0)) + c * Q ** eq * T ** et
        ratios = []
        for idx, ea in enumerate(As):
            s = idx
            expected = (-1) ** s * cs[s] * norm if s < len(cs) else None
            ratios.append(None if not expected else rows[ea] / expected)
        results.append({"Q": str(Q), "T": str(T), "ratios": [str(r) for r in ratios], "n_rows": len(As), "n_hooks": len(cs)})
    # 比が行によらないか、点ごとに
    same = [len(set(r["ratios"])) == 1 for r in results]
    # 比が単項式 ±Q^i T^j か：2 点から i, j を解いて 3 点目で検算
    mono = None
    try:
        r0, r1, r2 = [Fraction(r["ratios"][0]) for r in results]
        import math
        (Q0, T0), (Q1, T1), (Q2, T2) = points
        # |r| = |c| Q^i T^j：log で解く（近似）→ 整数に丸めて厳密に検算
        import itertools
        found = None
        for i in range(-60, 61):
            for j in range(-60, 61):
                c0 = r0 / (Q0 ** i * T0 ** j)
                if c0 == r1 / (Q1 ** i * T1 ** j) == r2 / (Q2 ** i * T2 ** j):
                    found = (i, j, str(c0))
                    break
            if found:
                break
        mono = found
    except Exception as e:  # noqa: BLE001
        mono = f"error {e}"
    log(f"{K} [2,1]：N の a の行 {len(As)}、hook の数 {len(cs)}、比が行によらない {same}、比の単項式 (i_Q, j_T, 定数) {mono}")
    return {"knot": K, "rows_equal_ratio": same, "monomial": mono, "points": results}


def highbox(m, n, Nb):
    """T(m, n) の |λ| = Nb の全ての λ：rCS_λ を GN（bpsgen.torus_rcs_factory）から補間し、N_λ の符号を数える。"""
    import numpy as np
    from rcsknots.exact import reconstruct_points
    from rcsknots.extract import A_OFFSET, _arr
    from rcsknots.field import V
    from rcsknots.bpsgen import torus_rcs_factory
    fac = torus_rcs_factory(m, n)
    out = {}
    for lam in partitions(Nb):
        t0 = time.time()

        def ev(xs, ys, p, lam=lam):
            k = len(xs)
            Qv, Tv = V(xs, p), V(ys, p)
            lo, cs = fac(lam, Qv, Tv)
            arr = np.zeros((k, 2 * A_OFFSET + 1), dtype=np.int64)
            for i, c in enumerate(cs):
                arr[:, lo + i + A_OFFSET] = _arr(c, k, p) % p
            return arr
        raw = reconstruct_points(ev, seed=7, log=lambda s: None)
        red = {(d - A_OFFSET, eq, et): c for (d, eq, et), c in raw.items()}
        Nl = N_lambda(red, lam)
        sm, rs, se = sign_summary(Nl), sign_summary(bold(red)), series_summary(Nl, lam)
        out[str(lam)] = {"N": sm, "red_neg": rs["neg"], "red_terms": rs["terms"], "series": se}
        log(f"T({m},{n}) λ={str(lam):16s} reduced 負 {rs['neg']:4d}/{rs['terms']:5d} → N 負 {sm['neg']:4d}/{sm['terms']:5d} 冪級数の負 𝐪→0 {se['neg_to0']:4d} 鏡像 {se['neg_toinf']:4d} ({time.time() - t0:.0f}s)")
    return out


# ---- 仮説 H_u″：P 正規化（Macdonald の P_λ）の hook の Schur 係数は冪級数として非負 ---------------------

def macdonald_denominators_bold(lam):
    """Π_□(1 − q^{a(□)} t^{l(□)+1}) の各因子の単項式（太字）：q^a t^{l+1} = 𝐪^{2a+2l+2} 𝐭^{2a}。"""
    return [(2 * arm(lam, i, j) + 2 * leg(lam, i, j) + 2, 2 * arm(lam, i, j)) for i, j in boxes(lam)]


def series_signs(N, lam, order=40):
    """N_λ の a の各行を Π_□(1 − q^{a(□)} t^{l(□)+1}) で割った冪級数の符号。
    "to0"：𝐪 → 0（𝐭 固定）の向き。"toinf"：(𝐪, 𝐭) ↦ (1/𝐪, 1/𝐭)（鏡像の向き）。後者では 1/(1 − 1/m) = −m/(1 − m) と書き直す
    （全体の符号 (−1)^{因子の数}、単項式倍は符号に無関係）。返り値 {行: {"to0": (負, 正), "toinf": (負, 正)}}。order は行の 𝐪 の幅に足す次数。"""
    import flint
    ctx = flint.fmpz_mpoly_ctx.get(("q", "t"), "lex")
    dens = macdonald_denominators_bold(lam)
    rows = {}
    for (a, q, t), c in N.items():
        rows.setdefault(a, {})[(q, t)] = c

    def expand(d, flip):
        pts = {((-k, -t) if flip else (k, t)): c for (k, t), c in d.items()}
        lo = (min(k[0] for k in pts), min(k[1] for k in pts))
        P = ctx.from_dict({(k[0] - lo[0], k[1] - lo[1]): v for k, v in pts.items()})
        top = (max(k[0] for k in pts) - lo[0]) + order
        sign = (-1) ** len(dens) if flip else 1
        for kq, kt in dens:
            acc, mk, m = ctx.from_dict({}), ctx.from_dict({(0, 0): 1}), ctx.from_dict({(kq, kt): 1})
            for _ in range(top // kq + 2):
                acc += P * mk
                mk = mk * m
            P = acc
        dd = P.to_dict()
        neg = sum(1 for e, c in dd.items() if int(e[0]) <= top and sign * c < 0)
        pos = sum(1 for e, c in dd.items() if int(e[0]) <= top and sign * c > 0)
        return neg, pos
    return {a: {"to0": expand(d, False), "toinf": expand(d, True)} for a, d in sorted(rows.items())}


def series_summary(N, lam, order=40):
    ss = series_signs(N, lam, order)
    n0 = sum(v["to0"][0] for v in ss.values()); ninf = sum(v["toinf"][0] for v in ss.values())
    rows_ok = sum(1 for v in ss.values() if v["to0"][0] == 0 or v["toinf"][0] == 0)
    return {"neg_to0": n0, "neg_toinf": ninf, "rows": len(ss), "rows_ok_either": rows_ok}


def series_check(knots, order=40):
    res = {}
    for K in knots:
        try:
            R = load_reduced(K, 3)
        except AssertionError as e:
            log(f"{K}: {e}")
            continue
        res[K] = {}
        for lam in sorted(R, key=lambda l: (sum(l), l)):
            sm = series_summary(N_lambda(R[lam], lam), lam, order)
            res[K][str(lam)] = sm
            log(f"{K:6s} λ={str(lam):10s} 冪級数の負：𝐪→0 {sm['neg_to0']:4d}、鏡像 {sm['neg_toinf']:4d}（行 {sm['rows']}、どちらかで非負の行 {sm['rows_ok_either']}）")
    return res


if __name__ == "__main__":
    cmd = sys.argv[1]
    log(f"=== dgr_hu {time.strftime('%Y-%m-%d %H:%M:%S')} {' '.join(sys.argv[1:])} ===")
    if cmd == "check":
        knots = sys.argv[2:] or list(TORUS) + ALGEBRAIC + NONALG
        res = check(knots)
        prev = json.loads(OUT.read_text()) if OUT.exists() else {}
        prev.update(res)
        OUT.write_text(json.dumps(prev, ensure_ascii=False, indent=1))
    elif cmd == "convention":
        for K in sys.argv[2:]:
            convention(K)
    elif cmd == "series":
        knots = sys.argv[2:] or list(TORUS) + ALGEBRAIC + NONALG
        res = series_check(knots)
        path = Path("data/dgr_lmov/hu_series.json")
        prev = json.loads(path.read_text()) if path.exists() else {}
        prev.update(res)
        path.write_text(json.dumps(prev, ensure_ascii=False, indent=1))
    elif cmd == "highbox":
        m, n, Nb = map(int, sys.argv[2:5])
        res = highbox(m, n, Nb)
        prev = json.loads(OUT.read_text()) if OUT.exists() else {}
        prev.setdefault(f"T{m}_{n}", {}).update(res)
        OUT.write_text(json.dumps(prev, ensure_ascii=False, indent=1))
