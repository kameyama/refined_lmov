"""N 箱の \\hat f_ρ の公式（scripts/dgr_box3._formula3 の一般化）：
  \\hat f_ρ = (Σ_i n_{ρ,i} R_i) / D_ρ、R_i = reduced の rCS の単項式（Adams 像 r_λ(x^d) を含む）。
導出は bps.f_formulas(N)（plethystic log、記号 C_λ, C1_d, C_λ_d）と M_matrix(N)⁻¹ から、unreduced の C の短い式を経由し、
C_λ(x^d) = unknot_λ(q^d, t^d) · r_λ(x^d) を代入して共通の分母にそろえる（多項式演算は flint）。

使い方（検査）：PYTHONPATH=src:scripts python scripts/lmov4_formula.py check3    3 箱で dgr_box3.formula3 と一致するか
              PYTHONPATH=src:scripts python scripts/lmov4_formula.py check4    T(2,3) の 4 箱で data/highbox の \\hat f と一致するか
"""
import re
import sys
from functools import lru_cache
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dgr_box3 import _to_flint  # noqa: E402
from v4_closed import L  # noqa: E402
from rcsknots.bps import M_matrix, Q, T, f_formulas  # noqa: E402
from rcsknots.partitions import arm, boxes, leg, partitions  # noqa: E402


def parse_symbol(name):
    """'C_21' → ((2,1), 1)、'C1_3' → ((1,), 3)、'C_11_2' → ((1,1), 2)。"""
    m = re.fullmatch(r"C1_(\d+)", name)
    if m:
        return (1,), int(m.group(1))
    m = re.fullmatch(r"C_(\d+)(?:_(\d+))?", name)
    lam = tuple(int(ch) for ch in m.group(1))
    return lam, int(m.group(2) or 1)


@lru_cache(maxsize=None)
def formula_n(N):
    """{ρ: ({R の単項式（記号ごとの指数の組）: 分子 n (L)}, 分母 D (L))} と、記号の並び [(λ, d), ...]。"""
    import flint
    Asym = sp.Symbol("A")

    def unknot(lam, d=1):
        k = sum(lam)
        out = T ** (d * k) * Asym ** (-d * k)
        for i, j in boxes(lam):
            out *= (T ** (2 * d * i) - Asym ** (2 * d) * Q ** (2 * d * j)) / (1 - Q ** (2 * d * arm(lam, i, j)) * T ** (2 * d * (leg(lam, i, j) + 1)))
        return out

    F = f_formulas(N)
    parts = partitions(N)
    syms = sorted({s for mu in parts for s in F[mu][0].free_symbols if s.name not in ("Q", "T")}, key=lambda s: s.name)
    keys = [parse_symbol(s.name) for s in syms]
    names = [Asym, Q, T] + syms
    ctx = flint.fmpz_mpoly_ctx.get(tuple(map(str, names)), "lex")
    gen = dict(zip(names, ctx.gens()))
    conv = lambda e: _to_flint(e, ctx, gen)[0]  # noqa: E731
    red = _to_flint(sp.Integer(1), ctx, gen)[1]
    unk = {}                                                   # 記号の指数の組 → unknot の積（分子, 分母）
    fq = [conv(F[mu][0]) for mu in parts]
    Minv = M_matrix(N).inv()

    def to_L(x):
        return L.from_dict({tuple(int(v) for v in e[:3]): int(c) for e, c in x.to_dict().items()})

    def unk_of(expo):
        if expo not in unk:
            u = sp.Integer(1)
            for (lam, d), k in zip(keys, expo):
                for _ in range(k):
                    u *= unknot(lam, d)
            unk[expo] = conv(u)
        return unk[expo]
    out = {}
    for i, rho in enumerate(parts):
        n, d = ctx.from_dict({}), ctx.from_dict({(0,) * len(names): 1})
        for j in range(len(parts)):                           # \\hat f_ρ = Σ_j (M⁻¹)_{ρj} f^q_j
            if Minv[i, j] == 0:
                continue
            mn, md = conv(sp.together(Minv[i, j]))
            a, b = red(mn * fq[j][0], md * fq[j][1])
            g = d.gcd(b)
            n, d = red(n * (b / g) + a * (d / g), d * (b / g))
        partsd = {}
        for e, c in n.to_dict().items():                     # 記号の単項式ごとの係数
            expo = tuple(int(v) for v in e[3:])
            partsd.setdefault(expo, {})[tuple(e[:3]) + (0,) * len(syms)] = c
        coeffs = {}
        for expo, v in partsd.items():
            un, ud = unk_of(expo)
            coeffs[expo] = red(ctx.from_dict(v) * un, d * ud)
        D = ctx.from_dict({(0,) * len(names): 1})
        for _, dk in coeffs.values():
            D = D * (dk / D.gcd(dk))
        out[rho] = ({expo: to_L(nk * (D / dk)) for expo, (nk, dk) in coeffs.items()}, to_L(D))
    return out, keys


def inputs_n(small, N):
    """small = {λ: reduced の rCS（小文字 {(e_A, e_Q, e_T): c}）}（|λ| ≤ N）→ 記号 (λ, d) → L（Adams 像は指数を d 倍）。"""
    _, keys = formula_n(N)
    out = {}
    for lam, d in keys:
        P = small.get(lam)
        if P is None:
            out[(lam, d)] = None
            continue
        out[(lam, d)] = L.from_dict({(d * i, d * j, d * k): c for (i, j, k), c in P.items()})
    return out


def hat_n(small, N):
    """全ての rCS（|λ| ≤ N）が既知のときの \\hat f_ρ（|ρ| = N、小文字、a 次数のずらしの前）。"""
    fo, keys = formula_n(N)
    R = inputs_n(small, N)
    out = {}
    for rho, (ns, den) in fo.items():
        num = None
        for expo, n in ns.items():
            term = n
            for (lam, d), k in zip(keys, expo):
                for _ in range(k):
                    term = term * R[(lam, d)]
            num = term if num is None else num + term
        out[rho] = (num / den).to_dict()
    return out


if __name__ == "__main__":
    import time
    cmd = sys.argv[1]
    if cmd == "check3":
        from dgr_box3 import KEYS, formula3
        t0 = time.time()
        fo, keys = formula_n(3)
        print("3 箱の記号", keys, f"({time.time() - t0:.1f}s)")
        f3 = formula3()
        # 単項式の対応：KEYS = [r3, r21, r111, r2r1, r11r1, r1^3, r1r1d2, r1d3]
        name_of = {}
        for expo in next(iter(fo.values()))[0]:
            parts_ = [f"{''.join(map(str, lam))}" + (f"d{d}" if d > 1 else "") for (lam, d), k in zip(keys, expo) for _ in range(k)]
            name_of[expo] = parts_
        ok = True
        for rho in fo:
            ns, D = fo[rho]
            ns3, D3 = f3[rho]
            for expo, n in ns.items():
                tag = sorted(name_of[expo])
                k3 = {"r3": ["3"], "r21": ["21"], "r111": ["111"], "r2r1": ["1", "2"], "r11r1": ["1", "11"], "r1^3": ["1", "1", "1"],
                      "r1r1d2": ["1", "1d2"], "r1d3": ["1d3"]}
                idx = next(i for i, kname in enumerate(KEYS) if sorted(k3[kname]) == tag)
                if not (n * D3 == ns3[idx] * D):
                    ok = False
                    print("不一致", rho, tag)
        print("formula3 と一致:", ok)
    elif cmd == "check4":
        import json
        import lzma
        import numpy as np
        from rcsknots.exact import reconstruct_points
        from rcsknots.extract import A_OFFSET, _arr
        from rcsknots.field import V
        from rcsknots.bpsgen import torus_rcs_factory
        from framing_check import load_reduced
        t0 = time.time()
        fo, keys = formula_n(4)
        print("4 箱の記号", keys, f"導出 {time.time() - t0:.1f}s、分母", {str(r): len(D.to_dict()) for r, (_, D) in fo.items()})
        small = load_reduced("T2_3", 3)
        fac = torus_rcs_factory(2, 3)
        for lam in partitions(4):
            def ev(xs, ys, p, lam=lam):
                k = len(xs)
                Qv, Tv = V(xs, p), V(ys, p)
                lo, cs = fac(lam, Qv, Tv)
                arr = np.zeros((k, 2 * A_OFFSET + 1), dtype=np.int64)
                for i, c in enumerate(cs):
                    arr[:, lo + i + A_OFFSET] = _arr(c, k, p) % p
                return arr
            raw = reconstruct_points(ev, seed=7, log=lambda s: None)
            small[lam] = {(d - A_OFFSET, eq, et): c for (d, eq, et), c in raw.items()}
        hats = hat_n(small, 4)
        ref = json.loads(lzma.open("/Users/masaya/projects/ai_research/refined_cs/data/highbox/T2_3_N4.json.xz").read())
        for rho, h in hats.items():
            key = "".join(map(str, rho))
            stored = {tuple(map(int, k.split(","))): v for k, v in ref["rho"][key]["hat_f"].items()}
            print(f"ρ={rho}: 項数 {len(h)}（保存 {len(stored)}）、一致 {h == stored}")
        print(f"({time.time() - t0:.0f}s)")
