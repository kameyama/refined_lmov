"""3 箱（docs/dgr_lmov/research_plan.md の D7・D10 の 3 箱への拡張）：\\hat f_ρ（|ρ| = 3）の公式と、rCS_{[2,1]} を未知数とする構成。

公式（scripts/v4_closed.py の formula の 3 箱版）：
  \\hat f_ρ = (Σ_i n_{ρ,i} R_i) / D_ρ、R = (r3, r21, r111, r2·r1, r11·r1, r1³, r1·r1(x²), r1(x³))、
  r は reduced の rCS（KN の小文字 (A, Q, T)、framing 0）、r1(x^d) は指数を d 倍したもの（Adams）。
  n, D は bps.f_formulas(3) と M_matrix(3) から sympy で一度だけ導く（各項の係数を個別に約分してから共通分母にそろえる）。
"""
import sys
from functools import lru_cache
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from t4_colored import bold_to_small  # noqa: E402
from v4_closed import L  # noqa: E402
from rcsknots.bps import C, C1, M_matrix, Q, T, f_formulas, to_bold  # noqa: E402
from rcsknots.partitions import arm, boxes, leg, partitions  # noqa: E402

RHOS3 = [(3,), (2, 1), (1, 1, 1)]
KEYS = ["r3", "r21", "r111", "r2r1", "r11r1", "r1^3", "r1r1d2", "r1d3"]


CACHE = Path(__file__).resolve().parents[1] / "data/dgr_lmov/formula3.json"   # git で管理する（31 KB）


@lru_cache(maxsize=None)
def formula3():
    """導出に約 15 分かかるので data/dgr_lmov/formula3.json に保存して再利用する。"""
    import json
    if CACHE.exists():
        raw = json.loads(CACHE.read_text())
        dec = lambda d: L.from_dict({tuple(map(int, k.split(","))): v for k, v in d.items()})  # noqa: E731
        return {tuple(map(int, r.split(","))): ({int(k): dec(v) for k, v in ns.items()}, dec(den)) for r, (ns, den) in raw.items()}
    out = _formula3()
    enc = lambda l_: {",".join(map(str, k)): v for k, v in l_.to_dict().items()}  # noqa: E731
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps({",".join(map(str, r)): ({str(k): enc(v) for k, v in ns.items()}, enc(den)) for r, (ns, den) in out.items()}))
    return out


def _formula3():
    Asym = sp.Symbol("A")

    def unknot(lam, d=1):
        k = sum(lam)
        out = T ** (d * k) * Asym ** (-d * k)
        for i, j in boxes(lam):
            out *= (T ** (2 * d * i) - Asym ** (2 * d) * Q ** (2 * d * j)) / (1 - Q ** (2 * d * arm(lam, i, j)) * T ** (2 * d * (leg(lam, i, j) + 1)))
        return out

    F = f_formulas(3)
    parts3 = partitions(3)
    fq = sp.Matrix([F[mu][0] for mu in parts3])
    hat = M_matrix(3).inv() * fq
    syms = [C((3,)), C((2, 1)), C((1, 1, 1)), C((2,)), C((1, 1)), C((1,)), C1(2), C1(3)]
    mon_key = {(1, 0, 0, 0, 0, 0, 0, 0): 0, (0, 1, 0, 0, 0, 0, 0, 0): 1, (0, 0, 1, 0, 0, 0, 0, 0): 2,
               (0, 0, 0, 1, 0, 1, 0, 0): 3, (0, 0, 0, 0, 1, 1, 0, 0): 4, (0, 0, 0, 0, 0, 3, 0, 0): 5,
               (0, 0, 0, 0, 0, 1, 1, 0): 6, (0, 0, 0, 0, 0, 0, 0, 1): 7}
    unk = [unknot((3,)), unknot((2, 1)), unknot((1, 1, 1)), unknot((2,)) * unknot((1,)), unknot((1, 1)) * unknot((1,)),
           unknot((1,)) ** 3, unknot((1,)) * unknot((1,), 2), unknot((1,), 3)]
    out = {}
    for rho, h in zip(parts3, hat):
        h = sp.together(h)
        P = sp.Poly(sp.numer(h), *syms)
        den0 = sp.denom(h)
        parts = {}
        for mon, c in P.terms():
            key = mon_key[mon]
            parts[key] = sp.cancel(sp.together(c * unk[key] / den0))
        dens = [sp.denom(parts[k]) for k in parts]
        D = sp.lcm_list(dens) if hasattr(sp, "lcm_list") else sp.lcm(dens)
        coeffs = {k: sp.cancel(parts[k] * D) for k in parts}

        def to_L(e):
            P_ = sp.Poly(sp.expand(e * (Asym * Q * T) ** 80), Asym, Q, T)
            return L.from_dict({(i - 80, j - 80, k - 80): int(c) for (i, j, k), c in P_.terms()})
        out[rho] = ({k: to_L(c) for k, c in coeffs.items()}, to_L(D))
    return out


def inputs_R(small):
    """small = {λ: 小文字の reduced の rCS（{(e_A, e_Q, e_T): c}）}（|λ| ≤ 3）→ R の各項（L）。r21 が無ければ None。"""
    Lp = {lam: L.from_dict(p) for lam, p in small.items()}
    r1 = Lp[(1,)]
    adams = lambda d: L.from_dict({(d * i, d * j, d * k): c for (i, j, k), c in small[(1,)].items()})  # noqa: E731
    return {"r3": Lp[(3,)], "r21": Lp.get((2, 1)), "r111": Lp[(1, 1, 1)], "r2r1": Lp[(2,)] * r1, "r11r1": Lp[(1, 1)] * r1,
            "r1^3": r1 * r1 * r1, "r1r1d2": r1 * adams(2), "r1d3": adams(3)}


def hat3(small):
    """全ての rCS が既知のときの \\hat f_ρ（|ρ| = 3、小文字、a 次数のずらしの前）。"""
    R = inputs_R(small)
    out = {}
    for rho, (ns, den) in formula3().items():
        num = None
        for k, n in ns.items():
            term = n * R[KEYS[k]]
            num = term if num is None else num + term
        out[rho] = (num / den).to_dict()
    return out


# ---- rCS_{[2,1]} を未知数とする構成（整数計画） ----

def _mul(A, B):
    out = {}
    for k1, v1 in A.items():
        for k2, v2 in B.items():
            k = tuple(k1[i] + k2[i] for i in range(len(k1)))
            out[k] = out.get(k, 0) + v1 * v2
    return {k: v for k, v in out.items() if v}


def build21(small, H21u, window_t, sigmas, conds=("U", "T", "H", "E", "U2", "L")):
    """small：[2,1] 以外の reduced の rCS（小文字、|λ| ≤ 3）。H21u：非 refined の [2,1]（太字 {(e_𝐚, e_𝐪): c}）。
    未知数 x_m = rCS_{[2,1]} の太字の係数（符号は自由）、h_{ρ,e} ≥ 0 = σ_ρ \\hat f_ρ の太字の係数（ρ = [3], [2,1], [1,1,1]）。"""
    from dgr_d7_construct import Lin
    from rcsknots.dgr import specialize
    bold = lambda P: to_bold(P, 0 if all(k[0] % 2 == 0 for k in P) else 1)  # noqa: E731
    r1, r11 = bold(small[(1,)]), bold(small[(1, 1)])
    lin = Lin()
    As = sorted({a for a, _ in H21u}); Qs = sorted({q for _, q in H21u})
    mons = [(a, q, t) for a in range(As[0], As[-1] + 1, 2) for q in range(Qs[0], Qs[-1] + 1, 2) for t in range(*window_t)]
    mset = set(mons)
    for m in mons:
        lin.var(("x", m))

    def special(n, m, target):
        img = {}
        for mm in mons:
            for k, v in specialize({mm: 1}, n, m).items():
                img.setdefault(k, {})[("x", mm)] = v
        for k in set(img) | set(target):
            lin.eq(img.get(k, {}), target.get(k, 0))
    if "U" in conds:
        for a in range(As[0], As[-1] + 1, 2):
            for q in range(Qs[0], Qs[-1] + 1, 2):
                lin.eq({("x", (a, q, t)): (-1) ** (t % 2) for t in range(*window_t)}, H21u.get((a, q), 0))
    if "T" in conds:
        for m in mons:
            s = (m[0], -m[1], m[2] - m[1])
            if s not in mset:
                lin.eq({("x", m): 1}, 0)
            elif s > m:
                lin.eq({("x", m): 1, ("x", s): -1}, 0)
    if "H" in conds:
        (cq, ct), cs = next(iter(specialize(r1, 1, 0).items()))
        A = specialize(r1, 1, 1)
        special(1, 1, {(3 * eq, et + eq + 2 * ct): cs * cs * v for (eq, et), v in A.items()})     # c²·A(𝐪³𝐭, 𝐭)
    if "E" in conds:
        prod = _mul(r1, r11)
        target = {}
        for (a, q, t), v in prod.items():
            target[(a, t - q)] = target.get((a, t - q), 0) + v
        img = {}
        for m in mons:
            img.setdefault((m[0], m[2] - m[1]), {})[("x", m)] = 1
        for k in set(img) | set(target):
            lin.eq(img.get(k, {}), target.get(k, 0))
    if "U2" in conds:
        special(2, 0, _mul(specialize(r1, 2, 0), specialize(r11, 2, 0)))
    if "L" in conds:
        R = inputs_R({**small, (2, 1): {(0, 0, 0): 0}})
        for idx, (rho, (ns, den)) in enumerate(formula3().items()):
            const = None
            for k, n in ns.items():
                if KEYS[k] == "r21":
                    continue
                term = n * R[KEYS[k]]
                const = term if const is None else const + term
            const = const.to_dict()
            lin_part = {}
            for m in mons:
                v = (ns[1] * L.from_dict(bold_to_small({m: 1}))).to_dict()
                for k, w in v.items():
                    lin_part.setdefault(k, {})[("x", m)] = w
            dd = den.to_dict()
            keys = set(lin_part) | set(const)
            lo = [min(k[i] for k in keys) - min(k[i] for k in dd) for i in range(3)]
            hi = [max(k[i] for k in keys) - max(k[i] for k in dd) for i in range(3)]
            par = {k[0] % 2 for k in keys}
            assert len(par) == 1 and len({k[0] % 2 for k in dd}) == 1, "a の指数の偶奇がそろわない"
            hpar = (par.pop() - next(iter(dd))[0]) % 2       # \hat f の A の偶奇 = 分子 − 分母（分母は A³ を含み奇数）
            shift = hpar                                     # 太字にするための a 次数のずらし（A の単位）
            hm = [(a, q, t) for a in range(lo[0], hi[0] + 1) for q in range(lo[1], hi[1] + 1) for t in range(lo[2], hi[2] + 1)
                  if a % 2 == hpar]
            sgn = {e: next(iter(to_bold({e: 1}, shift).values())) for e in hm}
            conv = {}
            for e in hm:
                for k, w in dd.items():
                    key = (e[0] + k[0], e[1] + k[1], e[2] + k[2])
                    conv.setdefault(key, {})[("h", rho, e)] = w * sgn[e] * sigmas[idx]
            for key in set(conv) | keys:
                row = dict(conv.get(key, {}))
                for n_, w in lin_part.get(key, {}).items():
                    row[n_] = row.get(n_, 0) - w
                lin.eq(row, const.get(key, 0))
    return lin, mons
