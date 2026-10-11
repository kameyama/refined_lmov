"""3 箱（docs/dgr_lmov/research_plan.md の D7・D10 の 3 箱への拡張）：\\hat f_ρ（|ρ| = 3）の公式と、rCS_{[2,1]} を未知数とする構成。

公式（scripts/v4_closed.py の formula の 3 箱版）：
  \\hat f_ρ = (Σ_i n_{ρ,i} R_i) / D_ρ、R = (r3, r21, r111, r2·r1, r11·r1, r1³, r1·r1(x²), r1(x³))、
  r は reduced の rCS（KN の小文字 (A, Q, T)、framing 0）、r1(x^d) は指数を d 倍したもの（Adams）。
  n, D は bps.f_formulas(3) と M_matrix(3) から導く（_formula3：unreduced の C の短い式を経由し、多項式演算は flint）。
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


@lru_cache(maxsize=None)
def formula3():
    return _formula3()


def _to_flint(e, ctx, gen):
    """sympy の有理式 → flint の (分子, 分母)。足し算・掛け算のたびに gcd で約分する（sympy の展開・約分を使わない）。"""
    one = ctx.from_dict({(0,) * ctx.nvars(): 1})

    def red(n, d):
        g = n.gcd(d)
        n, d = n / g, d / g
        return (-n, -d) if d.leading_coefficient() < 0 else (n, d)

    def conv(e):
        if e.is_Symbol:
            return gen[e], one
        if e.is_Rational:
            return one * int(e.p), one * int(e.q)
        if e.is_Add or e.is_Mul:
            n, d = conv(e.args[0])
            for a in e.args[1:]:
                n2, d2 = conv(a)
                if e.is_Mul:
                    n, d = red(n * n2, d * d2)
                else:
                    g = d.gcd(d2)
                    n, d = red(n * (d2 / g) + n2 * (d / g), d * (d2 / g))
            return n, d
        if e.is_Pow and e.exp.is_Integer:
            n, d = conv(e.base)
            k = int(e.exp)
            return (n ** k, d ** k) if k >= 0 else red(d ** -k, n ** -k)
        raise TypeError(e)
    return conv(e), red


def _formula3():
    """\\hat f_ρ = Σ_i k_{ρ,i}(Q, T)·（unreduced の C の単項式 i）（ユーザーのノートブック new derivation of link rf_qt と同じ形）。
    k は f^q = M(t) \\hat f を M(t)⁻¹ で解いたものの係数で、(1 − q^a t^b) 型の因子の短い比になる。
    C = （自明な結び目の値）× r を代入して共通分母にそろえる。多項式演算は flint（約 1 秒）。"""
    import flint
    Asym = sp.Symbol("A")

    def unknot(lam, d=1):
        k = sum(lam)
        out = T ** (d * k) * Asym ** (-d * k)
        for i, j in boxes(lam):
            out *= (T ** (2 * d * i) - Asym ** (2 * d) * Q ** (2 * d * j)) / (1 - Q ** (2 * d * arm(lam, i, j)) * T ** (2 * d * (leg(lam, i, j) + 1)))
        return out

    syms = [C((3,)), C((2, 1)), C((1, 1, 1)), C((2,)), C((1, 1)), C((1,)), C1(2), C1(3)]
    names = [Asym, Q, T] + syms
    ctx = flint.fmpz_mpoly_ctx.get(tuple(map(str, names)), "lex")
    gen = dict(zip(names, ctx.gens()))
    conv = lambda e: _to_flint(e, ctx, gen)[0]  # noqa: E731
    red = _to_flint(sp.Integer(1), ctx, gen)[1]
    mon_key = {(1, 0, 0, 0, 0, 0, 0, 0): 0, (0, 1, 0, 0, 0, 0, 0, 0): 1, (0, 0, 1, 0, 0, 0, 0, 0): 2,
               (0, 0, 0, 1, 0, 1, 0, 0): 3, (0, 0, 0, 0, 1, 1, 0, 0): 4, (0, 0, 0, 0, 0, 3, 0, 0): 5,
               (0, 0, 0, 0, 0, 1, 1, 0): 6, (0, 0, 0, 0, 0, 0, 0, 1): 7}
    unk = [conv(u) for u in [unknot((3,)), unknot((2, 1)), unknot((1, 1, 1)), unknot((2,)) * unknot((1,)),
                             unknot((1, 1)) * unknot((1,)), unknot((1,)) ** 3, unknot((1,)) * unknot((1,), 2), unknot((1,), 3)]]
    F = f_formulas(3)
    parts3 = partitions(3)
    fq = [conv(F[mu][0]) for mu in parts3]
    Minv = M_matrix(3).inv()

    def to_L(x):
        return L.from_dict({tuple(int(v) for v in e[:3]): int(c) for e, c in x.to_dict().items()})
    out = {}
    for i, rho in enumerate(parts3):
        n, d = ctx.from_dict({}), ctx.from_dict({(0,) * len(names): 1})
        for j in range(3):                                   # \\hat f_ρ = Σ_j (M⁻¹)_{ρj} f^q_j
            mn, md = conv(sp.together(Minv[i, j]))
            a, b = red(mn * fq[j][0], md * fq[j][1])
            g = d.gcd(b)
            n, d = red(n * (b / g) + a * (d / g), d * (b / g))
        parts = {}
        for e, c in n.to_dict().items():                     # C の単項式ごとの係数 k_{ρ,i} = parts[i] / d
            parts.setdefault(mon_key[tuple(int(v) for v in e[3:])], {})[tuple(e[:3]) + (0,) * len(syms)] = c
        coeffs = {k: red(ctx.from_dict(v) * unk[k][0], d * unk[k][1]) for k, v in parts.items()}
        D = ctx.from_dict({(0,) * len(names): 1})
        for _, dk in coeffs.values():
            D = D * (dk / D.gcd(dk))
        out[rho] = ({k: to_L(nk * (D / dk)) for k, (nk, dk) in coeffs.items()}, to_L(D))
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
