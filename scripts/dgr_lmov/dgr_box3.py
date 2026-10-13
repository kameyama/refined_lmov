"""3 箱（docs/dgr_lmov/research_plan.md の D7・D10 の 3 箱への拡張）：\\hat f_ρ（|ρ| = 3）の公式と、rCS_{[2,1]} を未知数とする構成。

公式（scripts/framed_knot/v4_closed.py の formula の 3 箱版）：
  \\hat f_ρ = (Σ_i n_{ρ,i} R_i) / D_ρ、R = (r3, r21, r111, r2·r1, r11·r1, r1³, r1·r1(x²), r1(x³))、
  r は reduced の rCS（KN の小文字 (A, Q, T)、framing 0）、r1(x^d) は指数を d 倍したもの（Adams）。
  n, D は bps.f_formulas(3) と M_matrix(3) から導く（_formula3：unreduced の C の短い式を経由し、多項式演算は flint）。
"""
import sys
from functools import lru_cache
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "framed_knot"))
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

def t_floor(r1_small, nbox=3):
    """【仮説 B】\\hat f_ρ（|ρ| = nbox）の spin 表示で、𝐚 の行ごとの 𝐭 の最低次の下限 {e_𝐚: e_𝐭}。
    [1] の線 P(e_𝐚) = min{e_𝐭 : \\hat f_{[1]} の項} の nbox 乗（min-plus の畳み込み）を (e_𝐚, e_𝐭) → (e_𝐚 − 2, e_𝐭 + nbox − 2) と平行移動したもの。
    保存済みの 22 結び目（data/framed_knot/bps）のうち 4_1 と DAHA 入力の K5b, K5c 以外の全てで、2・3 箱の全ての ρ について等号で成り立つ（report.md）。"""
    hf1 = (L.from_dict({(1, 0, 0): 1, (-1, 0, 0): -1}) * L.from_dict(r1_small)).to_dict()
    b1 = to_bold(hf1, 0 if all(k[0] % 2 == 0 for k in hf1) else 1)
    P = {}
    for ea, _, et in b1:
        P[ea] = min(P.get(ea, 10 ** 9), et)
    out = dict(P)
    for _ in range(nbox - 1):
        out = {b: min(f1 + f2 for b1_, f1 in out.items() for b2_, f2 in P.items() if b1_ + b2_ == b)
               for b in {b1_ + b2_ for b1_ in out for b2_ in P}}
    return {b - 2: f + nbox - 2 for b, f in out.items()}


def _mul(A, B):
    out = {}
    for k1, v1 in A.items():
        for k2, v2 in B.items():
            k = tuple(k1[i] + k2[i] for i in range(len(k1)))
            out[k] = out.get(k, 0) + v1 * v2
    return {k: v for k, v in out.items() if v}


def rcs_window21(r1_bold, r2_bold):
    """【予想 R】rCS_{[2,1]}（太字）の 𝐚 の行ごとの 𝐭 の最低次と最高次 {e_𝐚: (lo, hi)}。
    lo は floor(rCS_{[2]}) と floor(rCS_{[1]}) の min-plus の畳み込み、hi は最高次の max-plus の畳み込み。
    一般に「rCS_λ の 𝐭 の最低次・最高次 = 行ごとの対称表現 [λ_i] のものの畳み込み」（Jacobi–Trudi の対角項）。
    保存済みの 13〜16 結び目（トーラス、4_1、代数的ケーブル）の [2,1] で等号で成り立つ（report.md）。σ で移した列の版は条件 T から従う。"""
    def ends(b):
        o = {}
        for ea, _, et in b:
            lo, hi = o.get(ea, (10 ** 9, -10 ** 9))
            o[ea] = (min(lo, et), max(hi, et))
        return o
    e1, e2 = ends(r1_bold), ends(r2_bold)
    out = {}
    for a1, (l1, h1) in e2.items():
        for a2, (l2, h2) in e1.items():
            lo, hi = out.get(a1 + a2, (10 ** 9, -10 ** 9))
            out[a1 + a2] = (min(lo, l1 + l2), max(hi, h1 + h2))
    return out


def build21(small, H21u, window_t, sigmas, conds=("U", "T", "H", "E", "U2", "L")):
    """small：[2,1] 以外の reduced の rCS（小文字、|λ| ≤ 3）。H21u：非 refined の [2,1]（太字 {(e_𝐚, e_𝐪): c}）。
    未知数 x_m = rCS_{[2,1]} の太字の係数（符号は自由）、h_{ρ,e} ≥ 0 = σ_ρ \\hat f_ρ の太字の係数（ρ = [3], [2,1], [1,1,1]）。"""
    from dgr_d7_construct import Lin
    from rcsknots.dgr_lmov.dgr import specialize
    bold = lambda P: to_bold(P, 0 if all(k[0] % 2 == 0 for k in P) else 1)  # noqa: E731
    r1, r11 = bold(small[(1,)]), bold(small[(1, 1)])
    lin = Lin()
    As = sorted({a for a, _ in H21u}); Qs = sorted({q for _, q in H21u})
    mons = [(a, q, t) for a in range(As[0], As[-1] + 1, 2) for q in range(Qs[0], Qs[-1] + 1, 2) for t in range(*window_t)]
    mset = set(mons)
    for m in mons:
        lin.var(("x", m))
    if "R" in conds:                                         # 𝐭 の最低次・最高次の線の外の係数は 0
        win = rcs_window21(r1, bold(small[(2,)]))
        for m in mons:
            if m[0] not in win or not win[m[0]][0] <= m[2] <= win[m[0]][1]:
                lin.eq({("x", m): 1}, 0)

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
            if "B" in conds:                                 # 下限より下の \\hat f の項は 0（h を置かない）
                floor = t_floor(small[(1,)])
                ok = lambda e: (lambda k: k[0] in floor and k[2] >= floor[k[0]])(next(iter(to_bold({e: 1}, shift))))  # noqa: E731
                hm = [e for e in hm if ok(e)]
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


# ---- 解の一意性と選択（整数計画） ----

def _matrix(lin, extra_cols=0):
    import scipy.sparse as sps
    A = sps.lil_matrix((len(lin.rows), len(lin.col) + extra_cols))
    for i, row in enumerate(lin.rows):
        for j, v in row.items():
            A[i, j] = v
    return A


def other_solution(lin, mons, sol, cvec, bound, xbound=400):
    """x ≠ sol の整数解を探す：c·x ≥ c·sol + 1 または c·x ≤ c·sol − 1（bound(c·sol) が下限・上限を返す）を足した実行可能性問題。"""
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp
    n = len(lin.col)
    lb = np.array([-xbound if nm[0] == "x" else 0 for nm in lin.col], float)
    ub = np.array([xbound if nm[0] == "x" else np.inf for nm in lin.col], float)
    cv = sum(cvec[lin.col[("x", m)]] * v for m, v in sol.items())
    rhs = np.array(lin.rhs, float)
    cons = [LinearConstraint(_matrix(lin).tocsr(), rhs, rhs), LinearConstraint(cvec.reshape(1, -1), *bound(cv))]
    return milp(np.zeros(n), constraints=cons, bounds=Bounds(lb, ub), integrality=np.ones(n), options={"time_limit": 3000})


def l1_min(lin, mons, xbound=400):
    """【仮説】rCS_{[2,1]} の係数の絶対値の和が最小の解（打ち消しの最小）。x = p − n、p, n ≥ 0 として Σ(p + n) を最小化する。"""
    import numpy as np
    import scipy.sparse as sps
    from scipy.optimize import Bounds, LinearConstraint, milp
    n = len(lin.col); xs = [lin.col[("x", m)] for m in mons]; k = len(xs)
    A = sps.vstack([_matrix(lin, 2 * k).tocsr(),
                    sps.hstack([sps.coo_matrix(([1] * k, (range(k), xs)), shape=(k, n)),
                                -sps.identity(k), sps.identity(k)])]).tocsr()
    rhs = np.array(list(lin.rhs) + [0] * k, float)
    lb = np.array([-xbound if nm[0] == "x" else 0 for nm in lin.col] + [0] * (2 * k), float)
    ub = np.array([xbound if nm[0] == "x" else np.inf for nm in lin.col] + [xbound] * (2 * k), float)
    obj = np.array([0.0] * n + [1.0] * (2 * k))
    rr = milp(obj, constraints=[LinearConstraint(A, rhs, rhs)], bounds=Bounds(lb, ub), integrality=np.ones(n + 2 * k),
              options={"time_limit": 3000})
    if rr.x is None:
        return None
    return {m: int(round(rr.x[lin.col[("x", m)]])) for m in mons if round(rr.x[lin.col[("x", m)]]) != 0}


# ---- 自由度 Φ·g を未知数にした小さい定式化 ----

def phi21():
    """【導出】rCS_{[2,1]} の 2 つの解の差 Δ を必ず割り切る多項式 Φ（太字 {(e_𝐚, e_𝐪, e_𝐭): c}）。
    U（t = q）から (1+𝐭)、H・U2 と転置から (1+u_{1|1})(1+u_{2|0})(1+u_{0|2})、
    \\hat f の Laurent 性から k_{ρ,21}·（自明な結び目の [2,1]）の分母 (q−1)(q²t−1)(qt²−1)(t−1)
    = (𝐪²𝐭²−1)(𝐪⁶𝐭⁴−1)(𝐪⁶𝐭²−1)(𝐪²−1)（EGP の因子 𝐪𝐭−1 を含む）。因子は互いに素なので積が割り切る。"""
    import flint
    ctx = flint.fmpz_mpoly_ctx.get(("a", "q", "t"), "lex")
    a, q, t = ctx.gens()
    P = (1 + t) * (q**2 * t**2 - 1) * (q**6 * t**4 - 1) * (q**6 * t**2 - 1) * (q**2 - 1) \
        * (1 + a**2 * t**3) * (q**4 + a**2 * t) * (1 + a**2 * q**4 * t**5)          # (q⁴ + 𝐚²𝐭) = 𝐪⁴(1 + u_{2|0})
    return {tuple(int(v) for v in e): int(c) for e, c in P.to_dict().items()}


def solve21_g(small, xp, qrange, sigmas, conds=("B", "R"), cvec_seed=0, objectives=("l1",), check_unique=True):
    """x = xp + Φ·g（xp：条件の等式を全て満たす rCS_{[2,1]}（太字）、例えば最初の実行可能解か正解）として g を未知数にする。
    条件：R（𝐭 の両端の線の外は 0）、B（\\hat f の 𝐭 の下限より下は 0）、転置で不変、\\hat f の spin 表示の非負性（σ_ρ を掛けて）。
    qrange：𝐪 の範囲（非 refined の [2,1] の 𝐪 の範囲。build21 の未知数の範囲と同じ）。
    objectives：【仮説】選択の規則（辞書式に最小化）。"l1"＝rCS の係数の絶対値の和、"terms"＝rCS の項数、"bps"＝\\hat f の spin 表示の係数の和。
    返り値 {"unique": 別の g があるか, "best": 選んだ x, "values": 各目的の最適値, "n_g": 変数の数}。"""
    import numpy as np
    import scipy.sparse as sps
    from scipy.optimize import Bounds, LinearConstraint, milp
    bold = lambda P: to_bold(P, 0 if all(k[0] % 2 == 0 for k in P) else 1)  # noqa: E731
    Phi = phi21()
    win = rcs_window21(bold(small[(1,)]), bold(small[(2,)]))
    box = [(min(win), max(win)), qrange, (min(lo for lo, _ in win.values()), max(hi for _, hi in win.values()))]
    ext = [(min(k[i] for k in Phi), max(k[i] for k in Phi)) for i in range(3)]
    G = [(a, q, t) for a in range(box[0][0] - ext[0][0], box[0][1] - ext[0][1] + 1, 2)
         for q in range(box[1][0] - ext[1][0], box[1][1] - ext[1][1] + 1)
         for t in range(box[2][0] - ext[2][0], box[2][1] - ext[2][1] + 1)]
    gi = {e: i for i, e in enumerate(G)}
    if not G:                                                  # Φ が窓に入らない：自由度なし
        return {"n_g": 0, "n_eq": 0, "n_ineq": 0, "feasible": True, "unique": True, "best": dict(xp)}
    # x の各係数 = xp + Σ_e g_e Φ[m − e]
    xrow = {}
    for e in G:
        for f, c in Phi.items():
            m = (e[0] + f[0], e[1] + f[1], e[2] + f[2])
            xrow.setdefault(m, {})[gi[e]] = c
    eqs, ineqs = [], []                                        # (行 {列: 係数}, 右辺)：行·g = 右辺、行·g ≥ 右辺
    if "R" in conds:
        for m, row in xrow.items():
            if m[0] not in win or not win[m[0]][0] <= m[2] <= win[m[0]][1]:
                eqs.append((row, -xp.get(m, 0)))
    for m, row in xrow.items():                               # 転置で不変
        s = (m[0], -m[1], m[2] - m[1])
        if s > m:
            r2 = dict(row)
            for j, v in xrow.get(s, {}).items():
                r2[j] = r2.get(j, 0) - v
            r2 = {j: v for j, v in r2.items() if v}
            if r2:
                eqs.append((r2, xp.get(s, 0) - xp.get(m, 0)))
    # \hat f：hf(x) = hf(xp) + Σ_e g_e · s_e·shift(Q_ρ)、Q_ρ = n_{ρ,21}·small(Φ)/D_ρ
    hp = hat3({**small, (2, 1): bold_to_small(xp)})
    phis = L.from_dict(bold_to_small(Phi))
    floor = t_floor(small[(1,)]) if "B" in conds else None
    for idx, (rho, (ns, den)) in enumerate(formula3().items()):
        Qr = (ns[1] * phis / den).to_dict()
        shift = 0 if all(k[0] % 2 == 0 for k in hp[rho]) else 1
        hb = to_bold(hp[rho], shift)
        rows = {}
        for e in G:
            (se, ce), = bold_to_small({e: 1}).items()
            for k, v in Qr.items():
                (kb, cb), = to_bold({(k[0] + se[0], k[1] + se[1], k[2] + se[2]): v * ce}, shift).items()
                rows.setdefault(kb, {})[gi[e]] = rows.get(kb, {}).get(gi[e], 0) + cb
        for kb in set(rows) | set(hb):
            row = {j: v for j, v in rows.get(kb, {}).items() if v}
            if floor is not None and (kb[0] not in floor or kb[2] < floor[kb[0]]):
                eqs.append((row, -hb.get(kb, 0)))
            elif row:
                ineqs.append(({j: sigmas[idx] * v for j, v in row.items()}, -sigmas[idx] * hb.get(kb, 0)))
    n = len(G)

    def mat(rows, ncol):
        A = sps.lil_matrix((len(rows), ncol))
        for i, (row, _) in enumerate(rows):
            for j, v in row.items():
                A[i, j] = v
        return A.tocsr()
    rhs_e = np.array([r for _, r in eqs], float); rhs_i = np.array([r for _, r in ineqs], float)

    def feasible(extra=None, obj=None, ncol=n, extra_cons=()):
        cons = [LinearConstraint(mat(eqs, ncol), rhs_e, rhs_e), LinearConstraint(mat(ineqs, ncol), rhs_i, np.inf), *extra_cons]
        if extra is not None:
            cons.append(extra)
        lb = np.full(ncol, -np.inf); ub = np.full(ncol, np.inf)
        return milp(np.zeros(ncol) if obj is None else obj, constraints=cons, bounds=Bounds(lb, ub), integrality=np.ones(ncol),
                    options={"time_limit": 3000})
    out = {"n_g": n, "n_eq": len(eqs), "n_ineq": len(ineqs)}
    base = feasible()
    out["feasible"] = base.x is not None
    if base.x is None:
        return out
    if check_unique:
        rng = np.random.default_rng(cvec_seed)
        c = rng.integers(1, 10 ** 6, size=n).astype(float)
        g0 = np.round(base.x); v0 = c @ g0
        st = [feasible(LinearConstraint(c.reshape(1, -1), lo, hi)).status for lo, hi in ((v0 + 1, np.inf), (-np.inf, v0 - 1))]
        out["unique"] = all(s == 2 for s in st)
    if objectives:
        # 変数 [g (n), p (k), m (k), z (k)]：x_m = xp_m + Σ g Φ = p − m、|x_m| ≤ M z_m（z は項の有無、0/1）
        ms = sorted(set(xrow) | set(xp))
        k = len(ms)
        N = n + 3 * k
        X = sps.lil_matrix((k, N))
        Z = sps.lil_matrix((k, N))
        for i, m in enumerate(ms):
            for j, v in xrow.get(m, {}).items():
                X[i, j] = v
            X[i, n + i] = -1; X[i, n + k + i] = 1
            Z[i, n + i] = 1; Z[i, n + k + i] = 1; Z[i, n + 2 * k + i] = -400
        rhs_x = np.array([-xp.get(m, 0) for m in ms], float)
        cons = [LinearConstraint(X.tocsr(), rhs_x, rhs_x), LinearConstraint(Z.tocsr(), -np.inf, 0),
                LinearConstraint(sps.hstack([mat(eqs, n), sps.csr_matrix((len(eqs), 3 * k))]).tocsr(), rhs_e, rhs_e),
                LinearConstraint(sps.hstack([mat(ineqs, n), sps.csr_matrix((len(ineqs), 3 * k))]).tocsr(), rhs_i, np.inf)]
        lb = np.array([-np.inf] * n + [0.0] * (3 * k)); ub = np.array([np.inf] * (n + 2 * k) + [1.0] * k)
        objs = {"l1": np.array([0.0] * n + [1.0] * (2 * k) + [0.0] * k),
                "terms": np.array([0.0] * (n + 2 * k) + [1.0] * k),
                "bps": np.concatenate([np.asarray(mat(ineqs, n).sum(axis=0)).ravel(), np.zeros(3 * k)])}
        out["values"] = {}
        for name in objectives:                               # 辞書式：前の目的の最適値を制約にして次を最小化
            rr = milp(objs[name], constraints=cons, bounds=Bounds(lb, ub), integrality=np.ones(N), options={"time_limit": 3000})
            if rr.x is None:
                out["status"] = f"{name}: {rr.message}"
                break
            out["values"][name] = round(float(objs[name] @ rr.x))
            cons.append(LinearConstraint(objs[name].reshape(1, -1), -np.inf, out["values"][name] + 1e-6))
        if rr.x is not None:
            g = np.round(rr.x[:n])
            x = dict(xp)
            for e, ge in zip(G, g):
                if ge:
                    for f, c in Phi.items():
                        m = (e[0] + f[0], e[1] + f[1], e[2] + f[2])
                        x[m] = x.get(m, 0) + int(ge) * c
            out["best"] = {m: v for m, v in x.items() if v}
    return out
