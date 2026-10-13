"""refined LMOV 不変量の既約表示・種数表示と、DGR 型の構造を調べる道具（docs/dgr_lmov/report.md §11 の辞書）。

多項式は {(e_𝐚, e_𝐪, e_𝐭): c}（太字、DGR の変数）か {(e_A, e_Q, e_T): c}（KN の小文字の平方根 (A, Q, T) = (a, q, t)^{1/2}）。
行の座標 (β, k, f)：𝐚^{2β} w^k 𝐭^f、w = 𝐪²𝐭（SU(2)_L の重みの変数）、𝐭（U(1)_R の電荷の変数）。
  e_𝐚 = 2β, e_𝐪 = 2k, e_𝐭 = f + k。
W：(a, q, t) ↦ (a, 1/t, 1/q) ＝ 太字で (𝐚, w, 𝐭) ↦ (𝐚, 1/w, 𝐭) ＝ DGR の対称性 (e_𝐚, e_𝐪, e_𝐭) ↦ (e_𝐚, −e_𝐪, e_𝐭 − e_𝐪)。
既約表示：\\hat f = Σ N^{(s)}_{β,j,f} 𝐚^{2β} χ_j(w) 𝐭^f、χ_j(w) = w^{2j} + w^{2j−2} + … + w^{−2j}（j は半整数、2j = k ≥ 0）。
  各 (β, f) 行が w ↔ 1/w 対称なら N^{(s)}_{β,j,f} = c_{2j} − c_{2j+2}。非負 ⇔ 行が偶奇ごとに |k| について単峰。
種数表示：KN の G 展開（rcsknots.bps.g_expand）。G = w + 1/w − 𝐭 − 1/𝐭。
"""
import json
import lzma
from collections import defaultdict
from pathlib import Path

from ..bps import g_expand, to_bold


# ---- 入出力 ----------------------------------------------------------------------------------------

def _parse(poly):
    return {tuple(map(int, k.split(","))): v for k, v in poly.items()}


def load_bps(path):
    """data/framed_knot/bps/{knot}_{rho}.json(.xz) → (meta, 小文字の \\hat f, 太字の \\hat f（a 次数のずらしと全体の符号を適用済み）)。"""
    path = Path(path)
    opener = lzma.open if path.suffix == ".xz" else open
    with opener(path, "rt") as fp:
        d = json.load(fp)
    lower = _parse(d["hat_f"])
    bold = bold_from_lower(lower, d["a_shift_A_units"], d.get("overall_sign_applied", 1))
    meta = {k: d[k] for k in ("knot", "rho", "a_shift_A_units", "overall_sign_applied") if k in d}
    return meta, lower, bold


def load_highbox(path):
    """data/framed_knot/highbox/T{m}_{n}_N{N}.json.xz → {ρ（文字列）: (小文字の \\hat f, 太字の \\hat f)}。全体の符号は判定の sign。"""
    with lzma.open(path, "rt") as fp:
        d = json.load(fp)
    out = {}
    for rho, v in d["rho"].items():
        lower = _parse(v["hat_f"])
        out[rho] = (lower, bold_from_lower(lower, d["a_shift_A_units"], v["judge"]["sign"]))
    return d, out


def bold_from_lower(lower, a_shift, sign=1):
    b = to_bold(lower, a_shift)
    return {k: sign * v for k, v in b.items()}


# ---- 行・対称性 ------------------------------------------------------------------------------------

def rows(bold):
    """{(β, f): {k: c}}。𝐪 の指数が奇数なら ValueError。"""
    out = defaultdict(dict)
    for (ea, eq, et), c in bold.items():
        if eq % 2 or ea % 2:
            raise ValueError(f"odd exponent in {(ea, eq, et)}")
        k = eq // 2
        out[(ea // 2, et - k)][k] = c
    return dict(out)


def from_rows(rws):
    return {(2 * b, 2 * k, f + k): c for (b, f), r in rws.items() for k, c in r.items() if c}


def W(bold):
    """DGR の対称性 (e_𝐚, e_𝐪, e_𝐭) ↦ (e_𝐚, −e_𝐪, e_𝐭 − e_𝐪)（= 小文字の (a, q, t) ↦ (a, 1/t, 1/q)）。"""
    return {(ea, -eq, et - eq): c for (ea, eq, et), c in bold.items()}


def is_W_invariant(bold):
    return W(bold) == bold


def W_lower(lower):
    """小文字 (A, Q, T) での W：(e_A, e_Q, e_T) ↦ (e_A, −e_T, −e_Q)。"""
    return {(a, -t, -q): c for (a, q, t), c in lower.items()}


# ---- 既約表示 ------------------------------------------------------------------------------------

def spin(bold):
    """{(β, 2j, f): N^{(s)}}。行が w ↔ 1/w 対称でなければ ValueError。"""
    out = {}
    for (b, f), r in rows(bold).items():
        if any(r.get(-k) != c for k, c in r.items()):
            raise ValueError(f"row (beta, f) = {(b, f)} is not symmetric under w <-> 1/w")
        K = max(abs(k) for k in r)
        for k in range(0, K + 1):
            v = r.get(k, 0) - r.get(k + 2, 0)
            if v:
                out[(b, k, f)] = v
    return out


def from_spin(sp):
    """既約表示 → 太字。χ_j(w) = Σ_{k ≡ 2j (2), |k| ≤ 2j} w^k。"""
    rws = defaultdict(lambda: defaultdict(int))
    for (b, j2, f), n in sp.items():
        for k in range(-j2, j2 + 1, 2):
            rws[(b, f)][k] += n
    return from_rows({key: {k: c for k, c in r.items() if c} for key, r in rws.items()})


def spin_count(bold):
    """(既約表示の項数, 負の項数)。"""
    s = spin(bold)
    return len(s), sum(1 for v in s.values() if v < 0)


def unimodal_rows(bold):
    """各 (β, f) 行が、偶奇ごとに |k| について単峰（非負・非増加）かどうか。既約表示の非負性と同値。"""
    try:
        return all(v >= 0 for v in spin(bold).values())
    except ValueError:
        return False


# ---- 次数 ------------------------------------------------------------------------------------------

def delta_values(bold):
    """δ = f − 2β（= e_𝐭 − e_𝐪/2 − e_𝐚）の値の集合（昇順）。"""
    return sorted({et - eq // 2 - ea for (ea, eq, et) in bold})


def genus(bold):
    """種数表示（KN の G 展開）{(g, β, 2J_r): \\hat N}。閉じなければ ValueError。"""
    return g_expand(bold)


# ---- 特殊化 ----------------------------------------------------------------------------------------

def specialize_lower(lower, N, M=0):
    """小文字で a = t^N q^M（A → T^N Q^M）。{(e_Q, e_T): c}。"""
    out = defaultdict(int)
    for (ea, eq, et), c in lower.items():
        out[(eq + M * ea, et + N * ea)] += c
    return {k: v for k, v in out.items() if v}


def specialize_bold(bold, m, n, sign=1):
    """太字で 𝐚² → sign·𝐪^m 𝐭^n。{(e_𝐪, e_𝐭): c}。DGR の d_N の消去点は (m, n, sign) = (2N, −1, −1)、すなわち a = t^N。"""
    out = defaultdict(int)
    for (ea, eq, et), c in bold.items():
        b = ea // 2
        out[(eq + m * b, et + n * b)] += c * (sign ** (b % 2) if sign < 0 else 1)
    return {k: v for k, v in out.items() if v}


def is_monomial(poly):
    return len(poly) == 1


def a_derivative_at_one(lower):
    """小文字で d/da \\hat f |_{a=1}（A = a^{1/2}、d/da = (1/2A) d/dA）。{(e_Q, e_T): c}（有理数になりうるので 2 倍した整数を返す）。"""
    out = defaultdict(int)
    for (ea, eq, et), c in lower.items():
        out[(eq, et)] += ea * c
    return {k: v for k, v in out.items() if v}


def count_terms(poly):
    return len(poly)
