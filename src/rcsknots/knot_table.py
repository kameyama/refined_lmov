"""結び目表の reduced HOMFLY-PT ホモロジー（yui-kr、Nakagane–Sano arXiv:2111.00388）から [1] の refined LMOV 不変量を作る。

yui-kr の (i, j, k) は Poincaré 多項式 t^{(k−j)/2} q^i a^j（KnotInfo のブレイドのカイラリティ）。
DGR の太字の変数では (e_𝐚, e_𝐪, e_𝐭) = (j, i, (j−k)/2)。正のブレイドの三葉結び目が、第 1 期の T(2,3) の
superpolynomial 𝐚²𝐪⁻² + 𝐚²𝐪²𝐭² + 𝐚⁴𝐭³ になる（data/framed_knot/superpoly/T2_3.json）。
[1] では \\hat f_{[1]} = (1 + 𝐚²𝐭)/(𝐚²𝐭) · P（scripts/framed_knot/superpoly_lmov.py で T(2,3) と代数的なケーブルについて確認した因子）。
"""
import json
from pathlib import Path

from rcsknots.bps import g_expand

DATA = Path(__file__).resolve().parents[2] / "data/framed_knot/reference/yui_kr"


def load_homology():
    """{結び目名: {(i, j, k): 階数}}"""
    raw = json.loads((DATA / "homology.json").read_text())
    return {name: {(i, j, k): r for i, j, k, r in v} for name, v in raw.items()}


def superpoly(hom):
    """yui-kr の (i, j, k) → DGR の太字の Poincaré 多項式 {(e_𝐚, e_𝐪, e_𝐭): 階数}。"""
    return {(j, i, (j - k) // 2): r for (i, j, k), r in hom.items()}


def mirror(P):
    """DGR の鏡像 P(K*)(𝐚, 𝐪, 𝐭) = P(K)(𝐚⁻¹, 𝐪⁻¹, 𝐭⁻¹)。"""
    return {(-a, -q, -t): c for (a, q, t), c in P.items()}


def homfly(P):
    """𝐭 = −1 で得る HOMFLY-PT {(e_𝐚, e_𝐪): 係数}。"""
    out = {}
    for (a, q, t), c in P.items():
        out[(a, q)] = out.get((a, q), 0) + (-1) ** (t % 2) * c
    return {k: v for k, v in out.items() if v}


def delta_gradings(P):
    """δ 次数 e_𝐚 + e_𝐪/2 − e_𝐭 の値の集合（thin ⇔ 1 個）。yui-kr の i + j + k の −1/2 倍に当たる。"""
    return sorted({a + q / 2 - t for a, q, t in P})


def hat_f1(P):
    """[1] の \\hat f = (1 + 𝐚²𝐭)/(𝐚²𝐭) · P（太字）。"""
    out = {}
    for (a, q, t), c in P.items():
        for k in ((a - 2, q, t - 1), (a, q, t)):
            out[k] = out.get(k, 0) + c
    return {k: v for k, v in out.items() if v}


def judge1(P):
    """[1] の判定。G 展開の前の係数 \\hat{𝐍} と後の \\hat N、(g, β) ごとの 2J_r の偶奇。"""
    f = hat_f1(P)
    res = {"before": f, "before_neg": sum(v < 0 for v in f.values())}
    try:
        N = g_expand(f)
    except ValueError as e:
        res.update(closes=False, error=str(e))
        return res
    par = {}
    for g, b, F in N:
        par.setdefault((g, b), set()).add(F % 2)
    res.update(closes=True, after=N, after_neg=sum(v < 0 for v in N.values()),
               parity_aligned=all(len(s) == 1 for s in par.values()))
    return res


def parity_uniform(P):
    """[1] の (g, β) ごとの 2J_r の偶奇がそろう条件：𝐚 の次数が連続する範囲（2 ずつ）ごとに、e_𝐭 + e_𝐚/2 の偶奇が一定。

    \\hat f の 𝐚^{2β} の行には P の隣り合う 2 行しか入らず、G^g 𝐭^F の F は e_𝐭 − g と同じ偶奇（conventions.md §11）。
    """
    parity = {}
    for a, q, t in P:
        parity.setdefault(a, set()).add((t + a // 2) % 2)
    if any(len(v) > 1 for v in parity.values()):
        return False
    degs = sorted(parity)
    return all(parity[x] == parity[y] for x, y in zip(degs, degs[1:]) if y - x == 2)
