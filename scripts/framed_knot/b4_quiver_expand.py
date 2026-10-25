"""B4：Stošić（arXiv:2402.03066）の quiver 型の式 (specred2) の q-Pochhammer を展開して通常の quiver（水準 1 の節点）に直し、
KRSS（arXiv:1707.04017）の予想 (quad-C)（knots-quivers-long-arXiv.tex 740-746 行）で S^r の四重次数の Poincaré 多項式を作る。

展開の手順（Stošić 付録 A、2951-3010 行の 5_1 の例と同じ q-二項展開。順序は節点番号の順、(op1) 型）：
- (a^{-1}q;q^{-1})_{Σ r_p d_p} = Π_p (a^{-1}q^{1-o_p};q^{-1})_{d_p}（o_p は前の r 型の節点の d の和）。
- (a^{-1}q^{-r};q^{-1})_{Σ l_p d_p} = Π_p (a^{-1}q^{-r-o'_p};q^{-1})_{d_p}（o'_p は前の l 型の節点の d の和）。r と l の両方をもつ節点では、
  r で分けた 2 つの子に順に (·;q^{-1})_{E}(·q^{-E};q^{-1})_{F} と割り当てる。
- (z;q^{-1})_d = Σ_{e+f=d} (−z)^e q^{−e(e−1)/2 − ef} (q;q)_d/((q;q)_e (q;q)_f)。
展開後の q の指数は子の重複度の 2 次式で、そこから C（対角 = t 次数）を読む。
(quad-C)：P_r(a,Q,t_r,t_c) = Σ_{Σd=r} (t_c²;t_c²)_r/Π(t_c²;t_c²)_{d_i} a^{Σ a_i d_i} Q^{Σ q_i d_i} t_r^{Σ t_i d_i} t_c^{Σ C_ij d_i d_j}
（a, Q は DGR の太字の次数。(a_i, q_i, t_i) は [1] の生成元、t_i = C_ii）。
使い方：python scripts/framed_knot/b4_quiver_expand.py（5_1 で Stošić の C^{(1)} を再現することを確かめ、8_19 で展開する）
"""
import json
from itertools import permutations
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[2]


def expand(rec):
    """Stošić の (C, n, a, q, l, r)（水準 1 のみ）→ 子の節点の一覧 [(a_i, q_i, C の行)]（Stošić の単位：a = 𝐚²、q = 𝐪²）。"""
    C, n, av, qv, lv, rv = (rec[x] for x in ("C", "n", "a", "q", "l", "r"))
    assert set(n) == {1}
    k = len(n)
    kids, dsym = [], []          # kids[i] = (親 p, 'E'/'F'/'G'... の印, 記号)
    expo = sp.Integer(0)
    aexp = sp.Integer(0)
    sign = sp.Integer(0)        # (−1) の指数
    # 子の記号
    parts = []
    for p in range(k):
        if rv[p] and lv[p]:
            names = ["rE_lG", "rE_lK", "rF_lH", "rF_lK"]   # E = G + E', F = H + F'
        elif rv[p] or lv[p]:
            names = ["take", "keep"]
        else:
            names = ["only"]
        syms = [sp.Symbol(f"e{p}_{s}") for s in names]
        parts.append(syms)
        kids += [(p, s) for s in names]
        dsym += syms
    d = [sum(s) for s in parts]
    rtot = sum(dsym)
    expo += sum(qv[p] * d[p] for p in range(k)) + sp.Rational(1, 2) * (
        sum(C[p][s] * d[p] * d[s] for p in range(k) for s in range(k)) - sum(C[p][p] * d[p] for p in range(k)))
    aexp += sum(av[p] * d[p] for p in range(k))
    sign += sum(C[p][p] * d[p] for p in range(k))

    def qbin_split(zq_exp, e, f):
        """(z;q^{-1})_{e+f} の項 (−z)^e q^{-e(e-1)/2 - ef}、z = a^{-1} q^{zq_exp}。"""
        return -e, e, zq_exp * e - e * (e - 1) / 2 - e * f

    o = sp.Integer(0)
    for p in range(k):
        if not rv[p]:
            continue
        S = parts[p]
        E, F = (S[0] + S[1], S[2] + S[3]) if lv[p] else (S[0], S[1])
        da, ds, dq = qbin_split(1 - o, E, F)
        aexp, sign, expo = aexp + da, sign + ds, expo + dq
        o += d[p]
    o = sp.Integer(0)
    for p in range(k):
        if not lv[p]:
            continue
        S = parts[p]
        if rv[p]:
            for (G, K), extra in (((S[0], S[1]), 0), ((S[2], S[3]), S[0] + S[1])):
                da, ds, dq = qbin_split(-rtot - o - extra, G, K)
                aexp, sign, expo = aexp + da, sign + ds, expo + dq
        else:
            da, ds, dq = qbin_split(-rtot - o, S[0], S[1])
            aexp, sign, expo = aexp + da, sign + ds, expo + dq
        o += d[p]
    expo = sp.expand(expo)
    m = len(dsym)
    Cn = [[int(sp.diff(expo, dsym[i], dsym[j])) for j in range(m)] for i in range(m)]
    zero = {s: 0 for s in dsym}
    nodes = []
    for i in range(m):
        lin = sp.diff(expo, dsym[i]).subs(zero)
        qi = lin + sp.Rational(Cn[i][i], 2)
        ai = sp.diff(aexp, dsym[i])
        si = sp.diff(sign, dsym[i])
        assert (si - Cn[i][i]) % 2 == 0, (kids[i], si, Cn[i][i])
        nodes.append((int(ai), int(qi), Cn[i][i]))
    return kids, nodes, Cn


def compositions(r, m):
    if m == 1:
        yield (r,)
        return
    for x in range(r + 1):
        for rest in compositions(r - x, m - 1):
            yield (x,) + rest


def quiver_poly(nodes, Cn, r):
    """展開した quiver（Stošić の (kq2) の形、水準 1）から P_r(a, q)（Stošić 単位）。"""
    a, q = sp.symbols("a q")
    m = len(nodes)
    tot = sp.Integer(0)
    qp = lambda k: sp.prod([1 - q ** (i + 1) for i in range(k)])  # noqa: E731
    for d in compositions(r, m):
        quad = sum(Cn[i][j] * d[i] * d[j] for i in range(m) for j in range(m))
        diag = sum(Cn[i][i] * d[i] for i in range(m))
        tot += sp.cancel(qp(r) / sp.prod([qp(x) for x in d])) * (-1) ** diag * a ** sum(nodes[i][0] * d[i] for i in range(m)) \
            * q ** sum(nodes[i][1] * d[i] for i in range(m)) * q ** sp.Rational(quad - diag, 2)
    return sp.expand(tot)


def permuted(rec, perm):
    """圧縮した節点の順番を入れ替える（展開の順序の選び方）。"""
    out = {"C": [[rec["C"][i][j] for j in perm] for i in perm]}
    for x in "naqlr":
        out[x] = [rec[x][i] for i in perm]
    return out


def quadC(nodes, Cn, r):
    """(quad-C) の四重次数 {(a, Q, t_r, t_c): 係数}（a, Q は太字の次数）。"""
    tc = sp.Symbol("tc")
    m = len(nodes)
    out = {}

    def poch(k):
        return sp.prod([1 - tc ** (2 * (i + 1)) for i in range(k)])
    for d in compositions(r, m):
        ratio = sp.Poly(sp.cancel(poch(r) / sp.prod([poch(x) for x in d])), tc)
        A = sum(2 * nodes[i][0] * d[i] for i in range(m))
        Qd = sum(2 * nodes[i][1] * d[i] for i in range(m))
        R = sum(nodes[i][2] * d[i] for i in range(m))
        base = sum(Cn[i][j] * d[i] * d[j] for i in range(m) for j in range(m))
        for (e,), c in ratio.terms():
            key = (A, Qd, R, base + e)
            out[key] = out.get(key, 0) + int(c)
    return {k: v for k, v in out.items() if v}


def check_knot(name, rec, rs=(2, 3), eha=None, log=print):
    """展開 → [1] の照合 → P_r(a,q) の照合 → (quad-C) を t_c 次数に直して EHA（あれば）と照合。結果の dict を返す。"""
    import sys
    sys.path.insert(0, str(ROOT / "scripts/framed_knot"))
    sys.path.insert(0, str(ROOT / "src"))
    from b4_819_quad import compare, load_eha, mirror11, quad_to_tc
    from b4_stosic_quiver import P_r
    from rcsknots.knot_table import load_homology, superpoly
    from t4_colored import bold_to_small
    kids, nodes, Cn = expand(rec)
    gens = {}
    for (ai, qi, ti) in nodes:
        gens[(2 * ai, 2 * qi, ti)] = gens.get((2 * ai, 2 * qi, ti), 0) + 1
    res = {"gens_match_yui": gens == superpoly(load_homology()[name])}
    for r in rs:
        res[f"P{r}_match"] = sp.expand(quiver_poly(nodes, Cn, r) - P_r(rec, r)) == 0
        quad = quadC(nodes, Cn, r)
        res[f"quad{r}"] = quad
        tcP = quad_to_tc(quad)
        res[f"dim{r}"] = sum(tcP.values())
        if eha:
            for lab, P in ((f"[{r}]", tcP), ("[" + ",".join(["1"] * r) + "]", mirror11(tcP))):
                ref = eha.get(lab)
                if ref:
                    res[f"{lab} vs EHA"] = compare(bold_to_small(P), load_eha(ref))
    return res


def main():
    # 5_1：Stošić 付録 A の C^{(1)}（3017-3023 行）を再現するか
    rec51 = {"C": [[0, 1, 3], [1, 3, 4], [3, 4, 5]], "n": [1, 1, 1], "a": [2, 3, 3], "q": [-2, -1, 1], "l": [0, 0, 0], "r": [0, 1, 1]}
    kids, nodes, Cn = expand(rec51)
    ok = Cn == [[0, 1, 1, 3, 3], [1, 2, 2, 3, 4], [1, 2, 3, 3, 4], [3, 3, 3, 4, 4], [3, 4, 4, 4, 5]]
    print("5_1：展開した C が Stošić の C^{(1)}（3017-3023 行）に一致：", ok)
    data = json.loads((ROOT / "data/framed_knot/reference/stosic_2402.03066/quivers.json").read_text())["knots"]
    rec = data["8_19"]
    kids, nodes, Cn = expand(rec)
    print("8_19 の子の節点（Stošić 単位 (a, q, t)）", list(zip(kids, nodes)))
    print("8_19 の C（節点番号の順で展開）", Cn)
    import sys
    sys.path.insert(0, str(ROOT / "scripts/framed_knot"))
    from b4_819_quad import ggs_s2_T34, quad_terms
    eha = {"[2]": "T3_4_2", "[1,1]": "T3_4_11", "[3]": "T3_4_3", "[1,1,1]": "T3_4_111"}
    res = check_knot("8_19", rec, eha=eha)
    print("8_19（節点番号の順）：[1] の生成元が yui-kr と一致", res["gens_match_yui"], "、P_2, P_3 が Stošić の式と一致", res["P2_match"], res["P3_match"],
          "、次元", res["dim2"], res["dim3"])
    print("  (quad-C) の S² が GGS の四重次数（2058-2072 行）と全項一致：", res["quad2"] == quad_terms(ggs_s2_T34()))
    for k, v in res.items():
        if "EHA" in k:
            print(f"  {k}：{v}")
    # 展開の順序（圧縮した節点の並べ方）への依存
    summary = {}
    for perm in permutations(range(len(rec["n"]))):
        r2 = check_knot("8_19", permuted(rec, perm), rs=(2, 3), eha=eha)
        key = tuple(("一致" in str(r2[k])) for k in sorted(r2) if "EHA" in k)
        summary.setdefault(key, []).append(perm)
    print("  順序 120 通りの集計（", sorted(k for k in res if "EHA" in k), "の一致）：", {k: len(v) for k, v in summary.items()})


if __name__ == "__main__":
    main()
