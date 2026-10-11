"""S1：11 交点以下の全素結び目 801 個（yui-kr の HOMFLY-PT ホモロジー）の [1] で、既約表示の非負性・δ の台・特殊化を調べる。

[1] では \\hat f_[1] = (1 + 𝐚^{−2}𝐭^{−1}) P（太字）。特殊化は太字で 𝐚² → −𝐪^{2N}𝐭^{−1}（a = t^N、DGR の d_N の対が打ち消す点）。
出力：data/lmovstruct/knot_table.json、data/lmovstruct/knot_table.md、logs/lmovstruct_knot_table.log。
"""
import json
import time
from collections import Counter
from pathlib import Path

from rcsknots.knot_table import hat_f1, load_homology, superpoly
from rcsknots.lmovstruct import delta_values, is_W_invariant, specialize_bold, spin, spin_count

OUT = Path("data/lmovstruct")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    hom = load_homology()
    res = {}
    for name, h in hom.items():
        P = superpoly(h)
        f = hat_f1(P)
        assert is_W_invariant(f)
        n_terms, n_neg = spin_count(f)
        sP = spin(P)
        r = {"terms": len(f), "spin_terms": n_terms, "spin_neg": n_neg, "P_spin_neg": sum(1 for v in sP.values() if v < 0),
             "delta_P": delta_values(P), "delta_f": delta_values(f),
             "spec": {f"a=t^{N}": len(specialize_bold(f, 2 * N, -1, -1)) for N in range(-3, 4)},
             "P_at_a_eq_t_monomial": len(specialize_bold(P, 2, -1, -1)) == 1,
             "P_at_a_eq_1_terms": len(specialize_bold(P, 0, -1, -1)),
             "P_at_a_eq_t2_terms": len(specialize_bold(P, 4, -1, -1)), "dim": sum(P.values())}
        # a = t^2 の特殊化（sl(2) の消去点）で残る項数と、符号が揃うか
        s2 = specialize_bold(P, 4, -1, -1)
        r["P_at_a_eq_t2_signs_uniform"] = len({v > 0 for v in s2.values()}) <= 1
        res[name] = r
    (OUT / "knot_table.json").write_text(json.dumps(res, ensure_ascii=False, indent=0))
    n = len(res)
    spin_bad = [k for k, r in res.items() if r["spin_neg"]]
    P_spin_bad = [k for k, r in res.items() if r["P_spin_neg"]]
    mono = sum(r["P_at_a_eq_t_monomial"] for r in res.values())
    width = Counter(len(r["delta_P"]) for r in res.values())
    uni = sum(r["P_at_a_eq_t2_signs_uniform"] for r in res.values())
    lines = [f"# 801 結び目の [1]（scripts/lmovstruct_knot_table.py、{time.strftime('%Y-%m-%d')}）", "",
             f"- 結び目の数：{n}",
             f"- \\hat f_[1] の既約表示に負の項がある結び目：{len(spin_bad)}（{spin_bad[:10]}）",
             f"- P 自身の既約表示に負の項がある結び目：{len(P_spin_bad)}（{P_spin_bad[:10]}）",
             f"- P(a = t)（𝐚² = −𝐪²/𝐭、d_1 の消去点）が単項式：{mono} / {n}",
             f"- P の δ = f − 2β の値の数の分布：{dict(sorted(width.items()))}",
             f"- P(a = t²)（d_2 の消去点）で残る項の符号が一様：{uni} / {n}",
             "", "| 結び目 | 階数 | \\hat f の項数 | 既約項数 | 既約負 | δ(P) | P(a=t) 単項式 | P(a=1) 項数 | P(a=t²) 項数 | 符号一様 |", "|---|---|---|---|---|---|---|---|---|---|"]
    for k, r in res.items():
        lines.append(f"| {k} | {r['dim']} | {r['terms']} | {r['spin_terms']} | {r['spin_neg']} | {r['delta_P']} | {r['P_at_a_eq_t_monomial']} | {r['P_at_a_eq_1_terms']} | {r['P_at_a_eq_t2_terms']} | {r['P_at_a_eq_t2_signs_uniform']} |")
    (OUT / "knot_table.md").write_text("\n".join(lines) + "\n")
    with open("logs/lmovstruct_knot_table.log", "a") as fp:
        fp.write("\n".join(lines[:8]) + f"\n{time.time() - t0:.1f} 秒\n")
    print("\n".join(lines[:8]))


if __name__ == "__main__":
    main()
