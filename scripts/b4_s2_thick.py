"""B4：Stošić の quiver を展開し KRSS (quad-C) で作った S²（t_c 次数）から、thick な非トーラス結び目 10_139, 10_152 の 2 箱を判定する。

手順は scripts/b4_quiver_expand.py（8_19 で EHA の [2], [1,1], [3], [1,1,1] と全項一致した変換）と同じ。
- 10_152 の水準 3 の節点は Σ n_i d_i = 2 に入らないので r = 2 では落とす（P_2 が元の式と一致することで確かめる）。
- 対照：10_124 = T(3,5) を EHA の T3_5_2, T3_5_11 と照合する。
- 展開の順序への依存：圧縮した節点の並べ方をいくつか変えて S² が変わらないかを見る。
- Euler 標数：t = q 極限を、(2, s)-ケーブルの HOMFLY-PT（data/reference/two_cable_homfly.txt、scripts/two_cable_homfly.sage）から
  scripts/colored_from_cables.py で逆に解いた [2], [1,1] 彩色 HOMFLY-PT と比べる。
- 判定：scripts/t4_colored.py の judge2（\\hat f_{[2]}, \\hat f_{[1,1]} の G 展開の前後の負の項の数）。[1] は yui-kr の H^{[1]}。
出力：data/b4/<結び目>_s2_tc.json（t_c 次数の [2], [1,1]、DGR の太字 {"e_a,e_q,e_t": 係数}）。
使い方：python scripts/b4_s2_thick.py
"""
import json
import random
import sys
import time
from itertools import permutations
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))
from b4_819_quad import mirror11, quad_to_tc  # noqa: E402
from b4_quiver_expand import check_knot, expand, permuted, quadC  # noqa: E402
from b4_stosic_quiver import P_r  # noqa: E402
from rcsknots.knot_table import judge1, load_homology, superpoly  # noqa: E402
from t4_colored import judge2  # noqa: E402


def euler_vs_cables(name, P2, P11):
    """t_c 次数の [2], [1,1] の t = q 極限が、2 本ケーブルから逆に解いた彩色 HOMFLY-PT（colored_from_cables.colored）と一致するか。"""
    from colored_from_cables import colored
    from rcsknots.homfly import Q, a
    from t4_colored import bold_to_small
    col = colored(name)
    out = {}
    for mu, P in (((2,), P2), ((1, 1), P11)):
        e = sum(c * a ** sp.Rational(A, 2) * Q ** (q + t) for (A, q, t), c in bold_to_small(P).items())
        out[mu] = sp.expand(e - sp.expand(col[mu])) == 0
    return out


def level1(rec):
    keep = [i for i, n in enumerate(rec["n"]) if n == 1]
    return permuted(rec, keep)


def main():
    data = json.loads((ROOT / "data/reference/stosic_2402.03066/quivers.json").read_text())["knots"]
    H = load_homology()
    out_dir = ROOT / "data/b4"
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, eha in (("10_124", {"[2]": "T3_5_2", "[1,1]": "T3_5_11"}), ("10_139", None), ("10_152", None)):
        t0 = time.time()
        rec = data[name]
        rec1 = level1(rec)
        res = check_knot(name, rec1, rs=(2,), eha=eha)
        p2_full = sp.expand(P_r(rec, 2) - P_r(rec1, 2)) == 0
        P1 = superpoly(H[name])
        dim1 = sum(P1.values())
        print(f"{name}：[1] の生成元が yui-kr と一致 {res['gens_match_yui']}、P_2 が Stošić の式と一致 {res['P2_match']}"
              f"（水準 1 以外を落としても P_2 は不変 {p2_full}）、dim H^[2] = {res['dim2']}（(dim H^[1])² = {dim1 ** 2}）", flush=True)
        for k, v in res.items():
            if "EHA" in k:
                print(f"  {k}：{v}")
        P2 = quad_to_tc(res["quad2"])
        P11 = mirror11(P2)
        # 順序への依存（無作為に 10 通り）
        perms = list(permutations(range(len(rec1["n"]))))
        random.seed(0)
        same = 0
        for perm in random.sample(perms, min(10, len(perms))):
            kids, nodes, Cn = expand(permuted(rec1, perm))
            same += quad_to_tc(quadC(nodes, Cn, 2)) == P2
        print(f"  展開の順序を変えた 10 通りで S²（t_c 次数）が同じ：{same}/10")
        j1 = judge1(P1)
        if name != "10_124":
            print(f"  t = q 極限がケーブルから求めた [2], [1,1] 彩色 HOMFLY-PT と一致：{euler_vs_cables(name, P2, P11)}")
        print(f"  [1]：G 展開の前 負 {j1['before_neg']}/{len(j1['before'])}、後 負 {j1.get('after_neg')}/{len(j1.get('after', {}))}")
        print(f"  2 箱：{judge2(P1, P2, P11)}  （{time.time() - t0:.0f} s）", flush=True)
        (out_dir / f"{name}_s2_tc.json").write_text(json.dumps({
            "source": "Stošić arXiv:2402.03066 の quiver を展開（scripts/b4_quiver_expand.py）、KRSS arXiv:1707.04017 (quad-C)、"
                      "t_c 次数：P(𝐚,𝐪,𝐭) = P̃(a=𝐚, Q=𝐪, t_r=𝐪⁻¹, t_c=𝐭𝐪)。[1,1] はミラー (e_a,e_q,e_t) ↦ (e_a,−e_q,e_t−e_q)。予想に基づく",
            "variables": "DGR の太字 (e_a, e_q, e_t)",
            "[2]": {",".join(map(str, k)): v for k, v in sorted(P2.items())},
            "[1,1]": {",".join(map(str, k)): v for k, v in sorted(P11.items())}}, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
