r"""二重ツイスト結び目 K_{m,n} の refined の reduced 不変量を、本研究の data と同じ形式（KN の小文字 (A, Q, T) の指数、framing 0、
(N1)+(N2) の規格化）で data/refined_kauffman/torus/{名前}_{彩色}.json に書く。使い方：python scripts/rk_k5_41_data.py [名前 m n]（既定は 4_1 −1 1）。
名前と (m, n) の対応は KNTZ 1902.02275 の図 1（Cyclotomic-Rectangular.tex :521-531 で確認）：3_1 = K_{1,1}, 4_1 = K_{−1,1}, 5_2 = K_{2,1},
6_1 = K_{−2,1}, 7_2 = K_{3,1}, 7_4 = K_{2,2}, 8_1 = K_{−3,1}, 8_3 = K_{−2,2}。
  長方形の彩色 r^s：KNTZ 1902.02275 の refined 版（scripts/x1_check.triply → t4_colored.bold_to_small）。この経路は trefoil（K_{1,1}）で
    GN の data（data/torus/T2_3_*）と [1], [2], [1,1], [2,2], [3], [1,1,1] の全係数が一致する（2026-10-19 に確認）。
  [2,1]：KN 1703.05408 のノートブックの候補（ミラー対称性などを課して求めたもの、確立していない）。4_1_21.json の meta に明記する。"""
import sys, json, time
from pathlib import Path
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
from x1_check import triply, bold_dict
from t4_colored import bold_to_small
from rcsknots.laurent import load_table

OUT = Path("data/refined_kauffman/torus")


def write(name, poly, source, knot="4_1"):
    meta = {"knot": knot, "color": name, "source": source, "reduced": True, "framing": "0 (normalization (N1)+(N2))",
            "variables": "A=a^(1/2), Q=q^(1/2), T=t^(1/2) of 1703.05408; exponents of (A,Q,T)",
            "poly": {",".join(map(str, k)): v for k, v in sorted(poly.items())}}
    (OUT / f"{knot}_{name}.json").write_text(json.dumps(meta, indent=1))


if __name__ == "__main__":
    knot, m, n = (sys.argv[1], int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else ("4_1", -1, 1)
    nb = load_table("data/reference/notebook_torus.txt")
    import os
    colors = [(1,), (2,), (1, 1), (3,), (1, 1, 1), (2, 2), (2, 2, 2), (1, 1, 1, 1), (3, 3)]
    if os.environ.get("SMALL"):        # N = 2 の SO 側の検査と c の決定に要る彩色だけ
        colors = colors[:6]
    for lam in colors:
        t0 = time.time()
        name = "".join(map(str, lam))
        poly = bold_to_small(bold_dict(triply(lam[0], len(lam), m, n)))
        ok = (knot, name) not in nb or nb[(knot, name)] == poly
        write(name, poly, f"KNTZ 1902.02275 refined (scripts/x1_check.triply), double twist K_{{{m},{n}}}", knot)
        print(f"{knot} [{name}]: {len(poly)} 項、KN のノートブックと一致 {ok}（{time.time()-t0:.1f}s）", flush=True)
    if knot == "4_1":
        write("21", nb[("4_1", "21")], "KN 1703.05408 notebook: candidate fixed by mirror symmetry etc. (NOT established; contains negative coefficients)")
        print("4_1 [21]: KN の候補を書いた（確立していない）", flush=True)
