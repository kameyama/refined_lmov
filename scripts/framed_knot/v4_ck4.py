"""第 4 期 V4：sl_2 正値性（v = 𝐪²𝐭 について行ごとに [k]_v の正係数の和）を、ツイスト結び目以外の結び目で調べる。

入力：data/framed_knot/b4/<結び目>_s2_tc.json（Stošić の quiver と KRSS の予想式による t_c 次数の [2], [1,1]、scripts/framed_knot/b4_s2_thick.py）と
yui-kr の H^{[1]}（rcsknots.knot_table）。10_124 = T(3,5)、10_139, 10_152 は thick な非トーラス結び目。
あわせて [1] の superpolynomial（yui-kr、ツイスト結び目は GNSSS）の sl_2 正値性も調べる。
使い方：python scripts/framed_knot/v4_ck4.py
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/framed_knot"))
from v4_closed import RHOS, hat_from_inputs, load_inputs  # noqa: E402
from v4_ck2 import sl2_report  # noqa: E402
from rcsknots.knot_table import load_homology, superpoly  # noqa: E402


def main():
    H = load_homology()
    for name in ("10_124", "10_139", "10_152"):
        d = json.loads((ROOT / f"data/framed_knot/b4/{name}_s2_tc.json").read_text())
        P2 = {tuple(map(int, k.split(","))): v for k, v in d["[2]"].items()}
        P1 = superpoly(H[name])
        hat = hat_from_inputs(P1, P2)
        print(f"{name}: " + "；".join(f"{rho} (σ 不変, 正, 負) {sl2_report(hat[rho])}" for rho in RHOS), flush=True)
    print("[1] の superpolynomial（yui-kr）:", {n: sl2_report(superpoly(H[n])) for n in ("5_2", "6_1", "8_19", "10_124", "10_139", "10_152")}, flush=True)
    cache = load_inputs()
    print("[1] の superpolynomial（GNSSS のツイスト結び目 K_p）:", {p: sl2_report(cache[p][0]) for p in (1, 2, 3, 5, -1, -2, -3, -5)}, flush=True)


if __name__ == "__main__":
    main()
