"""open_problems B2：GNSSS の (2s−1,1,2t−1)-pretzel K_{s,t} を yui-kr の結び目表（11 交点以下）と照合し、thin/thick と [1] の一致を調べる。

GNSSS の [1]（t_c 次数 = DGR）の 𝐭 = −1 の HOMFLY-PT を、yui-kr の各結び目とその鏡像の HOMFLY-PT と比べて結び目を特定する。
特定できたものについて、yui-kr のホモロジーの δ 次数の個数（thin ⇔ 1）と、GNSSS の [1] が実際のホモロジーに一致するかを出力する。
2 箱の判定（logs/framed_knot/t3_pretzel.log の「次元が二乗か」）と並べる。結果は logs/framed_knot/b2_pretzel_thickness.log。
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from t3_pretzel import tc_graded  # noqa: E402
from rcsknots.knot_table import delta_gradings, homfly, load_homology, mirror, superpoly  # noqa: E402

squares = {}
for line in open("logs/framed_knot/t3_pretzel.log"):
    m = re.match(r"K_\{(-?\d+),(-?\d+)\}: rank (\d+), (\d+)", line)
    if m:
        squares[(int(m[1]), int(m[2]))] = int(m[4]) == int(m[3]) ** 2

H = {name: superpoly(h) for name, h in load_homology().items()}
by_homfly = {}
for name, P in H.items():
    for chir, Q in (("", P), ("*", mirror(P))):
        by_homfly.setdefault(tuple(sorted(homfly(Q).items())), []).append((name + chir, Q))

for (s, t), sq in sorted(squares.items()):
    P1 = tc_graded(s, t, 1)
    hits = by_homfly.get(tuple(sorted(homfly(P1).items())), [])
    if not hits:
        print(f"K_{{{s},{t}}}（2 箱：次元{'二乗' if sq else '二乗でない'}）：11 交点以下に該当なし", flush=True)
        continue
    for name, Q in hits:
        thick = len(delta_gradings(Q)) > 1
        print(f"K_{{{s},{t}}}（2 箱：次元{'二乗' if sq else '二乗でない'}）= {name}：{'thick' if thick else 'thin'}、"
              f"GNSSS の [1] が実際のホモロジーと{'一致' if Q == P1 else '不一致'}（階数 {sum(P1.values())} / {sum(Q.values())}）", flush=True)
