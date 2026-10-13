# 結び目表（data/framed_knot/reference/yui_kr/braid.csv、KnotInfo のブレイド）の HOMFLY-PT を Sage で計算する（ステージ T0 の Euler 標数の照合用）。
# 実行：DOT_SAGE=$TMPDIR/dotsage /opt/homebrew/Caskroom/miniforge/base/envs/sage/bin/sage scripts/framed_knot/knot_table_homfly.sage
# 出力：data/framed_knot/reference/yui_kr/sage_homfly.txt（名前<TAB>Sage の HOMFLY（変数 L, M、P(unknot)=1））
import csv

rows = list(csv.DictReader(open("data/framed_knot/reference/yui_kr/braid.csv")))
with open("data/framed_knot/reference/yui_kr/sage_homfly.txt", "w") as out:
    for r in rows:
        word = sage_eval(r["braid"])
        n = max(abs(g) for g in word) + 1
        h = Link(BraidGroup(n)(word)).homfly_polynomial()
        out.write(f"{r['name']}\t{h}\n")
print(len(rows))
