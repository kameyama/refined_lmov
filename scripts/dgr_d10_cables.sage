# D10（docs/dgr_lmov/research_plan.md）：yui-kr のブレイド（KnotInfo）から、結び目の HOMFLY-PT と (2, s)-ケーブル（s = 1, 3）の HOMFLY-PT を計算する。
# ケーブルの作り方は scripts/two_cable_homfly.sage と同じ。結び目自身の値は、[1]（yui-kr）と向きがそろっているかの確認に使う。
# 実行：DOT_SAGE=$TMPDIR/dotsage /opt/homebrew/Caskroom/miniforge/base/envs/sage/bin/sage scripts/dgr_d10_cables.sage 5_2 6_2 ...
# 出力：data/dgr_lmov/d10/homfly.txt（名前<TAB>s（結び目自身は "knot"）<TAB>Sage の HOMFLY（変数 L, M））に追記。既にある（名前, s）は飛ばす。
import csv, os, sys, time

def cable_word(word, nstr, p2, q2):
    w = sum(1 if g > 0 else -1 for g in word)
    out = []
    for g in word:
        i, s = abs(g), (1 if g > 0 else -1)
        block = []
        for r in range(p2):
            for c in range(p2):
                block.append(p2 * (i - 1) + p2 - r + c)
        out += [s * x for x in block]
    twist = q2 - p2 * w
    gens = list(range(1, p2))
    out += [(1 if twist > 0 else -1) * x for _ in range(abs(twist)) for x in gens]
    return out, nstr * p2

path = "data/dgr_lmov/d10/homfly.txt"
have = {tuple(l.split("\t")[:2]) for l in open(path)} if os.path.exists(path) else set()
braids = {r["name"]: sage_eval(r["braid"]) for r in csv.DictReader(open("data/reference/yui_kr/braid.csv"))}
with open(path, "a") as out:
    for name in sys.argv[1:]:
        word = braids[name]
        n = max(abs(g) for g in word) + 1
        for s in ("knot", 1, 3):
            if (name, str(s)) in have:
                continue
            t0 = time.time()
            if s == "knot":
                cw, m = word, n
            else:
                cw, m = cable_word(word, n, 2, s)
            h = Link(BraidGroup(m)(cw)).homfly_polynomial()
            out.write(f"{name}\t{s}\t{h}\n")
            out.flush()
            print(name, s, len(cw), f"{time.time() - t0:.1f}s", flush=True)
