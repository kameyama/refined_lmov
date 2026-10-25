# 結び目の (2, s)-ケーブル（s は Seifert framing）の HOMFLY-PT をブレイドから Sage で計算する。
# ケーブル公式を逆に解いて [2], [1,1] 彩色 HOMFLY-PT を求めるため（ステージ T4、scripts/framed_knot/colored_from_cables.py）。
# 実行：DOT_SAGE=$TMPDIR/dotsage /opt/homebrew/Caskroom/miniforge/base/envs/sage/bin/sage scripts/framed_knot/two_cable_homfly.sage 3_1 4_1 9_42
# 出力：data/framed_knot/reference/two_cable_homfly.txt（名前<TAB>s<TAB>Sage の HOMFLY（変数 L, M））に追記
import csv, sys, time

def cable_word(word, nstr, p2, q2):
    """scripts/framed_knot/braid_homfly.sage の cable_word と同じ（各ストランドを p2 本にし、束 1 に framing 補正 (σ1…σ_{p2−1})^{q2 − p2·w}）。"""
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

braids = {r["name"]: sage_eval(r["braid"]) for r in csv.DictReader(open("data/framed_knot/reference/yui_kr/braid.csv"))}
with open("data/framed_knot/reference/two_cable_homfly.txt", "a") as out:
    for name in sys.argv[1:]:
        word = braids[name]
        n = max(abs(g) for g in word) + 1
        for s in (1, 3):
            t0 = time.time()
            cw, m = cable_word(word, n, 2, s)
            h = Link(BraidGroup(m)(cw)).homfly_polynomial()
            out.write(f"{name}\t{s}\t{h}\n")
            out.flush()
            print(name, s, len(cw), f"{time.time() - t0:.1f}s", flush=True)
