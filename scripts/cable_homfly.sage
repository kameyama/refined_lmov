# 結び目の (p, s)-ケーブル（s は Seifert framing、p と互いに素）の HOMFLY-PT をブレイドから Sage で計算する。
# ケーブル公式を逆に解いて |λ| = p の彩色 HOMFLY-PT を求めるため（第 3 期 U2、scripts/colored_from_cables.py）。
# 実行：DOT_SAGE=$TMPDIR/dotsage /opt/homebrew/Caskroom/miniforge/base/envs/sage/bin/sage scripts/cable_homfly.sage P S1,S2,... 名前 ...
# 出力：data/reference/cable_homfly_p{P}.txt（名前<TAB>s<TAB>Sage の HOMFLY（変数 L, M））に追記
import csv, sys, time

def cable_word(word, nstr, p2, q2):
    """scripts/braid_homfly.sage の cable_word と同じ（各ストランドを p2 本にし、束 1 に framing 補正 (σ1…σ_{p2−1})^{q2 − p2·w}）。"""
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

p = int(sys.argv[1])
slopes = [int(x) for x in sys.argv[2].split(",")]
braids = {r["name"]: sage_eval(r["braid"]) for r in csv.DictReader(open("data/reference/yui_kr/braid.csv"))}
with open(f"data/reference/cable_homfly_p{p}.txt", "a") as out:
    for name in sys.argv[3:]:
        word = braids[name]
        n = max(abs(g) for g in word) + 1
        for s in slopes:
            t0 = time.time()
            cw, m = cable_word(word, n, p, s)
            h = Link(BraidGroup(m)(cw)).homfly_polynomial()
            out.write(f"{name}\t{s}\t{h}\n")
            out.flush()
            print(name, s, len(cw), f"{time.time() - t0:.1f}s", flush=True)
