# 別研究 [2,1] の X0：ダブルツイスト結び目 K_{m,n}（KNTZ 1902.02275 の図 1）のブレイドと、ケーブルの HOMFLY-PT を Sage で計算する。
# K_{m,n} は 2m 交点の水平なツイストと 2n 交点の垂直なツイストを持つ図式。PD コードを作り Link.braid()（Vogel）でブレイドにする。
# 向き（どちらのツイストを正とするか）とミラーは、[1] の HOMFLY-PT を KNTZ の式と比べて決める（scripts/framed_knot/x0_check.py、tests/framed_knot/test_x0_double_twist.py）。
#
# ブレイド：ツイスト結び目 K_{p,1} は短い語の規則 twist_word（p > 0 は 5_2 → 7_2 → 9_2 の規則、p < 0 は 6_1 → 8_1 → 10_1 の規則。
#   どちらも yui-kr の表（data/framed_knot/reference/yui_kr/braid.csv）の語と |p| ≤ 4 で一致）、表にあるダブルツイスト結び目は表の語、それ以外は PD → Vogel。
# 実行：DOT_SAGE=$TMPDIR/dotsage /opt/homebrew/Caskroom/miniforge/base/envs/sage/bin/sage scripts/framed_knot/x0_double_twist.sage 種類 m,n ...
#   種類 braid ：data/framed_knot/reference/x0_double_twist.txt（名前<TAB>本数<TAB>ブレイド語<TAB>[1] の HOMFLY（Sage の L, M））に追記
#   種類 p2    ：2 本ケーブル（s = 1, 3, 5）→ data/framed_knot/reference/x0_cable_p2.txt
#   種類 p3    ：3 本ケーブル（s = 1, 2, 4, 5）→ data/framed_knot/reference/x0_cable_p3.txt
#   種類 p3w   ：3 本ケーブルで傾き 3w ± 1, 3w ± 2（w はブレイドの writhe。framing 補正のねじれが小さく、libhomfly の 64 ビットの係数があふれにくい）→ x0_cable_p3.txt
#   種類 p22   ：2 段の 2 本ケーブル C_{2,s'}(C_{2,s}K)（(s, s') ∈ {1,3,5}×{1,3}）→ data/framed_knot/reference/x0_cable_p22.txt
#   出力の行：名前<TAB>傾き<TAB>Sage の HOMFLY（変数 L, M）。名前は D{m}_{n}（負は m を "m" で書く、例 Dm1_1 = K_{−1,1}）。
# 既に出力にある（名前, 傾き）は飛ばす（再開できる）。
import os, sys, time

def cable_word(word, nstr, p2, q2):
    """scripts/framed_knot/cable_homfly.sage の cable_word と同じ（各ストランドを p2 本にし、束 1 に framing 補正 (σ1…σ_{p2−1})^{q2 − p2·w}）。"""
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


def double_twist_pd(m, n):
    """K_{m,n} の PD コード（Sage の規約：入ってくる下の辺から反時計回り）。

    交点のポートは反時計回りに 0 = NE, 1 = NW, 2 = SW, 3 = SE。糸は (0,2) と (1,3)。
    水平なツイスト（2|m| 交点）を上に、垂直なツイスト（2|n| 交点）を下に置き、
    水平の下の 2 端を垂直の上の 2 端に、垂直の下の 2 端を左右の外側を回して水平の上の 2 端につなぐ。
    under[c] は交点 c で下を通る糸のポートの組。m > 0 の水平は (1,3) が下、n > 0 の垂直は (0,2) が下（どちらも全交点で同じ）。
    """
    H, V = list(range(2 * abs(m))), list(range(2 * abs(m), 2 * abs(m) + 2 * abs(n)))
    under = {}
    for c in H:
        under[c] = (1, 3) if m > 0 else (0, 2)
    for c in V:
        under[c] = (0, 2) if n > 0 else (1, 3)
    link = {}

    def join(x, y):
        link[x], link[y] = y, x
    for i in range(len(H) - 1):
        join((H[i], 0), (H[i + 1], 1))
        join((H[i], 3), (H[i + 1], 2))
    for i in range(len(V) - 1):
        join((V[i], 2), (V[i + 1], 1))
        join((V[i], 3), (V[i + 1], 0))
    join((H[0], 2), (V[0], 1))
    join((H[-1], 3), (V[0], 0))
    join((V[-1], 2), (H[0], 1))
    join((V[-1], 3), (H[-1], 0))
    # 一周たどって辺に番号を付ける
    label, incoming = {}, {}
    cur, e = (H[0], 1), 1
    for _ in range(2 * len(link)):
        c, p = cur
        incoming[(c, p)] = True
        out = (c, (p + 2) % 4)
        label[out] = e
        nxt = link[out]
        label[nxt] = e
        e += 1
        cur = nxt
        if cur == (H[0], 1):
            break
    assert len(label) == 4 * (len(H) + len(V)), "結び目でない（成分が 2 つ以上）"
    pd = []
    for c in H + V:
        u = under[c][0] if incoming.get((c, under[c][0])) else under[c][1]
        pd.append([label[(c, (u + k) % 4)] for k in range(4)])
    return pd


def twist_word(p):
    """K_{p,1} のブレイド語（p ≠ 0）。"""
    if p == 1:
        return [1, 1, 1]
    if p == -1:
        return [1, -2, 1, -2]
    if p > 0:
        w = [1, 1, 1, 2, -1, 2]
        for k in range(2, p):
            w += [k + 1, -k, k + 1]
        return w
    w = [1, 1, 2, -1]
    for k in range(2, -p):
        w += [k, k + 1, -k]
    w += [p - 1, -p, p - 1]
    return [-g for g in w]   # 6_1 → 8_1 → 10_1 の表の語は K_{p,1} の鏡像なので反転する（[1] で確認）


TABLE = {(2, 2): "7_4", (-2, 2): "8_3", (3, 2): "9_5", (-3, 2): "10_3"}
MIRROR = {(-3, 2)}   # 表の 10_3 の語は K_{−3,2} の鏡像（[1] で確認）


def braid_of(m, n):
    if n == 1:
        w = twist_word(m)
        return w, max(abs(g) for g in w) + 1, "rule"
    if (m, n) in TABLE:
        import csv
        rows = {r["name"]: sage_eval(r["braid"]) for r in csv.DictReader(open("data/framed_knot/reference/yui_kr/braid.csv"))}
        w = rows[TABLE[(m, n)]]
        w = [-g for g in w] if (m, n) in MIRROR else w
        return w, max(abs(g) for g in w) + 1, "table " + TABLE[(m, n)]
    b = Link(double_twist_pd(m, n)).braid()
    return list(b.Tietze()), b.strands(), "vogel"


def name_of(m, n):
    f = lambda x: str(x) if x >= 0 else "m" + str(-x)
    return f"D{f(m)}_{f(n)}"


def load_braids(path="data/framed_knot/reference/x0_double_twist.txt"):
    out = {}
    if os.path.exists(path):
        for line in open(path):
            name, k, word, _ = line.rstrip("\n").split("\t")
            out[name] = ([int(x) for x in word.split()], int(k))
    return out


def done(path):
    if not os.path.exists(path):
        return set()
    return {tuple(line.split("\t")[:2]) for line in open(path)}


kind = sys.argv[1]
pairs = [tuple(int(x) for x in a.split(",")) for a in sys.argv[2:]]
os.makedirs("data/framed_knot/reference", exist_ok=True)
if kind == "braid":
    path = "data/framed_knot/reference/x0_double_twist.txt"
    have = load_braids(path)
    with open(path, "a") as out:
        for m, n in pairs:
            nm = name_of(m, n)
            if nm in have:
                continue
            word, k, src = braid_of(m, n)
            h = Link(BraidGroup(k)(word)).homfly_polynomial()
            out.write(f"{nm}\t{k}\t{' '.join(map(str, word))}\t{h}\n")
            out.flush()
            print(nm, k, len(word), src, flush=True)
else:
    braids = load_braids()
    path = f"data/framed_knot/reference/x0_cable_{'p3' if kind == 'p3w' else kind}.txt"
    have = done(path)
    with open(path, "a") as out:
        for m, n in pairs:
            nm = name_of(m, n)
            word, k = braids[nm]
            if kind == "p3w":
                w0 = 3 * sum(1 if g > 0 else -1 for g in word)
                jobs = [(s, str(s)) for s in (w0 - 2, w0 - 1, w0 + 1, w0 + 2)]
            elif kind == "p22":
                jobs = [((s, s2), f"{s},{s2}") for s in (1, 3, 5) for s2 in (1, 3)]
            else:
                jobs = [(s, str(s)) for s in ((1, 3, 5) if kind == "p2" else (1, 2, 4, 5))]
            for s, key in jobs:
                if (nm, key) in have:
                    continue
                t0 = time.time()
                if kind == "p22":
                    w1, k1 = cable_word(word, k, 2, s[0])
                    cw, kk = cable_word(w1, k1, 2, s[1])
                else:
                    cw, kk = cable_word(word, k, 2 if kind == "p2" else 3, s)
                # libhomfly は係数を 64 ビット整数で持ち、あふれても警告しない。2^62 以上は「あふれの疑い」として読む側で除く

                h = Link(BraidGroup(kk)(cw)).homfly_polynomial()
                out.write(f"{nm}\t{key}\t{h}\n")
                out.flush()
                print(nm, key, kk, len(cw), f"{time.time() - t0:.1f}s", flush=True)
