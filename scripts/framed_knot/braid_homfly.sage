# ケーブルのブレイドから HOMFLY-PT を計算する（Stage 4 の独立検証）。
# 実行：DOT_SAGE=$TMPDIR/dotsage /opt/homebrew/Caskroom/miniforge/base/envs/sage/bin/sage scripts/framed_knot/braid_homfly.sage
# 出力：data/framed_knot/reference/braid_homfly.txt（名前<TAB>ストランド数<TAB>ブレイド語<TAB>Sage の HOMFLY（変数 L, M、P(unknot)=1））
#
# companion T(p1, q1) のブレイド (σ1…σ_{p1−1})^{q1}（writhe w = q1(p1−1)）の各ストランドを p2 本の平行線にし、
# 1 つの束に (σ1…σ_{p2−1})^{q2 − p2·w} を足すと、(p2, q2)-ケーブル（q2 は companion の Seifert framing）になる。

import os

def companion_word(p1, q1):
    return [i for _ in range(q1) for i in range(1, p1)]

def cable_word(word, nstr, p2, q2):
    """Cable a braid word on nstr strands by p2 parallels, then add the framing correction on bundle 1."""
    w = sum(1 if g > 0 else -1 for g in word)
    out = []
    for g in word:
        i, s = abs(g), (1 if g > 0 else -1)
        # bundle i (strands p2(i−1)+1 … p2 i) crosses bundle i+1
        block = []
        for r in range(p2):
            for c in range(p2):
                block.append(p2 * (i - 1) + p2 - r + c)
        out += [s * x for x in block]
    twist = q2 - p2 * w
    gens = list(range(1, p2))
    out += [(1 if twist > 0 else -1) * x for _ in range(abs(twist)) for x in gens]
    return out, nstr * p2

def homfly(word, n):
    if not word:
        return 1
    return Link(BraidGroup(n)(word)).homfly_polynomial()

cases = []
# framing-correction checks: unknot drawn as σ1 (w = 1) on 2 strands, cabled → T(p2, q2)
for q2 in (1, 3, 5, -3):
    w, n = cable_word([1], 2, 2, q2)
    cases.append((f"check_kinkunknot_2_{q2}", n, w))
    cases.append((f"check_T2_{q2}", 2, [1 if q2 > 0 else -1] * abs(q2)))
w, n = cable_word([1, 1, 1], 2, 1, 3)        # p2 = 1 returns the companion
cases.append(("check_p1_trefoil", n, w))
cases.append(("T2_3", 2, [1, 1, 1]))
# the knots of the research plan: (p2, q2)-cables of T(p1, q1)
for name, (p1, q1, p2, q2) in {"K1": (2, 3, 2, 13), "K2": (2, 3, 2, 1), "K3": (2, 3, 2, 15), "K4": (2, 3, 2, 3),
                               "K5a": (2, 3, 2, 5), "K5b": (2, 3, 2, 7), "K5c": (2, 3, 2, -1),
                               "K6a": (2, 5, 2, 21), "K6b": (2, 5, 2, 1)}.items():
    w, n = cable_word(companion_word(p1, q1), p1, p2, q2)
    cases.append((name, n, w))

rows = []
for name, n, w in cases:
    rows.append(f"{name}\t{n}\t{' '.join(map(str, w))}\t{homfly(w, n)}")
    print(name, n, len(w))
os.makedirs("data/framed_knot/reference", exist_ok=True)
open("data/framed_knot/reference/braid_homfly.txt", "w").write("\n".join(rows) + "\n")
