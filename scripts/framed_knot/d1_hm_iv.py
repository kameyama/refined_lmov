"""open_problems D1 (IV) の一般の証明の補助計算（Hogancamp–Mellit arXiv:1909.00418 の漸化式 def:the polys）。

記号：P(v,w) は HM の多項式。トーラス結び目 T(m,n)（互いに素）では reduced の Poincaré 級数は R = P(10^{m−1},10^{n−1})/(1+a)。
W(a^i q^j t^k) = 2j − i。示したいこと (IV_HM)：R の各項で W ≤ 2δ、等号は a⁰q^δ の 1 項（係数 1）だけ。
証明の筋（報告参照）：係数はすべて非負なので打ち消しがなく、R の q 次数の最大値 D と最大値を与える項は
「漸化式の分岐 (5) の選び方（路）」の上の最大で決まる。補題 M：D(1v,1w) ≤ D(0v,0w) なら、q 次数が最大の路は
(5) で常に q の枝を選ぶ一本道（貪欲路）に限られ、その路では 1 の個数 l が常に 1 なので、項は (1+a) q^D t^{−D}。
ここで確かめること：
  M1  (10^{m−1},10^{n−1}) から到達する状態はすべて ''/全 0 の状態を経由しない（底は (∅,∅) だけ）
  M2  補題 M：到達する状態 (v0,w0) すべてで D(1v,1w) ≤ D(0v,0w)
  M3  貪欲路の q の枝の数 = δ、終点は (1,1)
  M4  漸化式で R を厳密に計算し、(IV_HM) が成り立つ
使い方：PYTHONPATH=src python scripts/framed_knot/d1_hm_iv.py（ログ logs/framed_knot/d1_hm_iv.log）。
"""
import sys
import time
from functools import lru_cache
from math import gcd

sys.setrecursionlimit(1000000)


def step(v, w):
    """漸化式の 1 段：[(係数の単項式 (di,dj,dk), 次の状態), ...]。底なら None。"""
    x, y, v0, w0 = v[-1], w[-1], v[:-1], w[:-1]
    l = v0.count("1")
    if x == "1" and y == "1":
        return [((0, 0, l), (v0, w0)), ((1, 0, 0), (v0, w0))]
    if x == "0" and y == "1":
        return [((0, 0, 0), (v0, "1" + w0))]
    if x == "1" and y == "0":
        return [((0, 0, 0), ("1" + v0, w0))]
    return [((0, 0, -l), ("1" + v0, "1" + w0)), ((0, 1, -l), ("0" + v0, "0" + w0))]


def is_bad(v, w):
    return (v == "") != (w == "") or (v and set(v) == {"0"} and set(w) == {"0"})


@lru_cache(None)
def poly(v, w):
    if v == "" and w == "":
        return {(0, 0, 0): 1}
    assert not is_bad(v, w), (v, w)
    out = {}
    for (di, dj, dk), st in step(v, w):
        for (i, j, k), c in poly(*st).items():
            key = (i + di, j + dj, k + dk)
            out[key] = out.get(key, 0) + c
    return out


@lru_cache(None)
def Dmax(v, w):
    if v == "" and w == "":
        return 0
    return max(dj + Dmax(*st) for (_, dj, _), st in step(v, w))


def reachable(m, n):
    seen, stack = set(), [("1" + "0" * (m - 1), "1" + "0" * (n - 1))]
    while stack:
        s = stack.pop()
        if s in seen or s == ("", ""):
            continue
        seen.add(s)
        stack += [st for _, st in step(*s)]
    return seen


def greedy(m, n):
    v, w, c = "1" + "0" * (m - 1), "1" + "0" * (n - 1), 0
    while not (v[-1] == "1" and w[-1] == "1"):
        nxt = step(v, w)
        _, (v, w) = nxt[-1]          # (5) では q の枝（2 番目）、それ以外は唯一の枝
        c += len(nxt) == 2
    return c, (v, w)


def main():
    log = open("logs/framed_knot/d1_hm_iv.log", "w")

    def out(s):
        print(s, flush=True)
        log.write(s + "\n")

    pairs = [(m, n) for m in range(2, 14) for n in range(2, 14) if gcd(m, n) == 1]
    m1 = m2 = True
    nst = 0
    for m, n in pairs:
        R = reachable(m, n)
        nst += len(R)
        m1 &= not any(is_bad(*s) for s in R)
        for v, w in R:
            if v[-1] == "0" and w[-1] == "0":
                m2 &= Dmax("1" + v[:-1], "1" + w[:-1]) <= Dmax("0" + v[:-1], "0" + w[:-1])
    out(f"M1/M2：2 ≤ m,n ≤ 13 の互いに素な {len(pairs)} 組、到達状態 {nst}：M1 {m1}、M2 {m2}")
    m3 = all(greedy(m, n) == ((m - 1) * (n - 1) // 2, ("1", "1"))
             for m in range(1, 41) for n in range(1, 41) if gcd(m, n) == 1)
    out(f"M3：1 ≤ m,n ≤ 40 の互いに素な組すべてで貪欲路の q の枝の数 = δ、終点 (1,1)：{m3}")
    for m, n in pairs:
        if m >= n:
            continue
        t0 = time.time()
        P = poly("1" + "0" * (m - 1), "1" + "0" * (n - 1))
        # 直接割り算：P = (1+a)R ⇔ R_i = P_i − R_{i−1}（a の低い方から）
        R = {}
        for key in sorted(P, key=lambda k: k[0]):
            i, j, k = key
            c = P[key] - R.get((i - 1, j, k), 0)
            if c:
                R[key] = c
        check = {}
        for (i, j, k), c in R.items():
            for kk in ((i, j, k), (i + 1, j, k)):
                check[kk] = check.get(kk, 0) + c
        delta = (m - 1) * (n - 1) // 2
        top = {k: c for k, c in R.items() if 2 * k[1] - k[0] >= 2 * delta}
        ok = check == P and all(c > 0 for c in R.values()) and top == {(0, delta, -delta): 1}
        out(f"M4：T({m},{n}) δ={delta}、項数 {len(R)}、(IV_HM) {ok}（{time.time() - t0:.1f} 秒）")
        poly.cache_clear()


if __name__ == "__main__":
    main()
