"""open_problems D1 (IV)：Hogancamp–Mellit（arXiv:1909.00418）の漸化式（Def. 1.1 = def:the polys、Thm. 1.2/main thm）で
トーラス結び目の reduced HOMFLY-PT ホモロジーの Poincaré 級数を計算し、
  H1  本研究の [1] の値（data/framed_knot/torus/*_1.json、KN の変数）との辞書：rCS_[1](a,q,t) = a^δ (qt)^{δ/2} · P_HM(−a, 1/q, t)
  H2  (IV) の HM 版：P_HM の各項 a^i q^j t^k について 2j − i ≤ 2δ、等号は (i,j) = (0,δ) の 1 項だけ（係数 1）
      （辞書の下で KN の h = x + 2y は 2δ + i − 2j。よって H2 ⇔「最低 𝐭 次数 0 の項は最低 a 次数の 1 項だけ」）
を確かめる。reduced は P_HM(0^m,0^n)·(1−q)/(1+a)（HM Remark：H_KR(K) = H^red ⊗ R[x,θ]）。
使い方：PYTHONPATH=src python scripts/framed_knot/d1_hm.py（ログ logs/framed_knot/d1_hm.log）。
"""
import glob
import json
import sys
from functools import lru_cache
from math import gcd

import sympy as sp

q, t, a = sp.symbols("q t a")
sys.setrecursionlimit(10000)


@lru_cache(None)
def P(v, w):
    """HM Def. (def:the polys)。v, w は '0'/'1' の文字列。"""
    if v == "":
        return ((1 + a) / (1 - q)) ** len(w)
    if w == "":
        return ((1 + a) / (1 - q)) ** len(v)
    if set(v) == {"0"} and set(w) == {"0"}:
        return sp.cancel(P("1" + v[1:], "1" + w[1:]) / (1 - q))      # 式 (eq:allzeroes)
    x, y, v0, w0 = v[-1], w[-1], v[:-1], w[:-1]
    l = v0.count("1")
    if x == "1" and y == "1":
        return sp.expand((t ** l + a) * P(v0, w0))
    if x == "0" and y == "1":
        return P(v0, "1" + w0)
    if x == "1" and y == "0":
        return P("1" + v0, w0)
    return sp.expand(t ** (-l) * P("1" + v0, "1" + w0) + q * t ** (-l) * P("0" + v0, "0" + w0))


def reduced(m, n):
    red = sp.expand(sp.cancel(P("0" * m, "0" * n) * (1 - q) / (1 + a)))
    out = {}
    for term in sp.Add.make_args(red):
        c, mono = term.as_coeff_Mul()
        d = mono.as_powers_dict()
        out[(int(d.get(a, 0)), int(d.get(q, 0)), int(d.get(t, 0)))] = int(c)
    return out


def main():
    log = open("logs/framed_knot/d1_hm.log", "w")

    def out(s):
        print(s, flush=True)
        log.write(s + "\n")

    out("H1：本研究の [1] の値との辞書")
    for f in sorted(glob.glob("data/framed_knot/torus/*_1.json")):
        d = json.load(open(f))
        m, n = d["torus"]
        delta = (m - 1) * (n - 1) // 2
        kn = {tuple(map(int, k.split(","))): v for k, v in d["poly"].items()}   # (A, Q, T) = (a, q, t)^{1/2} の指数
        img = {}
        for (i, j, k), c in reduced(m, n).items():
            key = (2 * delta + 2 * i, delta - 2 * j, delta + 2 * k)
            img[key] = img.get(key, 0) + (-1) ** i * c
        out(f"  T({m},{n})：全係数一致 {img == kn}")
    out("H2：2j − i ≤ 2δ、等号は (0,δ) の 1 項（係数 1）")
    for m in range(2, 7):
        for n in range(m + 1, 10):
            if gcd(m, n) != 1 or (m - 1) * (n - 1) > 40:
                continue
            delta = (m - 1) * (n - 1) // 2
            H = reduced(m, n)
            top = {k: c for k, c in H.items() if 2 * k[1] - k[0] >= 2 * delta}
            ok = len(top) == 1 and all(k[0] == 0 and k[1] == delta and c == 1 for k, c in top.items())
            out(f"  T({m},{n})：項数 {len(H)}、係数はすべて正 {all(c > 0 for c in H.values())}、H2 {ok} {top}")


if __name__ == "__main__":
    main()
