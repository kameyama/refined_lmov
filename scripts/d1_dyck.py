"""open_problems D1 (IV)：Mellit "Homology of torus knots"（arXiv:1704.07630, khro.tex:573 Thm (dyck paths)）の
Dyck 路の式と本研究の [1] の値・HM の漸化式の照合。

Mellit の式：𝒫_{m,n} = (a(qt)^{−1/2})^δ (1/(1−t)) Σ_π t^{area(π)} q^{h₊(π)} Π_{v∈v*(π)} (1 − a q^{−k(v)})。
  π：(m,n)-Dyck 路（y = (n/m)x の上）、area：路と対角線の間の正方形の数、h₊：水平の歩み s₁（左）と垂直の歩み s₂ の組で
  対角線に平行な直線が両方と交わるものの数、v*(π)：外側の頂点（N の次が E）から対角線に最も遠いものを除いた集合、
  k(v)：v を通る対角線に平行な直線が交わる垂直の歩みの数（v で終わる歩みも数える）。
本研究の規約への辞書（照合で確かめるもの）：
  rCS_[1](a,q,t) = Σ_π Σ_{S⊆v*(π)} (−1)^{|S|} a^{δ+|S|} q^{δ/2−area(π)} t^{h₊(π)−δ/2+Σ_{v∈S}(1−k(v))}
  ＝ t^{−δ} · [(1 − t_M)𝒫_{m,n}]（a_M = a t, q_M = t, t_M = 1/q）。
  このとき h = (a の次数) + 2(q の次数) = 2(δ − area(π)) + |S| ≥ 0。等号は area = δ（π = N^n E^m、v* = ∅）かつ S = ∅ の 1 項だけ。
照合：
  D1  保存済みのトーラス結び目 9 個の [1]（data/torus/*_1.json）と全係数一致
  D2  新しいトーラス結び目で GN の正規化前の値 × M_□（conventions.md §13）と全係数一致
  D3  HM の漸化式（scripts/d1_hm_iv.py）と、HM の (a,q,t) の指数 = (|S|, area, |S| + h₊ − Σk − δ) で一致（m,n ≤ 9）
使い方：PYTHONPATH=src:scripts python scripts/d1_dyck.py（ログ logs/d1_dyck.log）。
"""
import glob
import json
import time
from fractions import Fraction as F
from itertools import product
from math import gcd

from d1_hm_iv import poly as hm_poly


def paths(m, n):
    def rec(x, y):
        if (x, y) == (m, n):
            yield ""
            return
        if y < n:
            for r in rec(x, y + 1):
                yield "N" + r
        if x < m and m * y >= n * (x + 1):
            for r in rec(x + 1, y):
                yield "E" + r
    return list(rec(0, 0))


def stats(p, m, n):
    pts = [(0, 0)]
    for s in p:
        x, y = pts[-1]
        pts.append((x + 1, y) if s == "E" else (x, y + 1))
    hor = [pts[i] for i, s in enumerate(p) if s == "E"]
    ver = [pts[i] for i, s in enumerate(p) if s == "N"]

    def c(x, y):
        return F(y) - F(n * x, m)
    area = sum(y + (n * (x + 1)) // -m for x, y in hor)          # y − ⌈n(x+1)/m⌉
    hp = sum(1 for x1, y1 in hor for x2, y2 in ver
             if x1 + 1 <= x2 and max(c(x1 + 1, y1), c(x2, y2)) < min(c(x1, y1), c(x2, y2 + 1)))
    outer = [pts[i] for i in range(1, len(p)) if p[i - 1] == "N" and p[i] == "E"]
    far = max(outer, key=lambda v: c(*v))
    ks = [sum(1 for x2, y2 in ver if c(x2, y2) <= c(*v) <= c(x2, y2 + 1)) for v in outer if v != far]
    return area, hp, ks


def terms(m, n):
    """[(area, h₊, |S|, Σ_S k)] with multiplicity."""
    out = []
    for p in paths(m, n):
        area, hp, ks = stats(p, m, n)
        for S in product([0, 1], repeat=len(ks)):
            out.append((area, hp, sum(S), sum(k for k, s in zip(ks, S) if s)))
    return out


def kn_form(m, n):
    """{(2·e_a, 2·e_q, 2·e_t): 係数}（A = a^{1/2} などの指数、data/ と同じ）。"""
    delta = (m - 1) * (n - 1) // 2
    out = {}
    for area, hp, s, sk in terms(m, n):
        key = (2 * (delta + s), delta - 2 * area, 2 * (hp - sk + s) - delta)
        out[key] = out.get(key, 0) + (-1) ** s
    return {k: v for k, v in out.items() if v}


def main():
    log = open("logs/d1_dyck.log", "w")

    def out(s):
        print(s, flush=True)
        log.write(s + "\n")

    for f in sorted(glob.glob("data/torus/*_1.json")):
        d = json.load(open(f))
        m, n = d["torus"]
        kn = {tuple(map(int, k.split(","))): v for k, v in d["poly"].items()}
        K = kn_form(m, n)
        delta = (m - 1) * (n - 1) // 2
        hmin = min(k[0] // 2 + k[1] for k in K)       # h = e_a + 2 e_q（半指数では A/2 + Q）
        bottom = {k: v for k, v in K.items() if k[0] // 2 + k[1] == hmin}
        out(f"D1 T({m},{n})：全係数一致 {K == kn}、最小の h = {hmin}、その項 {bottom}（期待 {{({2 * delta}, {-delta}, {-delta}): 1}}）")

    from rcsknots.exact import reconstruct
    for m, n in [(2, 13), (3, 8), (4, 7), (5, 6), (3, 10), (5, 7)]:
        t0 = time.time()
        raw = reconstruct(lambda G: G.torus_reduced_apoly((1,), m, -n), log=lambda s: None)
        delta = (m - 1) * (n - 1) // 2
        sign, ea, eq = (-1) ** (m - 1), F(delta), F(delta, 2) + m - 1      # M_□ = (−q/t)^{m−1}(a(q/t)^{1/2})^δ
        norm = {}
        for (a_, q_, t_), c in raw.items():
            key = (2 * (a_ + ea), 2 * (q_ + eq), 2 * (t_ - eq))
            norm[tuple(int(x) for x in key)] = sign * c
        out(f"D2 T({m},{n})：GN × M_□ = Dyck 路の式 {norm == kn_form(m, n)}（{time.time() - t0:.1f} 秒）")

    ok = True
    for m in range(2, 10):
        for n in range(2, 10):
            if gcd(m, n) != 1:
                continue
            delta = (m - 1) * (n - 1) // 2
            D = {}
            for area, hp, s, sk in terms(m, n):
                key = (s, area, s + hp - sk - delta)
                D[key] = D.get(key, 0) + 1
            P = hm_poly("1" + "0" * (m - 1), "1" + "0" * (n - 1))
            R = {}
            for key in sorted(P, key=lambda k: k[0]):
                i, j, k = key
                c = P[key] - R.get((i - 1, j, k), 0)
                if c:
                    R[key] = c
            ok &= R == D
            hm_poly.cache_clear()
    out(f"D3 2 ≤ m,n ≤ 9：Mellit の式と HM の漸化式が全係数一致 {ok}")


if __name__ == "__main__":
    main()
