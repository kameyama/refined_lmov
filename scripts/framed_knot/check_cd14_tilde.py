"""(N1)(N2) の [1] と CD14 の tilde 正規化（1408.4348 Conj. CONCONJ (ii)、cable.tex:1516-1531）の対応を、CD14 の明示例で確かめる。
本研究の太字の P と、CD14 の DAHA superpolynomial（DAHA の変数、a_DAHA = −a_KN）を太字に直したものの比を求める。
使い方：PYTHONPATH=src python scripts/framed_knot/check_cd14_tilde.py > logs/framed_knot/check_cd14_tilde.log
"""
import json, sympy as sp
from pathlib import Path
from rcsknots.daha_compare import daha_to_kn, q, t, a
from rcsknots.bps import to_bold

def stored(k, c):
    return {tuple(map(int, x.split(","))): v for x, v in json.loads(Path(f"data/framed_knot/cables/{k}_{c}.json").read_text())["poly"].items()}

def sym_ok(P):   # DGR の対称性 (e_a, e_q, e_t) ↔ (e_a, −e_q, e_t − e_q) と 𝐪 の指数が偶数
    return all(P.get((ea, -eq, et - eq)) == c for (ea, eq, et), c in P.items()), all(eq % 2 == 0 for _, eq, _ in P)

def ratio(A, B):
    ka, kb = max(A), max(B)
    sh = tuple(x - y for x, y in zip(ka, kb)); s = A[ka] // B[kb]
    return (s, sh) if {tuple(k[i] + sh[i] for i in range(3)): s * v for k, v in B.items()} == A else None

for line in open("data/framed_knot/reference/daha_examples.txt"):
    name, knot, color, expr = line.rstrip("\n").split("\t")
    if color != "1":
        continue
    ours = to_bold(stored(knot, "1"), 0)
    th = to_bold(daha_to_kn(sp.sympify(expr, locals={"q": q, "t": t, "a": a})), 0)
    amin = min(k[0] for k in th)
    low = {k: v for k, v in th.items() if k[0] == amin}
    print(f"{name:7s} {knot:4s} 本研究 P：対称性・𝐪偶数 {sym_ok(ours)}、最低𝐭 {min(k[2] for k in ours)}、負 {sum(v<0 for v in ours.values())} ｜ "
          f"CD14 tilde：負 {sum(v<0 for v in th.values())}、最低𝐚の行に 𝐪⁰𝐭⁰ の項 {low.get((amin,0,0))}、最低𝐚の行の負 {sum(v<0 for v in low.values())} ｜ 比 本研究/CD14 = {ratio(ours, th)}")
