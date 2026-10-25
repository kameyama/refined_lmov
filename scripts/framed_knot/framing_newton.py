"""framing を記号として扱う：\\hat f^{(x)}_ρ の x についての指数和とニュートン形で、全ての framing x での正値性を有限個の判定に落とす。

構造【導出】：|ρ| = n の \\hat f_ρ は、入力 rCS_λ（|λ| ≤ n）の積と Adams 作用の写しの、結び目によらない係数による一次結合。
framing x は各入力に framing[λ]^x を掛ける（conventions.md §17）。太字の変数では framing[□] = 𝐭、
framing[□](q², t²) = −𝐭²、framing[□](q³, t³) = 𝐭³、framing[2] = 𝐪²𝐭⁴、framing[1,1] = 𝐪⁻²𝐭²、
framing[3] = 𝐪⁶𝐭⁹、framing[2,1] = 𝐭⁵、framing[1,1,1] = 𝐪⁻⁶𝐭³。各項の因子は 𝐭^{nx} × r^x で、
  n = 2：r ∈ {𝐪²𝐭², 𝐪⁻², 1, −1}
  n = 3：r ∈ {𝐪⁶𝐭⁶, 𝐪⁻⁶, 𝐪²𝐭², 𝐪⁻², 𝐭², 1, −1}
よって S(x) = 𝐭^{−nx} \\hat f^{(x)} は x の指数和。x = ε + 2y（ε = 0, 1）とすれば根は r² で、全て係数 1 の単項式
（n = 2：{1, 𝐪⁴𝐭⁴, 𝐪⁻⁴}、n = 3：{1, 𝐪¹²𝐭¹², 𝐪⁻¹², 𝐪⁴𝐭⁴, 𝐪⁻⁴, 𝐭⁴}）。y ≥ 0 の向きは根そのもの、y ≤ 0 の向きは逆数を使う。
ニュートン形 S(ε + 2δy) = Σ_k c_k h_{y−k}(根_0..根_k)（δ = ±1）。h は単項式の完全斉次対称式で非負係数なので、
c_k が全て同じ符号なら、その (ε, δ) の全ての y で \\hat f^{(x)} は G 展開の前で一様な符号【導出】。
σ（G 展開が閉じる条件の対称性）は根の集合を保つので、根の個数の点で σ 不変なら全ての y で σ 不変【導出】。
すると命題 1 の補題により \\hat N も非負（REPORT_phase4.md §5.1）。
使い方：python scripts/framed_knot/framing_newton.py 結び目 N [起点のずらしの上限]（値は data/framed_knot/framing/ に保存して再利用）
"""
import json
import sys
import time
from itertools import permutations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from framing_check import framed, hat_f_nonzero, load_reduced, zero_rhos  # noqa: E402
from rcsknots.bps import g_expand, to_bold  # noqa: E402
from rcsknots.extract import A_SHIFT  # noqa: E402
from rcsknots.partitions import partitions  # noqa: E402
from v4_closed import L, mono  # noqa: E402

ROOTS = {2: [(0, 0, 0), (0, 4, 4), (0, -4, 0)],
         3: [(0, 0, 0), (0, 12, 12), (0, -12, 0), (0, 4, 4), (0, -4, 0), (0, 0, 4)]}   # r² の (e_𝐚, e_𝐪, e_𝐭)
OUT = Path("data/framed_knot/framing")


def S_values(name, n, x, reduced):
    """{ρ: S(x) = 𝐭^{−nx} \\hat f^{(x)}_ρ（太字、nontorus の a 次数のずらし込み）}。data/framed_knot/framing/ に保存。"""
    path = OUT / f"{name.replace(':', '')}_{n}_{x}.json"
    if path.exists():
        d = json.loads(path.read_text())
        return {tuple(json.loads(k)): {tuple(map(int, e.split(","))): c for e, c in v.items()} for k, v in d.items()}
    R = framed(reduced, x)
    zeros = zero_rhos(R, n)
    hats = hat_f_nonzero(R, n, [r for r in partitions(n) if r not in zeros])
    out = {}
    for rho in partitions(n):
        if rho in zeros:
            out[rho] = {}
            continue
        assert hats[rho] is not None, (name, n, x, rho)
        bold = to_bold(hats[rho], A_SHIFT["nontorus"][n])
        out[rho] = {(ea, eq, et - n * x): c for (ea, eq, et), c in bold.items()}
    OUT.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({json.dumps(list(r)): {",".join(map(str, e)): c for e, c in v.items()} for r, v in out.items()}))
    return out


def newton(seq, order):
    cur = [L.from_dict(s) for s in seq]
    cs = []
    for k in range(len(order)):
        cs.append(cur[0].to_dict())
        cur = [cur[i + 1] - mono(order[k]) * cur[i] for i in range(len(cur) - 1)]
    return cs


def uniform_sign(cs):
    vals = [v for c in cs for v in c.values() if v]
    if all(v > 0 for v in vals):
        return 1
    if all(v < 0 for v in vals):
        return -1
    return 0


def direct_sign(P):
    """直接の判定：G 展開の前の全体の符号（一様でなければ 0）と、その符号で G 展開した後に負が無いか。"""
    if not P:
        return ("0", True)
    sg = uniform_sign([P])
    if not sg:
        return (0, None)
    try:
        N = g_expand({k: sg * v for k, v in P.items()})
    except ValueError:
        return (sg, "閉じない")
    return (sg, all(v >= 0 for v in N.values()))


def sigma(P):
    return {(a, -q, t - q): c for (a, q, t), c in P.items()}


def run(name, n, max_shift=3, log=print):
    """(ε, 向き) ごとに、起点 x0 = ε + 2δb（b = 0..max_shift）を順に試し、c_k が一様な符号になる最初の起点を探す。
    起点より 0 に近い x（同じ偶奇）は直接の計算で確かめる必要がある（framing_batch の結果、または S_values で判定）。"""
    reduced = load_reduced(name, n)
    K = len(ROOTS[n])
    result = {}
    for eps in (0, 1):
        for delta in (1, -1):
            roots = [r if delta == 1 else tuple(-e for e in r) for r in ROOTS[n]]
            for rho in partitions(n):
                found = None
                for b in range(max_shift + 1):
                    t0 = time.time()
                    xs = [eps + delta * 2 * (b + y) for y in range(K)]
                    vals = {x: S_values(name, n, x, reduced) for x in xs}
                    seq = [vals[x][rho] for x in xs]
                    if all(not s for s in seq):
                        found = (xs[0], None, 0, True)
                        break
                    closed = all(sigma(s) == s for s in seq) and all(e[1] % 2 == 0 for s in seq for e in s)
                    for perm in permutations(range(K)):
                        sg = uniform_sign(newton(seq, [roots[i] for i in perm]))
                        if sg:
                            found = (xs[0], perm, sg, closed)
                            break
                    if found:
                        break
                gap = []
                if found:                                               # 起点と 0 の間（同じ偶奇）を直接判定する
                    for x in range(eps, found[0], 2 * delta):
                        gap.append((x, direct_sign(S_values(name, n, x, reduced)[rho])))
                found = found and found + (gap,)
                result[(eps, delta, rho)] = found
                msg = (f"起点 x0 = {found[0]}、並び {found[1]}、符号 {found[2]}、σ 不変 {found[3]}、間の直接判定 {found[4]}" if found and found[1] is not None
                       else ("恒等的に 0" if found else f"起点を {max_shift} 回ずらしても一様な符号の並び無し"))
                log(f"{name} |ρ|={n} ε={eps} 向き {delta:+d} ρ={list(rho)}：{msg}（{time.time() - t0:.1f}s）")
    return result


if __name__ == "__main__":
    run(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 3, log=lambda s: print(s, flush=True))
