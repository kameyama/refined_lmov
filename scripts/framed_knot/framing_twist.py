"""ツイスト結び目 K_p（GNSSS の予想式）の 2 箱：全ての (p, framing x) での正値性を、2 次元のニュートン形の有限個の判定に落とす。

構造【導出】：S(p, x) = 𝐭^{−2x} \\hat f^{(x)}_ρ(K_p) は p について根 m ∈ {1, X, X², Y, Ȳ}（p < 0 は逆数、REPORT_phase4.md §5.1）、
x = ε + 2δy について根 r ∈ {1, 𝐪⁴𝐭⁴, 𝐪⁻⁴}（δ = −1 は逆数、scripts/framed_knot/framing_newton.py）の積 m^{|p|} r^{y} の一次結合。
テンソル積のニュートン形 S = Σ_{k,l} c_{kl} h_{|p|−1−k}(m の節点) h_{y−y0−l}(r の節点) の係数 c_{kl} を |p| = 1..5、y = y0..y0+2 から作り、
全ての c_{kl} が同じ符号なら、その (p の符号, ε, δ) の |p| ≥ 1、y ≥ y0 の全てで G 展開の前が一様な符号【導出】。
σ 不変性も同じ格子の点から全体に従う（積の根の張る空間で、格子の点での値が全体を決める）【導出】。
起点 y0 と 0 の間の x は、p についての 1 次元のニュートン形（V4 と同じ）で別に確かめる。
入力の作り方と 2 箱の公式は scripts/framed_knot/v4_closed.py（hat_from_inputs）と同じで、framing の単項式を入力に掛ける。
使い方：python scripts/framed_knot/framing_twist.py [YSHIFT_MAX]　または　python scripts/framed_knot/framing_twist.py direct（|p| ≤ 6、|x| ≤ 14 の直接の判定）
"""
import sys
import time
from itertools import permutations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from framing_newton import ROOTS, newton, sigma, uniform_sign  # noqa: E402
from t4_colored import bold_to_small  # noqa: E402
from v4_closed import L, RHOS, formula, inputs, load_inputs, mir, roots_for, ORDER  # noqa: E402
from rcsknots.bps import to_bold  # noqa: E402

FRAME = {"r1": (-1, (0, 1, -1)), "r2": (1, (0, 4, -2)), "r11": (1, (0, 2, -4))}   # framing[λ] の (符号, (A, Q, T) の指数)


def times_mono(d, sign, e, x):
    return {(a + x * e[0], q + x * e[1], t + x * e[2]): (sign ** x) * c for (a, q, t), c in d.items()}


def hat_framed(P1, P2, x):
    """framing x の \\hat f_[2], \\hat f_[1,1]（太字、𝐭^{−2x} 倍、{(e_𝐚, e_𝐪, e_𝐭): c}）。"""
    r1 = times_mono(bold_to_small(P1), *FRAME["r1"], x)
    r2 = times_mono(bold_to_small(P2), *FRAME["r2"], x)
    r11 = times_mono(bold_to_small(mir(P2)), *FRAME["r11"], x)
    R = [L.from_dict(r2), L.from_dict(r11), L.from_dict(r1) * L.from_dict(r1),
         L.from_dict({(2 * i, 2 * j, 2 * k): c for (i, j, k), c in r1.items()})]
    out = {}
    for rho, (ns, den) in formula().items():
        num = ns[0] * R[0] + ns[1] * R[1] + ns[2] * R[2] + ns[3] * R[3]
        out[rho] = {(a, q, t - 2 * x): c for (a, q, t), c in to_bold((num / den).to_dict(), 0).items()}
    return out


def grid_ok(cache, sp, eps, delta, y0, log):
    """格子 |p| = 1..5、y = y0..y0+2 のテンソル積ニュートン形。ρ ごとに (y の並び, 符号, σ 不変) か None。"""
    xs = [eps + delta * 2 * (y0 + y) for y in range(3)]
    vals = {(n, x): hat_framed(*inputs(sp * n, cache), x) for n in range(1, 6) for x in xs}
    proots = roots_for(sp)
    rroots = [r if delta == 1 else tuple(-e for e in r) for r in ROOTS[2]]
    out = {}
    for rho in RHOS:
        closed = all(sigma(v[rho]) == v[rho] and all(e[1] % 2 == 0 for e in v[rho]) for v in vals.values())
        ck = {x: newton([vals[(n, x)][rho] for n in range(1, 6)], proots) for x in xs}        # p についての c_k(x)
        res = None
        for perm in permutations(range(3)):
            order = [rroots[i] for i in perm]
            ckl = [c for k in range(5) for c in newton([ck[x][k] for x in xs], order)]
            sg = uniform_sign(ckl)
            if sg:
                res = (perm, sg, closed)
                break
        out[rho] = res
    return out, xs


def gap_ok(cache, sp, x):
    """framing x を固定し、p についての 1 次元のニュートン形（V4 の並び）で全ての |p| ≥ 1 を判定。"""
    vals = [hat_framed(*inputs(sp * n, cache), x) for n in range(1, 6)]
    return {rho: (uniform_sign(newton([v[rho] for v in vals], roots_for(sp))),
                  all(sigma(v[rho]) == v[rho] for v in vals)) for rho in RHOS}


def main(yshift_max=4, log=print):
    cache = load_inputs()
    signs = set()
    for sp in (1, -1):
        for eps in (0, 1):
            for delta in (1, -1):
                for y0 in range(yshift_max + 1):
                    t0 = time.time()
                    res, xs = grid_ok(cache, sp, eps, delta, y0, log)
                    if all(res[r] for r in RHOS):
                        break
                log(f"p の符号 {sp:+d}、ε={eps}、向き {delta:+d}：起点 x0 = {xs[0]}（x = {xs}）"
                    + "、".join(f" ρ={list(r)} {'並び ' + str(res[r][0]) + ' 符号 ' + str(res[r][1]) + ' σ 不変 ' + str(res[r][2]) if res[r] else '一様な符号の並び無し'}" for r in RHOS)
                    + f"（{time.time() - t0:.1f}s）")
                signs |= {res[r][1] for r in RHOS if res[r]}
                for x in range(eps, xs[0], 2 * delta):                    # 起点と 0 の間（同じ偶奇）
                    g = gap_ok(cache, sp, x)
                    log(f"   間の framing x = {x}："
                        + "、".join(f"ρ={list(r)} 符号 {g[r][0]} σ 不変 {g[r][1]}" for r in RHOS))
                    signs |= {g[r][0] for r in RHOS}
    log(f"現れた全体の符号 {sorted(signs)}（0 は一様でない）")


def direct(pmax=6, xmax=14, log=print):
    """直接の判定：|p| ≤ pmax、|x| ≤ xmax の全てで G 展開の前後とも一様な符号か。"""
    from rcsknots.bps import g_expand
    cache = load_inputs()
    bad = tot = 0
    for p in [s * n for n in range(1, pmax + 1) for s in (1, -1)]:
        for x in range(-xmax, xmax + 1):
            for rho, P in hat_framed(*inputs(p, cache), x).items():
                tot += 1
                sg = uniform_sign([P])
                try:
                    ok = sg != 0 and all(v >= 0 for v in g_expand({k: sg * v for k, v in P.items()}).values())
                except ValueError:
                    ok = False
                if not ok:
                    bad += 1
                    log(f"   負または閉じない：p = {p}, x = {x}, ρ = {list(rho)}")
    log(f"直接の判定：|p| ≤ {pmax}、|x| ≤ {xmax} の {tot} 件中、一様でないもの {bad} 件")


if __name__ == "__main__":
    if sys.argv[1:2] == ["direct"]:
        direct(log=lambda s: print(s, flush=True))
    else:
        main(int(sys.argv[1]) if len(sys.argv) > 1 else 4, log=lambda s: print(s, flush=True))
