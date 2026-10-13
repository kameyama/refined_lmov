"""別研究 [2,1] の X0：ダブルツイスト結び目 K_{m,n} の t = q の彩色 HOMFLY-PT を、文献の式とケーブルの逆解きで全係数照合する。

入力：data/framed_knot/reference/x0_cable_p2.txt（s = 1, 3, 5）、x0_cable_p3.txt（s = 1, 2, 4, 5 と、|p| ≥ 5 では 3w ± 1, 3w ± 2）、x0_cable_p22.txt（(s, s') ∈ {1,3,5}×{1,3}）
      （scripts/framed_knot/x0_double_twist.sage）。逆解きは colored_from_cables.colored_fast（解くのに使わない傾きでデータとの一致も確かめる）。
照合：
  [2], [1,1], [3], [1,1,1], [2,2], [4], [1,1,1,1]：KNTZ 1902.02275 の長方形の式（rcsknots.double_twist.H_rect）
  [2,1]：Morozov 1612.00422 の (H21de)（KNTZ と同じ (m, n)）と、1903.00259 の 𝓑^{[2,1]}（Z を F^{(1)}F^{(−1)} で割って読む。conventions.md §14）
変数は A² = a, q = Q（rcsknots.double_twist.to_kn）。
出力：data/framed_knot/x0/colored_t_eq_q.txt（名前<TAB>λ<TAB>P_λ（reduced、framing 0、KN の (a, Q)））、ログは標準出力。
使い方：python scripts/framed_knot/x0_check.py [p2] [p3] [p22]
"""
import os
import sys
import time

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from colored_from_cables import SUSPECT, load, colored_fast  # noqa: E402
from rcsknots.double_twist import H_rect, H21_1612, H21_1903, to_kn  # noqa: E402

OUT = "data/framed_knot/x0/colored_t_eq_q.txt"


def mn(name):
    return tuple(int(x.replace("m", "-")) for x in name[1:].split("_"))


def reference(lam, m, n):
    """λ の文献の値（KN の (a, Q)）。[2,1] は (1612 の値, 1903 の値)。"""
    if lam == (2, 1):
        return to_kn(H21_1612(m, n)), to_kn(H21_1903(m, n))
    if len(set(lam)) == 1:
        return (to_kn(H_rect(lam[0], len(lam), m, n)),)
    return None


def main(kinds):
    os.makedirs("data/framed_knot/x0", exist_ok=True)
    saved = set()
    if os.path.exists(OUT):
        saved = {tuple(line.split("\t")[:2]) for line in open(OUT)}
    p_of = {"p2": 2, "p3": 3, "p22": 4}
    bad = 0
    with open(OUT, "a") as out:
        for kind in kinds:
            data = load(f"data/framed_knot/reference/x0_cable_{kind}.txt")
            print(f"{kind}：あふれの疑いで除いた（名前, 傾き）{SUSPECT}", flush=True)
            SUSPECT.clear()
            names = sorted({n for (n, s) in data}, key=lambda x: (abs(mn(x)[1]), abs(mn(x)[0]), mn(x)))
            for name in names:
                m, n = mn(name)
                t0 = time.time()
                col, checked = colored_fast(name, data, p_of[kind])
                t1 = time.time() - t0
                msg = []
                for lam, val in col.items():
                    if (name, str(lam)) not in saved:
                        out.write(f"{name}\t{lam}\t{val}\n")
                        out.flush()
                    t0 = time.time()
                    ref = reference(lam, m, n)
                    if ref is None:
                        msg.append(f"{lam}:—")
                        continue
                    ok = [sp.expand(r - val) == 0 for r in ref]
                    bad += not all(ok)
                    msg.append(f"{lam}:{'/'.join('一致' if x else '不一致' for x in ok)}({len(sp.Add.make_args(val))}項,{time.time() - t0:.1f}s)")
                print(f"K_{{{m},{n}}} {kind} 逆解き {t1:.1f}s（検証の傾き {checked}）", " ".join(msg), flush=True)
    print("不一致の数", bad)


if __name__ == "__main__":
    main(sys.argv[1:] or ["p2", "p22", "p3"])
