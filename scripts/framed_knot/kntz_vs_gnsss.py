"""2 箱の入力の基礎を KNTZ（1902.02275）に置くための照合：KNTZ の refined の式（src/rcsknots/double_twist_refined.poincare）と
GNSSS（1512.07883）の四重次数の式（scripts/framed_knot/t3_pretzel.quadruple）が、ダブルツイスト（pretzel）結び目 K_{m,n} の [r], [1^r]（r = 1, 2）で一致するか。

両方の式は m, n のそれぞれについて根 r + 1 個の指数和（符号の象限ごと）なので、各象限で |m|, |n| ∈ {1, 2, 3} の 3 × 3 の格子で一致すれば、
全ての (m, n) で一致する。[1^r] は GGS のミラー対称性 t_r ↔ t_c を介して比べる（tests/framed_knot/test_x1_refined.py と同じ）。
使い方：PYTHONPATH=src:scripts/framed_knot:scripts/dgr_lmov python scripts/framed_knot/kntz_vs_gnsss.py > logs/framed_knot/kntz_vs_gnsss.log
"""
import time

import sympy as sp

from rcsknots.double_twist_refined import poincare, tc, tr
from t3_pretzel import quadruple


def main():
    bad, n = [], 0
    for sm in (1, -1):
        for sn in (1, -1):
            for am in (1, 2, 3):
                for an in (1, 2, 3):
                    m, k = sm * am, sn * an
                    t0 = time.time()
                    for r in (1, 2):
                        g = sp.expand(sp.cancel(quadruple(m, k, r)))
                        ok1 = sp.expand(poincare(r, 1, m, k) - g) == 0
                        ok2 = sp.expand(poincare(1, r, m, k) - g.subs({tr: tc, tc: tr}, simultaneous=True)) == 0
                        n += 2
                        if not (ok1 and ok2):
                            bad.append((m, k, r, ok1, ok2))
                    print(f"K_{{{m},{k}}}：{'一致' if not bad or bad[-1][:2] != (m, k) else '不一致'}（{time.time() - t0:.1f}s）", flush=True)
    print(f"照合 {n} 件、不一致 {len(bad)} 件 {bad}")


if __name__ == "__main__":
    main()
