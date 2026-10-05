"""ステージ T3：(2s−1, 1, 2t−1)-pretzel 結び目 K_{s,t} の 2 箱の refined LMOV 不変量。

入力：Gukov–Nawata–Saberi–Stošić–Sułkowski（arXiv:1512.07883）の四重次数の superpolynomial (superpoly-pretzel), (twisting3)
（TeX references/1512.07883/src/BPS-SS.tex の 3186-3205 行）。予想に基づく式である。s = 1 はツイスト結び目。
t_c 次数の三重次数：P(a, q, t) = \\tilde P(a, Q = q, t_r = 1/q, t_c = t q)（t_r 次数は t_r = t/q, t_c = q と原典に明記。
各生成元の q 次数は Q − t_r + t_c なので、同じ規則で t_c を残す）。a, q は DGR の太字。
使い方：python scripts/t3_pretzel.py S_MIN S_MAX T_MIN T_MAX
"""
import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from t4_colored import judge2  # noqa: E402

a, Q, tr, tc, q, t = sp.symbols("a Q t_r t_c q t")


def qp(x, k, base):
    out = 1
    for i in range(k):
        out *= 1 - x * base ** i
    return out


def qbinom(n, k, base):
    return sp.cancel(qp(base, n, base) / (qp(base, k, base) * qp(base, n - k, base)))


def varpi(p, i):
    s = 0
    for j in range(i + 1):
        s += (-1) ** j * (a * tr * tc) ** (2 * p * j) * tc ** ((2 * p + 1) * j * (j - 1)) \
            * (1 - a ** 2 * tr ** 2 * tc ** (4 * j)) / qp(a ** 2 * tr ** 2 * tc ** (2 * j), i + 1, tc ** 2) * qbinom(i, j, tc ** 2)
    return a ** i * tc ** sp.Rational(i * (i - 1), 2) * s


def quadruple(s_, t_, r):
    sgn = (sp.sign(s_) + sp.sign(t_)) / 2
    total = 0
    for i in range(r + 1):
        total += varpi(s_, i) * varpi(t_, i) * (-a ** 2 * tc ** (2 * (r - 1))) ** (-i) \
            * qp(-a ** 2 * Q ** -2 * tr * tc, i, tc ** 2) * qp(-a ** 2 * Q ** 2 * tr ** 3 * tc ** (1 + 2 * r), i, tc ** 2) \
            * qbinom(r, i, tc ** 2)
    return sp.factor(sp.cancel(sp.together((-tr * tc) ** (-r * sgn) * total)))


def tc_graded(s_, t_, r):
    """t_c 次数の Poincaré 多項式 {(e_𝐚, e_𝐪, e_𝐭): c}。"""
    e = sp.expand(sp.cancel(quadruple(s_, t_, r).subs({Q: q, tr: 1 / q, tc: t * q}, simultaneous=True)) * (a * q * t) ** 300)
    return {(i - 300, j - 300, k - 300): int(c) for (i, j, k), c in sp.Poly(e, a, q, t).terms()}


def run(s_, t_):
    P1, P2 = tc_graded(s_, t_, 1), tc_graded(s_, t_, 2)
    P11 = {(ea, -eq, et - eq): c for (ea, eq, et), c in P2.items()}      # ミラー対称性（DGR 型、conventions.md §12）
    return f"K_{{{s_},{t_}}}: rank {sum(P1.values())}, {sum(P2.values())}: {judge2(P1, P2, P11)}"


if __name__ == "__main__":
    s0, s1, t0, t1 = map(int, sys.argv[1:5])
    for s_ in range(s0, s1 + 1):
        for t_ in range(max(t0, s_), t1 + 1):         # K_{s,t} = K_{t,s}
            if s_ and t_:
                print(run(s_, t_), flush=True)
