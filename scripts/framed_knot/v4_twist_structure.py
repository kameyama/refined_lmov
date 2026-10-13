"""第 4 期 V4：ツイスト結び目 K_p の \\hat f_[2], \\hat f_[1,1] の p についての構造。

p > 0 で \\hat f(K_p) = Σ_j A_j m_j^p（m_j は NRZS の式のねじれの因子 (a t² q^{ℓ−1})^ℓ から来る単項式）と書けるかを調べる。
太字の \\hat f を p = 1..8 で計算し、m_j を候補から選んで p = 1..4 で A_j を解き、p = 5..8 で確かめる。
"""
import sys
from itertools import combinations
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from t3_pretzel import tc_graded  # noqa: E402
from t4_colored import bold_to_small  # noqa: E402
from rcsknots.bps import to_bold  # noqa: E402
from rcsknots.extract import A_SHIFT, hat_f_exact  # noqa: E402

a, q, t = sp.symbols("a q t")    # 太字 (𝐚, 𝐪, 𝐭)


def mir(P):
    return {(x, -y, z - y): c for (x, y, z), c in P.items()}


def hat_bold(p):
    P1, P2 = tc_graded(1, p, 1), tc_graded(1, p, 2)
    red = {(1,): bold_to_small(P1), (2,): bold_to_small(P2), (1, 1): bold_to_small(mir(P2))}
    hats = hat_f_exact(red, 2, log=lambda s: None)
    return {rho: sum(c * a ** i * q ** j * t ** k for (i, j, k), c in to_bold(h, A_SHIFT["nontorus"][2]).items())
            for rho, h in hats.items()}


def main(sign=1):
    ps = [sign * k for k in range(1, 13)]
    import pickle
    cache = Path(f"data/cache/v4_hat_sign{sign}.pkl")
    F = pickle.load(open(cache, "rb")) if cache.exists() else {}
    for p in ps:
        if p not in F:
            F[p] = hat_bold(p)
            pickle.dump(F, open(cache, "wb"))          # p ごとに保存（再開用）
    # 候補の比 m（太字の単項式）：連続する p の \\hat f の「最高次」の比から推定する代わりに、小さい単項式の集合を総当たり
    cands = [a ** i * q ** j * t ** k for i in range(0, 9, 2) for j in range(-6, 7, 2) for k in range(-4, 13)]
    for rho in [(2,), (1, 1)]:
        seq = [sp.expand(F[p][rho]) for p in ps]
        # 線形漸化式（定数係数＝ p によらない有理関数）の次数を、ハンケル行列式で調べる（評価点で数値的に）
        import random
        rnd = random.Random(0)
        vals = {a: sp.Rational(rnd.randint(2, 9), rnd.randint(2, 9)), q: sp.Rational(rnd.randint(2, 9), rnd.randint(2, 9)),
                t: sp.Rational(rnd.randint(2, 9), rnd.randint(2, 9))}
        s = [x.subs(vals) for x in seq]
        for order in range(1, 7):
            if 2 * order >= len(s):
                break
            Hm = sp.Matrix([[s[i + j] for j in range(order + 1)] for i in range(order + 1)])
            print(f"ρ={rho} p の符号 {sign}: 漸化式の次数 {order} のハンケル行列式 = {'0' if Hm.det() == 0 else '≠0'}", flush=True)
    return F


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
