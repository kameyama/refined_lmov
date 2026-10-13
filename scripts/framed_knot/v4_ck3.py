"""第 4 期 V4：c_2 の障害の確認。τ = 𝐭 → 0 の形式的な展開で c_2 = T1 + T2、
T1 = F_{X²}/((μ_{X²}−μ_Y)(μ_{X²}−μ_Ȳ))、T2 = ∂_σ(F_Y/(μ_Y−μ_{X²}))（∂_σ f = (f − σf)/(μ_Y − μ_Ȳ)）とし、
各級数の τ の最低次から 14 次分の係数の符号を数える（(𝐚, v) の単項式ごと）。
p > 0 では T1 ≤ 0、T2 ≥ 0 で、c_2 ≤ 0 は T1 が T2 を項ごとに上回ることに当たる（打ち消しが必須）。
使い方：python scripts/framed_knot/v4_ck3.py
"""
import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from v4_ck import ORDER, RHOS, ROOTS, a, amplitudes, input_blocks, parts, residues, tau, to_avt, v  # noqa: E402

N = 14


def series_signs(T):
    nu, de = sp.fraction(sp.cancel(T))
    dp, np_ = sp.Poly(sp.expand(de), tau), sp.Poly(sp.expand(nu), tau)
    dl, nl = min(m[0] for m in dp.monoms()), min(m[0] for m in np_.monoms())
    d = [sp.expand(dp.coeff_monomial(tau**(dl + i))) for i in range(N)]
    n = [sp.expand(np_.coeff_monomial(tau**(nl + i))) for i in range(N)]
    s, pos, neg = [], 0, 0
    for i in range(N):
        s.append(sp.cancel((n[i] - sum(s[j] * d[i - j] for j in range(i))) / d[0]))
        nn, dd = sp.fraction(sp.together(s[-1]))
        (_, sd), = sp.Poly(dd, a, v).terms()
        for _, co in sp.Poly(sp.expand(nn), a, v).terms():
            pos, neg = (pos + 1, neg) if co * sd > 0 else (pos, neg + 1)
    return pos, neg


def main():
    AJ = amplitudes(parts(), input_blocks())
    for sign in (1, -1):
        for rho in RHOS:
            F = residues(AJ[rho], sign)
            mu = {r: to_avt(ROOTS[r]**sign) for r in ORDER}
            FY, FX2 = to_avt(F["Y"]), to_avt(F["X2"])
            T1 = sp.cancel(FX2 / ((mu["X2"] - mu["Y"]) * (mu["X2"] - mu["Yb"])))
            EY = sp.cancel(FY / (mu["Y"] - mu["X2"]))
            T2 = sp.cancel((EY - EY.subs(v, 1 / v)) / (mu["Y"] - mu["Yb"]))
            print(f"p の符号 {sign} ρ={rho}: T1 の級数の係数（正, 負）{series_signs(T1)}、T2 {series_signs(T2)}", flush=True)


if __name__ == "__main__":
    main()
