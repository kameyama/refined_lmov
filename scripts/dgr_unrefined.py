"""段階 D5（docs/dgr_lmov/research_plan.md）：非 refined（t = q、𝐭 = −1）の商と、その genus 表示の符号の規則。

R_ρ（refined の商）を 𝐭 = −1 にすると、d_{1|1} が d_{0|0} に重なるので |ρ| ≥ 2 では (1 − 𝐚²) でもう一度割れる。
Ru_ρ = R_ρ|_{𝐭=−1} / (1 − 𝐚²)（|ρ| = 1 は R_{[1]}|_{𝐭=−1} = HOMFLY-PT）を、非 refined の genus 表示
Ru = Σ n_{g,β} 𝐚^{2β} z^{2g}、z² = (𝐪 − 𝐪⁻¹)²（𝐭 = −1 で G = −z²）に展開し、n_{g,β} の符号を調べる。
また refined の genus 表示の 𝐭 の偶奇が (g, β) ごとにそろうか（KN の 2J_r の偶奇）と、そろうときに符号が (−1)^F で決まるかを見る。
使い方：PYTHONPATH=src python scripts/dgr_unrefined.py 結び目...（twist:p はツイスト結び目 K_p の 2 箱）
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rcsknots.bps import g_expand, to_bold  # noqa: E402
from rcsknots.dgr import quotient  # noqa: E402
from rcsknots.extract import A_SHIFT  # noqa: E402
from rcsknots.lmov import RHOS, load_lmov  # noqa: E402


def unrefined_genus(N):
    """refined の genus 表示 {(g, β, F): c} → 非 refined の Σ n_{g,β} 𝐚^{2β} z^{2g}（G = −z²、𝐭 = −1）。"""
    out = {}
    for (g, b, F), c in N.items():
        k = (int(g), int(b))
        out[k] = out.get(k, 0) + int(c) * (-1) ** ((int(F) + int(g)) % 2)
    return {k: v for k, v in out.items() if v}


def divide_1ma2(n):
    """Σ n_{g,β} 𝐚^{2β} z^{2g} を (1 − 𝐚²) で割る（genus ごとに β の多項式の割り算）。"""
    out = {}
    for g in {k[0] for k in n}:
        row = {b: c for (gg, b), c in n.items() if gg == g}
        lo, hi = min(row), max(row)
        acc = 0
        for b in range(lo, hi):
            acc += row.get(b, 0)           # 商の係数 = 累積和
            if acc:
                out[(g, b)] = acc
        if acc + row.get(hi, 0) != 0:
            raise ValueError("not divisible by 1 - a^2")
    return out


def main():
    for K in sys.argv[1:]:
        if K.startswith("twist:"):
            from framing_twist import hat_framed
            from v4_closed import inputs
            hats = hat_framed(*inputs(int(K[6:])), 0)
        else:
            F = load_lmov(K)
            hats = {rho: to_bold(F[("hat", rho)], A_SHIFT["nontorus"][sum(rho)]) for rho in RHOS}
        res = []
        for rho, b in hats.items():
            R = quotient(b, sum(rho))
            N = g_expand(R)
            par = {}
            for (g, b, Fj) in N:
                par.setdefault((g, b), set()).add(Fj % 2)
            aligned = all(len(s) == 1 for s in par.values())
            n = unrefined_genus(N)
            if sum(rho) > 1:
                n = divide_1ma2(n)
            # 符号の規則：sign(n_{g,β}) = ε (−1)^{jg + kβ} となる (j, k)。割り算の前（n0）と後（n）
            n0 = unrefined_genus(N)
            rule = lambda d, jk: len({(-1) ** ((jk[0] * g + jk[1] * b) % 2) * (1 if c > 0 else -1) for (g, b), c in d.items()}) == 1  # noqa: E731
            JK = [(0, 0), (0, 1), (1, 0), (1, 1)]
            sg = [jk for jk in JK if rule(n0, jk)], [jk for jk in JK if rule(n, jk)]
            cancel = sum(abs(c) for c in N.values()) - sum(abs(c) for c in unrefined_genus(N).values())
            res.append(f"{''.join(map(str, rho))}: 偶奇一致 {aligned} 打ち消し {cancel} 符号 (−1)^(jg+kβ) の (j,k)：前 {sg[0]} 後 {sg[1]}")
        print(K, " | ".join(res), flush=True)


if __name__ == "__main__":
    main()
