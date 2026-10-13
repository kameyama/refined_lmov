"""段階 D2（docs/dgr_lmov/research_plan.md）：d_{1|1} での \\hat f_ρ の普遍式（2 箱）と、最大 genus の c_ρ。

hook 性質（docs/dgr_lmov/report.md §5）を仮定すると、2 箱の入力の d_{1|1} での値は A = rCS_{[1]}|_{d11} だけで決まる。
  [1] → A(𝐪, 𝐭)、[2] → c·A(𝐪²𝐭, 𝐭)、[1,1] → c·A(𝐪², 𝐭)、Adams の [1](q², t²) → A(𝐪², −𝐭²)（c = (−𝐭)^{−Σ}）。
\\hat f_ρ は入力の（係数が a, q, t の有理式の）多項式なので、d_{1|1} での値は入力の d_{1|1} での値だけで決まる。
そこで 𝐚 によらない持ち上げ P1 = A、P2 = c·A(𝐪²𝐭, 𝐭) を 2 箱の公式（scripts/framed_knot/v4_closed.py）に入れ、結果を d_{1|1} で特殊化する。
（持ち上げた \\hat f 自体には意味がなく、d_{1|1} での値だけを使う。）

検証：実際の結び目の \\hat f_ρ|_{d11} と全係数で一致すること（ツイスト結び目、T(2,5) など）。
観察：genus g の乱数の A（σ 不変、A = Σ_{j≤g} α_j(𝐭) G^j）で、\\hat f_ρ|_{d11} の最大 genus が何で決まるか。
使い方：PYTHONPATH=src python scripts/dgr_lmov/dgr_d11_universal.py
"""
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "framed_knot"))
from t4_colored import bold_to_small  # noqa: E402
from v4_closed import L, formula, hat_from_inputs, inputs, load_inputs, mir  # noqa: E402
from rcsknots.bps import to_bold  # noqa: E402
from rcsknots.bps import g_expand  # noqa: E402
from rcsknots.dgr_lmov.dgr import specialize  # noqa: E402

G = {(2, 1): 1, (0, 1): -1, (0, -1): -1, (-2, -1): 1}       # (e_𝐪, e_𝐭) の G = 𝐪²𝐭 − 𝐭 − 𝐭⁻¹ + 𝐪⁻²𝐭⁻¹


def mul(u, v):
    out = {}
    for (a1, b1), x in u.items():
        for (a2, b2), y in v.items():
            k = (a1 + a2, b1 + b2)
            out[k] = out.get(k, 0) + x * y
    return {k: w for k, w in out.items() if w}


def from_genus(alphas):
    """A = Σ_j α_j(𝐭) G^j、α_j = {e_𝐭: c}。{(e_𝐪, e_𝐭): c}。"""
    A, Gj = {}, {(0, 0): 1}
    for al in alphas:
        for k, v in mul(Gj, {(0, e): c for e, c in al.items()}).items():
            A[k] = A.get(k, 0) + v
        Gj = mul(Gj, G)
    return {k: v for k, v in A.items() if v}


def lift(A, shift_t=0):
    """{(e_𝐪, e_𝐭)} → 𝐚 によらない太字の多項式。"""
    return {(0, eq, et + shift_t): c for (eq, et), c in A.items()}


def hook(A, n, arm):
    """A(𝐪^n 𝐭^{arm}, 𝐭)。"""
    return {(n * eq, et + arm * eq): c for (eq, et), c in A.items()}


def universal(A, Sigma=0):
    """hook 性質から作る \\hat f_ρ|_{d11}（ρ = [2], [1,1]）。"""
    P1 = lift(A)
    P2 = lift(hook(A, 2, 1), -Sigma)                      # c = (−𝐭)^{−Σ}：符号は (−1)^Σ
    if Sigma % 2:
        P2 = {k: -v for k, v in P2.items()}
    return {rho: specialize(b, 1, 1) for rho, b in hat_from_inputs(P1, P2).items()}


def universal_qdeg(A, Sigma=0):
    """割り算をせず、分子と分母を別々に d_{1|1} で特殊化し、𝐪 の最高次数の差 / 2（= 最大 genus）を返す（ρ = [2], [1,1]）。"""
    P1, P2 = lift(A), lift(hook(A, 2, 1), -Sigma)
    if Sigma % 2:
        P2 = {k: -v for k, v in P2.items()}
    r1, r2, r11 = (L.from_dict(bold_to_small(P)) for P in (P1, P2, mir(P2)))
    r1d = L.from_dict({(2 * i, 2 * j, 2 * k): c for (i, j, k), c in bold_to_small(P1).items()})
    R = [r2, r11, r1 * r1, r1d]
    out = {}
    for rho, (ns, den) in formula().items():
        num = ns[0] * R[0] + ns[1] * R[1] + ns[2] * R[2] + ns[3] * R[3]
        sn, sd = specialize(to_bold(num.to_dict(), 0), 1, 1), specialize(to_bold(den.to_dict(), 0), 1, 1)
        hi = lambda S: max(k[0] for k in S)  # noqa: E731
        lo = lambda S: min(k[0] for k in S)  # noqa: E731
        out[rho] = ((hi(sn) - hi(sd)) / 2, (lo(sn) - lo(sd)) / 2) if sn else None
    return out


def top_genus(S):
    """σ 不変な {(e_𝐪, e_𝐭)} の genus 表示の最大 genus。"""
    if not S:
        return None
    return max(k[0] for k in g_expand({(0, eq, et): c for (eq, et), c in S.items()}))


def main():
    cache = load_inputs()
    # 1. 実データとの照合
    for p in (-3, -2, -1, 1, 2, 3, 5):
        P1, P2 = inputs(p, cache)
        A = specialize(P1, 1, 1)
        ((_, ct), cs), = specialize(P1, 1, 0).items()
        Sigma = -ct
        real = {rho: specialize(b, 1, 1) for rho, b in hat_from_inputs(P1, P2).items()}
        uni = universal(A, Sigma)
        print(f"K_{p}: Σ = {Sigma}  普遍式と一致 {all(real[r] == uni[r] for r in real)}  "
              f"最大 genus（\\hat f） {[top_genus(real[r]) for r in real]}  次数の差 {universal_qdeg(A) if Sigma == 0 else '-'}  A の genus {top_genus(A)}", flush=True)
    # 2. 乱数の A
    random.seed(1)
    for g in range(0, 6):
        for trial in range(3):
            alphas = [{e: random.randint(-3, 3) for e in range(-3, 4)} for _ in range(g + 1)]
            alphas[-1][0] = random.choice([1, 2, -1])
            A = from_genus(alphas)
            if top_genus(A) != g:
                continue
            deg = universal_qdeg(A)
            print(f"乱数 A（genus {g}、試行 {trial}）：\\hat f|_{{d11}} の (𝐪 最高次/2, 最低次/2) [2] {deg[(2,)]}, [1,1] {deg[(1, 1)]}", flush=True)


if __name__ == "__main__":
    main()
