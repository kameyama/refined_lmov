"""H1（docs/dgr_lmov/research_plan.md の D3、docs/dgr_lmov/open_problems.md1）：(S_{1|0}), (S_{0|1}) の導出に使う補題の厳密な検査。

GN の構成（src/rcsknots/gn.py）：\\bar{rCS}_λ(T_{n,m}) = ε(P_λ[p_k ↦ P̃_{kn,km}]·1)、ε は tilde 側で f̃[p_k ↦ 1 − a^k]。
q, t に有理数を入れ、全て Fraction で厳密に計算する。検査する内容：
  (a) 核：a = u で ev(H̃_μ) = H̃_μ[p_k ↦ 1 − u^k] が 0 になる μ
      - u = t：ℓ(μ) ≥ 2
      - u = q⁻¹：μ₁ ≥ 2
      - u = t/q：μ ⊇ (2,2)
  (b) 両立：P̃_{kn,km} は核を保ち、商（1 次元ずつ）にスカラー s_k(j) で作用する。閉じた形は
      - u = t：s = γ^k q^{(k²−k)nm/2 + kmj}、γ = (−t)^{n−1} q^{(n−1)(m−1)/2}
      - u = q⁻¹：s = γ'^k t^{−(k²−k)nm/2 − kmj}、γ' はここで決める（表示する）
  (c) 最終式（reduced）：rCS_{[r]}|_{a=t} = γ^r q^{nm r(r−1)/2}、rCS_{[1^r]}|_{a=q^{-1}} = γ''^r t^{−nm r(r−1)/2}（framing 補正の前）。
      二次の部分はちょうど (N2') の T_λ^{nm}。ケーブルでは F = n₂m₂ + n₂²n₁m₁（= p₂q₂）の二次の部分を除くと c^r 型。
使い方：PYTHONPATH=src python scripts/dgr_lmov/dgr_h1_lemma.py
"""
from fractions import Fraction as Fr

from rcsknots.gn import GN
from rcsknots.partitions import partitions

q, t = Fr(3, 7), Fr(5, 11)


def ev(G, vec, u):
    s = Fr(0)
    for rho, c in G.to_p(vec).items():
        w = c
        for k in rho:
            w *= 1 - u ** k
        s += w
    return s


def check_kernel(G, u, rule, nmax=5):
    return all((ev(G, {mu: Fr(1)}, u) == 0) == rule(mu) for n in range(1, nmax + 1) for mu in partitions(n))


def scalars(G, u, n, m, k, jmax=3):
    """[s_k(j)]（核を保たなければ None）。"""
    out = []
    for j in range(jmax + 1):
        ratios, ok = set(), True
        for mu in (partitions(j) if j else [()]):
            v = {mu: Fr(1)}
            w0, w1 = ev(G, v, u), ev(G, G.op(k, n, m, v), u)
            if w0 == 0:
                ok &= w1 == 0
            else:
                ratios.add(w1 / w0)
        if not ok or len(ratios) != 1:
            return None
        out.append(ratios.pop())
    return out


def main():
    Gt, Gq = GN(q, t, t, lambda x: x), GN(q, t, 1 / q, lambda x: x)
    print("(a) 核 u = t：ℓ(μ) ≥ 2", check_kernel(Gt, t, lambda mu: len(mu) >= 2))
    print("(a) 核 u = 1/q：μ₁ ≥ 2", check_kernel(Gt, 1 / q, lambda mu: mu[0] >= 2))
    print("(a) 核 u = t/q：μ ⊇ (2,2)", check_kernel(Gt, t / q, lambda mu: len(mu) >= 2 and mu[1] >= 2))
    slopes = [(1, 0), (1, 1), (1, 2), (1, -1), (2, 1), (2, 3), (3, 2), (2, -3), (3, 4)]
    for n, m in slopes:
        gam = (-t) ** (n - 1) * q ** Fr((n - 1) * (m - 1), 2)
        for k in (1, 2):
            s = scalars(Gt, t, n, m, k)
            pred = [gam ** k * q ** Fr((k * k - k) * n * m, 2) * q ** (k * m * j) for j in range(4)]
            s2 = scalars(Gt, 1 / q, n, m, k)
            # u = 1/q：j 依存は t^{−kmj}、k 依存は t^{−(k²−k)nm/2} と予想し、γ' = s2[0] / t^{-(k²-k)nm/2} の k 乗性を見る
            ok2 = s2 is not None and all(s2[j] == s2[0] * t ** (-k * m * j) for j in range(4))
            g2 = s2[0] / t ** Fr(-(k * k - k) * n * m, 2) if s2 else None
            print(f"(b) (n,m)=({n},{m}) k={k}: u=t 両立・閉じた形 {s == pred}   u=1/q 両立・j 依存 {ok2}  γ'_k = {g2}")
    # (c) reduced（自明な結び目の値で割る）で比べる。u = 1/q の商のスカラーは c_k(u)·(γ'' U')^k、c_k(u) = (1 − u^k)/(1 − t^k)
    for n, m in [(2, 1), (2, 3), (3, 2), (2, -3)]:
        gam = (-t) ** (n - 1) * q ** Fr((n - 1) * (m - 1), 2)
        g2 = scalars(Gt, 1 / q, n, m, 1)[0] / ((1 - q ** -1) / (1 - t))
        for r in (1, 2, 3):
            row, col = (r,), tuple([1] * r)
            red = Gt.torus(row, n, m) / Gt.torus(row, 1, 0)
            red2 = Gq.torus(col, n, m) / Gq.torus(col, 1, 0)
            print(f"(c) T({n},{m}) r={r}: a=t 行 reduced = γ^r q^(nm r(r−1)/2) {red == gam ** r * q ** Fr(n * m * r * (r - 1), 2)}  "
                  f"a=1/q 列 reduced = γ''^r t^(−nm r(r−1)/2) {red2 == g2 ** r * t ** Fr(-n * m * r * (r - 1), 2)}")
    # ケーブル：inner = (n2, m2)、outer = (n1, m1)。framing F = n2 m2 + n2² n1 m1 の二次の部分を除いた reduced が c^r 型か
    for inner, outer in [((2, 1), (2, 3)), ((2, 3), (2, 1)), ((3, 1), (2, 3)), ((2, 1), (3, 2))]:
        (n2, m2), (n1, m1) = inner, outer
        F = n2 * m2 + n2 * n2 * n1 * m1
        vals = [Gt.cable((r,), inner, outer) / Gt.torus((r,), 1, 0) / q ** Fr(F * r * (r - 1), 2) for r in (1, 2, 3)]
        vals2 = [Gq.cable(tuple([1] * r), inner, outer) / Gq.torus(tuple([1] * r), 1, 0) / t ** Fr(-F * r * (r - 1), 2) for r in (1, 2, 3)]
        print(f"(c) ケーブル inner={inner} outer={outer}（F = {F}）：a=t で c^r 型 {vals[1] == vals[0] ** 2 and vals[2] == vals[0] ** 3}  "
              f"a=1/q で c^r 型 {vals2[1] == vals2[0] ** 2 and vals2[2] == vals2[0] ** 3}")

if __name__ == "__main__":
    main()
