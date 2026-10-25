r"""案 B：Γ = ⟨n, k⟩ の平面曲線特異点（x^n = y^k、リンクは T(n,k)）の Hilbert scheme のセル分解（ORS 1201.2115 定理 unnesteddims）。
  半群イデアル 𝔧 ⊂ Γ（余長 l = |Γ∖𝔧| < ∞）ごとにセル A^{N(𝔧)}、N(𝔧) = Σ_i |Γ_{>γ_i}∖𝔧| − Σ_i |Γ_{>σ_i}∖𝔧|（γ_i：生成元、σ_i：隣り合う生成元の間の凸角）。
  曲線 x^n = y^k では、高さが 0 と n の間の列が k 個以下の階段だけが現れる（stable 極限 T(n,∞) では制限なし）。
  P_alg = (a/q)^{μ−1} Σ_𝔧 q^{2l} t^{2N(𝔧)} Π_{j=1}^{r(𝔧)} (1 + t^{2j−1} a²)（r = 生成元の数）を ORS の系（T(2,2k+1)）と照合する。"""
import sys, itertools
import sympy as sp

a, q, t = sp.symbols("a q t")


def semigroup(n, k, bound):
    return sorted({i * n + j * k for i in range(bound // n + 2) for j in range(bound // k + 2) if i * n + j * k <= bound})


def ideals(n, k, lmax):
    """余長 ≤ lmax の半群イデアル（Γ の有限部分集合の補集合として列挙）。"""
    c = n * k  # 導手より大きい範囲を見る
    bound = c + lmax * max(n, k) + n * k
    G = semigroup(n, k, bound)
    Gset = set(G)
    # 𝔧 = Γ ∖ D、D は有限で「Γ-イデアルの補集合」：g ∈ D, h ∈ Γ, g − h ∈ Γ ⇒ g − h ∈ D（下に閉じている）
    small = [g for g in G if g <= c + lmax * max(n, k)]
    out = []
    def rec(D, cand):
        out.append(frozenset(D))
        if len(D) >= lmax:
            return
        for idx, g in enumerate(cand):
            # g を加えられるのは、g より小さい Γ の元で g − h ∈ Γ となるものがすべて D にあるとき
            if all((g - h) in D for h in G if 0 < h <= g and (g - h) in Gset):
                rec(D | {g}, cand[idx + 1:])
    rec(frozenset(), small)
    return list({D for D in out})


def generators(J, G, n, k):
    """イデアル 𝔧 = Γ∖D の極小生成元。"""
    Jl = [g for g in G if g not in J]  # placeholder
    return None


def cell_data(D, n, k, G):
    """D = Γ∖𝔧。生成元 γ（𝔧 の元で、γ − h ∈ 𝔧 となる h ∈ Γ_{>0} が無いもの）、σ、N(𝔧)、r。"""
    Gset = set(G)
    J = [g for g in G if g not in D]
    Jset = set(J)
    gens = [g for g in J if not any((g - h) in Jset for h in G if 0 < h <= g)]
    # 列の高さの制限：曲線 x^n = y^k では、Γ の元 m = a k + b n（0 ≤ a < n）；γ を a（行）で並べる
    r = len(gens)
    # 生成元を階段を上る順（a の昇順 = 行番号）に並べる：γ = a k + b n、0 ≤ a < n
    def ab(m):
        for aa in range(n):
            if (m - aa * k) % n == 0 and m - aa * k >= 0:
                return aa, (m - aa * k) // n
        raise ValueError
    gens_sorted = sorted(gens, key=lambda g: ab(g)[0])            # 階段を上る順 = 行 a の昇順
    m = len(gens_sorted)
    N = 0
    for i, g in enumerate(gens_sorted):
        a1, b1 = ab(g)
        if i + 1 < m:
            a2, _ = ab(gens_sorted[i + 1])
        else:
            a2 = ab(gens_sorted[0])[0] + n                             # 巡回：a_{m+1} = a_1 + n（ORS の例 σ_3 = 28）
        sigma = a2 * k + b1 * n
        N += sum(1 for h in D if h > g) - sum(1 for h in D if h > sigma)
    return gens_sorted, N, r


def P_alg(n, k, lmax, stable=False):
    G = semigroup(n, k, 4 * n * k + lmax * max(n, k))
    mu = (n - 1) * (k - 1)
    tot = 0
    for D in ideals(n, k, lmax):
        gens, N, r = cell_data(D, n, k, G)
        # 階段の列の高さの制限（x^n = y^k）：高さが 0 < h < n の列の数 ≤ k
        if not stable:
            # 列 b（b = 0,1,…）の高さ = #{a : a k + b n ∈ D}
            heights = {}
            for h in D:
                for aa in range(n):
                    if (h - aa * k) % n == 0 and h - aa * k >= 0:
                        bb = (h - aa * k) // n
                        heights[bb] = heights.get(bb, 0) + 1
                        break
            if sum(1 for v in heights.values() if 0 < v < n) > k:
                continue
        l = len(D)
        tot += q**(2 * l) * t**(2 * N) * sp.prod([1 + t**(2 * j - 1) * a**2 for j in range(1, r + 1)])
    return sp.expand((a / q)**(mu - 1) * tot)


if __name__ == "__main__":
    n, k, lmax = 2, 3, 8
    Pa = P_alg(n, k, lmax)                                                 # (a/q)^{μ−1} を含む
    # 目標：DGR の unreduced superpolynomial（ORS の規約：(at + a^{-1})/(q^{-1} − q) × 𝒫、𝒫(T(2,3)) = a²q^{-2} + a⁴t³ + a²q²t²）
    target = sp.expand(sp.series((a * t + 1 / a) * q / (1 - q**2) * (a**2 * q**-2 + a**4 * t**3 + a**2 * q**2 * t**2), q, 0, 2 * lmax).removeO())
    diff = sp.expand((Pa - target) * q)
    low = sum(term for term in sp.Add.make_args(diff) if sp.Poly(term, q).degree() <= 2 * lmax - 2)
    print("T(2,3)：セルの和（l ≤ 8）− DGR の unreduced（q^{2l−1} ≤ q^{13}）=", sp.factor(low), flush=True)
    print("セルの和の低次：", sp.expand(sum(term for term in sp.Add.make_args(sp.expand(Pa * q)) if sp.Poly(term, q).degree() <= 4) / q), flush=True)
    print("目標の低次：", sp.expand(sum(term for term in sp.Add.make_args(sp.expand(target * q)) if sp.Poly(term, q).degree() <= 4) / q), flush=True)
    G = semigroup(n, k, 60)
    print("セル（l, 生成元, N, r）：", sorted((len(D), cell_data(D, n, k, G)[0], cell_data(D, n, k, G)[1], cell_data(D, n, k, G)[2]) for D in ideals(n, k, 5)), flush=True)
