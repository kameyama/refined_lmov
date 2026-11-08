"""γ = τ_+^m τ_−^k（トーラス結び目 T(mk+1, k)）のデータのタグ・結び目名・r·s・δ（progress.md 続き 102）。
k = 2 が従来の T(2,2m+1)（タグ m{m}、δ = m）。k ≠ 2 はタグ m{m}k{k}（例：T(3,4) は m = 1, k = 3）。
m が文字列（例 "T3_5"）なら daha_dn.GAMMA_WORDS の語（τ_± の一般の語、続き 103）で、タグ・結び目名はその文字列、k は無視する。"""
from rcsknots.refined_kauffman.daha_dn import gamma_factors, gamma_first_column

# data（KN の規格化）の行列式表現のずれ δ：U(2) で R_(1,1)(a = t) = (t/q)^δ（conventions.md「一様な規格化」）。
# トーラス結び目では種数 (r − 1)(s − 1)/2 に等しい：T(2,2m+1) の δ = m と、T(3,4) = 3、T(3,5) = 4、T(3,7) = 6、T(4,5) = 6 を data で確認（test_rk_torus_tag、続き 102）。


def first_column(m, k=2):
    return gamma_first_column(gamma_factors(m, k))


def tag(m, k=2):
    if isinstance(m, str):
        return m
    return f"m{m}" if k == 2 else f"m{m}k{k}"


def knot_name(m, k=2):
    if isinstance(m, str):
        return m
    return f"T2_{2 * m + 1}" if k == 2 else f"T{k}_{m * k + 1}"


def rs(m, k=2):
    """T(r, s) の r·s（framing の単項式 T_λ^{−rs}。k = 2 では 2(2m + 1)）。"""
    r, s = first_column(m, k)
    return r * s


def delta(m, k=2):
    r, s = first_column(m, k)
    return (r - 1) * (s - 1) // 2
