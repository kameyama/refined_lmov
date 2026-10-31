"""γ = τ_+^m τ_−^k（トーラス結び目 T(mk+1, k)）のデータのタグ・結び目名・r·s・δ（progress.md 続き 102）。
k = 2 が従来の T(2,2m+1)（タグ m{m}、δ = m）。k ≠ 2 はタグ m{m}k{k}（例：T(3,4) は m = 1, k = 3）。"""

# data（KN の規格化）の行列式表現のずれ：U(2) で R_(1,1)(a = t) = (t/q)^δ（conventions.md「一様な規格化」。T(2,2m+1) は δ = m）
DELTA = {(1, 3): 3}


def tag(m, k=2):
    return f"m{m}" if k == 2 else f"m{m}k{k}"


def knot_name(m, k=2):
    return f"T2_{2 * m + 1}" if k == 2 else f"T{k}_{m * k + 1}"


def rs(m, k=2):
    """T(r, s) の r·s（framing の単項式 T_λ^{−rs}。k = 2 では 2(2m + 1)）。"""
    return (m * k + 1) * k


def delta(m, k=2):
    return m if k == 2 else DELTA[(m, k)]
