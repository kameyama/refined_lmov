r"""refined Kauffman（第 3 の研究）：refined unknot の閉じた式。

Aganagic–Schaeffer 1202.4456 の SO(2N) refined Chern–Simons の unknot
  dim_λ = S_{0λ}/S_{00} = M^{D_N}_λ(t^ρ)（D_N 型 Macdonald 多項式の主特殊化）
と計量 g_λ = <M_λ, M_λ>（AS (B.4)）を、U = t^N（u = t^{N−1} = U/t、Q = t^{2N−1} = U²/t）の有理式として与える。
一般 t の評価式（Macdonald/Cherednik）：
  P_λ(t^ρ) = t^{-<λ,ρ>} Π_{α>0} Π_{k=0}^{<λ,α^∨>-1} (1 − q^k t^{<ρ,α^∨>+1})/(1 − q^k t^{<ρ,α^∨>})
を D_N（正ルート e_i ± e_j、ρ_i = N − i）で ℓ(λ) ≤ N−1 の安定域で書き直したもの（scripts/rk_so_unknot_check.py で N = 3..6 の積と一致を確認）。

U(M) の composite 表現 [λ, μ]_M（ヤング図：λ の右に μ_1 列を足し、下に μ の補集合）の A 型 Macdonald 次元は (q,t)-hook の式
  P_ν(1, t, …, t^{M−1}) = t^{n(ν)} Π_{s∈ν} (1 − q^{a'(s)} t^{M − l'(s)})/(1 − q^{a(s)} t^{l(s)+1})
を対称な正規化（x_i = t^{(M+1)/2 − i}）に直し、t^M を A（a_U）に置き換えて与える。
"""
import sympy as sp

q, t, U, A = sp.symbols("q t U A", positive=True)


def _tpow(e_const, e_N):
    """t^{e_const + e_N·N} を U = t^N で表す。"""
    return t**e_const * U**e_N


def dim_DN(lam):
    """SO(2N) refined unknot dim_λ（AS）。λ は分割（ℓ(λ) ≤ N−1 の安定域）。U = t^N の有理式。"""
    lam = tuple(x for x in lam if x > 0)
    l = len(lam)
    res = sp.Integer(1)
    # 前因子 t^{-Σ λ_i (N − i)}
    for i, li in enumerate(lam, start=1):
        res *= _tpow(li*i, -li)
    for i, li in enumerate(lam, start=1):
        for k in range(li):
            # e_i − e_j (j = l+1..N): <ρ,α> = j−i → 伸縮して (1−q^k t^{N−i+1})/(1−q^k t^{l+1−i})
            res *= (1 - q**k*_tpow(-i+1, 1))/(1 - q**k*t**(l+1-i))
            # e_i + e_j (j = l+1..N): <ρ,α> = 2N−i−j ∈ [N−i, 2N−i−l−1] → (1−q^k t^{2N−i−l})/(1−q^k t^{N−i})
            res *= (1 - q**k*_tpow(-i-l, 2))/(1 - q**k*_tpow(-i, 1))
    for i in range(1, l+1):
        for j in range(i+1, l+1):
            for k in range(lam[i-1] - lam[j-1]):
                res *= (1 - q**k*t**(j-i+1))/(1 - q**k*t**(j-i))
            for k in range(lam[i-1] + lam[j-1]):
                res *= (1 - q**k*_tpow(-i-j+1, 2))/(1 - q**k*_tpow(-i-j, 2))
    return sp.factor(res)


def metric_DN(lam):
    """D_N の Macdonald 計量 g_λ（AS (B.4)、λ_i 整数の場合）。U = t^N の有理式。"""
    lam = tuple(x for x in lam if x > 0)
    lamT = [sum(1 for x in lam if x > j) for j in range(lam[0])] if lam else []
    res = sp.Integer(1)
    for i, li in enumerate(lam, start=1):
        for j in range(1, li+1):
            num = ((1 - _tpow(-2*i, 2)*q**(2*j-1)) * (1 - _tpow(-i+1, 1)*q**(j-1)) * (1 + _tpow(-i, 1)*q**(j-1)) * (1 - _tpow(-i-1, 1)*q**j))
            den = ((1 - _tpow(-2*i-1, 2)*q**(2*j)) * (1 - _tpow(-2*i-1, 2)*q**(2*j-1)) * (1 - _tpow(-i, 1)*q**j))
            res *= num/den
            res *= (1 - t**(lamT[j-1]-i)*q**(li-j+1))/(1 - t**(lamT[j-1]-i+1)*q**(li-j))
            res *= (1 - _tpow(-lamT[j-1]-i-1, 2)*q**(li+j))/(1 - _tpow(-lamT[j-1]-i, 2)*q**(li+j-1))
    return sp.factor(res)


def composite_diagram(lam, mu, M):
    """U(M) の composite 表現 [λ, μ]_M のヤング図（行の長さのリスト）。"""
    lam = [x for x in lam if x > 0]; mu = [x for x in mu if x > 0]
    assert len(lam) + len(mu) <= M
    m1 = mu[0] if mu else 0
    rows = [li + m1 for li in lam] + [m1]*(M - len(lam) - len(mu)) + [m1 - mu[-1-k] for k in range(len(mu))]
    return [r for r in rows if r > 0]


def dim_A_sym(nu, M):
    """U(M) の Macdonald 次元（対称な正規化）：t^{-(M-1)|ν|/2} P_ν(1,t,…,t^{M−1})。数値の M。"""
    nuT = [sum(1 for x in nu if x > j) for j in range(nu[0])] if nu else []
    res = t**(-sp.Rational(M-1, 2)*sum(nu) + sum(i*nu[i] for i in range(len(nu))))
    for i, row in enumerate(nu):
        for j in range(row):
            arm = row - j - 1; leg = nuT[j] - i - 1
            res *= (1 - q**j*t**(M-i))/(1 - q**arm*t**(leg+1))
    return res


def dim_composite(lam, mu):
    """composite [λ, μ] の対称正規化の Macdonald 次元を A = t^M（M 十分大）の有理式として与える。
    (q,t)-hook の式で、高さ M − ℓ(λ) − ℓ(μ)、幅 μ_1 の長方形の部分の積は列 j と列 μ_1−1−j を組にすると伸縮し、
      Π_{j<μ_1} Π_{s<μ^T_{j+1}} (1 − q^j A t^{−ℓ(λ)−s}) / (1 − q^j t^{ℓ(μ)−s})
    になる。残りの箱（上の ℓ(λ) 行と下の ℓ(μ) 行）は M に依らないか t^M → A で置き換えられる。"""
    lam = [x for x in lam if x > 0]; mu = [x for x in mu if x > 0]
    ll, lm = len(lam), len(mu)
    m1 = mu[0] if mu else 0
    muT = [sum(1 for x in mu if x > j) for j in range(m1)]
    M = ll + lm + 4  # 前因子の指数の決定に使う（結果は M に依らない）
    res = sp.Integer(1)
    # 長方形の部分
    for j in range(m1):
        for s_ in range(muT[j]):
            res *= (1 - q**j * A * t**(-ll - s_)) / (1 - q**j * t**(lm - s_))
    # 列 j（j < μ_1）の高さは M − c_j、c_j = #{k : μ_1 − μ_k ≤ j} = μ^T_{μ_1 − j}
    # 上の ℓ(λ) 行
    for i in range(ll):
        row = lam[i] + m1
        for j in range(row):
            arm = row - j - 1
            res *= (1 - q**j * A * t**(-i))
            if j < m1:
                res /= (1 - q**arm * A * t**(-muT[m1 - 1 - j] - i))   # leg+1 = M − c_j − i、c_j = μ^T_{μ_1 − j}
            else:
                lamT_j = sum(1 for x in lam if x > j - m1)
                res /= (1 - q**arm * t**(lamT_j - i))              # leg+1 = λ^T_{j−μ_1+1} − i
    # 下の ℓ(μ) 行：行 k（k = 0..ℓ(μ)−1、上から）、長さ μ_1 − μ_{ℓ(μ)−k}、位置 i = M − ℓ(μ) + k
    for k in range(lm):
        row = m1 - mu[lm - 1 - k]
        for j in range(row):
            arm = row - j - 1
            res *= (1 - q**j * t**(lm - k))                      # t^{M−i} = t^{ℓ(μ) − k}
            res /= (1 - q**arm * t**(lm - k - muT[m1 - 1 - j]))  # leg+1 = (M − c_j) − (M − ℓ(μ) + k)
    # 前因子 t^{n(ν) − (M−1)|ν|/2}：M の一次式（二次の項は消える）
    def expo(MM):
        nn = composite_diagram(lam, mu, MM)
        return sp.Rational(sum(i*nn[i] for i in range(len(nn)))) - sp.Rational(MM-1, 2)*sum(nn)
    e0, e1, e2 = expo(M), expo(M+1), expo(M+2)
    c1 = e1 - e0
    assert e2 - e1 == c1, "前因子の指数が M の一次式でない"
    c0 = e0 - c1*M
    return sp.factor(res * t**c0 * A**c1)


# ---- Koornwinder 多項式の主特殊化（van Diejen 1996、Sahi 1999；Rains–Warnaar 1506.02755 (6.?) の引用の形）----
def _C0(lam, z):
    return sp.prod([1 - z*q**(j-1)*t**(1-i) for i, row in enumerate(lam, 1) for j in range(1, row+1)])


def _Cm(lam, z):
    lamT = [sum(1 for x in lam if x > j) for j in range(lam[0])] if lam else []
    return sp.prod([1 - z*q**(row-j)*t**(lamT[j-1]-i) for i, row in enumerate(lam, 1) for j in range(1, row+1)])


def _Cp(lam, z):
    lamT = [sum(1 for x in lam if x > j) for j in range(lam[0])] if lam else []
    return sp.prod([1 - z*q**(row+j-1)*t**(2-lamT[j-1]-i) for i, row in enumerate(lam, 1) for j in range(1, row+1)])


def koornwinder_eval(lam, tbar, U=U):
    """K_λ(t_0(1, t, …, t^{n−1}); q, t; t_0, t_1, t_2, t_3)
       = t^{n(λ)} (t_0 t^{n−1})^{−|λ|} (t^n, t_0t_1t^{n−1}, t_0t_2t^{n−1}, t_0t_3t^{n−1}; q,t)_λ / [C^−_λ(t) C^+_λ(t_0t_1t_2t_3 t^{2n−2}/q)]
    を U = t^n の式として返す（t^{n−1} = U/t、t^{2n−2} = U²/t²）。tbar = (t_0, t_1, t_2, t_3)。"""
    lam = tuple(x for x in lam if x > 0)
    t0, t1, t2, t3 = tbar
    n_lam = sum(i*row for i, row in enumerate(lam))
    tn1 = U/t
    res = t**n_lam * (t0*tn1)**(-sum(lam))
    for a in (U, t0*t1*tn1, t0*t2*tn1, t0*t3*tn1):
        res *= _C0(lam, a)
    res /= _Cm(lam, t) * _Cp(lam, t0*t1*t2*t3*tn1**2/q)
    return sp.factor(res)
