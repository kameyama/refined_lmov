r"""refined Kauffman：refined unknot の a 依存性（極）の確認【計算】。

(1) AS 1202.4456 の SO(2N) refined CS の unknot dim_{q,t}^{D_N}(λ) = P^{D_N}_λ(t^ρ) を、一般 t の評価式
    P_λ(t^ρ) = t^{-<λ,ρ>} Π_{α>0} Π_{i=0}^{<λ,α^∨>-1} (1 - q^i t^{<ρ,α^∨>+1})/(1 - q^i t^{<ρ,α^∨>})
    で N = 4..7 について求め、u = t^{N-1}（a^{1/2} = u t^{1/2}, a = t^{2N-1}）の有理式として分母の因子を見る。
(2) U(N) の composite 表現 ([1],[1]) = 随伴表現（ヤング図 (2,1^{N-2})）の A 型 Macdonald 次元（(q,t)-hook の式）を a = t^N の有理式として見る。
結果（2026-10-14）：どちらも t ≠ q で分母に (1 − q t^{N-1}) 型の因子が残り、a の Laurent 多項式にならない（t = q では約分される）。
"""
import sympy as sp
q, t, u, a = sp.symbols('q t u a')


def eval_DN(lam, N):
    lam = list(lam) + [0]*(N-len(lam))
    rho = [N-i for i in range(1, N+1)]
    res = t**(-sum(l*r for l, r in zip(lam, rho)))
    for i in range(N):
        for j in range(i+1, N):
            for sgn in (+1, -1):
                la = lam[i] + sgn*lam[j]; ra = rho[i] + sgn*rho[j]
                for k in range(la):
                    res *= (1 - q**k * t**(ra+1))/(1 - q**k * t**ra)
    return res


def in_u(expr, N):
    """t^{N-1} の冪を u に置き換える（t の冪を t^{m(N-1)+r} と分解）。"""
    e = sp.factor(expr)
    num, den = sp.fraction(e)
    def conv(poly):
        poly = sp.Poly(sp.expand(poly), t)
        out = 0
        for (k,), c in poly.as_dict().items():
            m, r = divmod(k, N-1)
            out += c * u**m * t**r
        return sp.expand(out)
    # 分母・分子が t の多項式でない場合（負冪）は t^{big} を掛けて処理
    return sp.factor(conv(num)) / sp.factor(conv(den)) if (num.is_polynomial(t) and den.is_polynomial(t)) else e


def macdonald_dim_A(lam, N):
    """A 型：P_λ(1,t,…,t^{N-1}) = t^{n(λ)} Π_s (1 - q^{a'} t^{N-l'})/(1 - q^{a} t^{l+1})。"""
    lamT = [sum(1 for x in lam if x > j) for j in range(lam[0])]
    res = t**sum(i*lam[i] for i in range(len(lam)))
    for i, row in enumerate(lam):
        for j in range(row):
            arm = row - j - 1; leg = lamT[j] - i - 1
            res *= (1 - q**j * t**(N - i))/(1 - q**arm * t**(leg+1))
    return res


if __name__ == "__main__":
    print("== (1) D_N refined unknot（AS）、u = t^{N-1}")
    for lam in ((1,), (2,), (1, 1), (3,), (2, 1), (1, 1, 1)):
        for N in (5, 6):
            e = sp.factor(eval_DN(lam, N))
            den = sp.factor_list(sp.fraction(e)[1])
            print(f"  λ={lam} N={N}: 分母の因子 = {[f for f, m in den[1]]}")
    print("== (2) U(N) composite ([1],[1]) = (2,1^{N-2}) の A 型 Macdonald 次元")
    for N in (4, 5, 6):
        lam = [2] + [1]*(N-2)
        e = sp.factor(macdonald_dim_A(lam, N))
        print(f"  N={N}: {e}")
        print(f"        t=q: {sp.factor(e.subs(q, t))}")
