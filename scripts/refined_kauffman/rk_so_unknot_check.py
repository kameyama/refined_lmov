"""refined Kauffman：refined SO(2N) unknot の検算。

D_N 型 Macdonald 多項式の主特殊化 P_λ(t^ρ) の一般 t の評価式
   P_λ(t^ρ) = t^{-<λ,ρ>} Π_{α>0} Π_{i=0}^{<λ,α^∨>-1} (1 - q^i t^{<ρ,α^∨>+1})/(1 - q^i t^{<ρ,α^∨>})
を (i) t = q^β で Aganagic–Schaeffer 1202.4456 (5.5) の積と数値比較し、
(ii) 2018 年のノートブック check kauf の refined SO unknot（skew Macdonald の和、a = t^{2N-1}）と記号比較する。
結果（2026-10-14）：[1] は一致、[2]・[1,1] は (q − t) に比例する差が残る（docs/refined_kauffman/past_attempt.md）。"""
import sympy as sp
q, t, u = sp.symbols('q t u', positive=True)

def eval_generic(lam, N):
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

def dim_AS(lam, N, beta):
    lam = list(lam) + [0]*(N-len(lam))
    rho = [N-i for i in range(1, N+1)]
    num = sp.Integer(1); den = sp.Integer(1)
    for i in range(N):
        for j in range(i+1, N):
            for sgn in (+1, -1):
                la = lam[i] + sgn*lam[j]; ra = rho[i] + sgn*rho[j]
                for m in range(beta):
                    num *= q**sp.Rational(la+m, 2)*t**sp.Rational(ra, 2) - q**(-sp.Rational(la+m, 2))*t**(-sp.Rational(ra, 2))
                    den *= q**sp.Rational(m, 2)*t**sp.Rational(ra, 2) - q**(-sp.Rational(m, 2))*t**(-sp.Rational(ra, 2))
    return num/den

qv = sp.Rational(3, 7)
for N in (3, 4):
    for beta in (2, 3):
        for lam in ((1,), (2,), (1, 1)):
            A = dim_AS(lam, N, beta).subs({q: qv, t: qv**beta})
            G = eval_generic(lam, N).subs({q: qv, t: qv**beta})
            print(N, beta, lam, "AS==generic:", sp.simplify(A - G) == 0)

# 記号比較：u = t^{N-1}、a = u^2 t（a^{1/2} = t^{N-1/2}）
def eval_generic_u(lam, N):
    # N を記号的に扱えないので、N=6 で計算して t^{N-1} を u に置き換える（λ の長さ ≤ 2 なら N ≥ 4 で構造は同じ）
    e = sp.together(eval_generic(lam, N))
    return e
sa = u*sp.sqrt(t); a = u**2*t
st = sp.sqrt(t)
rkauf1 = 1 + (sa - 1/sa)/(st - 1/st)
rkauf2 = ((sa-1/sa)*(1+q)*(-1+t))/((st-1/st)*(-1+q*t)) + ((((sa-1/sa)**2*(1+q)*(-1+t))/((st-1/st)**2)) + (((a-1/a)*(-1+q)*(1+t))/(t-1/t)))/(-2+2*q*t)
rkauf11 = (sa-1/sa)/(st-1/st) + sp.Rational(1,2)*(((sa-1/sa)**2)/((st-1/st)**2) - ((a-1/a)/(t-1/t)))
for N in (5, 6):
    for lam, nb in (((1,), rkauf1), ((2,), rkauf2), ((1,1), rkauf11)):
        G = eval_generic(lam, N)
        # u = t^{N-1} を代入した nb と比べる
        d = sp.simplify(G - nb.subs(u, t**(N-1)))
        print("N=", N, lam, "generic - notebook =", sp.factor(d))
