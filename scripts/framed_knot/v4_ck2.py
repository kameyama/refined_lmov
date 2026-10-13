"""第 4 期 V4 の残り（続き）：c_k の sl_2 構造と、[1]・[r] への拡張の確認。

v = 𝐪²𝐭, τ = 𝐭 の座標（scripts/framed_knot/v4_ck.py）。KN の G = (𝐪 − 𝐪⁻¹)(𝐪𝐭 − 𝐪⁻¹𝐭⁻¹) = (v + v⁻¹) − (τ + τ⁻¹) なので
v は精密化 GV の左の変数 q_L、τ = 𝐭 は右の変数 q_R（𝐭^{2J_r}）に当たる。σ：v ↦ 1/v は SU(2)_L のワイル群。

1. sl_2 正値性（行ごとに [k]_v の正係数の和）を、\\hat f(K_p)（|p| ≤ 10）と pretzel K_{s,t}（小さい s, t）で調べる。
2. 導出の確認：\\hat f(K_p) = Σ_{k≤2} c_k h_{n−k}(1, X, X²) + Σ_m h_{n−3−m}(1, X, X²) ∂_σ(F_Y μ_Y^m)
   （ここで ∂_σ f = (f − σf)/(μ_Y − μ_Ȳ)。第 2 項は c_3 と同じ理由で sl_2 正）を |p| = 6, 7 で確かめる。
3. c_2 の障害：τ → 0 の展開で c_2 = T1 + T2（T1 = F_{X²}/((X²−Y)(X²−Ȳ)) ≤ 0、T2 = ∂_σ(F_Y/(Y−X²)) ≥ 0）となり打ち消しが要ることを、
   最低次の係数で示す。c_2 の sl_2 係数 f_k を D_2 = N[X², Y] の反対称部分の上からの交代和として書けることを確かめる。
4. [1]：\\hat f_{[1]} ∝ (1+𝐚²𝐭) P_[1] のニュートン形（並び 1, X）で c_1^{[1]} が単項式倍の Z = (1+𝐚²𝐭)x になることを確かめる。
5. [3]：GNSSS の一番上のブロック B^{(r)}_r が積になること（r ≤ 5、ϖ(1, r) が単項式）と、π_{ρ,[1³]} B̄_3 の分母因子が
   既知の根との差の「−」の半分で消えるかを調べる（[2,1] の入力が無いので \\hat f_{[3]} 自体は作れない）。
使い方：python scripts/framed_knot/v4_ck2.py
"""
import sys
import time
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from v4_closed import RHOS, hat_direct, hat_from_inputs, load_inputs  # noqa: E402
from v4_ck import (ORDER, ROOTS, SUBS, V, W, a, amplitudes, avt_terms, ck, input_blocks, parts, q, qint,  # noqa: E402
                   residues, t, tau, to_avt, to_dict, v, varpi1)
from t3_pretzel import tc_graded  # noqa: E402


def from_dict(d):
    return sum(c * a**i * q**j * t**k for (i, j, k), c in d.items())


def sl2_report(expr):
    """(σ 不変か, sl_2 係数の正の数, 負の数)。行は (𝐚 の指数, τ の指数)。expr は sympy 式または太字の dict。"""
    d = expr if isinstance(expr, dict) else to_dict(expr)
    rows = {}
    for (al, be, ga), x in d.items():
        assert be % 2 == 0
        rows.setdefault((al, ga - be // 2), {})[be // 2] = x
    sym = all(r.get(-n, 0) == x for r in rows.values() for n, x in r.items())
    fs = [r.get(k - 1, 0) - r.get(k + 1, 0) for r in rows.values() for k in {abs(n) + 1 for n in r}]
    return sym, sum(1 for f in fs if f > 0), sum(1 for f in fs if f < 0)


def hcomp(n, xs):
    if n < 0:
        return sp.Integer(0)
    return sp.expand(sum(sp.prod([xs[i]**e for i, e in enumerate(es)]) for es in compositions(n, len(xs))))


def compositions(n, k):
    if k == 1:
        yield (n,)
        return
    for i in range(n + 1):
        for rest in compositions(n - i, k - 1):
            yield (i,) + rest


def part1():
    print("== 1. sl_2 正値性：ツイスト結び目 \\hat f(K_p) ==", flush=True)
    cache = load_inputs()
    for p in [s * n for n in range(1, 11) for s in (1, -1)]:
        H = hat_direct(p, cache)
        print(f"p={p}: " + "；".join(f"{rho} {sl2_report(H[rho])}" for rho in RHOS), flush=True)
    print("== 1'. sl_2 正値性：pretzel K_{s,t}（入力は GNSSS の予想式） ==", flush=True)
    for s_, t_ in [(2, 2), (2, 3), (3, 3), (-2, 2), (-2, 3), (-3, 2), (-2, -2), (-1, 3), (-3, 4)]:
        H = hat_from_inputs(tc_graded(s_, t_, 1), tc_graded(s_, t_, 2))
        print(f"K_{{{s_},{t_}}}: " + "；".join(f"{rho} {sl2_report(H[rho])}" for rho in RHOS), flush=True)


def part2_3(AJ):
    print("== 2. 分解 \\hat f = (一重項の部分) + Σ_m h ∂_σ(F_Y μ_Y^m) の確認、3. c_2 の障害 ==", flush=True)
    cache = load_inputs()
    for sign in (1, -1):
        for rho in RHOS:
            F = residues(AJ[rho], sign)
            cs = ck(AJ[rho], sign)
            mu = {r: ROOTS[r]**sign for r in ORDER}
            sing = [mu["1"], mu["X"], mu["X2"]]
            ok = True
            for n in (6, 7):
                tot = sum(cs[k] * hcomp(n - 1 - k, sing[:k + 1]) for k in range(3))
                for m in range(0, n - 3):
                    g = F["Y"] * mu["Y"]**m
                    tot += hcomp(n - 4 - m, sing) * (g - g.subs({q: 1 / (q * t)}, simultaneous=True)) / (mu["Y"] - mu["Yb"])
                ok &= to_dict(sp.cancel(tot)) == hat_direct(sign * n, cache)[rho]
            # 3. τ → 0 の最低次での T1, T2
            FY, FX2 = to_avt(F["Y"]), to_avt(F["X2"])
            m_ = {r: to_avt(mu[r]) for r in ORDER}
            T1 = sp.cancel(FX2 / ((m_["X2"] - m_["Y"]) * (m_["X2"] - m_["Yb"])))
            EY = sp.cancel(FY / (m_["Y"] - m_["X2"]))
            T2 = sp.cancel((EY - EY.subs(v, 1 / v)) / (m_["Y"] - m_["Yb"]))
            low = []
            for T in (T1, T2):
                nu, de = sp.fraction(sp.cancel(T))
                lo = lambda e: sp.Poly(sp.expand(e), tau).terms()[-1]   # noqa: E731  τ の最低次
                (dn,), cn = lo(nu)
                (dd,), cd = lo(de)
                low.append(sp.factor(cn / cd * tau**(dn - dd)))
            # c_2 の sl_2 係数 = D_2 の反対称部分の上からの交代和
            D2 = sp.cancel((FY - FX2) / (m_["Y"] - m_["X2"]))
            lam = sp.cancel((m_["Y"] - m_["Yb"]) / (v**2 - v**-2))
            dt = {}
            for (al, n, ga), x in avt_terms(D2).items():
                dt[(al, n, ga)] = dt.get((al, n, ga), 0) + x
                dt[(al, -n, ga)] = dt.get((al, -n, ga), 0) - x
            # c_2 λ (v + v⁻¹)(v − v⁻¹) = Σ_{n>0} d̃_n (vⁿ − v⁻ⁿ)、c_2 = Σ f_k [k]_v ⇒ d̃_n = λ(f_{n−1} + f_{n+1})
            rec = 0
            rows = {}
            for (al, n, ga), x in dt.items():
                if n > 0 and x != 0:
                    rows.setdefault((al, ga), {})[n] = x
            for (al, ga), r in rows.items():
                f = {}
                for n in sorted(r, reverse=True):
                    f[n - 1] = r[n] - f.get(n + 1, 0)
                for k, x in f.items():
                    if k > 0:
                        rec += x * a**al * tau**ga * qint(k)
            okf = sp.expand(sp.cancel(rec / lam - to_avt(cs[2]))) == 0
            print(f"p の符号 {sign} ρ={rho}: 分解が直接計算と一致 {ok}；τ→0 の最低次 T1 = {low[0]}、T2 = {low[1]}；"
                  f"c_2 の sl_2 係数 = D_2 の交代和の式が一致 {okf}", flush=True)


def part4():
    print("== 4. [1]：ニュートン形（並び 1, X）の c_1^{[1]} ==", flush=True)
    cache = load_inputs()
    Z = a**2 * t * (1 + a**2 * t) * (a**2 * t + q**2) * (1 + a**2 * q**2 * t**3) / q**2
    for sign in (1, -1):
        s = {n: from_dict(cache[sign * n][0]) * (1 + a**2 * t) for n in (1, 2, 3, 4)}
        Xs = ROOTS["X"]**sign
        c0, c1 = s[1], sp.expand(s[2] - s[1])
        ok = all(sp.expand(s[n] - (c0 + c1 * sum(Xs**i for i in range(n - 1)))) == 0 for n in (3, 4))
        print(f"p の符号 {sign}: c_1^{{[1]}} / Z = {sp.factor(sp.cancel(c1 / Z))}（単項式なら予想どおり）、"
              f"ニュートン形が |p| = 3, 4 で成立 {ok}、c_1^{{[1]}} の sl_2 {sl2_report(c1)}", flush=True)


def part5():
    print("== 5. [r]：一番上のブロック B^{(r)}_r の形 ==", flush=True)
    for r in range(1, 6):
        w1 = sp.factor(sp.cancel(varpi1(r).subs(SUBS)))
        top = sp.factor(sp.cancel((varpi1(r) * V(r, r) * W(r, r) * (-sp.Symbol("t_r") * sp.Symbol("t_c"))**(-r)).subs(SUBS)))
        print(f"r={r}: ϖ(1, r) = {w1}；B^{{(r)}}_r = {top}", flush=True)
    print("== 5'. [3]：π_{ρ,[1³]} B̄_3 の分母と、一番上の根 Ȳ_3 と既知の根との差 ==", flush=True)
    from rcsknots.bps import C, M_matrix, Q as QS, T as TS, f_formulas
    from rcsknots.partitions import arm, boxes, leg
    A = sp.Symbol("A")
    lam = (1, 1, 1)
    u = TS**3 * A**-3
    for i, j in boxes(lam):
        u *= (TS**(2 * i) - A**2 * QS**(2 * j)) / (1 - QS**(2 * arm(lam, i, j)) * TS**(2 * (leg(lam, i, j) + 1)))
    F3 = f_formulas(3)
    hat = M_matrix(3).inv() * sp.Matrix([F3[(3,)][0], F3[(2, 1)][0], F3[(1, 1, 1)][0]])
    Ab = sp.Symbol("Ab")
    top = sp.cancel((varpi1(3) * V(3, 3) * W(3, 3) * (-sp.Symbol("t_r") * sp.Symbol("t_c"))**(-3)).subs(SUBS))
    topb = sp.cancel(top.subs({q: 1 / (q * t)}, simultaneous=True))
    M3 = a**6 * q**12 * t**18
    Y3b = sp.cancel(M3.subs({q: 1 / (q * t)}, simultaneous=True))
    X, Y, Yb = ROOTS["X"], ROOTS["Y"], ROOTS["Yb"]
    known = {"1": 1, "X": X, "X²": X**2, "X³": X**3, "Y": Y, "Ȳ": Yb, "XY": X * Y, "XȲ": X * Yb, "Y_3": M3}
    for rho, h in zip([(3,), (2, 1), (1, 1, 1)], hat):
        pi = sp.cancel(sp.diff(sp.together(h), C(lam)) * u * A)          # |ρ| = 3 は a^{+1/2} のずらし
        pib = sp.factor(sp.cancel(sp.expand(pi.subs(A, sp.sqrt(Ab))).subs({Ab: -a**2 * t, QS: -q * t, TS: q})))
        prod = sp.factor(sp.cancel(pib * topb * Y3b))
        den = sp.factor(sp.denom(sp.cancel(prod)))
        print(f"ρ={rho}: π_{{ρ,[1³]}} = {pib}", flush=True)
        print(f"   π B̄_3 Ȳ_3 = {prod}", flush=True)
        dfac = [f for f, _ in sp.factor_list(den)[1]]
        match = {}
        for name, m in known.items():
            for f, _ in sp.factor_list(sp.numer(sp.cancel(Y3b - m)))[1]:
                for g in dfac:
                    if sp.cancel(f / g).is_number:
                        match.setdefault(str(g), []).append(name)
        print(f"   分母の因子と、それを「−」の半分として含む差 (Ȳ_3 − μ)：{match}；対応の無い因子 {[str(g) for g in dfac if str(g) not in match]}", flush=True)

def main():
    t0 = time.time()
    part1()
    PI = parts()
    AJ = amplitudes(PI, input_blocks())
    part2_3(AJ)
    part4()
    part5()
    import resource
    print(f"{time.time() - t0:.0f} 秒、最大メモリ {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20:.0f} MB", flush=True)


if __name__ == "__main__":
    main()
