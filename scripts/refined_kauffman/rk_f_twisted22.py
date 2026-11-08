r"""凍結した状態の t の塔を、被覆の固定セクターの状態への Ĩ の作用の分け方として調べる（theory §12.7、progress.md 続き 69）。
  Hopf (4,4)・結び目の 2 箱（theory §9、rk_b_twisted_multiplet.py）と同じ整数計画：
    固定セクターの指数 P = s · B_unt · Pk（B_unt は捻らない境界の因子）、Pk = Σ_m c_m m（単項式 m = r^a s^b ν^c）。
    Ĩ の捻ったトレース E = s · B_tw · Σ_m e_m m、|e_m| ≤ |c_m|、e_m ≡ c_m (mod 2)（Ĩ で入れ替わる対は 2 つずつ）。
    境界の因子（§11.3）：brane n・image n の Ĩ 不変な曲線で B_unt = z_q^{2n−1}、B_tw = z_q^{n−1}(z_q^+)^n。
  条件 (ii)（t = q）：E|_{t=q} = −(1 + q) ψ²[F]|_{t=q}（ψ²：r, s, ν → 2 乗）。F は電荷が半分の凍結した状態：
    n = 1（brane–image (1),(1)）：F = g_(1) + 1（theory §9 と同じ。+1 は unknot の定数、§11.15）、
    n = 2（brane 2・image 2）：F = g_(2) + g_(1,1)（2 箱の p₁² 方向。unknot の p₁² は 0 なので定数は無い）。
  条件 (i)（t = −1、同じ電荷の SO 側の g）は、n = 1 では 2(1 + t) g_(2)、n = 2 では SO の 4 箱が要り、データが無いので課さない。
  ν のブロックごとに解く（どちらの条件も ν を混ぜない）。両端のブロックを主に見る。
  さらに refined の ψ² 型の候補 E = −(1 + q) ψ²[F]（refined のまま、t² の刻みの塔をもつ）が |e| ≤ |c|・偶奇を満たすかを単項式ごとに調べる。"""
import sys, io, contextlib, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_b_twisted_multiplet as TM
import rk_f_frozen_basis as FB

r, s, nu = TM.r, TM.s, TM.nu
t = s**2
zq, zqp = r - 1 / r, r + 1 / r
psi2 = lambda e: e.subs({r: r**2, s: s**2, nu: nu**2}, simultaneous=True)


def sectors():
    """n = 1, 2 の (P, F)。P はどちらも冪和の規格化込み（n = 2 は 4P）。"""
    import rk_b_knot_fixed_sector as K
    g = K.load_g(1)
    d = pickle.load(open("data/refined_kauffman/daha/largeN_cover_fixed22_m1.pkl", "rb"))
    P22 = sp.sympify(d["P"], locals={"r": r, "s": s, "nu": nu})
    x11 = sp.sympify(d["x11"], locals={"r": r, "s": s, "nu": nu})
    return {1: (sp.cancel((1 - t) * x11), sp.cancel(g[(1,)] + 1)),
            2: (sp.cancel(4 * P22), sp.cancel(g[(2,)] + g[(1, 1)]))}


def block(e, k):
    return FB.nu_blocks(e).get(k, sp.Integer(0)) * nu**k


def run(n, ks):
    P, F = sectors()[n]
    Bu, Bt = zq**(2 * n - 1), zq**(n - 1) * zqp**n
    Pk = sp.cancel(P / (s * Bu))
    assert sp.Poly(sp.fraction(Pk)[1], r, s, nu).is_monomial, "P / (s B_unt) が Laurent でない"
    Etgt = sp.cancel(-(1 + r**2) * psi2(F))
    print(f"== n = {n}（brane {n}・image {n}）：Pk の ν のブロック {sorted(FB.nu_blocks(Pk))}", flush=True)
    out = {}
    for k in ks:
        Pkb = block(Pk, k)
        mult = TM.monomials(Pkb)
        Eq = sp.cancel(block(Etgt, k).subs(s, r))
        res = TM.solve(mult, lambda key: s * Bt * TM.twisted(key, "mono"), [(Eq, {s: r})], uniq=True)
        msg = f"解あり、一意でない変数 {len(res['free'])}/{res['n']}" if res["feasible"] else f"解なし（変数 {res['n']}）"
        # refined の ψ² 型の候補
        cand = sp.cancel(block(Etgt, k) / (s * Bt))
        lau = sp.Poly(sp.fraction(cand)[1], r, s, nu).is_monomial
        cmsg = "s B_tw で割り切れない"
        if lau:
            ce = TM.monomials(cand)
            bad_sup = [m for m in ce if m not in mult]
            bad = [m for m in ce if m in mult and (abs(ce[m]) > abs(mult[m]) or (ce[m] - mult[m]) % 2)]
            miss = [m for m in mult if m not in ce and mult[m] % 2]
            cmsg = (f"候補の単項式 {len(ce)}、Pk の台の外 {len(bad_sup)}、|e| > |c| か偶奇の違反 {len(bad)}、"
                    f"候補に無い単項式で c が奇数 {len(miss)}")
        print(f"  ν^{k}：Pk の単項式 {len(mult)}；条件 (ii) だけの整数計画：{msg}；refined の ψ² 型の候補：{cmsg}", flush=True)
        out[k] = res
    return out


if __name__ == "__main__":
    run(1, [-12, -10, -2, 0])
    run(2, [-24, -22, -4, -2])
