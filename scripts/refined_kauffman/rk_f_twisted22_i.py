r"""brane n・image n の固定セクターの捻ったトレース E を、条件 (i)（t = −1）と (ii)（t = q）の整数計画で決める（theory §12.7、progress.md 続き 72）。
  設定は rk_f_twisted22 と同じ（P = s B_unt Pk、E = s B_tw Σ_m e_m m、|e_m| ≤ |c_m|、e_m ≡ c_m mod 2、ν のブロックごと）。
  条件 (i)：E|_{t=−1} = [2^n (1 + t) f^{SO}_{p_2^n}]_{t=−1}。
    不変な状態（brane の表現 λ、電荷 m、符号 ε）は f^{cov} に m s_λ(Y) s_λ(Y′)、g に ½ε m ψ²[s_λ](Y)/(1 + t) を与える。
    n = 1：E = Σεm = [2(1 + t) g]_{p_2}。n = 2：p_1(Y)²p_1(Y′)² の係数が ¼m（4P の 4）、ψ²[s_λ] の p_2² の係数が ½ なので E = 4(1 + t) g_{p_2²}。
    f^{cov} の分母は (1 − t) だけ（2 箱・3 箱で確認）なので、t = −1 では g の代わりに SO の f を使える。
    n = 1 は 2(1 + t) g_(2)（Schur、rk_b_twisted_multiplet と同じ）、n = 2 は 4(1 + t) f_{(2,2)}（rk_c_4box_so_f）。
  条件 (ii)：E|_{t=q} = −(1 + q) ψ²[F]|_{t=q}（F は rk_f_twisted22.sectors の電荷が半分の凍結した状態）。
  検算：(i) と (ii) の目標は q = t = −1（r = s = i）で一致しなければならない。"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_b_twisted_multiplet as TM
import rk_f_twisted22 as T
import rk_f_frozen_basis as FB

r, s, nu = T.r, T.s, T.nu
t, I = s**2, sp.I
OUT = "data/refined_kauffman/daha"


def target_i(n):
    if n == 1:
        import rk_b_knot_fixed_sector as K
        return sp.cancel(2 * (1 + t) * K.load_g(1)[(2,)])
    u = sp.Symbol("u", positive=True)
    f22 = sp.sympify(pickle.load(open(f"{OUT}/largeN_4box_so_m1.pkl", "rb"))["f"]["(2, 2)"], locals={"r": r, "s": s, "u": u})
    return sp.cancel(4 * (1 + t) * f22.subs(u, nu * s))


def run(n):
    P, F = T.sectors()[n]
    Bu, Bt = T.zq**(2 * n - 1), T.zq**(n - 1) * T.zqp**n
    Pk = sp.cancel(P / (s * Bu))
    Ri = target_i(n)
    Eq = sp.cancel(-(1 + r**2) * T.psi2(F))
    blocks = sorted(FB.nu_blocks(Pk))
    extra = sorted(set(FB.nu_blocks(Ri)) - set(blocks))
    vanish = all(sp.expand(sp.cancel(T.block(Ri, k).subs(s, I))) == 0 for k in extra)
    ok = sp.simplify(sp.cancel(Ri.subs(s, I) - Eq.subs(s, r)).subs(r, I)) == 0
    print(f"== n = {n}：Pk の ν のブロック {blocks}、(i) の目標だけにあるブロック {extra}（t = −1 で 0：{vanish}）、"
          f"(i) と (ii) の目標が r = s = i で一致 {ok}", flush=True)
    E, out = sp.Integer(0), {}
    for k in blocks:
        mult = TM.monomials(T.block(Pk, k))
        basis = lambda key: s * Bt * TM.twisted(key, "mono")
        conds = [(T.block(Ri, k), {s: I}), (T.block(Eq, k), {s: r})]
        res = TM.solve(mult, basis, conds, uniq=True)
        if res["feasible"]:
            Ek = sp.expand(s * Bt * sum(v * TM.twisted(key, "mono") for key, v in res["e"].items()))
            E += Ek
            same = sum(1 for key in mult if res["e"][key] == mult[key])
            msg = f"解あり、一意でない変数 {len(res['free'])}/{res['n']}（e = c が {same}）"
        else:
            msg = f"解なし（変数 {res['n']}）"
        print(f"  ν^{k}：{msg}", flush=True)
        out[k] = res
    return E, out


if __name__ == "__main__":
    ns = [int(a) for a in sys.argv[1:]] or [1, 2]
    for n in ns:
        E, out = run(n)
        if all(v["feasible"] and not v["free"] for v in out.values()):
            pickle.dump({"E": str(E), "e": {k: {str(kk): vv for kk, vv in v["e"].items()} for k, v in out.items()}},
                        open(f"{OUT}/twisted_E_n{n}_m1.pkl", "wb"))
            print(f"   E は一意。保存：{OUT}/twisted_E_n{n}_m1.pkl", flush=True)
