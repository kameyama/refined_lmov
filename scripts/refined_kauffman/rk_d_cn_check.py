r"""C_N の DAHA（src/rcsknots/daha_cn.py）の検証（progress.md 続き 50、案 D の第一歩）。
  (1) 関係式：2 次（t_i）、組紐（短い根どうしは 3 項、α_{N−1}・α_N は 4 項）、交差 T_i X_{ε_i} T_i = X_{s_i ε_i}、Y の可換性。
  (2) u = 1 で D_N の DAHA-Jones（daha_dn）の [1] に一致するか（N = 2, 3、T(3,2), T(5,2)）。
  (3) Cherednik 1111.6195 (super-c-2-3)：ℋ^C_{3,2}(ω_1; q,t,u,a) = 1 + qt + a(qt − qu) + a²(−qu + q²u − q²tu) + a³(−q²tu + q²u²)、
      a = −t^{m−1} で C_m の JD̃ に一致（tilde 正規化：単項式倍を除いて比べる）。m = 2, 3, 4。"""
import sys, random
sys.path.insert(0, "src")
import sympy as sp
from rcsknots.refined_kauffman.daha_cn import CN, daha_jones_omega1, _add, ONE
from rcsknots.refined_kauffman.flint_rf3 import R, S, U, conv

q, t, u, a = sp.symbols("q t u a")
r_, s_, u_ = sp.symbols("r s u", positive=True)


def relations(N):
    random.seed(0)
    D = CN(N)
    f = {}
    for _ in range(4):
        _add(f, {tuple(random.randint(-2, 2) for _ in range(N)): conv(random.randint(1, 5))})
    eq = lambda A, B: not _add(dict(A), B, conv(-1))
    ok = {}
    for i in range(N):
        th = D.tpar[i]
        ok[f"quad{i}"] = not _add(_add(D.T(i, D.T(i, f)), D.T(i, f), -(th - 1 / th)), f, conv(-1))
    if N >= 3:
        ok["braid3"] = eq(D.T(0, D.T(1, D.T(0, f))), D.T(1, D.T(0, D.T(1, f))))
    i, j = N - 2, N - 1
    ok["braid4"] = eq(D.T(i, D.T(j, D.T(i, D.T(j, f)))), D.T(j, D.T(i, D.T(j, D.T(i, f)))))
    for i in range(N):
        ok[f"cross{i}"] = eq(D.T(i, D.Xmul(D.e(i), D.T(i, f))), D.Xmul(D.reflect(i, D.e(i)), f))
    Y2 = lambda g: D.T(0, D.Y1(D.T(0, g, inv=True)), inv=True)
    ok["Ycomm"] = eq(D.Y1(Y2(f)), Y2(D.Y1(f)))
    print(f"C_{N} の関係式：{ok}", flush=True)
    return all(ok.values())


def jd_sympy(N, m):
    Z, X1 = daha_jones_omega1(N, m)
    return sp.cancel((Z / X1).to_sympy().subs({r_: sp.sqrt(q), s_: sp.sqrt(t), u_: sp.sqrt(u)}))


def vs_DN(N, m):
    from rcsknots.refined_kauffman.daha_dn import daha_jones
    Zd, Kd = daha_jones(N, (1,), m)
    rd, sd = sp.symbols("r s", positive=True)
    jd_D = sp.cancel((Zd / Kd).to_sympy().subs({rd: sp.sqrt(q), sd: sp.sqrt(t)}))
    jd_C = sp.cancel(jd_sympy(N, m).subs(u, 1))
    ratio = sp.cancel(jd_C / jd_D)
    print(f"N = {N}, T({2*m+1},2)：JD^C(u=1)/JD^D = {sp.factor(ratio)}", flush=True)
    return ratio


def vs_cherednik(m_rank):
    H = 1 + q * t + a * (q * t - q * u) + a**2 * (-q * u + q**2 * u - q**2 * t * u) + a**3 * (-q**2 * t * u + q**2 * u**2)
    Hm = sp.expand(H.subs(a, -t**(m_rank - 1)))
    jd = jd_sympy(m_rank, 1)
    res = {}
    for lab, sub in (("そのまま", {}), ("q→1/q,t→1/t,u→1/u", {q: 1 / q, t: 1 / t, u: 1 / u})):
        ratio = sp.cancel(jd.subs(sub, simultaneous=True) / Hm)
        n_, d_ = sp.fraction(ratio)
        res[lab] = (sp.Poly(n_, q, t, u).is_monomial and sp.Poly(d_, q, t, u).is_monomial, sp.factor(ratio) if len(str(ratio)) < 80 else "…")
    print(f"C_{m_rank}, T(3,2)：Cherednik の ℋ^C(a = −t^{m_rank-1}) との比：{res}", flush=True)
    return res


if __name__ == "__main__":
    for N in (2, 3):
        relations(N)
    for N in (2, 3):
        for m in (1, 2):
            vs_DN(N, m)
    for N in (2, 3, 4):
        vs_cherednik(N)
