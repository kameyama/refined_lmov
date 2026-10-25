r"""C_N の 2 箱（daha_cn.daha_jones_two_box）の検証（progress.md 続き 52）。
  (1) P_{2ε_1}, P_{ε_1+ε_2} が Macdonald 作用素の固有関数（固有値は 2ε_1 > ε_1+ε_2 > 0 の三角性で決まる）。
  (2) u = 1 で D_N の DAHA-Jones の (2), (1,1) に一致（N = 2, 3；D_N の K_λ は Koornwinder の保存データ）。
  (3) Cherednik (ℋ^C_{3,2}(ω_2)) に a = −t^{m−1} で単項式を除いて一致（C_3, C_4）。"""
import sys, time
sys.path.insert(0, "src")
import sympy as sp
from rcsknots.refined_kauffman.daha_cn import CN, two_box_P, P_poly, daha_jones_two_box, _add, ZERO
from rcsknots.refined_kauffman.flint_rf3 import R, S, U, conv

q, t, u, a = sp.symbols("q t u a")
r_, s_, u_ = sp.symbols("r s u", positive=True)
to_qtu = lambda e: sp.cancel(sp.sympify(str(e.to_sympy()), locals={"r": r_, "s": s_, "u": u_}).subs({r_: sp.sqrt(q), s_: sp.sqrt(t), u_: sp.sqrt(u)}))


def eigen_check(N):
    D = CN(N); coef = two_box_P(D)
    ok = {}
    for lam, key in (((2,), "2"), ((1, 1), "11")):
        P = P_poly(D, lam, coef)
        MP = D.macdonald_op(P)
        diff = _add(dict(MP), P, -coef["eig"][key])
        ok[lam] = not diff
    print(f"C_{N}：P_(2), P_(1,1) が固有関数 {ok}", flush=True)
    return all(ok.values())


def vs_DN(N, m):
    from rcsknots.refined_kauffman.daha_dn import daha_jones
    rd, sd = sp.symbols("r s", positive=True)
    out, _ = daha_jones_two_box(N, m)
    res = {}
    for lam in ((2,), (1, 1)):
        Zd, Kd = daha_jones(N, lam, m)
        jd_D = sp.cancel((Zd / Kd).to_sympy().subs({rd: sp.sqrt(q), sd: sp.sqrt(t)}))
        Zc, Pc = out[lam]
        jd_C = sp.cancel(to_qtu(Zc / Pc).subs(u, 1))
        res[lam] = sp.factor(sp.cancel(jd_C / jd_D))
    print(f"N = {N}, T({2*m+1},2)：JD^C(u=1)/JD^D = {res}", flush=True)
    return res


def vs_cherednik_omega2(N):
    H = (1 + q*t + q*t**2 + q**2*t**4 + a*(q + q*t + q**2*t**2 + q**2*t**3 - q*u - q*u/t - q**2*t*u - q**2*t**2*u)
         + a**2*(q**2*t - 3*q**2*u - q*u/t**2 + q**2*u/t**2 - q*u/t + q**3*u/t - 2*q**2*t*u + q**3*t*u - q**3*t**2*u - q**3*t**3*u + q**2*u**2/t)
         + a**3*(-q**2*u - q**2*u/t**2 - 2*q**2*u/t + q**3*u/t - 2*q**3*t*u - q**3*t**2*u + 2*q**3*u**2 + q**2*u**2/t**3 + 2*q**2*u**2/t**2 - q**3*u**2/t**2 + q**2*u**2/t + q**3*t*u**2)
         + a**4*(-q**3*u - q**3*u/t + q**3*u**2 - q**4*u**2 + q**2*u**2/t**3 - q**3*u**2/t**3 + q**3*u**2/t**2 + q**4*u**2/t**2 + 3*q**3*u**2/t - q**4*u**2/t + q**4*t*u**2 - q**3*u**3/t**3 - q**3*u**3/t**2)
         + a**5*(q**4*u**2 + q**3*u**2/t**3 + q**3*u**2/t**2 - q**4*u**2/t**2 - q**3*u**3/t**4 - q**3*u**3/t**3 + q**4*u**3/t**3 - q**4*u**3/t)
         + a**6*(q**4*u**2/t**2 - q**4*u**3/t**3 - q**5*u**3/t**3 - q**4*u**3/t**2 + q**5*u**3/t**2 + q**4*u**4/t**4))
    out, _ = daha_jones_two_box(N, 1)
    Z, P = out[(1, 1)]
    jd = to_qtu(Z / P)
    Hm = H.subs(a, -t**(N - 1))
    ratio = sp.cancel(jd / Hm)
    n_, d_ = sp.fraction(ratio)
    ok = sp.Poly(n_, q, t, u).is_monomial and sp.Poly(d_, q, t, u).is_monomial
    print(f"C_{N}, T(3,2), ω₂：Cherednik の ℋ^C(ω₂; a = −t^{N-1}) との比は単項式 {ok}：{sp.factor(ratio) if ok else '—'}", flush=True)
    return ok


if __name__ == "__main__":
    for N in (2, 3):
        eigen_check(N)
    for N in (2, 3):
        for m in (1, 2):
            vs_DN(N, m)
    for N in (3, 4):
        t0 = time.time(); vs_cherednik_omega2(N); print(f"   {time.time()-t0:.1f}s", flush=True)
