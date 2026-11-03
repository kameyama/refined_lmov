"""随伴の塔の母関数（続き 121、theory §12.51）：N_(2)∘σ − N_(1,1) = −(s³/r)(q−t)[R₁(ν)²κ_R(ν) + G_adj(ν)κ_G(ν) + 背景の差] の κ_R, κ_G を閉じた形で求める。
σ 共変性（§12.42）により、2 箱の行（b, c）と composite の G をすべて 0 にすると N_(2) − N_(1,1) は塔そのものになり、
そのときパイプライン（rk_f_n2_rows_sym.build）は級数を使わず ν の有理式で計算できる。a₂, a₄, gadj_j は記号。
実行：PYTHONPATH=src:scripts/refined_kauffman python scripts/refined_kauffman/rk_f_n2_tower.py"""
import sys, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_f_n2_rows_sym as RS
from rk_f_n2_rows_sym import r, s, q, t, nu, Psym, dimc, decomposition_flint, macdonald_P, P2, LAMS

a2, a4 = sp.symbols("a2 a4")
g6, g8, g10 = sp.symbols("gadj6 gadj8 gadj10")


def tower_rational():
    """塔 T(ν) = N_(2) − N_(1,1)（2 箱の行と G を 0 にしたもの）を ν の有理式で返す（Adams の項は λ によらず消える）。"""
    t0 = time.time()
    A = nu**2
    qt1 = (q - 1) * (t - 1)
    R1 = a2 * nu**-2 + a4 * nu**-4
    Gadj = g6 * nu**-6 + g8 * nu**-8 + g10 * nu**-10
    qn = (nu - 1 / nu) / (s - 1 / s)
    kappa = (1 - A) * (1 - q) / ((1 - t) * (1 - q * A / t))
    Kadj = qn * qn - kappa
    W1 = R1 * Psym((1,))
    Wadj = (R1**2 + qt1 * Gadj) * Kadj
    pair = {((), ()): sp.Integer(1), ((1,), ()): W1, ((), (1,)): W1, ((1,), (1,)): Wadj + kappa}
    for lam in P2:
        pair[(lam, ())] = pair[((), lam)] = sp.Integer(0)
    for alpha in P2:
        dc = dimc(alpha, (1,))
        d = sp.cancel((Psym(alpha) * Psym((1,)) - dc) / Psym((1,)))
        pair[(alpha, (1,))] = pair[((1,), alpha)] = d * W1
    for alpha in P2:
        for beta in P2:
            key = (alpha, beta) if (alpha, beta) != ((1, 1), (2,)) else ((2,), (1, 1))
            c1, c0, Kb = decomposition_flint(*key)
            pair[(alpha, beta)] = c1 * Wadj + c0
    PY = {(): {(): 1}}
    for n in (1, 2):
        for lam, v in macdonald_P(n, t, q, lambda c: sp.Rational(c.numerator, c.denominator)).items():
            PY[lam] = v
    T_ = lambda lam: tuple(sum(1 for x in lam if x > j) for j in range(lam[0])) if lam else ()
    Z = {}
    for (lam, mu), v in pair.items():
        if v == 0:
            continue
        sg = (-1)**(sum(lam) + sum(mu))
        for a_, x in PY[T_(lam)].items():
            for b_, y in PY[T_(mu)].items():
                Z[(a_, b_)] = Z.get((a_, b_), 0) + v * sg * x * y
    def mul2(A_, B_):
        out = {}
        for (a1, a2_), x in A_.items():
            for (b1, b2), y in B_.items():
                if sum(a1) + sum(b1) <= 2 and sum(a2_) + sum(b2) <= 2:
                    k = (tuple(sorted(a1 + b1, reverse=True)), tuple(sorted(a2_ + b2, reverse=True)))
                    out[k] = out.get(k, 0) + x * y
        return out
    X = {k: v for k, v in Z.items() if k != ((), ())}
    Lg, Xn = {}, dict(X)
    for n in range(1, 5):
        for k, v in Xn.items():
            Lg[k] = Lg.get(k, 0) + v * sp.Rational((-1)**(n + 1), n)
        if n < 4:
            Xn = mul2(Xn, X)
    zq = r - 1 / r
    phi = {(1, 1): zq, (2,): r + 1 / r}
    den = s * (t - nu**2) * (nu**2 - 1) * (nu**2 * q - 1)
    cross = Lg.get(((1, 1), (2,)), 0) + Lg.get(((2,), (1, 1)), 0)
    T = sp.cancel(2 * (1 - t) / (zq * phi[(1, 1)] * phi[(2,)]) * cross / den)      # N_(2) − N_(1,1)（CHI の差は交差項の ×2）
    print(f"  塔の有理式：{time.time() - t0:.0f} 秒", flush=True)
    return T


if __name__ == "__main__":
    T = tower_rational()
    pref = -(s**3 / r) * (q - t)
    # a₂² の係数（a₄ = 0, g = 0）→ κ_R(ν) ν^{−4}；gadj6 の係数 → κ_G(ν) ν^{−6}
    T2 = sp.cancel(T.subs({a4: 0, g6: 0, g8: 0, g10: 0}))
    kR = sp.cancel(sp.diff(T2, a2, 2) / 2 / pref * nu**4)
    bg = sp.cancel(T2.subs(a2, 0) / pref)
    kG = sp.cancel(sp.diff(T, g6) / pref * nu**6)
    # 検算：a₂a₄ の係数 = 2κ_R ν^{−6}、a₄² = κ_R ν^{−8}、gadj8 = κ_G ν^{−8}
    chk1 = sp.cancel(sp.diff(sp.diff(T, a2), a4) / pref - 2 * kR * nu**-6)
    chk2 = sp.cancel(sp.diff(T, a4, 2) / 2 / pref - kR * nu**-8)
    chk3 = sp.cancel(sp.diff(T, g8) / pref - kG * nu**-8)
    print(f"検算 a₂a₄：{chk1}、a₄²：{chk2}、gadj8：{chk3}", flush=True)
    x = sp.Symbol("x")
    for name, k in (("κ_R", kR), ("κ_G", kG), ("背景の差", bg)):
        kx = sp.factor(k.subs(nu, 1 / sp.sqrt(x)))
        print(f"{name}(ν) = {sp.factor(k)}", flush=True)
        print(f"   x = ν^-2 で：{kx}", flush=True)
        ser = sp.series(k.subs(nu, 1 / sp.sqrt(x)), x, 0, 5).removeO()
        print(f"   級数：{ {i: sp.factor(ser.coeff(x, i)) for i in range(5)} }", flush=True)
