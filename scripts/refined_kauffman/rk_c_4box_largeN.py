r"""SO の 4 箱（trefoil、m = 1）の大きな N の reduced 不変量 R_λ(U) = Σ_j c_j U^j を、rk_c_bounded_interp の一次方程式（仮定 A＋B、mode "alpha"）で
決め、使わない N すべてで全係数一致したものだけ data/refined_kauffman/daha/R_U_m1_{λ}.pkl（{j: c_j の文字列}）に保存する（progress.md 続き 71）。
3 箱では、同じ方法の解が従来の方法（N = 4〜13 の Vandermonde）の R_U_m1_{λ}.pkl と一致するかを確かめる。
(1,1,1,1) は N = 4〜12 で一意に決まる。N = 4〜11 では r^{−8} の部分が一意に決まらず、その部分を t = q（NRZ 1310.2240 の [1^4] 彩色の
Kauffman 多項式、rk_c_check_tq と同じ対応）から作り、N = 4〜11 の方程式と N = 12 を満たすかを確かめた（tq1111。N = 4〜12 の解と一致）。
使い方：python rk_c_4box_largeN.py [3box | tq1111]"""
import sys, os, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_c_bounded_interp as B

OUT = "data/refined_kauffman/daha"
r, s = sp.symbols("r s", positive=True)

# λ: (使う N の最小, 最大, 帯を取る N_s)
CFG4 = {(4,): (2, 8, 2), (2, 2): (2, 8, 4), (3, 1): (2, 9, 4), (2, 1, 1): (3, 11, 5), (1, 1, 1, 1): (4, 12, 5)}
CFG3 = {(3,): (2, 6, 2), (2, 1): (2, 8, 4), (1, 1, 1): (3, 9, 4)}


def to_sympy(c):
    return sp.expand(sum(sp.Rational(v.numerator, v.denominator) * r**a * s**b for (a, b), v in c.items()))


def _T(lam):
    T = {}
    for N in range(2, 40):
        t = B.terms_true(N, 1, lam)
        if t is not None:
            T[N] = t
    return T


def tq1111(al_free=-8):
    """(1,1,1,1)：r^{al_free} 以外を一次方程式で、r^{al_free} を t = q から決め、使う N すべてで照合する。"""
    from fractions import Fraction
    import rk_c_check_tq as C, rk_unrefined_check as rk
    lam = (1, 1, 1, 1); a, _, Ns = CFG4[lam]; b = 11              # N ≤ 11 では r^{−8} の部分が一意に決まらない
    T = _T(lam)
    _, d = B.slopes(T)
    used = [N for N in sorted(T) if a <= N <= b]
    info = []
    sol, st, _ = B.solve_bounded(T, used, *d[-1], "alpha", info, Ns, partial=True)
    free = [al for al, _, nv, rk_ in info if rk_ < nv]
    assert free == [al_free], free
    known = sum(to_sympy(cj) * C.U**j for j, cj in sol.items())
    known_tq = sp.expand(known.subs(s, r).subs(C.U, C.nu * r).subs(r, C.t))
    for sg in (1, -1):
        F = rk.nrz_trefoil(4, sg * C.nu, 1 / C.t, sp.Integer(-1))
        D = sp.expand(sp.cancel(C.mir(sp.cancel(F))) - known_tq)
        full = {j: dict(c) for j, c in sol.items()}
        for (j, k), v in sp.Poly(sp.expand(D * C.nu**40 * C.t**200), C.nu, C.t).terms():
            j, k = j - 40, k - 200
            full.setdefault(j, {})[(al_free, k - j - al_free)] = Fraction(int(v.p), int(v.q))
        ok = {N: B.predict(full, N) == T[N] for N in sorted(T)}
        js = sorted({j for j, c in full.items() for (a_, _) in c if a_ == al_free})
        print(f"[1^4] の符号 ν → {sg}ν：r^{al_free} の部分は U^{js[0]}〜U^{js[-1]}、N = {sorted(T)} の全データとの照合 {ok}", flush=True)
        if all(ok.values()):
            return full
    return None


def solve(lam, Nmin, Nmax, Ns):
    T = _T(lam)
    _, d = B.slopes(T)
    jmin, jmax = d[-1]
    used = [N for N in sorted(T) if Nmin <= N <= Nmax]
    sol, st, _ = B.solve_bounded(T, used, jmin, jmax, "alpha", Ns=Ns)
    checks = {N: B.predict(sol, N) == T[N] for N in T if N not in used} if sol is not None else {}
    return sol, st, used, checks


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "tq1111":
        full = tq1111()
        if full is not None:
            c = {j: to_sympy(cj) for j, cj in sorted(full.items()) if cj}
            path = os.path.join(OUT, "R_U_m1_1111.pkl")
            pickle.dump({j: str(v) for j, v in c.items()}, open(path, "wb"))
            print(f"   保存：{path}（U^{min(c)}〜U^{max(c)}。r^{{−8}} の部分は t = q から）", flush=True)
        sys.exit()
    three = len(sys.argv) > 1 and sys.argv[1] == "3box"
    for lam, (a, b, Ns) in (CFG3 if three else CFG4).items():
        sol, st, used, checks = solve(lam, a, b, Ns)
        name = "".join(map(str, lam))
        ok = sol is not None and checks and all(checks.values())
        print(f"λ={lam}：使う N = {used}（帯 N={Ns}）、{st}、使わない N での照合 {checks}", flush=True)
        if not ok:
            continue
        c = {j: to_sympy(cj) for j, cj in sorted(sol.items())}
        path = os.path.join(OUT, f"R_U_m1_{name}.pkl")
        if three:
            old = pickle.load(open(path, "rb"))
            same = set(old) | set(c) == set(c) and all(sp.expand(sp.sympify(old[j], locals={"r": r, "s": s}) - c.get(j, 0)) == 0 for j in old)
            print(f"   従来の R_U_m1_{name}.pkl と全係数一致：{same}", flush=True)
        else:
            pickle.dump({j: str(v) for j, v in c.items()}, open(path, "wb"))
            print(f"   保存：{path}（U^{min(c)}〜U^{max(c)}）", flush=True)
