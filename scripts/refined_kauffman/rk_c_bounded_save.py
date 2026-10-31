r"""rk_c_bounded_interp（r の次数ごとに U の冪の範囲を絞る補間、flint の約分）で SO の 2 箱の R_λ(U) を決め、
rk_c_read_largeN と同じ形式 {j: str(c_j)} で data/refined_kauffman/daha/R_U_m{m}_{λ}.pkl に保存する（progress.md 続き 94）。
rk_c_read_largeN（sympy の約分、Newton の補間）は m = 3 で 3 時間かかり、m = 4 では日の単位になる見込みなので、こちらを使う。
解が一意で、使わない N すべてで全係数一致したときだけ保存する。保存済みの R_U があれば全係数で比較して報告する（m = 3 で検証）。
使い方：python rk_c_bounded_save.py m λ... [--nmax=N]（例：4 1 2 11 --nmax=20）"""
import sys, os, pickle, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_c_bounded_interp as B

rr, ss = sp.symbols("r s", positive=True)


def to_sympy(c):
    return sp.Add(*[sp.Rational(v.numerator, v.denominator) * rr**a * ss**b for (a, b), v in c.items()])


def run(m, lam, Nmax_used=None):
    t0 = time.time()
    T = {N: t for N in range(2, 40) if (t := B.terms_true(N, m, lam)) is not None}
    rng, d = B.slopes(T)
    jmin, jmax = d[-1]
    assert d[-1] == d[-2], ("傾きが一定でない", d)
    Ns = sorted(T)
    print(f"m={m} λ={lam}：N = {Ns}、U^{jmin}〜U^{jmax}（{jmax - jmin + 1} 項）、約分 {time.time() - t0:.0f} 秒", flush=True)
    # 使う N の下限 Nmin と上限 Nmax を動かし、一意で使わない N での照合が全部通る最初の組を採る（小さい N は帯の仮定が破れることがある）
    for Nmin in range(Ns[0], Ns[0] + 8):
        for Nmax in range(Nmin + 1, (Nmax_used or Ns[-1]) + 1):
            used = [N for N in Ns if Nmin <= N <= Nmax]
            sol, st, nvar = B.solve_bounded(T, used, jmin, jmax, "all")
            if sol is None:
                if "一意でない" not in st:
                    print(f"   使う N = {Nmin}..{Nmax}：{st}", flush=True); break        # 帯の矛盾：Nmin を上げる
                continue
            ok = {N: B.predict(sol, N) == T[N] for N in Ns if N not in used}
            print(f"   使う N = {Nmin}..{Nmax}：{st}、未知数 {nvar}、使わない N での照合 {ok}", flush=True)
            if ok and all(ok.values()):
                c = {j: to_sympy(v) for j, v in sol.items()}
                path = os.path.join(B.OUT, f"R_U_m{m}_{''.join(map(str, lam))}.pkl")
                if os.path.exists(path):
                    old = pickle.load(open(path, "rb"))
                    same = set(old) == set(c) and all(sp.expand(sp.sympify(old[j], locals={"r": rr, "s": ss}) - c[j]) == 0 for j in c)
                    print(f"   保存済み {path} と全係数一致：{same}", flush=True)
                else:
                    pickle.dump({j: str(v) for j, v in c.items()}, open(path, "wb"))
                    print(f"   保存：{path}（{time.time() - t0:.0f} 秒）", flush=True)
                return
    print("   採用できる組が無い", flush=True)


if __name__ == "__main__":
    m = int(sys.argv[1])
    Nmax_used = None
    lams = []
    for a in sys.argv[2:]:
        if a.startswith("--nmax="):
            Nmax_used = int(a[7:])
        else:
            lams.append(tuple(int(x) for x in a))
    for lam in lams:
        run(m, lam, Nmax_used)
