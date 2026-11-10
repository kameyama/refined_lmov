r"""rk_c_bounded_interp（r の次数ごとに U の冪の範囲を絞る補間、flint の約分）で SO の 2 箱の R_λ(U) を決め、
rk_c_read_largeN と同じ形式 {j: str(c_j)} で data/refined_kauffman/daha/R_U_m{m}_{λ}.pkl に保存する（progress.md 続き 94）。
rk_c_read_largeN（sympy の約分、Newton の補間）は m = 3 で 3 時間かかり、m = 4 では日の単位になる見込みなので、こちらを使う。
解が一意で、使わない N すべてで全係数一致したときだけ保存する。保存済みの R_U があれば全係数で比較して報告する（m = 3 で検証）。
随伴 R_adj(A^{1/2})（rk_c_2box_largeN_m.adj_interp の Newton 補間は大きな N_A のノードが要る）も同じ方法で決める（--adj）：
  R_adj(N) = W_adj(N) / (qn² − κ) を flint で約分し、x = A^{1/2} = s^N の冪 j ごとの係数を r の次数ごとに解く（UNIT = 1、偶奇の制限なし）。
使い方：python rk_c_bounded_save.py m λ... [--nmax=N]（例：4 1 2 11 --nmax=20）、python rk_c_bounded_save.py m --adj [--nmax=N]、
  composite (α;1̄)：python rk_c_bounded_save.py m --comp=21（Rcomp_N*_m{m}_{α} から Rcomp_A_m{m}_{α}.pkl、続き 139）"""
import sys, os, pickle, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_c_bounded_interp as B

rr, ss = sp.symbols("r s", positive=True)


def to_sympy(c):
    return sp.Add(*[sp.Rational(v.numerator, v.denominator) * rr**a * ss**b for (a, b), v in c.items()])


def run(m, lam, Nmax_used=None, k=2):
    t0 = time.time()
    T = {N: t for N in range(2, 40) if (t := B.terms_true(N, m, lam, k)) is not None}
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
                path = os.path.join(B.OUT, f"R_U_{B.TG.tag(m, k)}_{''.join(map(str, lam))}.pkl")
                if os.path.exists(path):
                    old = pickle.load(open(path, "rb"))
                    same = all(sp.expand(sp.sympify(old.get(j, 0), locals={"r": rr, "s": ss}) - c.get(j, 0)) == 0 for j in set(old) | set(c))   # 係数 0 の j は保存の有無が違ってよい
                    print(f"   保存済み {path} と全係数一致：{same}", flush=True)
                else:
                    pickle.dump({j: str(v) for j, v in c.items()}, open(path, "wb"))
                    print(f"   保存：{path}（{time.time() - t0:.0f} 秒）", flush=True)
                return
    print("   採用できる組が無い", flush=True)


def terms_adj(N, m, k=2):
    """随伴 R_adj(N) = W / (qn² − κ) の項 {(α, β): Fraction}（flint で約分。data/cache にキャッシュ）。"""
    from fractions import Fraction
    from rcsknots.refined_kauffman.flint_rf import from_sympy, RF, R, S, conv
    path = os.path.join(B.OUT, f"Wadj_N{N}_{B.TG.tag(m, k)}.pkl")
    if not os.path.exists(path):
        return None
    cpath = os.path.join(B.CACHE, f"Radj_N{N}_{B.TG.tag(m, k)}.pkl")
    if os.path.exists(cpath) and os.path.getmtime(cpath) >= os.path.getmtime(path):
        return pickle.load(open(cpath, "rb"))
    W = from_sympy(pickle.load(open(path, "rb"))["W"])
    A2 = S**N
    qn = (A2 - 1 / A2) / (S - 1 / S)
    kappa = (1 - S**(2 * N)) * (1 - R * R) / ((1 - S * S) * (1 - R * R * S**(2 * N - 2)))
    Rad = W / (qn * qn - kappa)
    assert Rad.is_laurent(), (N, m)
    (dm, dc), = Rad.d.to_dict().items()
    out = {(e[0] - dm[0], e[1] - dm[1]): Fraction(int(c), int(dc)) for e, c in Rad.n.to_dict().items()}
    os.makedirs(B.CACHE, exist_ok=True)
    pickle.dump(out, open(cpath, "wb"))
    return out


def terms_comp(N, m, alpha):
    """composite (α;1̄) の簡約した値 R_{(α;1)}(N)（Rcomp_N{N}_m{m}_{α}.pkl、rk_c_fast_composite*.py）の項 {(α, β): Fraction}（続き 139）。"""
    from fractions import Fraction
    from rcsknots.refined_kauffman.flint_rf import from_sympy
    path = os.path.join(B.OUT, f"Rcomp_N{N}_m{m}_{''.join(map(str, alpha))}.pkl")
    if not os.path.exists(path):
        return None
    W = from_sympy(pickle.load(open(path, "rb"))["R"])
    assert W.is_laurent(), (N, m, alpha)
    (dm, dc), = W.d.to_dict().items()
    return {(e[0] - dm[0], e[1] - dm[1]): Fraction(int(c), int(dc)) for e, c in W.n.to_dict().items()}


def run_adj(m, Nmax_used=None, k=2, comp=None):
    """comp = α のときは随伴の代わりに composite (α;1̄) を A^{1/2} の冪で決め、Rcomp_A_m{m}_{α}.pkl に保存する（続き 139）。"""
    t0 = time.time()
    B.UNIT, B.PAR = 1, 1
    if comp is None:
        T = {N: t for N in range(2, 50) if (t := terms_adj(N, m, k)) is not None}
    else:
        T = {N: t for N in range(2, 50) if (t := terms_comp(N, m, comp)) is not None}
    rng, d = B.slopes(T)
    jmin, jmax = d[-1]
    assert d[-1] == d[-2], ("傾きが一定でない", d)
    Ns = sorted(T)
    print(f"{'随伴' if comp is None else 'composite ' + str(comp)} m={m}：N_A = {Ns}、A^{{1/2}} の冪 {jmin}〜{jmax}（{jmax - jmin + 1} 項）、約分 {time.time() - t0:.0f} 秒", flush=True)
    for Nmin in range(Ns[0], Ns[0] + 10):
        for Nmax in range(Nmin + 1, (Nmax_used or Ns[-1]) + 1):
            used = [N for N in Ns if Nmin <= N <= Nmax]
            sol, st, nvar = B.solve_bounded(T, used, jmin, jmax, "all")
            if sol is None:
                if "一意でない" not in st:
                    print(f"   使う N = {Nmin}..{Nmax}：{st}", flush=True); break
                continue
            ok = {N: B.predict(sol, N) == T[N] for N in Ns if N not in used}
            print(f"   使う N = {Nmin}..{Nmax}：{st}、未知数 {nvar}、使わない N での照合 {ok}", flush=True)
            if ok and all(ok.values()):
                c = {j: to_sympy(v) for j, v in sol.items()}
                path = os.path.join(B.OUT, f"Radj_A_{B.TG.tag(m, k)}.pkl" if comp is None else f"Rcomp_A_m{m}_{''.join(map(str, comp))}.pkl")
                if os.path.exists(path):
                    old = pickle.load(open(path, "rb"))
                    same = all(sp.expand(sp.sympify(old.get(j, 0), locals={"r": rr, "s": ss}) - c.get(j, 0)) == 0 for j in set(old) | set(c))   # 係数 0 の j は保存の有無が違ってよい
                    print(f"   保存済み {path} と全係数一致：{same}", flush=True)
                else:
                    pickle.dump({j: str(v) for j, v in c.items()}, open(path, "wb"))
                    print(f"   保存：{path}（{time.time() - t0:.0f} 秒）", flush=True)
                B.UNIT, B.PAR = 2, 2
                return
    B.UNIT, B.PAR = 2, 2
    print("   採用できる組が無い", flush=True)


def run_fixed(m, lam, mode, Nmin, Nmax, Ns, comp=False):
    """使う N（Nmin..Nmax）・mode・帯を取る N_s を指定して解き、一意で使わない N すべてで全係数一致したときだけ保存する（続き 139）。
    comp = True なら λ は composite (λ;1̄) の α（A^{1/2} の冪、Rcomp_A_m{m}_{λ}.pkl）。"""
    t0 = time.time()
    if comp:
        B.UNIT, B.PAR = 1, 1
        T = {N: t for N in range(2, 50) if (t := terms_comp(N, m, lam)) is not None}
        path = os.path.join(B.OUT, f"Rcomp_A_m{m}_{''.join(map(str, lam))}.pkl")
    else:
        T = {N: t for N in range(2, 40) if (t := B.terms_true(N, m, lam)) is not None}
        path = os.path.join(B.OUT, f"R_U_{B.TG.tag(m)}_{''.join(map(str, lam))}.pkl")
    rng, d = B.slopes(T)
    jmin, jmax = d[-1]
    assert d[-1] == d[-2], ("傾きが一定でない", d)
    Ns_all = sorted(T)
    used = [N for N in Ns_all if Nmin <= N <= Nmax]
    sol, st, nvar = B.solve_bounded(T, used, jmin, jmax, mode, Ns=Ns)
    ok = {N: B.predict(sol, N) == T[N] for N in Ns_all if N not in used} if sol is not None else {}
    print(f"m={m} {'composite ' if comp else ''}{lam}：冪 {jmin}〜{jmax}、mode {mode}、使う N = {Nmin}..{Nmax}、N_s = {Ns}：{st}、未知数 {nvar}、"
          f"使わない N での照合 {ok}", flush=True)
    B.UNIT, B.PAR = 2, 2
    if sol is None or not ok or not all(ok.values()):
        print("   保存しない", flush=True)
        return False
    c = {j: to_sympy(v) for j, v in sol.items()}
    if os.path.exists(path):
        old = pickle.load(open(path, "rb"))
        same = all(sp.expand(sp.sympify(old.get(j, 0), locals={"r": rr, "s": ss}) - c.get(j, 0)) == 0 for j in set(old) | set(c))
        print(f"   保存済み {path} と全係数一致：{same}", flush=True)
    else:
        pickle.dump({j: str(v) for j, v in c.items()}, open(path, "wb"))
        print(f"   保存：{path}（{time.time() - t0:.0f} 秒）", flush=True)
    return True


def run_adaptive(m, lam, Nmax, Nmin=2):
    """r の次数 α ごとに帯の仮定を選んで解く（続き 141）：mode alpha、次に alphaA の順に、帯を取る N₀ = Nmin, Nmin + 1, … を試し、
    N₀ 以上の使う N（Nmin..Nmax）で一意に解けた最初のものを採る。使わない N（Nmin 未満と Nmax より大）がすべて全係数一致したときだけ保存する。"""
    t0 = time.time()
    T = {N: t for N in range(2, 40) if (t := B.terms_true(N, m, lam)) is not None}
    path = os.path.join(B.OUT, f"R_U_{B.TG.tag(m)}_{''.join(map(str, lam))}.pkl")
    rng, d = B.slopes(T)
    jmin, jmax = d[-1]
    assert d[-1] == d[-2], ("傾きが一定でない", d)
    used = [N for N in sorted(T) if Nmin <= N <= Nmax]
    alphas = sorted({a for N in used for a, _ in T[N]})
    sol, choice, fail = {}, {}, []
    for al in alphas:
        got = None
        for mode in ("alpha", "alphaA"):
            for N0 in range(used[0], used[-1] - 1):
                s_, st, _ = B.solve_bounded(T, [N for N in used if N >= N0], jmin, jmax, mode, Ns=N0, alphas=[al])
                if s_ is not None:
                    got = s_; choice[al] = (mode, N0); break
            if got is not None:
                break
        if got is None:
            fail.append(al); continue
        for j, c in got.items():
            sol.setdefault(j, {}).update(c)
    ok = {N: B.predict(sol, N) == T[N] for N in sorted(T) if N not in used} if not fail else {}
    from collections import Counter
    print(f"{B.TG.tag(m)} {lam}：冪 {jmin}〜{jmax}、使う N = {Nmin}..{Nmax}、帯の選択 {dict(Counter(choice.values()))}、解けない α {fail}、"
          f"使わない N での照合 {ok}（{time.time() - t0:.0f} 秒）", flush=True)
    if fail or not ok or not all(ok.values()):
        print("   保存しない", flush=True)
        return False
    c = {j: to_sympy(v) for j, v in sol.items()}
    pickle.dump({j: str(v) for j, v in c.items()}, open(path, "wb"))
    print(f"   保存：{path}", flush=True)
    return True


if __name__ == "__main__":
    m = int(sys.argv[1]) if sys.argv[1].isdigit() else sys.argv[1]          # 文字列は GAMMA_WORDS の語（例 T3_5、続き 103）
    Nmax_used = None
    lams = []; adj = False; k = 2; comps = []
    for a in sys.argv[2:]:
        if a.startswith("--nmax="):
            Nmax_used = int(a[7:])
        elif a.startswith("--k="):            # γ = τ_+^m τ_−^k（T(mk+1, k)、続き 102）
            k = int(a[4:])
        elif a == "--adj":
            adj = True
        elif a.startswith("--comp="):         # composite (α;1̄)（例 --comp=21、続き 139）
            comps.append(tuple(int(x) for x in a[7:]))
        else:
            lams.append(tuple(int(x) for x in a))
    for lam in lams:
        run(m, lam, Nmax_used, k)
    if adj:
        run_adj(m, Nmax_used, k)
    for al in comps:
        run_adj(m, Nmax_used, k, comp=al)
