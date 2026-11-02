"""Ω_λ の構成（続き 118、theory §12.48）：共有の核も記号にした記号ブロック（n2_rows_ker_shared_m1_k14.pkl）から、
各 (λ, K) で「共有の核も G も含まない純粋な積の単項式」を取り出し、brane の単項式 × image の単項式の双線形形式として階数を調べる。
階数 1 なら N_λ^{prod} = c_K · Ω_K^b Ω_K^i で、Ω_K（値）を取り出して ν^{−12} の X_λ・Ŷ_λ と比べる。
実行：PYTHONPATH=src:scripts/refined_kauffman python scripts/refined_kauffman/rk_f_n2_omega.py"""
import sys, pickle, itertools, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_f_n2_rows_sym as RS
from rk_f_n2_rows_sym import r, s, q, t
import rk_f_n2_grows as G

PKL = "data/refined_kauffman/daha/n2_rows_ker_shared_m1_k14.pkl"


def load():
    d = pickle.load(open(PKL, "rb"))
    safe = lambda n: n.replace("_-", "_m")          # 負の準位の記号名 D11_-2b は読み戻せないので置き換える
    syms = {}
    for name, (side, base, v) in d["syms"].items():
        syms[safe(name)] = (side, base, sp.sympify(v, locals={"r": r, "s": s}))
    order = [safe(n) for n in d["order"]]
    gens = [sp.Symbol(n) for n in order]
    loc = {n: g for n, g in zip(order, gens)}; loc.update({"r": r, "s": s})
    N = {k: sp.Poly(sp.expand(sp.sympify(safe(v), locals=loc)), *gens, domain=RS.DOM) for k, v in d["N"].items()}
    return N, syms, order, gens


def classify(poly, syms, order):
    """単項式 → (brane の内容（基底名のタプル）, image の内容, 共有を含む, G を含む, 係数)。"""
    out = []
    for mon, c in poly.terms():
        cb, ci, sh, g = [], [], False, False
        for name, e in zip(order, mon):
            if not e:
                continue
            if name not in syms:
                g = True; continue
            side, base, v = syms[name]
            if v == 1:
                continue                      # 値 1 の記号（M0_0：空の分割の係数）は無視する
            if side == "s":
                sh = True
            elif side == "b":
                cb += [name] * e
            else:
                ci += [name] * e
        ce = sp.sympify(c).as_expr() if hasattr(c, "as_expr") else sp.sympify(c)
        out.append((tuple(sorted(cb)), tuple(sorted(ci)), sh, g, ce))
    merged = {}
    for cb, ci, sh, g, ce in out:
        if g:
            merged[("G", len(merged))] = (cb, ci, sh, g, ce); continue
        key = (cb, ci, sh)
        if key in merged:
            cb0, ci0, sh0, g0, ce0 = merged[key]; merged[key] = (cb, ci, sh, g, sp.cancel(ce0 + ce))
        else:
            merged[key] = (cb, ci, sh, g, ce)
    return list(merged.values())


def main():
    t0 = time.time()
    N, syms, order, gens = load()
    print(f"読み込み {time.time() - t0:.0f} 秒", flush=True)
    val = lambda names: sp.prod([syms[n][2] for n in names]) if names else sp.Integer(1)
    pts = [{r: sp.Rational(3, 2), s: sp.Rational(5, 7)}, {r: sp.Rational(7, 3), s: sp.Rational(2, 11)}]
    X = {lab: sp.cancel(G.X_of(lab).subs(RS.actual_rows())) for lab in ("2", "11")}
    for lab in ("11", "2"):
        for K in (10, 12, 14):
            if (lab, K) not in N:
                continue
            cls = classify(N[(lab, K)], syms, order)
            prod = [(cb, ci, ce) for cb, ci, sh, g, ce in cls if not sh and not g and (cb or ci)]
            n_sh = sum(1 for cb, ci, sh, g, ce in cls if sh and not g); n_g = sum(1 for cb, ci, sh, g, ce in cls if g); n_bg = sum(1 for cb, ci, sh, g, ce in cls if not cb and not ci and not g)
            bs = sorted({cb for cb, ci, ce in prod}); is_ = sorted({ci for cb, ci, ce in prod})
            M = sp.zeros(len(bs), len(is_))
            for cb, ci, ce in prod:
                M[bs.index(cb), is_.index(ci)] += ce
            rk = [M.subs(pt).rank() for pt in pts]
            print(f"== λ = {lab}、ν^-{K}：単項式 {len(cls)}（純粋な積 {len(prod)}、共有の核を含む {n_sh}、G {n_g}、記号なし {n_bg}）；積の部分の brane 基底 {len(bs)}・image 基底 {len(is_)}、双線形形式の階数 {rk}", flush=True)
            print(f"   brane の基底：{[tuple(syms[n][1] for n in cb) for cb in bs]}", flush=True)
            if rk == [1, 1]:
                # M = c · u uᵀ（対称と仮定）：u を第 1 行から
                i0 = 0
                u = [sp.cancel(M[i, i0] / M[i0, i0]) for i in range(len(bs))] if M[i0, i0] != 0 else None
                if u is not None:
                    c = sp.cancel(M[i0, i0])
                    Omega = sp.cancel(sum(ui * val(cb) for ui, cb in zip(u, bs)))
                    print(f"   N^prod = c·Ω^b Ω^i、c = {sp.factor(c)}、Ω = {sp.factor(Omega)}", flush=True)
                    print(f"   c·Ω² = {sp.factor(sp.cancel(c * Omega**2))}", flush=True)
                    if K == 12:
                        print(f"   Ω/X_λ = {sp.factor(sp.cancel(Omega / X[lab]))}；c·Ω² + (s³/r)X_λ² = {sp.cancel(c * Omega**2 + s**3 / r * X[lab]**2)}", flush=True)
            print(f"   （{time.time() - t0:.0f} 秒）", flush=True)


if __name__ == "__main__" and "--rank" not in sys.argv:
    main()


def bilinear(N, syms, order, lab, K, pts):
    """G を含まない単項式を brane の内容 × image の内容の双線形形式にする（共有の核の記号は値を掛けて係数に入れる）。"""
    poly = N[(lab, K)]
    M = {}
    for mon, c in poly.terms():
        cb, ci, g, shv = [], [], False, sp.Integer(1)
        for name, e in zip(order, mon):
            if not e:
                continue
            if name not in syms:
                g = True; break
            side, base, v = syms[name]
            if v == 1:
                continue
            if side == "s":
                shv *= v**e
            elif side == "b":
                cb += [name] * e
            else:
                ci += [name] * e
        if g:
            continue
        ce = sp.sympify(c).as_expr() if hasattr(c, "as_expr") else sp.sympify(c)
        key = (tuple(sorted(cb)), tuple(sorted(ci)))
        M[key] = M.get(key, 0) + ce * shv
    bs = sorted({k[0] for k in M}); is_ = sorted({k[1] for k in M})
    Mat = sp.zeros(len(bs), len(is_))
    for (cb, ci), v in M.items():
        Mat[bs.index(cb), is_.index(ci)] = sp.cancel(v)
    rk = [Mat.subs(pt).rank() for pt in pts]
    return Mat, bs, is_, rk


def part_rank():
    t0 = time.time()
    N, syms, order, gens = load()
    pts = [{r: sp.Rational(3, 2), s: sp.Rational(5, 7)}, {r: sp.Rational(7, 3), s: sp.Rational(2, 11)}]
    val = lambda names: sp.prod([syms[n][2] for n in names]) if names else sp.Integer(1)
    X = {lab: sp.cancel(G.X_of(lab).subs(RS.actual_rows())) for lab in ("2", "11")}
    for lab in ("11", "2"):
        for K in (10, 12, 14):
            Mat, bs, is_, rk = bilinear(N, syms, order, lab, K, pts)
            print(f"== λ = {lab}、ν^-{K}：G を含まない部分の双線形形式（共有の核は値）：brane 基底 {len(bs)}・image 基底 {len(is_)}、階数 {rk}（{time.time() - t0:.0f} 秒）", flush=True)
            # 対称性と、brane の内容ごとの値（行ベクトルの値）
            sym_ok = (set(bs) == set(is_)) and all(sp.cancel(Mat[i, j] - Mat[is_.index(bs[i]), bs.index(is_[j])]) == 0 for i in range(len(bs)) for j in range(len(is_)) if bs[i] in is_ and is_[j] in bs)
            print(f"   対称：{sym_ok}；brane 基底：{[tuple(syms[n][1] for n in cb) for cb in bs]}", flush=True)
            # 値の双線形形式：brane の値 × image の値を掛けた総和 = N の G を含まない部分
            tot = sp.cancel(sum(Mat[i, j] * val(bs[i]) * val(is_[j]) for i in range(len(bs)) for j in range(len(is_))))
            print(f"   G を含まない部分の値 = {sp.factor(tot)}"[:400], flush=True)
            if rk[0] <= 3:
                # 階数分解：M の列空間から Ω の候補（行の内容の値の線形結合）
                Mp = Mat.subs(pts[0])
                piv = Mp.rref()[1]
                print(f"   ピボット列（image 基底）：{[tuple(syms[n][1] for n in is_[j]) for j in piv]}", flush=True)


if __name__ == "__main__" and "--rank" in sys.argv:
    part_rank()
