"""composite [α;β] の大きな N の因子化（theory §12.36、続き 106）：
  R_{[α;β]}(ν) = R_α(ν) R_β(ν) + (q − 1)(t − 1) G_{αβ}(ν) が、随伴 [1;1]（§12.33 の A = R² + (q−1)(t−1)G）以外の composite でも成り立つか。
  R_α は被覆の彩色データ（鏡像、Rdata_nu）、R_{[α;β]} は GL の DAHA の補間（Rcomp_A_m1_{α}.pkl（[α;1]、ν の冪）、Rcomp22_A_m1_{α}_{β}.pkl（A = ν² の冪））。
  R_{[α;β]} の全体の単項式は最上行を積の最上行に合わせて決める（最上行が比例しなければ不成立）。
実行：PYTHONPATH=src:scripts/refined_kauffman python scripts/refined_kauffman/rk_f_comp_factor.py"""
import sys, os, json, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp

nu = sp.Symbol("nu")
r_, s_ = sp.symbols("r s", positive=True)          # r = q^{1/2}、s = t^{1/2}：以下はすべて r, s の Laurent 多項式
q, t = r_**2, s_**2
OUT = "data/refined_kauffman/daha"
SUB = {}


def R_data(c, knot="T2_3"):
    """被覆の彩色データの鏡像 Σ v ν^{−i} q^{−j/2} t^{−l/2}。"""
    d = json.load(open(f"data/framed_knot/torus/{knot}_{c}.json"))
    out = 0
    for key, v in d["poly"].items():
        i, j, l = map(int, key.split(","))
        out += v * nu**(-i) * r_**(-j) * s_**(-l)
    return sp.expand(out)


def R_comp(path, unit):
    d = pickle.load(open(path, "rb"))
    c = d["c"] if "c" in d else d
    return sp.expand(sum(sp.sympify(v, locals={"r": r_, "s": s_}).subs(SUB) * nu**(unit * int(j)) for j, v in c.items()))


def rows(e):
    """{ν の冪: 係数}（0 を除く）。"""
    P = sp.Poly(sp.expand(e * nu**200), nu)
    return {int(m[0]) - 200: sp.expand(c) for m, c in zip(P.monoms(), P.coeffs()) if c != 0}


def mono_ratio(a, b):
    """a/b が単項式ならそれを返す。"""
    rt = sp.cancel(a / b)
    n_, d_ = sp.fraction(rt)
    if sp.Poly(n_, r_, s_).is_monomial and sp.Poly(d_, r_, s_).is_monomial:
        return rt
    return None


def laurent_div(e, d):
    e = sp.expand(e)
    if e == 0:
        return sp.Integer(0)
    e2 = sp.expand(e * r_**120 * s_**120)
    Q, Rm = sp.div(sp.Poly(e2, r_, s_), sp.Poly(d, r_, s_))
    if not Rm.is_zero:
        return None
    return sp.expand(Q.as_expr() * r_**(-120) * s_**(-120))


def test_pair(name, Rab, Ra, Rb):
    ra, rb, rab = rows(Ra), rows(Rb), rows(Rab)
    prod = rows(sp.expand(Ra * Rb))
    top_ab, top_p = max(rab), max(prod)
    rt = mono_ratio(rab[top_ab], prod[top_p])
    print(f"== {name}：R_{{[α;β]}} の行 ν^{min(rab)}..ν^{top_ab}（{len(rab)} 行）、積の行 ν^{min(prod)}..ν^{top_p}（{len(prod)} 行）")
    if rt is None:
        print(f"   最上行は積の最上行に比例しない：{sp.factor(rab[top_ab])}  vs  {sp.factor(prod[top_p])}"[:300]); return
    shift = nu**(top_p - top_ab) / rt
    Rab2 = sp.expand(Rab * shift)
    diff = rows(sp.expand(Rab2 - Ra * Rb))
    print(f"   最上行は積の最上行 × {rt} × ν^{top_ab - top_p}：一致。差 R_{{[α;β]}} − R_α R_β の行：")
    allok = True
    for j in sorted(diff, reverse=True):
        g = laurent_div(diff[j], (q - 1) * (t - 1))
        ok = g is not None
        allok &= ok
        cs = sp.expand(g).as_coefficients_dict().values() if ok else []
        print(f"     ν^{j}：(q−1)(t−1) で割れる {ok}" + (f"、G の行 {len(cs)} 単項式 Σ|c| = {sum(abs(v) for v in cs)}：{sp.factor(g)}"[:200] if ok else f"：{sp.factor(diff[j])}"[:200]))
    print(f"   全行で (q−1)(t−1) で割れる：{allok}")


def G_rows(Rab, Ra, Rb):
    """R_{[α;β]} = R_α R_β + (q−1)(t−1) G の G を {ν の冪: 行} で返す（最上行で単項式を合わせる）。"""
    ra, rb, rab = rows(Ra), rows(Rb), rows(Rab)
    prod = rows(sp.expand(Ra * Rb))
    rt = mono_ratio(rab[max(rab)], prod[max(prod)])
    assert rt is not None
    Rab2 = sp.expand(Rab * nu**(max(prod) - max(rab)) / rt)
    return {j: laurent_div(c, (q - 1) * (t - 1)) for j, c in rows(sp.expand(Rab2 - Ra * Rb)).items()}, ra, rb


def g1_test(name, Rab, Ra, Rb):
    """G の最上行が ρ₀^α ρ₁^β + ρ₁^α ρ₀^β − ρ₁^α ρ₁^β × 単項式 か（[1;1] の g₁ = (2ρ₀ρ₁ − ρ₁²)/q の一般形）。"""
    G, ra, rb = G_rows(Rab, Ra, Rb)
    a0, b0 = max(ra), max(rb)
    X = sp.expand(ra[a0] * rb[b0 - 2] + ra[a0 - 2] * rb[b0] - ra[a0 - 2] * rb[b0 - 2])
    top = max(G)
    rt = mono_ratio(G[top], X)
    print(f"   {name}：G の最上行（ν^{top}）/ (ρ₀^αρ₁^β + ρ₁^αρ₀^β − ρ₁^αρ₁^β) = {rt}")
    return rt



if __name__ == "__main__":
    knot = "T2_3"
    w = s_ / r_                              # (t/q)^{δ/2}、δ = 1
    R1, R2, R11 = (sp.expand(w**k * R_data(c, knot)) for c, k in (("1", 1), ("2", 2), ("11", 2)))
    radj = pickle.load(open(f"{OUT}/Radj_A_m1.pkl", "rb"))
    Radj = sp.expand(sum(sp.sympify(v, locals={"r": r_, "s": s_}).subs(SUB) * nu**int(j) for j, v in radj.items()))
    test_pair("[1;1]（随伴）", Radj, R1, R1)
    test_pair("[2;1]", R_comp(f"{OUT}/Rcomp_A_m1_2.pkl", 1), R2, R1)
    test_pair("[1,1;1]", R_comp(f"{OUT}/Rcomp_A_m1_11.pkl", 1), R11, R1)
    test_pair("[2;2]", R_comp(f"{OUT}/Rcomp22_A_m1_2_2.pkl", 2), R2, R2)
    test_pair("[2;1,1]", R_comp(f"{OUT}/Rcomp22_A_m1_2_11.pkl", 2), R2, R11)
    test_pair("[1,1;1,1]", R_comp(f"{OUT}/Rcomp22_A_m1_11_11.pkl", 2), R11, R11)
    print("== g₁ の一般形：G の最上行 vs ρ₀^αρ₁^β + ρ₁^αρ₀^β − ρ₁^αρ₁^β")
    g1_test("[1;1]", Radj, R1, R1)
    g1_test("[2;1]", R_comp(f"{OUT}/Rcomp_A_m1_2.pkl", 1), R2, R1)
    g1_test("[1,1;1]", R_comp(f"{OUT}/Rcomp_A_m1_11.pkl", 1), R11, R1)
    g1_test("[2;2]", R_comp(f"{OUT}/Rcomp22_A_m1_2_2.pkl", 2), R2, R2)
    g1_test("[2;1,1]", R_comp(f"{OUT}/Rcomp22_A_m1_2_11.pkl", 2), R2, R11)
    g1_test("[1,1;1,1]", R_comp(f"{OUT}/Rcomp22_A_m1_11_11.pkl", 2), R11, R11)
