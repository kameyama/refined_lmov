"""n = 1 の N を [1] と随伴の a の行で書く母関数（続き 103、theory §12.33）。
x = ν^{−2}。[1] の鏡像 R_1 = q^{−δ/2} t^{−δ/2} x^δ R(x)、R(y) = Σ_j ρ_j y^j（ρ_0 は最上行、符号つき）。
随伴 R_adj = q^{−2δ} x^{2δ} A(x)、A(y) = Σ_j α_j y^j。κ の項は Sym^k V の系列。
【導出】N = (1 − t)[R_adj(q_n² − κ) + κ − w²R_1²P_1²]/D を x で展開すると
  (q − 1)(t − 1) Σ_k (n_k − Sym^k V) x^k = x^{2δ−1}[A(qx)/(1 − tx) − R(qx)²(1 − qx)/((1 − x)(1 − qtx))]。
N が Laurent なら A = R² + (q − 1)(t − 1) G（G は Laurent）で、
  Σ_k n_k x^k = [1 − x^{2δ} R(qx)²]/((1 − x)(1 − tx)(1 − qtx)) + x^{2δ−1} G(qx)/(1 − tx)。
移り目 n_{2δ} = Sym^{2δ}V − (ρ_0 − ρ_1)² は g_0 = 0、g_1 = (2ρ_0ρ_1 − ρ_1²)/q と同値。

実行：PYTHONPATH=src:scripts/refined_kauffman python scripts/refined_kauffman/rk_f_rows_gen.py [m [--k=K] | T3_5] ...
"""
import sys, json, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_torus_tag as TG

q, t, x, y = sp.symbols("q t x y")
r_, s_ = sp.symbols("r s", positive=True)
OUT = "data/refined_kauffman/daha"


def sym_V(k):
    return sum(q**a * t**(a + b) for a in range(k + 1) for b in range(k + 1 - a))


def rows_R1(m, k=2):
    """[1] の鏡像の行 {i: 係数}（ν^{−i} の係数、r, s の Laurent 多項式）。"""
    d = json.load(open(f"data/framed_knot/torus/{TG.knot_name(m, k)}_1.json"))
    P = {}
    for key, v in d["poly"].items():
        i, j, l = map(int, key.split(","))
        P[i] = P.get(i, 0) + v * r_**(-j) * s_**(-l)
    return {i: sp.expand(e) for i, e in P.items()}


def R_poly(m, k=2):
    """R(y) = Σ ρ_j y^j（q, t の Laurent 多項式）。R_1 = q^{−δ/2} t^{−δ/2} x^δ R(x)。"""
    de = TG.delta(m, k)
    P = rows_R1(m, k)
    assert min(P) == 2 * de, (min(P), de)
    return sp.expand(sum(e * r_**de * s_**de * y**((i - 2 * de) // 2) for i, e in P.items()).subs({r_: sp.sqrt(q), s_: sp.sqrt(t)}))


def A_poly(m, k=2):
    """A(y) = Σ α_j y^j。R_adj = q^{−2δ} x^{2δ} A(x)（Radj_A の鍵は ν の冪）。"""
    de = TG.delta(m, k)
    radj = pickle.load(open(f"{OUT}/Radj_A_{TG.tag(m, k)}.pkl", "rb"))
    A = 0
    for j, v in radj.items():
        e = sp.sympify(v, locals={"r": r_, "s": s_})
        if e == 0:
            continue
        assert j % 2 == 0 and -j // 2 - 2 * de >= 0, (j, de)
        A += e * y**(-j // 2 - 2 * de)
    return sp.expand(sp.expand(A * q**(2 * de)).subs({r_: sp.sqrt(q), s_: sp.sqrt(t)}))


def _coeffs(P, var):
    Pp = sp.Poly(P, var)
    return {int(mm[0]): c for mm, c in zip(Pp.monoms(), Pp.coeffs())}


def _laurent_div(c, d):
    """Laurent 多項式 c（q, t）を多項式 d で割り切る（割り切れなければ例外）。"""
    c = sp.expand(c)
    if c == 0:
        return sp.Integer(0)
    Pc = sp.Poly(c * q**60 * t**60, q, t)
    Q, Rm = sp.div(Pc, sp.Poly(d, q, t))
    assert Rm.is_zero, "割り切れない"
    return sp.expand(Q.as_expr() * q**(-60) * t**(-60))


def G_poly(m, k=2):
    """G = (A − R²)/((q − 1)(t − 1))。Laurent でなければ例外。"""
    A, R = A_poly(m, k), R_poly(m, k)
    diff = _coeffs(sp.expand(A - sp.expand(R**2)), y)
    G = sum(_laurent_div(c, (q - 1) * (t - 1)) * y**j for j, c in diff.items())
    return sp.expand(G), A, R


def blocks_from_rows(m, k=2, kmax=None):
    """母関数から n_k（k ≤ kmax）。"""
    de = TG.delta(m, k)
    G, A, R = G_poly(m, k)
    K = TG.rs(m, k) - 3
    kmax = K + 2 if kmax is None else kmax
    # 畳み込みで展開：n_k = Sym^k V − Σ_{a+b=k−2δ} [R²]_a q^a Sym^b V + Σ_{j+n=k−2δ+1} g_j q^j t^n
    R2 = _coeffs(sp.expand(R**2), y)
    g = _coeffs(G, y) if G != 0 else {}
    out = {}
    for kk in range(kmax + 1):
        e = sym_V(kk)
        e -= sum(R2[a] * q**a * sym_V(kk - 2 * de - a) for a in R2 if 0 <= kk - 2 * de - a)
        e += sum(g[j] * q**j * t**(kk - 2 * de + 1 - j) for j in g if 0 <= kk - 2 * de + 1 - j)
        out[kk] = sp.expand(e)
    return out, G, A, R


def direct_blocks(m, k=2):
    """cover_sectors の x セクターから直接 N = (1 − t)x/(s²(t − ν²)(ν² − 1)(q − 1)(ν²q − 1)) を作りブロックに分ける（SO の 2 箱は不要）。"""
    import rk_b_cover_split as C
    import rk_f_frozen_basis as FB
    import io, contextlib
    u = sp.Symbol("u", positive=True)
    nu = FB.nu
    with contextlib.redirect_stdout(io.StringIO()):
        sec = C.cover_sectors(m, k)
    tt, qq = s_**2, r_**2
    P = sp.cancel((1 - tt) * sp.cancel(sec["x"].to_sympy().subs(u, nu * s_)))
    N = sp.cancel(P / (s_**2 * (tt - nu**2) * (nu**2 - 1) * (qq - 1) * (nu**2 * qq - 1)))
    bl = FB.nu_blocks(N)
    return {(-b - 6) // 2: sp.expand((blk * r_**(4 + 2 * ((-b - 6) // 2))).subs({r_: sp.sqrt(q), s_: sp.sqrt(t)})) for b, blk in bl.items()}


def main(m, k=2, check=True):
    de, K = TG.delta(m, k), TG.rs(m, k) - 3
    nk, G, A, R = blocks_from_rows(m, k)
    print(f"== {TG.knot_name(m, k)}（δ = {de}、K = rs − 3 = {K}）")
    Pr = sp.Poly(R, y); Pa = sp.Poly(A, y); Pg = sp.Poly(G, y) if G != 0 else None
    print(f"  R(y) の行 ρ_j：", {int(mm[0]): sp.factor(c) for mm, c in zip(Pr.monoms(), Pr.coeffs())})
    print(f"  A(y) の次数 {Pa.degree()}、G(y) = (A − R²)/((q−1)(t−1)) の行 g_j：")
    if Pg is not None:
        for mm, c in sorted(zip(Pg.monoms(), Pg.coeffs())):
            cs = list(sp.expand(c).as_coefficients_dict().values()) if c != 0 else []
            print(f"    g_{mm[0]} = {sp.factor(c)}   （{len(cs)} 単項式、Σ|c| = {sum(abs(v) for v in cs)}）")
    rho = {int(mm[0]): c for mm, c in zip(Pr.monoms(), Pr.coeffs())}
    St = sp.expand(rho.get(0, 0) - rho.get(1, 0))
    print(f"  ρ_0 − ρ_1 = {sp.factor(St)}")
    print(f"  g_0 = 0：{sp.expand(G.coeff(y, 0)) == 0}、g_1 = (2ρ_0ρ_1 − ρ_1²)/q：{sp.expand(G.coeff(y, 1) - (2 * rho[0] * rho[1] - rho[1]**2) / q) == 0}")
    print(f"  移り目 n_{2*de} − Sym^{2*de}V + (ρ_0 − ρ_1)² = {sp.cancel(nk[2 * de] - sym_V(2 * de) + St**2)}")
    print(f"  n_k = 0（k > K）：{all(sp.cancel(nk[kk]) == 0 for kk in range(K + 1, K + 3))}、n_K = {sp.factor(nk[K])}")
    if check:
        try:
            dn = direct_blocks(m, k)
        except FileNotFoundError as ex:
            print(f"  直接の照合は省略（被覆の 2 箱データが無い：{ex.filename}）"); return nk, G, A, R
        ok = all(sp.cancel(nk.get(kk, 0) - dn.get(kk, 0)) == 0 for kk in set(nk) | set(dn))
        print(f"  cover_sectors から直接作った N の全ブロックと一致：{ok}")
    return nk, G, A, R


if __name__ == "__main__":
    args = sys.argv[1:] or ["1", "2", "3", "1 --k=3"]
    kk = next((int(a[4:]) for a in args if a.startswith("--k=")), 2)
    for a in args:
        if a.startswith("--"):
            continue
        m = int(a) if a.isdigit() else a
        main(m, kk)
