r"""iterated torus 結び目（ケーブル）の DAHA-Jones 多項式（D_N と GL_N の多項式表現、一般の t。progress.md 続き 140）。
Cherednik–Danilenko 1408.4348 定理 2.4（cable.tex:1152-1174）：任意の被約な根系と重み b に対し
  JD_{r⃗,s⃗}(b) = { γ̂_1( … γ̂_{ℓ−1}( (γ̂_ℓ(P_b)/P_b(q^{−ρ_k}))⇓ )⇓ … ) }_ev、H⇓ = H(1)（多項式表現で 1 に作用）。
γ_ℓ（最初に作用）が最も内側の pattern、γ_1 が companion。γ_i は第 1 列が (r_i, s_i) の SL(2,Z) の元の持ち上げ。
2 段のケーブル（companion の (p₁,q₁) と pattern の (p₂,q₂)）では、Newton 対 r = (p₁,p₂)、s = (q₁, q₂ − p₁q₁p₂)
（docs/framed_knot/conventions.md「ケーブル」）。代数的 ⇔ q₂ > p₁q₁p₂。
実装：pattern の段は旧版 daha_dn.gamma_K（D_N）/ gl_gamma_P（GL_N）で γ̂(P_b)·1 を多項式として作り、
companion の段は得た多項式 f(X) = Σ c_a X^a に f(γ̂(X))·1 = Σ c_a Π_i γ̂(X_{ε_i})^{a_i}·1 を当てる
（γ̂(X_{ε_i}) どうしは可換。指数の接頭辞ごとに途中の多項式を共有する）。"""
from .daha_dn import DN, GLN, gamma_K, gl_gamma_P, gl_P_poly, K_poly, _add, GAMMA_WORDS
from .flint_rf import ONE, conv

# pattern の語（第 1 列 (r, s)、左から順に掛けた τ_+^{a} τ_−^{b}）
GAMMA_WORDS.setdefault("P2_1", [(1, 1), (-1, 1)])               # τ_+ τ_− = [[2,1],[1,1]]：第 1 列 (2, 1)
GAMMA_WORDS.setdefault("P2_3", [(-1, 1), (1, 1), (-1, 1)])      # τ_− τ_+ τ_− = [[2,1],[3,2]]：第 1 列 (2, 3)

# 代数的なケーブル（trefoil の (2, q₂)-ケーブル）：名前 → (pattern の語, companion の (m, k))。companion は T(3,2)（m = 1, k = 2）
CABLES = {"K1": ("P2_1", (1, 2)),       # trefoil の (2,13)-ケーブル、Newton 対 s₂ = 13 − 12 = 1
          "K3": ("P2_3", (1, 2))}       # trefoil の (2,15)-ケーブル、s₂ = 3


def gamma_poly(D, f, m, k=2):
    """f(γ̂(X))·1（f は Laurent 多項式 {指数: 係数}、γ = τ_+^m τ_−^k または語の名前 m）。"""
    N = D.N
    words = [D.gamma_X(i, m, k) for i in range(N)]
    inv = [D.inverse_word(w) for w in words]
    one = {tuple([0] * N): ONE}
    cache = {(): one}

    def get(pre):
        # pre = 指数の接頭辞 (a_1, …, a_j)。γ̂(X_j)^{a_j} … γ̂(X_1)^{a_1}·1（可換なので順序は自由）
        if pre in cache:
            return cache[pre]
        *head, last = pre
        head = tuple(head)
        if last == 0:
            g = get(head)
        else:
            step = 1 if last > 0 else -1
            g = D.apply_word(words[len(pre) - 1] if step > 0 else inv[len(pre) - 1], get(head + (last - step,)))
        cache[pre] = g
        return g

    out = {}
    for a, c in sorted(f.items()):
        if c.is_zero():
            continue
        # 末尾の 0 を落として接頭辞を短くする
        aa = tuple(a)
        while aa and aa[-1] == 0:
            aa = aa[:-1]
        _add(out, get(aa), c)
    return out


def iter_jones_dn(N, lam, cable):
    """D_N：(Z, K)。Z = {γ̂_comp((γ̂_pat(K_λ)·1))·1}_ev（非 reduced、P_b(q^{−ρ}) で割る前）、K = K_λ(q^{−ρ})。"""
    pat, (m, k) = CABLES[cable]
    D = DN(N)
    f1 = gamma_K(D, lam, pat, 2)
    return D.evaluate(gamma_poly(D, f1, m, k)), D.evaluate(K_poly(D, lam))


def iter_jones_gl(N, lam, cable):
    """GL_N：(Z, K)。P_λ は macdonald_P。"""
    pat, (m, k) = CABLES[cable]
    D = GLN(N)
    f1 = gl_gamma_P(D, lam, pat, 2)
    return D.evaluate(gamma_poly(D, f1, m, k)), D.evaluate(gl_P_poly(D, lam))


def iter_jones_dn_fast(N, lam, cable, F=None):
    """iter_jones_dn と同じ (Z, K)。companion の段を高速版で：f₁ = γ̂_pat(K_λ)·1 は B_N 対称なので
    f₁ = Σ_ν c_ν H_ν（|ν| ≤ deg f₁）と解き、Z = Σ c_ν {γ̂_comp(H_ν)}（daha_fast.Fast.value）。F を渡すと λ の間で共有する。"""
    from .daha_fast import Fast, _with_rho2, solve_symmetric, h_basis
    pat, (m, k) = CABLES[cable]
    if F is None:
        F = Fast(_with_rho2(DN(N)), m, k)
    D = F.D
    f1 = gamma_K(D, lam, pat, 2)
    deg = max(sum(abs(x) for x in a) for a in f1)
    x = solve_symmetric(D, f1, h_basis(F, range(deg + 1)))
    Z = sum((c * F.value(nu) for nu, c in x.items() if not c.is_zero()), conv(0))
    return Z, D.evaluate(K_poly(D, lam))


# ---- companion の段を約分なしの版（daha_fast_lp）で（続き 141） ------------------------------------------
def _expand_dn(D, f1, m, k):
    """D_N：B_N 対称な f₁ を H_ν（|ν| ≤ deg f₁）で展開した係数 {ν: RF}（0 は除く）。"""
    from .daha_fast import Fast, solve_symmetric, h_basis
    deg = max(sum(abs(x) for x in a) for a in f1)
    x = solve_symmetric(D, f1, h_basis(Fast(D, m, k), range(deg + 1)))
    return {nu: c for nu, c in x.items() if not c.is_zero()}


def dn_basis_lp(N, nus, m=1, k=2, V=None, lowmem=False):
    """V[ν] = {γ̂(H_ν)}（D_N、γ = τ_+^m τ_−^k、RF）を約分なしの版で。V に無い ν だけ計算して足す。
    lowmem：ν ごとに作り直し、B^j·(対称な多項式) を保存しない（ν の間の共有は少ないので時間はほぼ同じ。続き 141）。"""
    from .daha_fast_lp import FastLP, DNLP
    V = {} if V is None else V
    todo = [nu for nu in nus if nu not in V]
    F = None
    for nu in todo:
        if F is None or lowmem:
            F = FastLP(DNLP(N), m, k, 1, final="functional")
            if lowmem:
                F._B = _B_nocache(F)
        V[nu] = F.value(nu)
    return V


def _B_nocache(F):
    def _B(nu, j):
        g = F.get_sym(nu)
        for _ in range(j):
            g = F.D.apply_word(F.word, g)
        return g
    return _B


def iter_jones_dn_lp(N, lam, cable, V=None):
    """iter_jones_dn と同じ (Z, K)。companion の段は基底の値 V[ν]（dn_basis_lp、λ・ケーブルの間で共有できる）の一次結合。"""
    from .daha_fast import _with_rho2
    pat, (m, k) = CABLES[cable]
    D = _with_rho2(DN(N))
    x = _expand_dn(D, gamma_K(D, lam, pat, 2), m, k)
    V = dn_basis_lp(N, list(x), m, k, V)
    Z = sum((c * V[nu] for nu, c in x.items()), conv(0))
    return Z, D.evaluate(K_poly(D, lam))


def iter_adj_gl(N, cable):
    """GL_N の随伴 (1;1)：{γ̂_comp((γ̂_pat(P_adj)·1))}（旧版。照合用）。"""
    from .daha_dn import gl_adj_gamma
    pat, (m, k) = CABLES[cable]
    D = GLN(N)
    return D.evaluate(gamma_poly(D, gl_adj_gamma(D, pat, 2), m, k))


def _expand_gl_adj(D, f1, m, k):
    """GL_N：f₁ = γ̂_pat(P_adj)·1（X の総次数 0、対称）を H_μ(X) H_ν(X^{−1})（|μ| = |ν|）で展開した係数 {(μ, ν): RF}。"""
    from .daha_fast import Fast, solve_symmetric, h_basis, poly_mul
    d = max(sum(x for x in a if x > 0) for a in f1)
    Hp = h_basis(Fast(D, m, k, 1), range(d + 1))
    Hm = h_basis(Fast(D, m, k, -1), range(d + 1))
    basis = {(mu, nu): poly_mul(Hp[mu], Hm[nu]) for mu in Hp for nu in Hm if sum(mu) == sum(nu)}
    x = solve_symmetric(D, f1, basis)
    return {p: c for p, c in x.items() if not c.is_zero()}


def gl_basis_lp(N, pairs, m=1, k=2, V=None):
    """V[(μ, ν)] = {γ̂(H_μ(X) H_ν(X^{−1}))}（GL_N、ρ = (N−1, …, 0) の評価、|μ| = |ν| なので換算なし）。無いものだけ計算する。"""
    from .daha_fast_lp import FastLP, GLNLP
    V = {} if V is None else V
    todo = [p for p in pairs if p not in V]
    if todo:
        D = GLNLP(N)
        Fm = FastLP(D, m, k, -1)
        for nu in sorted({b for _, b in todo}):
            Fp = FastLP(D, m, k, 1, final="functional")
            Fp.sym = {(): Fm.get_sym(nu)}
            for mu, b in todo:
                if b == nu:
                    assert sum(mu) == sum(nu)
                    V[(mu, nu)] = Fp.value(mu)
    return V


def iter_adj_gl_lp(N, cable, V=None):
    """iter_adj_gl と同じ値。companion の段は基底の値 V[(μ, ν)]（gl_basis_lp）の一次結合。"""
    from .daha_fast import _with_rho2
    from .daha_dn import gl_adj_gamma
    pat, (m, k) = CABLES[cable]
    D = _with_rho2(GLN(N))
    x = _expand_gl_adj(D, gl_adj_gamma(D, pat, 2), m, k)
    V = gl_basis_lp(N, list(x), m, k, V)
    return sum((c * V[p] for p, c in x.items()), conv(0))
