r"""DAHA-Jones 多項式の高速版（D_N と GL_N の多項式表現）。src/rcsknots/daha_dn.py と同じ値を返す。

daha_dn.py との違いは計算の順序だけで、作用素（T_i、π_1、Y_{ε_1}、X）は daha_dn.DN / GLN のものをそのまま使う。

【導出】使う恒等式（記号は daha_dn.py に従う。γ̂ = τ_+^m τ_−^k、H は X の対称 Laurent 多項式）。
1. φ（X_b ↦ Y_b^{−1}、Y_b ↦ X_b^{−1}、T_i ↦ T_i の反自己同型、{φ(H)} = {H}）と φ τ_± φ = τ_∓ から
   {γ̂(H(X))} = {γ̂'(H(Y^{−1}))}、γ̂' = τ_−^m τ_+^k。
   B := γ̂'(Y_{ε_1}^{−1}) = q^{k/2} Y^{−1}(q^{−m/2} X^{−1} Y^{−m})^k は X の次数を k だけ上げる（γ̂(X_{ε_1}) は 2m + 1）。
2. 軌道 W·β（β = ±ε_1、minuscule）の剰余類代表 w に沿って Π_σ(F) := Σ_w t^{σ l(w)/2} T_w^{(σ)} F
   （σ = +1 は T、σ = −1 は T^{−1}、(b, α_j) = σ の向きに進む）と置く。(b, α_j) = ±1 のとき T_j^{±1} X_b = t^{∓1/2} X_{s_j b}
   なので Π_σ(X_β) = Σ_{b∈Wβ} X_b。H_d := Π_σ(X_{dβ}) は対称（H_aff の球 Hecke 代数）で、多項式表現が忠実なので
   H_aff の中で H_d e = Π_σ X_{dβ} e（e は対称化の冪等元）。φ と γ̂' を当てると、対称な g に対して
   γ̂'(H_d(Y^{−1})) g = Π_σ(B_β^d g)、B_β = γ̂'(Y_β^{−1})（B_β g は W_1 = Stab(β) で対称）。
3. 評価 ev（X_b ↦ t^{−(b,ρ)}、ρ は正ルートの和の半分）は ev ∘ T_j^{±1} = t^{±1/2} ev を満たすので
   ev(Π_σ(F)) = W^1(t^σ) ev(F)（W^1(t) = Σ_w t^{l(w)}）、ev(Y^{−1} F) = t^{−l(u_1)/2} ev(π_1^{−1} F)。
したがって {γ̂(H_ν(X))}（H_ν = Π H_{ν_i}）は、対称な g を Π_σ(B^{ν_i} g) で積み上げ、最後の一つだけ ev(B^d g) で評価すればよい。
K_λ（D_N）・P_λ（GL_N）は H_ν の一次結合に直してから使う（係数は線形方程式で決める）。
GL_N の daha_dn.py は ρ = (N−1, …, 0) で評価しているので、γ̂(H)·1 が X の総次数 c で斉次なら値は t^{−c(N−1)/2} だけ違う（変換して返す）。"""
from .daha_dn import DN, GLN, K_poly, gl_P_poly, _add, gamma_factors, gamma_subst, gamma_first_column
from .flint_rf import R, S, ONE, conv


def orbit_chain(D, beta, sign):
    """W·β の剰余類代表を幅優先で：[(b, 親, j, 長さ)]。(b, α_j) = sign のとき s_j b を T_j^{sign} で作る。"""
    out = [(beta, None, None, 0)]
    seen = {beta}
    k = 0
    while k < len(out):
        b, _, _, l = out[k]
        k += 1
        for j, al in enumerate(D.alpha):
            if al is not None and D._dot(b, al) == sign:
                c = D.reflect(j, b)
                if c not in seen:
                    seen.add(c)
                    out.append((c, b, j, l + 1))
    return out


def symmetrize(D, F, chain, sign):
    """Π_σ(F) = Σ_w t^{σ l(w)/2} T_w^{(σ)} F。"""
    vals = {chain[0][0]: F}
    out = dict(F)
    for b, p, j, l in chain[1:]:
        vals[b] = D.T(j, vals[p], inv=sign < 0)
        _add(out, vals[b], S**(sign * l))
    return out


def evaluate(D, f):
    """ev(f)：X_b ↦ t^{−(b, ρ)}、ρ = 正ルートの和の半分（2ρ は整数）。"""
    tot = conv(0)
    for a, v in f.items():
        tot = tot + v * S**(-D._dot(a, D.rho2))
    return tot


def poly_mul(f, g):
    out = {}
    for a, x in f.items():
        for b, y in g.items():
            _add(out, {tuple(p + q for p, q in zip(a, b)): x * y})
    return out


def dominant(D, a):
    if isinstance(D, GLN):
        return a == tuple(sorted(a, reverse=True))
    return all(a[i] >= a[i + 1] for i in range(D.N - 2)) and a[D.N - 2] >= abs(a[D.N - 1])


def solve_symmetric(D, target, basis):
    """target = Σ_ν x_ν basis[ν]（対称 Laurent 多項式、優勢な単項式の係数で比べる）を解く。解が一つでなければ余分は 0。"""
    keys = sorted({a for f in [target, *basis.values()] for a in f if dominant(D, a)})
    names = list(basis)
    rows = [[basis[n].get(a, conv(0)) for n in names] + [target.get(a, conv(0))] for a in keys]
    piv, r = [], 0
    for c in range(len(names)):
        p = next((i for i in range(r, len(rows)) if not rows[i][c].is_zero()), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        inv = 1 / rows[r][c]
        rows[r] = [x * inv for x in rows[r]]
        for i in range(len(rows)):
            if i != r and not rows[i][c].is_zero():
                f = rows[i][c]
                rows[i] = [x - f * y for x, y in zip(rows[i], rows[r])]
        piv.append(c)
        r += 1
    assert all(row[-1].is_zero() for row in rows[r:]), "解が無い"
    x = {n: conv(0) for n in names}
    for i, c in enumerate(piv):
        x[names[c]] = rows[i][-1]
    return x


class Fast:
    """一つの (D, β, σ, m, k) について {γ̂(H_ν(X))} を共有しながら計算する。"""

    def __init__(self, D, m, k=2, sign=1):
        self.D, self.m, self.k, self.sign = D, m, k, sign
        e1 = D.e(0)
        self.beta = e1 if sign > 0 else tuple(-x for x in e1)
        self.chain = orbit_chain(D, self.beta, sign)
        self.W1 = sum((S**(2 * sign * l) for *_, l in self.chain), conv(0))
        # B = γ̂'(Y_{ε_1}^{−1})、γ' = φγφ（τ_+ ↔ τ_−）。γ = τ_+^m τ_−^k では q^{k/2} Y^{−1} (q^{−m/2} X^{−1} Y^{−m})^k。B_{−ε_1} = B^{−1}
        w = gamma_subst([(-sg, e) for sg, e in gamma_factors(m, k)], [("Y", -1)], e1)
        assert w[0][0] == "c" and w[1] == ("Y", -1), w[:2]          # value() は先頭の Y^{−1} を評価で消す
        self.word = w if sign > 0 else D.inverse_word(w)
        self.one = {tuple([0] * D.N): ONE}
        self.sym = {(): self.one}       # ν → γ̂'(H_ν(Y^{−1}))·1（対称）
        self.Bpow = {}                  # (ν, j) → B^j · sym[ν]
        self.values = {}                # ν → {γ̂(H_ν(X))}

    def H(self, d):
        return symmetrize(self.D, {tuple(d * x for x in self.beta): ONE}, self.chain, self.sign)

    def _B(self, nu, j):
        if j == 0:
            return self.get_sym(nu)
        key = (nu, j)
        if key not in self.Bpow:
            self.Bpow[key] = self.D.apply_word(self.word, self._B(nu, j - 1))
        return self.Bpow[key]

    def get_sym(self, nu):
        if nu not in self.sym:
            d, rest = nu[0], nu[1:]
            self.sym[nu] = symmetrize(self.D, self._B(rest, d), self.chain, self.sign)
        return self.sym[nu]

    def value(self, nu):
        """{γ̂(H_ν(X))} = W^1(t^σ) ev(B^{ν_1} γ̂'(H_{ν_2…})·1)（ν は降順）。σ = +1 なら最後の Y^{−1} は評価で消す。"""
        if not nu:
            return ONE
        if nu not in self.values:
            d, rest = nu[0], nu[1:]
            g = self._B(rest, d - 1)
            D = self.D
            if self.sign > 0:
                h = D.apply_word(self.word[2:], g)                       # word = [c, Y^{−1}, …]
                ev = self.word[0][1] * S**(-len(D.word_u1)) * evaluate(D, D.pi1(h, inv=True))
            else:
                ev = evaluate(D, D.apply_word(self.word, g))
            self.values[nu] = self.W1 * ev
        return self.values[nu]


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for p in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - p, p):
            yield (p,) + rest


def h_basis(F, sizes):
    """{ν: H_ν}（|ν| ∈ sizes）。"""
    D = F.D
    Hd = {}
    out = {}
    for n in sizes:
        for nu in partitions(n):
            f = F.one
            for d in nu:
                if d not in Hd:
                    Hd[d] = F.H(d)
                f = poly_mul(f, Hd[d])
            out[nu] = f
    return out


def _with_rho2(D):
    D.rho2 = tuple(2 * x for x in D.rho) if not isinstance(D, GLN) else tuple(D.N - 1 - 2 * i for i in range(D.N))
    return D


def daha_jones_fast(N, lam, m, k=2, F=None):
    """daha_dn.daha_jones と同じ (Z, K)。F（Fast(DN(N), m, k)）を渡すと λ の間で途中結果を共有する。"""
    if F is None:
        F = Fast(_with_rho2(DN(N)), m, k)
    D = F.D
    Kp = K_poly(D, lam)
    x = solve_symmetric(D, Kp, h_basis(F, range(sum(lam) + 1)))
    Z = sum((c * F.value(nu) for nu, c in x.items() if not c.is_zero()), conv(0))
    return Z, D.evaluate(Kp)


def gl_jones_fast(N, lam, m, k=2, F=None):
    """daha_dn の D.evaluate(gl_gamma_P(D, λ, m, k)) と同じ値（GL_N、ρ = (N−1, …, 0) の評価に直す）。"""
    if F is None:
        F = Fast(_with_rho2(GLN(N)), m, k)
    D = F.D
    x = solve_symmetric(D, gl_P_poly(D, lam), h_basis(F, [sum(lam)]))
    Z = sum((c * F.value(nu) for nu, c in x.items() if not c.is_zero()), conv(0))
    return Z * S**(-sum(lam) * gamma_first_column(gamma_factors(m, k))[0] * (N - 1))


def gl_adj_fast(N, m, k=2):
    """daha_dn の D.evaluate(gl_adj_gamma(D, m, k)) と同じ値。{γ̂(p_1(X) p_1(X^{−1}))} − κ。"""
    D = _with_rho2(GLN(N))
    Fp, Fm = Fast(D, m, k, 1), Fast(D, m, k, -1)
    g = Fm.get_sym((1,))                                   # γ̂'(p_1(Y))·1 = Π_−(B^{−1}·1)
    Fp.sym[(1,)] = g                                       # 以後 ν = (1, 1) は「p_1(Y^{−1}) を g に」
    v = Fp.value((1, 1))
    q, t = R**2, S**2
    kappa = (1 - t**N) * (1 - q) / ((1 - t) * (1 - q * t**(N - 1)))
    return v - kappa
