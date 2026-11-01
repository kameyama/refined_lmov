r"""D_N の DAHA の多項式表現と DAHA-Jones 多項式（Cherednik 1111.6195、jones.tex §1, §2.1）。

規約（jones.tex :509-842）：単純ルート α_i = ε_i − ε_{i+1}（i < N）、α_N = ε_{N−1} + ε_N、(α, α) = 2。ρ = (N−1, …, 1, 0)、t = q^k。
  Demazure–Lusztig：T_i = t^{1/2} s_i + (t^{1/2} − t^{−1/2})(X_{α_i} − 1)^{−1}(s_i − 1)（(T − t^{1/2})(T + t^{−1/2}) = 0、T_i(1) = t^{1/2}）。
  π_1（ω_1 = ε_1 は minuscule）：π_1(X_b) = q^{(ω_1, b)} X_{u_1^{−1}(b)}、u_1 = w_0 w_0^{ω_1}。Y_{ε_1} = π_1 T_{u_1}（l(u_1) = 2N − 2）。
  Y_{ε_{i+1}} = T_i^{−1} Y_{ε_i} T_i^{−1}、X_{ε_{i+1}} = T_i X_{ε_i} T_i。
  τ_+：X ↦ X、T_i ↦ T_i、Y_{ε_1} ↦ q^{−1/2} X_{ε_1} Y_{ε_1}。τ_−：Y ↦ Y、T_i ↦ T_i、X_{ε_1} ↦ q^{1/2} Y_{ε_1} X_{ε_1}（(ω_1, ω_1) = 1）。
  評価 {H} = (H·1)(q^{−ρ_k})、X_b(q^{−ρ_k}) = t^{−(b, ρ)}。JD_{r,s}(b) = {γ̂_{r,s}(P_b)} / P_b(q^{−ρ_k})、γ ∈ ⟨τ_±⟩ の第 1 列が (r, s)。
多項式は {指数の組: RF}（RF は rcsknots.refined_kauffman.flint_rf の (r, s) = (q^{1/2}, t^{1/2}) の有理関数）。"""
import itertools
from fractions import Fraction

from .flint_rf import RF, R, S, ONE, conv

ZERO = conv(0)
G_ = S - 1 / S          # t^{1/2} − t^{−1/2}


def _add(acc, poly, c=None):
    for a, v in poly.items():
        w = v if c is None else c * v
        if a in acc:
            acc[a] = acc[a] + w
            if acc[a].is_zero():
                del acc[a]
        elif not w.is_zero():
            acc[a] = w
    return acc


class DN:
    def __init__(self, N):
        self.N = N
        e = lambda i: tuple(1 if j == i else 0 for j in range(N))
        self.alpha = [tuple(x - y for x, y in zip(e(i), e(i + 1))) for i in range(N - 1)]
        self.alpha.append(tuple(x + y for x, y in zip(e(N - 2), e(N - 1))) if N >= 2 else None)
        self.rho = tuple(N - 1 - i for i in range(N))
        self.pos = [tuple(v) for v in self._positive_roots()]
        self.u1 = self._u1()
        self.u1inv = self._inv(self.u1)
        self.word_u1 = self._reduced_word(self.u1)
        assert len(self.word_u1) == 2 * N - 2

    # ---- 根系 ----
    def _positive_roots(self):
        N = self.N
        out = []
        for i in range(N):
            for j in range(i + 1, N):
                for sg in (-1, 1):
                    v = [0] * N
                    v[i], v[j] = 1, sg
                    out.append(v)
        return out

    @staticmethod
    def _dot(a, b):
        return sum(x * y for x, y in zip(a, b))

    def reflect(self, i, a):
        al = self.alpha[i]
        k = self._dot(a, al)
        return tuple(x - k * y for x, y in zip(a, al))

    def _mat_apply(self, M, a):
        return tuple(sum(M[i][j] * a[j] for j in range(self.N)) for i in range(self.N))

    def _matmul(self, A, B):
        n = self.N
        return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

    def _inv(self, M):
        return [list(row) for row in zip(*M)]       # 符号つき置換行列：逆 = 転置

    def _u1(self):
        N = self.N
        def w0(ids):
            # ε_i (i ∈ ids) の上の D_{|ids|} の最長元：|ids| が偶数なら −1、奇数なら最後の座標以外 −1
            M = [[1 if i == j else 0 for j in range(N)] for i in range(N)]
            for k, i in enumerate(ids):
                M[i][i] = -1
            if len(ids) % 2 == 1:
                M[ids[-1]][ids[-1]] = 1
            return M
        return self._matmul(w0(list(range(N))), w0(list(range(1, N))))

    def _length_and_descent(self, M):
        neg = [i for i in range(len(self.alpha)) if not self._is_pos(self._mat_apply(M, self.alpha[i]))]
        return neg

    def _is_pos(self, v):
        return tuple(v) in self.pos

    def _reduced_word(self, M):
        """M = s_{j_1} … s_{j_l}（簡約）の (j_1, …, j_l)。右の降下 w(α_i) < 0 を順に外す。"""
        word = []
        cur = [row[:] for row in M]
        while True:
            neg = self._length_and_descent(cur)
            if not neg:
                break
            i = neg[0]
            si = [[1 if r == c else 0 for c in range(self.N)] for r in range(self.N)]
            al = self.alpha[i]
            for r in range(self.N):
                for c in range(self.N):
                    si[r][c] -= al[r] * al[c]
            cur = self._matmul(cur, si)
            word.append(i)
        return list(reversed(word))

    # ---- 作用素 ----
    def s_op(self, i, f):
        out = {}
        for a, v in f.items():
            _add(out, {self.reflect(i, a): v})
        return out

    def divdiff(self, i, f):
        """(s_i f − f)/(X_{α_i} − 1)。"""
        al = self.alpha[i]
        out = {}
        for a, v in f.items():
            k = self._dot(a, al)
            if k > 0:
                for j in range(k):
                    _add(out, {tuple(x - (k - j) * y for x, y in zip(a, al)): -v})
            elif k < 0:
                for j in range(-k):
                    _add(out, {tuple(x + j * y for x, y in zip(a, al)): v})
        return out

    def T(self, i, f, inv=False):
        out = _add({}, self.s_op(i, f), S)
        _add(out, self.divdiff(i, f), G_)
        if inv:                                 # T^{−1} = T − (t^{1/2} − t^{−1/2})
            _add(out, f, -G_)
        return out

    def pi1(self, f, inv=False):
        out = {}
        if not inv:
            for a, v in f.items():
                _add(out, {self._mat_apply(self.u1inv, a): v * R**(2 * a[0])})
        else:
            # π_1^{-1}：X_c ↦ q^{−(ε_1, u_1 c)} X_{u_1 c}
            for a, v in f.items():
                b = self._mat_apply(self.u1, a)
                _add(out, {b: v * R**(-2 * b[0])})
        return out

    def Y1(self, f, inv=False):
        """Y_{ε_1} = π_1 T_{j_1} … T_{j_l}（右から作用）。逆は T^{−1} を逆順に、最後に π_1^{−1} を先に。"""
        if not inv:
            g = f
            for j in reversed(self.word_u1):
                g = self.T(j, g)
            return self.pi1(g)
        g = self.pi1(f, inv=True)
        for j in self.word_u1:
            g = self.T(j, g, inv=True)
        return g

    def Xmul(self, b, f, power=1):
        b = tuple(power * x for x in b)
        out = {}
        for a, v in f.items():
            _add(out, {tuple(x + y for x, y in zip(a, b)): v})
        return out

    def evaluate(self, f):
        """f(q^{−ρ_k})：X_b ↦ t^{−(b, ρ)}。"""
        tot = ZERO
        for a, v in f.items():
            tot = tot + v * S**(-2 * self._dot(a, self.rho))
        return tot

    # ---- 語（右から作用する演算の列） ----
    def apply_word(self, word, f):
        g = f
        for op in reversed(word):
            kind = op[0]
            if kind == "c":
                g = {a: op[1] * v for a, v in g.items()}
            elif kind == "X":
                g = self.Xmul(op[1], g, op[2])
            elif kind == "T":
                g = self.T(op[1], g, inv=op[2] < 0)
            elif kind == "Y":
                g = self.Y1(g, inv=op[1] < 0)
            else:
                raise ValueError(op)
        return g

    @staticmethod
    def inverse_word(word):
        out = []
        for op in reversed(word):
            kind = op[0]
            if kind == "c":
                out.append(("c", 1 / op[1]))
            elif kind == "X":
                out.append(("X", op[1], -op[2]))
            elif kind == "T":
                out.append(("T", op[1], -op[2]))
            elif kind == "Y":
                out.append(("Y", -op[1]))
        return out

    def e(self, i):
        return tuple(1 if j == i else 0 for j in range(self.N))

    def gamma_X(self, i, m, k):
        """γ̂(X_{ε_{i+1}})：γ = τ_+^m τ_−^k（m が文字列なら GAMMA_WORDS の語、続き 103）。i は 0 始まり。"""
        e1 = self.e(0)
        w = gamma_subst(gamma_factors(m, k), [("X", e1, 1)], e1)
        for j in range(i):                     # X_{ε_{j+2}} = T_{j+1} X_{ε_{j+1}} T_{j+1}
            w = [("T", j, 1)] + w + [("T", j, 1)]
        return w

    def gamma_X_old(self, i, m, k):
        """γ̂(X_{ε_{i+1}})、γ = τ_+^m τ_−^k：τ_−^k(X_{ε_1}) = q^{k/2} Y^k X、τ_+^m(Y) = q^{−m/2} X^m Y。i は 0 始まり。"""
        e1 = self.e(0)
        w = [("c", R**k)]
        for _ in range(k):
            w += [("c", R**(-m)), ("X", e1, m), ("Y", 1)]
        w.append(("X", e1, 1))
        for j in range(i):                     # X_{ε_{j+2}} = T_{j+1} X_{ε_{j+1}} T_{j+1}
            w = [("T", j, 1)] + w + [("T", j, 1)]
        return w

    def gamma_S(self, d, m, k, f):
        """γ̂(S_d) f、S_d = Σ_i (X_{ε_i}^d + X_{−ε_i}^d)。"""
        out = {}
        for i in range(self.N):
            w = self.gamma_X(i, m, k)
            wi = self.inverse_word(w)
            g1, g2 = f, f
            for _ in range(d):
                g1 = self.apply_word(w, g1)
                g2 = self.apply_word(wi, g2)
            _add(out, g1)
            _add(out, g2)
        return out


# ---- D_N の Macdonald 多項式（data/refined_kauffman/koornwinder_DN の K_λ、U = t^N） ----
def K_poly(D, lam, store="data/refined_kauffman/koornwinder_DN"):
    """K_λ を X の Laurent 多項式に（S_d = Σ_i (X_{ε_i}^d + X_{−ε_i}^d)、係数は U = t^N を代入）。"""
    import os, pickle
    import sympy as sp
    from .flint_rf import from_sympy
    if not lam:
        return {tuple([0] * D.N): ONE}
    q_, t_, U_ = sp.symbols("q t U")
    r_, s_ = sp.symbols("r s", positive=True)
    raw = pickle.load(open(os.path.join(store, "K_" + "_".join(map(str, lam)) + ".pkl"), "rb"))
    out = {}
    Sd = {}
    for mu, v in raw.items():
        c = from_sympy(sp.sympify(v, locals={"q": q_, "t": t_, "U": U_}).subs({q_: r_**2, t_: s_**2, U_: s_**(2 * D.N)}))
        term = {tuple([0] * D.N): ONE}
        for d in mu:
            if d not in Sd:
                Sd[d] = {}
                for i in range(D.N):
                    _add(Sd[d], {tuple(d if j == i else 0 for j in range(D.N)): ONE})
                    _add(Sd[d], {tuple(-d if j == i else 0 for j in range(D.N)): ONE})
            new = {}
            for a, x in term.items():
                for b, y in Sd[d].items():
                    _add(new, {tuple(p + w for p, w in zip(a, b)): x * y})
            term = new
        _add(out, term, c)
    return out


def macdonald_op(D, f):
    """Σ_i (Y_{ε_i} + Y_{ε_i}^{−1}) f。"""
    out = {}
    for i in range(D.N):
        g, gi = f, f
        for j in reversed(range(i)):            # Y_{ε_{i+1}} = T_i^{-1} … T_1^{-1} Y_1 T_1^{-1} … T_i^{-1}
            pass
        def Yi(h, inv):
            # Y_{ε_{i+1}} = T_i^{−1} Y_{ε_i} T_i^{−1}：Y_{ε_{i+1}} = (T_i^{−1}⋯T_1^{−1}) Y_1 (T_1^{−1}⋯T_i^{−1})、逆は T を正の冪で
            for j in reversed(range(i)):
                h = D.T(j, h, inv=not inv)
            h = D.Y1(h, inv=inv)
            for j in range(i):
                h = D.T(j, h, inv=not inv)
            return h
        _add(out, Yi(f, False))
        _add(out, Yi(f, True))
    return out


# ---- γ ∈ PSL(2, Z) の τ_± の語（続き 103） ----
# γ = τ_+^m τ_−^k は T(mk+1, k)。それ以外の結び目は名前で指定する：GAMMA_WORDS[name] = [(+1, a_1), (−1, b_1), …]（左から順に掛けた τ_+^{a_1} τ_−^{b_1} …）。
# T(r, s) は γ の第 1 列 (r, s)（τ_+ = [[1,1],[0,1]]、τ_− = [[1,0],[1,1]]）。
GAMMA_WORDS = {
    "T3_5": [(1, 1), (-1, 1), (1, 1), (-1, 1)],      # τ_+ τ_− τ_+ τ_− = [[5,3],[3,2]]：T(5,3) = T(3,5)
    "T3_10": [(1, 3), (-1, 3)],      # 第 1 列 (10, 3)（続き 104、語は冪 ≤ 4 の最短のもの）
    "T4_7": [(1, 1), (-1, 1), (1, 2), (-1, 1)],      # 第 1 列 (7, 4)（続き 104、語は冪 ≤ 4 の最短のもの）
    "T3_8": [(1, 2), (-1, 1), (1, 1), (-1, 1)],      # 第 1 列 (8, 3)（続き 104、語は冪 ≤ 4 の最短のもの）
}


def gamma_factors(m, k=2):
    """γ の因子 [(±1, 冪), …]。m が文字列なら GAMMA_WORDS。"""
    if isinstance(m, str):
        return list(GAMMA_WORDS[m])
    return [(1, m), (-1, k)]


def gamma_first_column(factors):
    """γ の第 1 列 (r, s)。"""
    a, b, c, d = 1, 0, 0, 1
    for sg, e in factors:
        if sg > 0:
            a, b, c, d = a, a * e + b, c, c * e + d          # × τ_+^e = [[1,e],[0,1]]
        else:
            a, b, c, d = a + b * e, b, c + d * e, d          # × τ_−^e = [[1,0],[e,1]]
    return a, c


def _inv_word(word):
    out = []
    for op in reversed(word):
        kind = op[0]
        if kind == "c":
            out.append(("c", 1 / op[1]))
        elif kind == "X":
            out.append(("X", op[1], -op[2]))
        elif kind == "T":
            out.append(("T", op[1], -op[2]))
        elif kind == "Y":
            out.append(("Y", -op[1]))
    return out


def gamma_subst(factors, word, e1):
    """γ̂(word)：τ_+^a は Y ↦ q^{−a/2} X^a Y、τ_−^b は X ↦ q^{b/2} Y^b X（X, Y は ε_1 のもの）。
    γ = g_1 g_2 … g_n は ĝ_1(ĝ_2(…ĝ_n(word)))。"""
    for sg, e in reversed(factors):
        if e == 0:
            continue
        if sg > 0:
            blk = [("c", R**(-e)), ("X", e1, e), ("Y", 1)]
            out = []
            for op in word:
                if op[0] == "Y":
                    out += (blk if op[1] > 0 else _inv_word(blk)) * abs(op[1])
                else:
                    out.append(op)
        else:
            blk = [("c", R**e)] + [("Y", 1)] * e + [("X", e1, 1)]
            out = []
            for op in word:
                if op[0] == "X":
                    assert op[1] == e1, op
                    out += (blk if op[2] > 0 else _inv_word(blk)) * abs(op[2])
                else:
                    out.append(op)
        word = out
    # 定数を先頭にまとめる
    c = ONE
    rest = []
    for op in word:
        if op[0] == "c":
            c = c * op[1]
        else:
            rest.append(op)
    return [("c", c)] + rest


def gamma_K(D, lam, m, k, store="data/refined_kauffman/koornwinder_DN"):
    """γ̂(K_λ)·1（γ = τ_+^m τ_−^k）。K_λ の S_μ 表示から γ̂(S_μ) = Π γ̂(S_{μ_i})（可換）として作用させる。"""
    import os, pickle
    import sympy as sp
    from .flint_rf import from_sympy
    q_, t_, U_ = sp.symbols("q t U")
    r_, s_ = sp.symbols("r s", positive=True)
    one = {tuple([0] * D.N): ONE}
    if not lam:
        return one
    raw = pickle.load(open(os.path.join(store, "K_" + "_".join(map(str, lam)) + ".pkl"), "rb"))
    out = {}
    for mu, v in raw.items():
        c = from_sympy(sp.sympify(v, locals={"q": q_, "t": t_, "U": U_}).subs({q_: r_**2, t_: s_**2, U_: s_**(2 * D.N)}))
        g = one
        for d in mu:
            g = D.gamma_S(d, m, k, g)
        _add(out, g, c)
    return out


def daha_jones(N, lam, m, k=2):
    """(Z, K(q^{−ρ}))：Z = {γ̂(K_λ)}（非 reduced）、γ = τ_+^m τ_−^k。T(2m+1, 2) は (m, 2)。"""
    D = DN(N)
    return D.evaluate(gamma_K(D, lam, m, k)), D.evaluate(K_poly(D, lam))



class GLN(DN):
    """GL_N（A_{N−1}、U(N)）の DAHA の多項式表現。α_i = ε_i − ε_{i+1}（i < N）、ρ = (N−1, …, 0)、
    ω_1 = ε_1（minuscule）、u_1 = w_0 w_0^{ω_1}（巡回置換、l(u_1) = N − 1）、Y_{ε_1} = π_1 T_{u_1}、(ε_1, ε_1) = 1。"""

    def __init__(self, N):
        self.N = N
        e = lambda i: tuple(1 if j == i else 0 for j in range(N))
        self.alpha = [tuple(x - y for x, y in zip(e(i), e(i + 1))) for i in range(N - 1)]
        self.rho = tuple(N - 1 - i for i in range(N))
        self.pos = [tuple(1 if k == i else (-1 if k == j else 0) for k in range(N)) for i in range(N) for j in range(i + 1, N)]
        rev = lambda ids: [[1 if (r in ids and c == ids[len(ids) - 1 - ids.index(r)]) or (r not in ids and r == c) else 0
                            for c in range(N)] for r in range(N)]
        self.u1 = self._matmul(rev(list(range(N))), rev(list(range(1, N))))
        self.u1inv = self._inv(self.u1)
        self.word_u1 = self._reduced_word(self.u1)
        assert len(self.word_u1) == N - 1

    def gamma_P(self, d, m, k, f, sign=1):
        """γ̂(p_d(X^{sign})) f、p_d(X) = Σ_i X_{ε_i}^d。"""
        out = {}
        for i in range(self.N):
            w = self.gamma_X(i, m, k)
            if sign < 0:
                w = self.inverse_word(w)
            g = f
            for _ in range(d):
                g = self.apply_word(w, g)
            _add(out, g)
        return out


def _gl_pi1(self, f, inv=False):
    """GL_N の π_1：X_b ↦ q^{(ω_{1*}, b)} X_{u_1^{−1} b}、ω_{1*} = ω_{N−1} ≡ −ε_N（π(X_i) = X_{i+1}、π(X_N) = q^{−1} X_1）。"""
    out = {}
    if not inv:
        for a, v in f.items():
            _add(out, {self._mat_apply(self.u1inv, a): v * R**(-2 * a[-1])})
    else:
        for a, v in f.items():
            b = self._mat_apply(self.u1, a)
            _add(out, {b: v * R**(2 * a[-1] if False else 2 * b[-1])})
    return out


GLN.pi1 = _gl_pi1


def gl_gamma_P(D, lam, m, k):
    """GL_N：γ̂(P_λ(X))·1（P_λ は macdonald_P の冪和表示、p_d(X) = Σ_i X_{ε_i}^d）。"""
    from ..macdonald import macdonald_P
    one = {tuple([0] * D.N): ONE}
    if not lam:
        return one
    P = macdonald_P(sum(lam), R**2, S**2, conv)[tuple(lam)]
    out = {}
    for rho, c in P.items():
        g = one
        for d in rho:
            g = D.gamma_P(d, m, k, g)
        _add(out, g, c)
    return out


def gl_P_poly(D, lam):
    from ..macdonald import macdonald_P
    one = {tuple([0] * D.N): ONE}
    if not lam:
        return one
    P = macdonald_P(sum(lam), R**2, S**2, conv)[tuple(lam)]
    pd = lambda d: {tuple(d if j == i else 0 for j in range(D.N)): ONE for i in range(D.N)}
    out = {}
    for rho, c in P.items():
        g = one
        for d in rho:
            new = {}
            for a, x in g.items():
                for b, y in pd(d).items():
                    _add(new, {tuple(p + w for p, w in zip(a, b)): x * y})
            g = new
        _add(out, g, c)
    return out


def gl_adj_gamma(D, m, k):
    """GL_N：γ̂(P_adj)·1、P_adj = e_1(x) e_1(1/x) − κ、κ = (1 − t^N)(1 − q)/((1 − t)(1 − q t^{N−1}))（scripts/rk_k5_trefoil_3box.py で確認）。"""
    one = {tuple([0] * D.N): ONE}
    g = D.gamma_P(1, m, k, D.gamma_P(1, m, k, one, sign=-1), sign=1)
    q, t = R**2, S**2
    kappa = (1 - t**D.N) * (1 - q) / ((1 - t) * (1 - q * t**(D.N - 1)))
    return _add(g, one, -kappa)


def gl_adj_poly(D):
    q, t = R**2, S**2
    kappa = (1 - t**D.N) * (1 - q) / ((1 - t) * (1 - q * t**(D.N - 1)))
    out = {}
    for i in range(D.N):
        for j in range(D.N):
            a = tuple((1 if l == i else 0) - (1 if l == j else 0) for l in range(D.N))
            _add(out, {a: ONE})
    return _add(out, {tuple([0] * D.N): ONE}, -kappa)
