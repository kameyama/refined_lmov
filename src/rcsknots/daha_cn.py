r"""C_N の DAHA の多項式表現と DAHA-Jones 多項式（Cherednik 1111.6195 §Kauffman、twisted 設定）。

規約（jones.tex §Kauffman、§1）：短い根 ±ε_i ± ε_j（(α,α) = 2、ν = 1）、長い根 ±2ε_i（(α,α) = 4、ν = 2）。
  単純根 α_i = ε_i − ε_{i+1}（i < N）、α_N = 2ε_N。反射 s_α(b) = b − (b, α^∨) α、α^∨ = 2α/(α,α)。
  パラメータ：t = t_sht = q^{k_sht}、u = t_lng = q^{ν_lng k_lng}。u = 1 で D_N の DAHA-Jones に一致（Cherednik (kauf)）。
  Avetisyan–Mkrtchyan 2205.12832 の refined Sp(2N)（k_α = y が一様、q_α = q^{(α,α)/2}）は u = t² に当たる。
  Demazure–Lusztig：T_i = t_i^{1/2} s_i + (t_i^{1/2} − t_i^{−1/2})(X_{α_i} − 1)^{−1}(s_i − 1)、t_i = t（i < N）、t_N = u。
  π_1（ω_1 = ε_1 は minuscule）：π_1(X_b) = q^{(ω_1, b)} X_{u_1^{−1} b}、u_1 = w_0 w_0^{ω_1} = diag(−1, 1, …, 1)、l(u_1) = 2N − 1。
  Y_{ε_1} = π_1 T_{u_1}。τ_−(X_{ε_1}) = q^{1/2} Y X、τ_+(Y_{ε_1}) = q^{−1/2} X Y。
  評価 {H} = (H·1)(q^{−ρ_k})、X_b ↦ t^{−(b, ρ_sht)} u^{−(b, ρ_lng)/2}、ρ_sht = (N−1, …, 0)、ρ_lng = (1, …, 1)。
  JD_{r,s}(ω_1) = {γ̂(X_{ε_1})}/X_{ε_1}(q^{−ρ_k})（E_{ω_1} = X_{ε_1}、Cherednik §Kauffman の注意）。
多項式は {指数の組: RF}（RF は rcsknots.flint_rf3 の 3 変数 (R, S, U) = (q^{1/2}, t^{1/2}, u^{1/2}) の有理関数）。"""
from .flint_rf3 import RF, R, S, U, conv

ONE = conv(1)
ZERO = conv(0)


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


class CN:
    def __init__(self, N):
        self.N = N
        e = lambda i: tuple(1 if j == i else 0 for j in range(N))
        self.alpha = [tuple(x - y for x, y in zip(e(i), e(i + 1))) for i in range(N - 1)] + [tuple(2 * x for x in e(N - 1))]
        self.alpha_sq = [2] * (N - 1) + [4]
        self.tpar = [S] * (N - 1) + [U]                 # t_i^{1/2}
        self.rho_sht = tuple(N - 1 - i for i in range(N))
        self.rho_lng = tuple(1 for _ in range(N))
        self.pos = set()
        for i in range(N):
            self.pos.add(tuple(2 if k == i else 0 for k in range(N)))
            for j in range(i + 1, N):
                for sg in (-1, 1):
                    self.pos.add(tuple(1 if k == i else (sg if k == j else 0) for k in range(N)))
        self.u1 = [[(-1 if (r == 0 and c == 0) else (1 if r == c else 0)) for c in range(N)] for r in range(N)]
        self.u1inv = self.u1
        self.word_u1 = self._reduced_word(self.u1)
        assert len(self.word_u1) == 2 * N - 1, self.word_u1

    @staticmethod
    def _dot(a, b):
        return sum(x * y for x, y in zip(a, b))

    def coroot_pairing(self, a, i):
        return 2 * self._dot(a, self.alpha[i]) // self.alpha_sq[i]

    def reflect(self, i, a):
        k = self.coroot_pairing(a, i)
        return tuple(x - k * y for x, y in zip(a, self.alpha[i]))

    def _mat_apply(self, M, a):
        return tuple(sum(M[i][j] * a[j] for j in range(self.N)) for i in range(self.N))

    def _reduced_word(self, M):
        word = []
        cur = [row[:] for row in M]
        while True:
            neg = [i for i in range(self.N) if tuple(self._mat_apply(cur, self.alpha[i])) not in self.pos]
            if not neg:
                break
            i = neg[0]
            al = self.alpha[i]
            si = [[(1 if r == c else 0) - 2 * al[r] * al[c] // self.alpha_sq[i] for c in range(self.N)] for r in range(self.N)]
            cur = [[sum(cur[r][k] * si[k][c] for k in range(self.N)) for c in range(self.N)] for r in range(self.N)]
            word.append(i)
        return list(reversed(word))

    def s_op(self, i, f):
        out = {}
        for a, v in f.items():
            _add(out, {self.reflect(i, a): v})
        return out

    def divdiff(self, i, f):
        """(s_i f − f)/(X_{α_i} − 1)。(a, α_i^∨) = k なら s_i X^a = X^{a − kα}、(X^{a−kα} − X^a)/(X^α − 1) = −Σ_{j=0}^{k−1} X^{a−(k−j)α}（k > 0）。"""
        al = self.alpha[i]
        out = {}
        for a, v in f.items():
            k = self.coroot_pairing(a, i)
            if k > 0:
                for j in range(k):
                    _add(out, {tuple(x - (k - j) * y for x, y in zip(a, al)): -v})
            elif k < 0:
                for j in range(-k):
                    _add(out, {tuple(x + j * y for x, y in zip(a, al)): v})
        return out

    def T(self, i, f, inv=False):
        th = self.tpar[i]
        out = _add({}, self.s_op(i, f), th)
        _add(out, self.divdiff(i, f), th - 1 / th)
        if inv:
            _add(out, f, -(th - 1 / th))
        return out

    def pi1(self, f, inv=False):
        out = {}
        if not inv:
            for a, v in f.items():
                _add(out, {self._mat_apply(self.u1inv, a): v * R**(2 * a[0])})
        else:
            for a, v in f.items():
                b = self._mat_apply(self.u1, a)
                _add(out, {b: v * R**(-2 * b[0])})
        return out

    def Y1(self, f, inv=False):
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
        tot = ZERO
        for a, v in f.items():
            tot = tot + v * S**(-2 * self._dot(a, self.rho_sht)) * U**(-self._dot(a, self.rho_lng))
        return tot

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

    def e(self, i):
        return tuple(1 if j == i else 0 for j in range(self.N))

    def gamma_X1(self, m, k):
        """γ̂(X_{ε_1})、γ = τ_+^m τ_−^k（daha_dn.gamma_X の i = 0）。"""
        e1 = self.e(0)
        w = [("c", R**k)]
        for _ in range(k):
            w += [("c", R**(-m)), ("X", e1, m), ("Y", 1)]
        w.append(("X", e1, 1))
        return w


def daha_jones_omega1(N, m, k=2):
    """(Z, X_{ε_1}(q^{−ρ_k}))：Z = {γ̂(X_{ε_1})}、γ = τ_+^m τ_−^k（T(2m+1, 2) は (m, 2)）。JD = Z / X_{ε_1}(q^{−ρ_k})。"""
    D = CN(N)
    one = {tuple([0] * N): ONE}
    Z = D.evaluate(D.apply_word(D.gamma_X1(m, k), one))
    return Z, D.evaluate({D.e(0): ONE})
