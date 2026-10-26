r"""daha_fast（GL_N）の係数を Laurent 多項式（flint_lp.LP、約分なし）にした版。作用素と計算の順序は daha_fast・daha_dn と同じ。
有理関数の係数（P_λ を H_ν = Π p_{ν_i} で展開した係数）は途中に入れず、基底の値
  V[μ; ν] = {γ̂(H_μ(X) H_ν(X^{−1}))}（ρ = (N−1, …, 0) の評価、総次数 c = (|μ| − |ν|)(mk + 1) の換算込み）
を Laurent 多項式で計算し、最後に RF（flint_rf）で一次結合する。
加えて、最後の一段（評価の直前の X^{−1}・Y^{−1} の作用）を、多項式を作らずに評価の汎関数として計算できる（final="functional"）：
  点 p での評価 ev_p は環準同型で、ev_p ∘ X_b = p^b ev_p、ev_p ∘ π_1^{−1} = ev_{p′}（π_1^{−1} は単項式の写像）、
  ev_p ∘ T_j^{±1} は ev_p と ev_{s_j p} の一次結合（Demazure–Lusztig の式）。係数は p^{α_j} を含む有理関数になるので RF で持つ。"""
from collections import defaultdict

from .daha_dn import GLN, _add
from .daha_fast import orbit_chain, partitions
from .flint_lp import LP, R, S, ONE, ZERO, _mono
from . import flint_rf as RFm

G_ = S - S ** -1          # t^{1/2} − t^{−1/2}


class GLNLP(GLN):
    """GL_N の多項式表現（係数 LP）。T と π_1 だけを LP の定数で書き直す。"""

    def __init__(self, N):
        super().__init__(N)
        self.rho2 = tuple(N - 1 - 2 * i for i in range(N))

    def T(self, i, f, inv=False):
        out = _add({}, self.s_op(i, f), S)
        _add(out, self.divdiff(i, f), G_)
        if inv:
            _add(out, f, -G_)
        return out

    def pi1(self, f, inv=False):
        out = {}
        if not inv:
            for a, v in f.items():
                _add(out, {self._mat_apply(self.u1inv, a): v * R ** (-2 * a[-1])})
        else:
            for a, v in f.items():
                b = self._mat_apply(self.u1, a)
                _add(out, {b: v * R ** (2 * b[-1])})
        return out


def symmetrize(D, F, chain, sign):
    vals = {chain[0][0]: F}
    out = dict(F)
    for b, p, j, l in chain[1:]:
        vals[b] = D.T(j, vals[p], inv=sign < 0)
        _add(out, vals[b], S ** (sign * l))
    return out


def evaluate(D, f):
    """ev(f)：X_b ↦ t^{−(b, ρ)}、ρ = 正ルートの和の半分（2ρ = rho2）。"""
    tot = ZERO
    for a, v in f.items():
        tot = tot + v * S ** (-D._dot(a, D.rho2))
    return tot


class FastLP:
    """daha_fast.Fast の LP 版（GL_N）。"""

    def __init__(self, D, m, k=2, sign=1, final="poly"):
        self.D, self.m, self.k, self.sign, self.final = D, m, k, sign, final
        e1 = D.e(0)
        self.beta = e1 if sign > 0 else tuple(-x for x in e1)
        self.chain = orbit_chain(D, self.beta, sign)
        self.W1 = ZERO
        for *_, l in self.chain:
            self.W1 = self.W1 + S ** (2 * sign * l)
        w = [("c", R ** k), ("Y", -1)]
        for _ in range(k):
            w += [("c", R ** (-m)), ("X", e1, -1)] + [("Y", -1)] * m
        self.word = w if sign > 0 else D.inverse_word(w)
        self.one = {tuple([0] * D.N): ONE}
        self.sym = {(): self.one}
        self.Bpow = {}
        self.values = {}

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
        """{γ̂(H_ν(X) ·)} を RF で返す（daha_fast.Fast.value と同じ）。"""
        if not nu:
            return evaluate(self.D, self.sym[()]).to_rf()
        if nu not in self.values:
            d, rest = nu[0], nu[1:]
            g = self._B(rest, d - 1)
            D = self.D
            assert self.sign > 0
            if self.final == "poly":
                h = D.apply_word(self.word[2:], g)
                ev = (self.word[0][1] * S ** (-len(D.word_u1)) * evaluate(D, D.pi1(h, inv=True))).to_rf()
            else:
                ev = self.word[0][1].to_rf() * (S ** (-len(D.word_u1))).to_rf() * functional_value(D, self.word[2:], g)
            self.values[nu] = self.W1.to_rf() * ev
        return self.values[nu]


# ---- 最後の一段：評価の汎関数 ----------------------------------------------------------------------
# 点 p は各座標 p_i = q^{a_i/2} t^{b_i/2}（= r^{a_i} s^{b_i}）で、(a, b) の組で表す。p^c = r^{Σ a_i c_i} s^{Σ b_i c_i}。
def _pt_mono(p, c):
    a = sum(x * y for x, y in zip(p[0], c)); b = sum(x * y for x, y in zip(p[1], c))
    return RFm.R ** a * RFm.S ** b


def _start_point(D):
    """ev ∘ π_1^{−1}：X_c ↦ ev(q^{(…)} X_{u_1 c})。π_1^{−1}(X_c) = r^{2(u_1 c)_N} X_{u_1 c}、ev(X_b) = s^{−(b, 2ρ)}。
    c ↦ r^{2(u_1 c)_N} s^{−(u_1 c, 2ρ)} は c について乗法的なので点 p_0 の評価：p_0 の第 j 座標 = (r^{2 (u_1)_{N,j}}, s^{−Σ_i 2ρ_i (u_1)_{i,j}})。"""
    N = D.N
    U = D.u1
    a = tuple(2 * U[N - 1][j] for j in range(N))
    b = tuple(-sum(D.rho2[i] * U[i][j] for i in range(N)) for j in range(N))
    return (a, b)


def _reflect_pt(D, j, p):
    """s_j p（座標 j と j+1 の入れ替え）。"""
    a, b = list(p[0]), list(p[1])
    a[j], a[j + 1] = a[j + 1], a[j]
    b[j], b[j + 1] = b[j + 1], b[j]
    return (tuple(a), tuple(b))


def _T_dual(D, j, phi, inv):
    """φ ↦ φ ∘ T_j^{±1}。ev_p(T_j f) = t^{1/2} f(s_j p) + G (f(s_j p) − f(p))/(p^{α_j} − 1)、T^{−1} = T − G。"""
    s, G = RFm.S, RFm.S - 1 / RFm.S
    out = defaultdict(lambda: RFm.conv(0))
    for p, c in phi.items():
        x = _pt_mono(p, D.alpha[j])
        h = G / (x - 1)
        sp_ = _reflect_pt(D, j, p)
        out[sp_] = out[sp_] + c * (s + h)
        out[p] = out[p] + c * (-h - (G if inv else 0))
    return {p: c for p, c in out.items() if not c.is_zero()}


def functional_value(D, word, g, stats=None):
    """ev(π_1^{−1}(word · g)) を、汎関数 φ = ev ∘ π_1^{−1} ∘ word を左から作って φ(g) として計算する（RF）。
    word は X・c・Y^{−1}（= T_{u_1}^{−1} π_1^{−1}）だけからなること。"""
    phi = {_start_point(D): RFm.conv(1)}
    for op in word:                         # 左から：φ ∘ op
        kind = op[0]
        if kind == "c":
            cf = op[1].to_rf()
            phi = {p: c * cf for p, c in phi.items()}
        elif kind == "X":
            b = tuple(op[2] * x for x in op[1])
            phi = {p: c * _pt_mono(p, b) for p, c in phi.items()}
        elif kind == "Y":
            assert op[1] == -1
            # daha_dn.Y1(inv=True) は π_1^{−1} を先に、T^{−1} を word_u1 = (j_1, …, j_l) の順に作用させる：
            # Y^{−1} = T_{j_l}^{−1} ⋯ T_{j_1}^{−1} π_1^{−1}。よって φ ∘ Y^{−1} は T_{j_l}^{−1} から順に双対をとり、最後に π_1^{−1}
            for j in reversed(D.word_u1):
                phi = _T_dual(D, j, phi, inv=True)
            phi = _pi1inv_dual(D, phi)
        else:
            raise ValueError(op)
        if stats is not None:
            stats.append(len(phi))
    import numpy as np
    vals = list(g.values())
    A = np.array(list(g.keys()), dtype=np.int64).reshape(len(vals), D.N)
    E0 = np.array([v.er for v in vals], dtype=np.int64); E1 = np.array([v.es for v in vals], dtype=np.int64)
    tot = RFm.conv(0)
    for p, c in phi.items():                # 点ごとに g(p) を求め、RF の係数は最後に 1 回だけ掛ける
        xr = E0 + A @ np.array(p[0], dtype=np.int64); xs = E1 + A @ np.array(p[1], dtype=np.int64)
        er, es = int(xr.min()), int(xs.min())
        polys = [v.n * _mono(int(x) - er, int(y) - es) for v, x, y in zip(vals, xr, xs)]
        tot = tot + c * LP(_pairwise_sum(polys), er, es).to_rf()
    return tot


def _pairwise_sum(ps):
    """fmpz_mpoly の和を二分木で（逐次の和より大きな多項式のコピーが少ない）。"""
    if not ps:
        return ZERO.n
    while len(ps) > 1:
        ps = [ps[i] + ps[i + 1] if i + 1 < len(ps) else ps[i] for i in range(0, len(ps), 2)]
    return ps[0]


def _pi1inv_dual(D, phi):
    """φ ↦ φ ∘ π_1^{−1}。ev_p(π_1^{−1} X_c) = r^{2(u_1 c)_N} p^{u_1 c}：点 p′ で p′^c = r^{2(u_1 c)_N} p^{u_1 c}。"""
    N = D.N
    U = D.u1
    out = {}
    for p, c in phi.items():
        a = tuple(2 * U[N - 1][j] + sum(p[0][i] * U[i][j] for i in range(N)) for j in range(N))
        b = tuple(sum(p[1][i] * U[i][j] for i in range(N)) for j in range(N))
        q_ = (a, b)
        out[q_] = out.get(q_, RFm.conv(0)) + c
    return out


# ---- 基底の値と composite ----------------------------------------------------------------------------
def basis_values(N, mus, nus, m=1, k=2, final="poly"):
    """V[μ; ν] = {γ̂(H_μ(X) H_ν(X^{−1}))}（RF、ρ = (N−1, …, 0) の評価に換算）。H_ν = Π p_{ν_i}。"""
    D = GLNLP(N)
    Fm = FastLP(D, m, k, -1)
    out = {}
    for nu in nus:
        g = Fm.get_sym(nu)
        Fp = FastLP(D, m, k, 1, final=final)
        Fp.sym = {(): g}
        for mu in mus:
            c = (sum(mu) - sum(nu)) * (m * k + 1)
            out[(mu, nu)] = Fp.value(mu) * RFm.S ** (-c * (N - 1))
    return out


def composite_raw(N, alpha, beta, V):
    """{γ̂(P_α(X) P_β(X^{−1}))} = Σ x^α_μ y^β_ν V[μ; ν]（x, y は P を H の基底で展開した RF の係数）。"""
    from .daha_dn import gl_P_poly
    from .daha_fast import solve_symmetric, h_basis, Fast, _with_rho2
    D = _with_rho2(GLN(N))
    Fp, Fm = Fast(D, 1, 2, 1), Fast(D, 1, 2, -1)
    inv = lambda f: {tuple(-x for x in a): v for a, v in f.items()}
    x = solve_symmetric(D, gl_P_poly(D, alpha), h_basis(Fp, [sum(alpha)]))
    y = solve_symmetric(D, inv(gl_P_poly(D, beta)), h_basis(Fm, [sum(beta)]))
    tot = RFm.conv(0)
    for mu, xc in x.items():
        for nu, yc in y.items():
            if not xc.is_zero() and not yc.is_zero():
                tot = tot + xc * yc * V[(mu, nu)]
    return tot
