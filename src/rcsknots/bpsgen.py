"""任意の箱数の \\hat f_ρ（評価点ごと）。bps.py / bpsnum.py の |μ| ≤ 3 の記号式を、箱数によらない数値の手順に置き換えたもの。

評価点 (A, Q, T)（A = a^{1/2}, Q = q^{1/2}, T = t^{1/2}、法 p の数値の配列）で次を計算する。
  Z(x) = 1 + Σ_λ \\bar{rCS}_λ g_λ P_λ(x; q, t)、\\bar{rCS}_λ = (unknot_λ) · rCS_λ（bpsnum.unknot_sym と同じ）
  Log Z = Σ_d μ(d)/d ψ_d log Z（ψ_d：p_k ↦ p_{dk}、(A, Q, T) ↦ d 乗）
  f^q の p 基底の係数 G_κ = (Q − 1/Q)(T/Q)^{|κ|} (Log Z)_κ（bps.f_formulas と同じ尺度）
  \\hat f：f^q_μ = Σ_ρ M_{μρ}(t) \\hat f_ρ の M(t) は B(t) = Σ_σ B_σ s_σ との Kronecker 積なので、類関数の基底で対角：
      G_κ = B(κ) H_κ、B(κ) = T^{−(n−1)} Π_i (1 − T^{2κ_i}) / (1 − T²)（hook の和の母関数、分母は 1 回）、
      \\hat f_ρ = Σ_κ χ_ρ(κ) H_κ。
  検査：f^{\\bar t} = M(q^{−1}) \\hat f（KN の構造、bps.py と同じ）。
rCS_λ（reduced、framing 0）は関数 rcs(lam, Qv, Tv, d) で与える（A の多項式の係数のリスト、a = A² の冪ではなく A の冪）。
"""
from fractions import Fraction
from functools import lru_cache

import numpy as np
import sympy as sp

from .characters import chi
from .field import V
from .macdonald import macdonald_P, macdonald_norm
from .partitions import partitions, conjugate, boxes, arm, leg, n_stat


@lru_cache(maxsize=None)
def _parts_upto(N):
    return [lam for k in range(1, N + 1) for lam in partitions(k)]


def _mul_series(X, Y, nmax):
    """p 基底の級数の積（次数 nmax まで）。X, Y：{ρ: V}。"""
    out = {}
    for r1, x in X.items():
        for r2, y in Y.items():
            if sum(r1) + sum(r2) <= nmax:
                key = tuple(sorted(r1 + r2, reverse=True))
                out[key] = out[key] + x * y if key in out else x * y
    return out


def _log_graded(Z, nmax, zero):
    """log Z（Z = 1 + Σ_{n≥1} Z_n、p 基底）を次数ごとの漸化式 n L_n = n Z_n − Σ_{k<n} k L_k Z_{n−k} で。"""
    byn = {n: {r: c for r, c in Z.items() if sum(r) == n} for n in range(1, nmax + 1)}
    L = {}
    for n in range(1, nmax + 1):
        acc = {r: c * n for r, c in byn[n].items()}
        for k in range(1, n):
            Lk = {r: c * k for r, c in L[k].items()}
            for key, v in _mul_series(Lk, byn[n - k], n).items():
                acc[key] = acc[key] - v if key in acc else -v
        L[n] = {r: c / n for r, c in acc.items()}
    return L


def unknot_bar(lam, A, Q, T):
    """(t/a)^{|λ|/2} Π (t^{l'} − a q^{a'})/(1 − q^a t^{l+1})（bpsnum.unknot_sym の数値版）。"""
    q, t, a = Q * Q, T * T, A * A
    out = T ** sum(lam) / A ** sum(lam)
    for i, j in boxes(lam):
        out = out * (t ** i - a * q ** j) / (1 - q ** arm(lam, i, j) * t ** (leg(lam, i, j) + 1))
    return out


def _horner_A(coeffs, A):
    """Σ_k coeffs[k] A^k（coeffs は V または整数）。"""
    out = A * 0
    for c in reversed(coeffs):
        out = out * A + c
    return out


class Points:
    """評価点 (A, Q, T) の配列と、その上の \\hat f の計算。

    rcs(lam, Qv, Tv) -> (lo, [V, ...])：reduced rCS_λ の A の冪 lo, lo+1, ... の係数（Q, T の点の上の V）。
    (Q, T) の点は長さ P、A の点は長さ Na。内部の配列は (Q, T) の点ごとに A の点を並べた長さ P·Na。
    """

    def __init__(self, rcs, Qv, Tv, Av):
        from .field import make_conv
        self.rcs, self.conv = rcs, make_conv(Qv.p, len(Qv.a))      # conv：(Q, T) の点の長さ P の定数
        self.P, self.Na = len(Qv.a), len(Av)
        self.Q0, self.T0 = Qv, Tv
        self.p = Qv.p
        self.Aarr = np.asarray(Av, dtype=np.int64)
        self._pts, self._rcs = {}, {}

    def _QT(self, d):
        if d not in self._pts:
            self._pts[d] = (self.Q0 ** d, self.T0 ** d)
        return self._pts[d]

    def reduced(self, lam, d):
        """rCS_λ を (Q^d, T^d) の点で（A の冪 lo からの係数）。λ と d ごとに保存して f^q と f^{\\bar t} で共有する。"""
        if (d, lam) not in self._rcs:
            self._rcs[(d, lam)] = self.rcs(lam, *self._QT(d))
        return self._rcs[(d, lam)]

    def _tile(self, v):
        return V(np.repeat(v.a, self.Na), self.p)

    def _A(self, d):
        return V(np.tile(np.array([pow(int(x), d, self.p) for x in self.Aarr], dtype=np.int64), self.P), self.p)

    def zseries(self, d, nmax, kind="q"):
        """ψ_d を掛ける前の Z（(A, Q, T) を d 乗した点での値）の次数 nmax まで。kind="t" は f^{\\bar t} の母関数。"""
        Qd, Td = self._QT(d)
        q, t = Qd * Qd, Td * Td
        conv = self.conv
        A = self._A(d)
        Qt, Tt = self._tile(Qd), self._tile(Td)
        Z = {}
        for n in range(1, nmax + 1):
            Pm = macdonald_P(n, t, q, conv) if kind == "t" else macdonald_P(n, q, t, conv)
            for lam in partitions(n):
                lo, cs = self.reduced(lam, d)
                red = _horner_A([self._tile(c) for c in cs], A) * A ** lo
                cbar = unknot_bar(lam, A, Qt, Tt) * red
                if kind == "t":
                    w = cbar * (-1) ** n
                    vec = Pm[conjugate(lam)]
                else:
                    w = cbar * self._tile((q / t) ** n / macdonald_norm(lam, q, t))
                    vec = Pm[lam]
                for rho, c in vec.items():
                    term = w * self._tile(c)
                    Z[rho] = Z[rho] + term if rho in Z else term
        return Z

    def log_series(self, nmax, kind="q"):
        """Log Z（plethystic）の次数 nmax まで、p 基底 {κ: 値}。"""
        zero = self._A(1) * 0
        out = {}
        for d in range(1, nmax + 1):
            mu = int(sp.mobius(d))
            if mu == 0:
                continue
            L = _log_graded(self.zseries(d, nmax // d, kind), nmax // d, zero)
            for n, Ln in L.items():
                for r, c in Ln.items():
                    key = tuple(d * k for k in r)
                    term = c * Fraction(mu, d)
                    out[key] = out[key] + term if key in out else term
        return out

    def hat_f(self, nmax, check_ft=True):
        """{n: {ρ: V}}（ρ ⊢ n ≤ nmax の \\hat f_ρ、a 次数のずらし無し）と、f^{\\bar t} の検査の結果。"""
        Q, T = self._tile(self.Q0), self._tile(self.T0)
        Lq = self.log_series(nmax, "q")
        Lt = self.log_series(nmax, "t") if check_ft else None
        out, ok = {}, True
        for n in range(1, nmax + 1):
            H = {}
            for kappa in partitions(n):
                G = Lq.get(kappa, Q * 0) * (Q - 1 / Q) * (T / Q) ** n
                H[kappa] = G / B_class(kappa, T)
                if check_ft:
                    Gt = Lt.get(kappa, Q * 0) * (1 / T - T)
                    if Gt != B_class(kappa, 1 / Q) * H[kappa]:
                        ok = False
            out[n] = {rho: sum((H[k] * chi(rho, k) for k in partitions(n)), Q * 0) for rho in partitions(n)}
        return out, ok


def B_class(kappa, T):
    """B(t) = Σ_σ B_σ(t) χ_σ の類 κ での値 = T^{−(n−1)} Π_i (1 − T^{2κ_i}) / (1 − T²)（T = t^{1/2}）。"""
    n = sum(kappa)
    out = T ** (-(n - 1)) / (1 - T * T)
    for k in kappa:
        out = out * (1 - T ** (2 * k))
    return out


class EpsTable:
    """ε(H̃_μ) = Π_{(r,c)∈μ} (1 − a q^c t_w^r)（a の多項式）の表を、箱数ごとに 1 つ前の箱数の表から作る。

    gn.GN.evaluate_apoly と同じ値。μ ごとに 1 回の 1 次式の掛け算で済み、要求された箱数の表だけを残す（メモリ節約）。
    """

    def __init__(self, G):
        self.G, self.level, self.cur, self.kept = G, 0, {(): [G.one]}, {0: {(): [G.one]}}
        self._x = {}

    def _chi(self, r, c):
        if (r, c) not in self._x:
            self._x[(r, c)] = self.G.q ** c * self.G.tw ** r * self.G.one
        return self._x[(r, c)]

    def table(self, L):
        if L in self.kept:
            return self.kept[L]
        if L < self.level:
            self.level, self.cur = 0, {(): [self.G.one]}
        while self.level < L:
            nxt = {}
            for mu in partitions(self.level + 1):
                r = len(mu) - 1                       # 最後の行の最後の箱を除いた形から
                c = mu[-1] - 1
                parent = mu[:-1] if mu[-1] == 1 else mu[:-1] + (mu[-1] - 1,)
                poly, x = self.cur[parent], self._chi(r, c)
                nxt[mu] = [poly[0]] + [poly[i] - poly[i - 1] * x for i in range(1, len(poly))] + [-(poly[-1] * x)]
            self.level, self.cur = self.level + 1, nxt
        self.kept[L] = self.cur
        return self.cur

    def evaluate(self, vec):
        out = []
        for mu, c in vec.items():
            poly = self.table(sum(mu))[mu]
            while len(out) < len(poly):
                out.append(self.G.zero)
            for i, y in enumerate(poly):
                out[i] = out[i] + c * y
        return out


# ---- トーラス結び目の rCS_λ（GN の値 × conventions.md §13 の閉じた正規化式） ---------------------

def torus_rcs_factory(m, n, workers_cache=None):
    """T(m, n)（GN の傾き (m, −n)）の reduced rCS_λ を (Q, T) の点で返す関数。

    正規化（conventions.md §13）：rCS_λ = M_□^{|λ|} T_λ^{mn} × (GN の reduced 値)、
    M_□ = (−q/t)^{m−1} (a (q/t)^{1/2})^δ、δ = (m−1)(n−1)/2、T_λ = q^{n(λ')} t^{−n(λ)}。
    同じ (Q, T) の点（オブジェクトの id）では、p_ρ の演算子の作用 v_ρ = Π_k P̃_{k(m,−n)}·1 を λ をまたいで共有する。
    """
    from .gn import GN
    from .field import make_conv
    delta = (m - 1) * (n - 1) // 2
    state = {}

    def rcs(lam, Qv, Tv):
        key = (id(Qv), id(Tv))
        if state.get("key") != key:
            k = len(Qv.a)
            G = GN(Qv * Qv, Tv * Tv, V(np.full(k, 3), Qv.p), make_conv(Qv.p, k))
            state.clear()
            state.update(key=key, G=G, v={(): {(): G.one}}, e={}, eps=EpsTable(G), keep=(Qv, Tv))
        G, vs, es = state["G"], state["v"], state["e"]

        def v_of(rho):
            if rho not in vs:
                vs[rho] = G.op(rho[0], m, -n, v_of(rho[1:]))
            return vs[rho]

        def e_of(rho):                                        # ε(v_ρ) の a の係数（λ をまたいで共有）
            if rho not in es:
                es[rho] = state["eps"].evaluate(v_of(rho))
            return es[rho]

        apoly = []
        for rho, c in G.P(sum(lam))[lam].items():
            for i, x in enumerate(e_of(rho)):
                if i < len(apoly):
                    apoly[i] = apoly[i] + c * x
                else:
                    apoly.append(c * x)
        red = G.reduce_apoly(lam, apoly)                      # a の冪の係数
        size = sum(lam)
        q, t = Qv * Qv, Tv * Tv
        mono = ((-1) ** (m - 1) * (q / t) ** (m - 1) * (Qv / Tv) ** delta) ** size
        mono = mono * (q ** n_stat(conjugate(lam)) * t ** (-n_stat(lam))) ** (m * n)
        coeffs = []
        for c in red:                                         # a^k = A^{2k}
            coeffs += [c * mono, Qv * 0]
        return 2 * delta * size, coeffs[:-1]

    return rcs
