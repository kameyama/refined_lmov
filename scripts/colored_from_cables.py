"""2 本ケーブルの HOMFLY-PT からケーブル公式を逆に解き、[2], [1,1] 彩色 HOMFLY-PT（reduced、framing 0）を求める（ステージ T4）。

入力：data/reference/two_cable_homfly.txt（scripts/two_cable_homfly.sage、s = 1, 3）。
W_□(C_{2,s}(K)) = θ_□^{2s}·Σ_μ c^{(2)}_{□μ} θ_μ^{−s/2} W_μ(K)（src/rcsknots/homfly.py の cable_unreduced、sign = −1、ステージ 4 で検証済み）。
変数は KN の (a, Q = q^{1/2}) で、Sage の値は tests/test_stage4_cable_homfly.py の sage_to_kn で変換する（DGR の向き）。
使い方：python scripts/colored_from_cables.py 3_1 4_1 9_42
"""
import re
import sys

import sympy as sp

from rcsknots.homfly import Q, a, adams_coeffs, theta, unknot

L, M, Abold = sp.symbols("L M A")


def sage_to_kn(expr):
    e = sp.expand(expr.subs({L: sp.I * Abold, M: -sp.I * (Q - 1 / Q)}, simultaneous=True))
    e = e.subs({Abold: 1 / Abold, Q: 1 / Q}, simultaneous=True)
    return sp.expand(e.subs(Abold, sp.sqrt(a)))


OVERFLOW = 2**62
SUSPECT = []   # load が除いた（名前, 傾き）


def load(path="data/reference/two_cable_homfly.txt", names=None):
    """ケーブルの HOMFLY-PT を読む。Sage（libhomfly）は係数を 64 ビット整数で持ち、あふれても警告しないので、
    係数（Sage の L, M での値）が 2^62 以上の行は「あふれの疑い」として除き SUSPECT に記録する（X0 の K_{−7,1} の傾き 5 で実際にあふれた）。names を与えるとその結び目だけを読む。"""
    out = {}
    for line in open(path):
        name, s, h = line.rstrip("\n").split("\t")
        if names is not None and name not in names:
            continue
        if max(int(x) for x in re.findall(r"\d+", h)) >= OVERFLOW:
            SUSPECT.append((name, s))
            continue
        s = tuple(int(x) for x in s.split(",")) if "," in s else int(s)
        out[(name, s)] = sage_to_kn(sp.sympify(h, locals={"L": L, "M": M}))
    return out


def colored(name, data=None):
    """{(2,): P_[2], (1,1): P_[1,1]}（reduced、framing 0、KN の (a, Q)）。"""
    data = data or load()
    c = adams_coeffs((1,), 2)
    mus = [(2,), (1, 1)]
    W = sp.symbols("W2 W11")
    eqs = []
    for s in (1, 3):
        lhs = data[(name, s)] * unknot((1,))
        rhs = theta((1,), 2 * s) * sum(c[mu] * theta(mu, sp.Rational(-s, 2)) * w for mu, w in zip(mus, W))
        eqs.append(sp.Eq(lhs, rhs))
    sol = sp.solve(eqs, W, dict=True)[0]
    return {mu: sp.factor(sp.simplify(sol[w] / unknot(mu))) for mu, w in zip(mus, W)}


def colored_p(name, p, path=None):
    """p 本ケーブル（data/reference/cable_homfly_p{p}.txt の全ての s）から |λ| = p の彩色 HOMFLY-PT（reduced、framing 0）。"""
    from rcsknots.partitions import partitions
    data = load(path or f"data/reference/cable_homfly_p{p}.txt")
    slopes = sorted(s for (n, s) in data if n == name)
    c = adams_coeffs((1,), p)
    mus = [mu for mu in partitions(p) if mu in c]
    W = sp.symbols(" ".join(f"W{i}" for i in range(len(mus))))
    W = W if isinstance(W, tuple) else (W,)
    eqs = []
    for s in slopes:
        lhs = data[(name, s)] * unknot((1,))
        rhs = theta((1,), p * s) * sum(c[mu] * theta(mu, sp.Rational(-s, p)) * w for mu, w in zip(mus, W))
        eqs.append(sp.Eq(lhs, rhs))
    sol = sp.solve(eqs, W, dict=True)[0]
    return {mu: sp.factor(sp.simplify(sol[w] / unknot(mu))) for mu, w in zip(mus, W)}


def _cable_coeffs(p, s):
    """W_□(ケーブル) = pre · Σ_μ coef[μ] W_μ(K) の (pre, coef)。s が整数なら (p, s)-ケーブル、(s, s') なら C_{2,s'}(C_{2,s}K)（p = 4）。"""
    if isinstance(s, int):
        c = adams_coeffs((1,), p)
        return theta((1,), p * s), {mu: cm * theta(mu, sp.Rational(-s, p)) for mu, cm in c.items()}
    s1, s2 = s
    coef = {}
    for nu, cn in adams_coeffs((1,), 2).items():
        for mu, cm in adams_coeffs(nu, 2).items():
            coef[mu] = coef.get(mu, 0) + cn * theta(nu, sp.Rational(-s2, 2)) * theta(nu, 2 * s1) * cm * theta(mu, sp.Rational(-s1, 2))
    return theta((1,), 2 * s2), coef


def colored_fast(name, data, p, solve=None):
    """|μ| = p の彩色 HOMFLY-PT（reduced、framing 0）を、係数行列の逆行列を先に作り cancel 1 回で求める（colored_p の高速版）。

    solve：解くのに使う傾き（既定は先頭から未知数の個数だけ）。残りの傾きでは予測値とデータの一致を確かめる（不一致なら AssertionError）。
    戻り値：({μ: P_μ}, 確かめた傾きの数)。
    """
    slopes = sorted(s for (n, s) in data if n == name)
    pre0, coef0 = _cable_coeffs(p, slopes[0])
    mus = [mu for mu in coef0 if coef0[mu] != 0]
    solve = list(solve or slopes[:len(mus)])
    rows, rhs = [], []
    for s in solve:
        pre, coef = _cable_coeffs(p, s)
        rows.append([coef.get(mu, 0) for mu in mus])
        rhs.append(sp.expand(data[(name, s)] * unknot((1,)) / pre))
    Minv = sp.Matrix(rows).inv().applyfunc(sp.simplify)
    out = {}
    for i, mu in enumerate(mus):
        w = sum(Minv[i, j] * rhs[j] for j in range(len(solve)))
        out[mu] = sp.expand(sp.cancel(sp.together(w / unknot(mu))))
    checked = 0
    for s in slopes:
        if s in solve:
            continue
        pre, coef = _cable_coeffs(p, s)
        pred = pre * sum(coef.get(mu, 0) * out[mu] * unknot(mu) for mu in mus) / unknot((1,))
        assert sp.expand(sp.cancel(sp.together(pred)) - data[(name, s)]) == 0, (name, s)
        checked += 1
    return out, checked


if __name__ == "__main__":
    data = load()
    for name in sys.argv[1:]:
        for mu, p in colored(name, data).items():
            print(name, mu, sp.expand(p))
