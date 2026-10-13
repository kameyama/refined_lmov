"""An. Morozov（1211.4596）の一次ずれ公式の検査：P_R = P_□^{|R|} + (ħ ν_{R^T} − ħ̄ ν_R) σ_□^{|R|−2} σ_2 + O(2)。

q = e^ħ, t = e^{ħ̄}。単項式の正規化（a^α q^β t^γ 倍）は一次の項を σ^{|R|} の倍数だけ変えるので、
各 R について「δ_R − (βħ + γħ̄)|R|… 」のような自由度を許して比較する：
  η_R, η̄_R を P_R の ħ, ħ̄ 一次の係数（a の多項式）とし、η_R − |R| σ^{|R|−1} η_□ などの差を作ると、
  公式は  η_R − |R|σ^{|R|−1}η_□ = ν_{R^T} σ^{|R|−2} σ_2 + c_R σ^{|R|},  η̄_R − |R|σ^{|R|−1}η̄_□ = −ν_R σ^{|R|−2} σ_2 + c̄_R σ^{|R|}
  （c_R, c̄_R は正規化に由来する定数）を予言する。σ_2 は R = [2] の ħ 部分から決める。
使い方：python scripts/framed_knot/morozov_check.py data/framed_knot/cables/K1  （{prefix}_{color}.json を読む）
"""
import json
import sys

import sympy as sp

A = sp.Symbol("A")


def load(prefix, color):
    d = json.load(open(f"{prefix}_{color}.json"))
    return {tuple(map(int, k.split(","))): v for k, v in d["poly"].items()}


def parts(poly):
    sigma = sum(c * A ** ea for (ea, eq, et), c in poly.items())
    eta = sum(sp.Rational(eq, 2) * c * A ** ea for (ea, eq, et), c in poly.items())
    etab = sum(sp.Rational(et, 2) * c * A ** ea for (ea, eq, et), c in poly.items())
    return sp.expand(sigma), sp.expand(eta), sp.expand(etab)


def nu(lam):
    return sum(i * x for i, x in enumerate(lam))


def main(prefix):
    s1, e1, eb1 = parts(load(prefix, "1"))
    res = {}
    data = {}
    for color, lam in [("2", (2,)), ("11", (1, 1)), ("21", (2, 1)), ("3", (3,)), ("111", (1, 1, 1))]:
        s, e, eb = parts(load(prefix, color))
        k = sum(lam)
        assert sp.expand(s - s1 ** k) == 0, "special polynomial property fails"
        data[lam] = (sp.expand(e - k * s1 ** (k - 1) * e1), sp.expand(eb - k * s1 ** (k - 1) * eb1), k)
    # σ_2 from [2]: δη_[2] = ν_{[1,1]} σ_2 + c σ^2 = σ_2 + c σ^2  →  σ_2 is defined modulo σ^2
    sigma2 = data[(2,)][0]
    report = []
    for lam, (de, deb, k) in data.items():
        lamT = tuple(sum(1 for x in lam if x > j) for j in range(lam[0]))
        for name, val, coeff in (("hbar", de, nu(lamT)), ("hbar_bar", deb, -nu(lam))):
            rem = sp.expand(val - coeff * s1 ** (k - 2) * sigma2)
            # allowed: (constant) σ^k, plus the ambiguity of σ_2 (mod σ^2) times coeff σ^{k-2}
            c = sp.Symbol("c")
            ok = sp.solve(sp.Poly(sp.expand(rem - c * s1 ** k), A).all_coeffs(), c) != []
            report.append((lam, name, ok))
    for r in report:
        print(*r)
    return report


if __name__ == "__main__":
    main(sys.argv[1])
