"""別研究 [2,1] の X2：[2,1] のチャネル分解に EGP（q = 1）と U(2)（a = t²）を課したときの整合性と、正則なチャネルの Z̃_X の決まり方。

仮定（チャネル分解、conventions.md §16）：
  rCS_{[2,1]}(K_{m,n}) = (T/Q)^{c (sgn m + sgn n)/2} [ Σ_{X ∈ {∅,[1],[1,1],[2],[2,1]}} Z̃_X τ_X(m) τ_X(n)/τ_X(1) + W(m, n) ]
  τ_X は KNTZ の refined の正則な行（X1）。W は反常チャネル X̃₂ の寄与（固有値 ẽ = A⁴(Q/T)^k、係数は未知）。
条件（どちらも全ての (m, n) で）：
  EGP：Q = 1 で rCS_{[2,1]} = rCS_{[1]}·rCS_{[1,1]}（KN (EGP)、conventions.md §10）
  U(2)：A² = T⁴ で rCS_{[2,1]} = rCS_{[1]}·rCS_{[1,1]}（REPORT_phase4.md §4）
  右辺の rCS_{[1]}, rCS_{[1,1]} は KNTZ の refined の長方形の式（X1 で GNSSS・EHA と照合済み）。
固有値の組 (e_Y^m, e_{Y'}^n) ごとに係数を比べる（特殊化で固有値が一致するものはまとめる）。
未知数 z_X = Z̃_X（その特殊化での値）について線形な連立方程式を作り、解けるか（整合性）と解を出す。
使い方：python scripts/x2_specializations.py
"""
import itertools
import sys

import sympy as sp

from rcsknots.double_twist_family import A, Q, T, U_row, eig, tre0, M_rect

REG = [(), (1,), (1, 1), (2,), (2, 1)]


def grouped(entries, sub):
    """{(e_m, e_n): 係数} を特殊化 sub の後の固有値の値でまとめる。entries は [(e_m, e_n, 係数)]。"""
    out = {}
    for em, en, c in entries:
        key = (sp.simplify(em.subs(sub)), sp.simplify(en.subs(sub)))
        out[key] = out.get(key, 0) + c
    return out


def rhs_entries():
    """rCS_{[1]}·rCS_{[1,1]} の指数和の項 [(固有値 m, 固有値 n, 係数)]（N は (T/Q)^{1+2}）。"""
    M1, M11 = M_rect(1, 1), M_rect(1, 2)
    out = []
    for ((a1, b1), c1), ((a2, b2), c2) in itertools.product(M1.items(), M11.items()):
        out.append((eig(a1) * eig(a2), eig(b1) * eig(b2), c1 * c2))
    return out


def lhs_entries(z):
    out = []
    for X in REG:
        U = U_row(X)
        for mu, u in U.items():
            for nu, v in U.items():
                out.append((eig(mu), eig(nu), z[X] * u * v / tre0(X)))
    return out


def solve_at(sub, label, extra=()):
    """特殊化 sub で LHS = RHS を z_X について解く。extra は W の項（[(固有値, 固有値, 係数)]）。"""
    z = {X: sp.Symbol(f"z_{''.join(map(str, X)) or '0'}") for X in REG}
    L = grouped(lhs_entries(z) + list(extra), sub)
    R = grouped(rhs_entries(), sub)
    eqs = []
    for key in set(L) | set(R):
        e = sp.cancel(sp.together((L.get(key, 0) - R.get(key, 0)).subs(sub)))
        if e != 0:
            eqs.append(e)
    sol = sp.solve(eqs, list(z.values()), dict=True)
    print(f"{label}：固有値の組 {len(set(L) | set(R))}、式 {len(eqs)}、未知数 {len(z)}、解 {len(sol)}", flush=True)
    for s in sol:
        for k, v in s.items():
            print(f"   {k} = {sp.factor(v)}", flush=True)
    return sol


if __name__ == "__main__":
    solve_at({Q: 1}, "EGP（Q = 1、W = 0 と仮定）")
    solve_at({A: T**2}, "U(2)（A = T²、W = 0 と仮定）")
    sys.stdout.flush()
