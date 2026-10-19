r"""結び目の 2 箱で、Ĩ で固定される状態の数 E を固定セクターの状態の部分和として決める（progress.md 続き 34。Hopf 版は scripts/rk_k4_fixed_sector.py）。
  g_(2) = ½E/(1+t) + C_(2)、g_(1,1) = −½E/(1+t) + C_(1,1)（全 d、E_(2) = −E_(1,1)）。
  固定セクター：被覆を brane（Y）と image（Y'）に分けたときの p_1(Y) p_1(Y') のセクター x（scripts/rk_b_cover_split.py）。Ĩ のトレースは p_2 = s_(2) − s_(1,1) に入る。
  P := (1 − t) x を固定セクターの状態の和とみて E = P − 2Q（Q は Ĩ の固有値 −1 の状態の和、0 ≤ Q_m/P_m ≤ 1 の整数）。条件：
   (i)  (1 + t) | (R − E)、R = 2(1 + t) g_(2)、
   (ii) t = q で E = −(1 + q)(g_(1)(q², ν²) + 1)（C_(2)|_{t=q} = g^{odd}_(2) − g^{odd}_(2)(unknot)、Mariño との辞書、続き 24）。
使い方：python scripts/rk_b_knot_fixed_sector.py m"""
import sys, pickle, io, contextlib
sys.path.insert(0, "src"); sys.path.insert(0, "scripts")
import sympy as sp
import rk_k4_fixed_sector as H
import rk_b_cover_split as C

r, s, nu = H.r, H.s, H.nu
u = sp.Symbol("u", positive=True)
t = s**2


def load_g(m):
    d = pickle.load(open(f"data/refined_kauffman/daha/largeN_2box_m{m}.pkl", "rb"))
    return {eval(k): sp.cancel(sp.sympify(v, locals={"r": r, "s": s, "u": u}).subs(u, nu * s)) for k, v in d["g"].items()}


def main(m, sign=1):
    """sign = +1：E_(2) = E（Ĩ のトレースが +p_2）、−1：E_(2) = −E（フェルミオン的な核で −p_2）。"""
    g = load_g(m)
    with contextlib.redirect_stdout(io.StringIO()):
        sec = C.cover_sectors(m)
    x = sp.cancel(sec["x"].to_sympy().subs(u, nu * s))
    P = sp.cancel((1 - t) * x)
    n_, d_ = sp.fraction(P)
    print(f"T(2,{2*m+1})：P = (1 − t)x の分母 {sp.factor(d_)}、項数 {len(sp.Add.make_args(sp.expand(n_)))}", flush=True)
    R = sp.cancel(2 * (1 + t) * sign * g[(2,)])
    g1x2 = g[(1,)].subs({r: r**2, s: s**2, nu: nu**2}, simultaneous=True)
    Eq = sp.cancel((-(1 + r**2) * sign * (g1x2 + 1)).subs(s, r))
    res = H.solve_generic(P, R, Eq)
    print(f"   状態の単項式 {res['n']} 個、解あり = {res['feasible']}", flush=True)
    if not res["feasible"]:
        return res
    print(f"   一意でない変数 {len(res['free'])} 個 {res['free'][:8]}", flush=True)
    E = res["E"]
    Cp = {k: sp.cancel(g[k] - sign * sg * E / (2 * (1 + t))) for k, sg in (((2,), 1), ((1, 1), -1))}
    Q = sp.expand((P - E) / 2)
    print(f"   C_(2), C_(1,1) は整数 Laurent：{[sp.Poly(sp.fraction(v)[1], r, s, nu).is_monomial for v in Cp.values()]}", flush=True)
    pickle.dump({"E": str(E), "P": str(P), "Q": str(Q), "C": {str(k): str(v) for k, v in Cp.items()}, "free": res["free"]},
                open(f"data/refined_kauffman/daha/knot_fixed_sector_m{m}_sign{sign:+d}.pkl", "wb"))
    return res


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]) if len(sys.argv) > 2 else 1)
