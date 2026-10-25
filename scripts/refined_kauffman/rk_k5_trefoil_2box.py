r"""N を固定した refined SO(2N) の trefoil の 2 箱（probe 2 個、Y の Schur 基底 |R| ≤ 2）。N = 2, 3。
  Z~(Y) = ⟨Π_{α,i}(1 − y_α x_i)(1 − y_α/x_i)⟩ = 1 − p_1 W_(1) + ½p_1²(W_(2) + A W_(1,1) + B W_(1) + C) − ½p_2(W_(2) + A' W_(1,1) + B' W_(1) + C')
  （S_1² = K_(2) + A K_(1,1) + B K_(1) + C、S_2 = K_(2) + A' K_(1,1) + B' K_(1) + C' の係数は kauffman_hopf2.koornwinder_small に U = t^N を代入）。
  F = Σ_d (1/d) Σ_R f_R(x^d) s_R(Y^d)：f_(1) = L_(1)、f_(2) = H_(1,1) + H_(2)、f_(1,1) = H_(1,1) − H_(2)、H_(1,1) = L_(1,1)、H_(2) = L_(2) − ½ f_(1)(x²)。
  双対 Cauchy なので s_(2)(Y) には e_2(x^±) = Λ²V、s_(1,1)(Y) には h_2(x^±) = Sym²V が付く。S_1² と S_2 の K_(2) の係数はどちらも 1 なので
  f_(2) は W_(2) を含まない（1 probe の計算の f_2 と同じもの）。新しいのは f_(1,1)（D_N の (2) = Sym²V − 1）。
  W_λ = m_λ R_λ K_λ(t^ρ)（R は data/torus・data/refined_kauffman/torus の reduced 不変量の鏡像 = AS の向き）。
  SO(6) = SU(4)：(1) = [1,1]、(1,1) = [2,1,1]、(2) = [2,2]。SO(4) = SU(2)²：(1) = [1]⊗[1]、(1,1) = [2]⊗∅ ⊕ ∅⊗[2]、(2) = [2]⊗[2]。
  framing の単項式 m_λ の決め方（docs/refined_kauffman/conventions.md「固定 N の framing」）：
    (i) N 依存の極の打ち消し（m_(1,1), m_(2) は極の上で決まる）、(ii) t = q で 1（data は t = q で framing 0 の HOMFLY に合わせてあり、
    Spin(6) = SU(4)・SO(4) = SU(2)² の framing 0 の不変量は一致する）。m_(1) は固定 N の極に関わらないので別に与える：
    K3' の値 (t/q)^{1/2}（大きな N の [1] の Laurent 性）か、hatted の値 t/q（2 箱の ĝ の Laurent 性、rk_k5_trefoil_2box_cover.py）。
  (i)+(ii) で m_(1,1) = (t/q)^{N−1}、m_(2) = (t/q)^2（走査で一意）。2026-10-18 の初版は (ii) を課さず、極の上で最も簡単な
  m_(1,1) = s^4（N=2）・s^16（N=3）を使った（そのとき m_(2) = s^8・r^4 s^20）。"""
import sys, json, time, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rcsknots.refined_kauffman.kauffman_hopf2 import koornwinder_small
from rcsknots.refined_kauffman.kauffman_hopf import qs, ts, Us
from rcsknots.refined_kauffman.kauffman_knot import r, s

W_ = s / r                          # (t/q)^{1/2}
M1 = W_                             # K3' の値
M1_HAT = W_**2                      # hatted の値（t/q）
M11 = {2: W_**2, 3: W_**4}          # (t/q)^{N−1}：極の上の値 s^4・s^16 に t = q で 1 の条件を課したもの
M2 = {2: W_**4, 3: W_**4}           # (t/q)^2：走査で一意（scan_m2）


KNOT = "T2_3"      # data の結び目（T2_5 なども可。lru_cache を使う側は同じプロセスで切り替えない）


def Rdata(c, Nsu):
    path = f"data/refined_kauffman/torus/{KNOT}_{c}.json"
    try:
        d = json.load(open(path))
    except FileNotFoundError:
        d = json.load(open(f"data/framed_knot/torus/{KNOT}_{c}.json"))
    A = s**Nsu
    return sp.expand(sum(v * A**(-int(k.split(",")[0])) * r**(-int(k.split(",")[1])) * s**(-int(k.split(",")[2])) for k, v in d["poly"].items()))


def coeffs(N):
    ks = koornwinder_small()
    sub = {qs: r**2, ts: s**2, Us: s**(2 * N)}
    cv = lambda v: sp.cancel(v.as_expr().subs(sub))
    return {name: {lam: cv(v) for lam, v in ks[name].items()} for name in ("S1sq", "S2")}, cv(ks["kappa"])


def unknot_K(N, lam):
    from rcsknots.refined_kauffman.kauffman_refined import koornwinder_eval, q as _q, t as _t, U as _U
    e = koornwinder_eval(lam, (1, -1, sp.sqrt(_q), -sp.sqrt(_q)))
    return sp.cancel(e.subs({_q: r**2, _t: s**2, _U: s**(2 * N)}))


def P_A(lam, n):
    """A 型 Macdonald P_λ(t^ρ)（n 変数、主特殊化）。"""
    from rk_k5_framing_scan import P_sym
    return P_sym(lam, n)


def free_energy_2(N, W):
    """W = {λ: W_λ}（λ ∈ {(1), (1,1), (2)}）→ {(1): f, (2): f, (1,1): f}。"""
    co, kappa = coeffs(N)
    assert kappa == 0
    Z1 = -W[(1,)]
    Z11 = sp.Rational(1, 2) * (W[(2,)] + co["S1sq"][(1, 1)] * W[(1, 1)] + co["S1sq"][(1,)] * W[(1,)] + co["S1sq"][()])
    Z2 = -sp.Rational(1, 2) * (W[(2,)] + co["S2"][(1, 1)] * W[(1, 1)] + co["S2"][(1,)] * W[(1,)] + co["S2"][()])
    L1, L11, L2 = Z1, Z11 - Z1**2 / 2, Z2
    f1 = sp.cancel(L1)
    H11, H2 = L11, L2 - f1.subs({r: r**2, s: s**2}, simultaneous=True) / 2
    return {(1,): f1, (2,): sp.cancel(H11 + H2), (1, 1): sp.cancel(H11 - H2)}


def trefoil_W(N, m2, m1=None, m11=None):
    m1 = M1 if m1 is None else m1
    m11 = M11[N] if m11 is None else m11
    if N == 3:
        return {(1,): m1 * Rdata("11", 4) * unknot_K(3, (1,)), (1, 1): m11 * Rdata("211", 4) * unknot_K(3, (1, 1)),
                (2,): m2 * Rdata("22", 4) * unknot_K(3, (2,))}
    return {(1,): m1 * Rdata("1", 2)**2 * unknot_K(2, (1,)), (1, 1): m11 * Rdata("2", 2) * unknot_K(2, (1, 1)),
            (2,): m2 * Rdata("2", 2)**2 * unknot_K(2, (2,))}


def laurent(e):
    return sp.Poly(sp.fraction(sp.cancel(e))[1], r, s).is_monomial


def scan_m2(N, s0s=(sp.Rational(3, 7), sp.Rational(5, 11)), arange=range(-12, 13), brange=range(-50, 51)):
    """f_(1,1) を Laurent にする m_(2) = ± r^α s^β（s を 2 つの有理数に固定して r の Laurent 性で篩う）。"""
    M = sp.Symbol("M")
    e = sp.expand(sp.cancel(free_energy_2(N, trefoil_W(N, M))[(1, 1)]))
    A1 = sp.cancel(sp.diff(e, M)); A0 = sp.cancel(e - M * A1)
    pre = {s0: (sp.cancel(A0.subs(s, s0)), sp.cancel(A1.subs(s, s0))) for s0 in s0s}
    cands = []
    for sg in (1, -1):
        for a in arange:
            for b in brange:
                if all(sp.Poly(sp.fraction(sp.cancel(a0 + sg * r**a * s0**b * a1))[1], r).is_monomial for s0, (a0, a1) in pre.items()):
                    cands.append((sg, a, b))
    return cands


def stats(e):
    """Laurent 性と整数性。分母の単項式の係数（例えば 2）で分子を割ってから係数を見る。"""
    n, d = sp.fraction(sp.cancel(e))
    D = sp.Poly(d, r, s)
    cs = [c / D.LC() for c in sp.Poly(sp.expand(n), r, s).coeffs()]
    return {"laurent": D.is_monomial, "integer": D.is_monomial and all(c.is_integer for c in cs), "terms": len(cs),
            "pos": sum(int(bool(c > 0)) for c in cs), "maxabs": max(abs(c) for c in cs)}


if __name__ == "__main__":
    from rcsknots.bps import M_matrix
    out = {}
    for N in (2, 3):
        t0 = time.time()
        # 1. 主特殊化の同一視
        if N == 3:
            ids = {(1,): P_A((1, 1), 4), (1, 1): P_A((2, 1, 1), 4), (2,): P_A((2, 2), 4)}
        else:
            P1, P2 = P_A((1,), 2), P_A((2,), 2)
            ids = {(1,): P1**2, (1, 1): 2 * P2, (2,): P2**2}
        print(f"N={N} K_λ(t^ρ) = A 型の積・和:", all(sp.cancel(unknot_K(N, lam) - v) == 0 for lam, v in ids.items()), flush=True)
        # 2. unknot
        fu = free_energy_2(N, {lam: unknot_K(N, lam) for lam in ((1,), (1, 1), (2,))})
        print(f"N={N} unknot: f_(2) = {fu[(2,)]}、f_(1,1) = {fu[(1, 1)]}", flush=True)
        # 3. m_(2) の走査
        cands = scan_m2(N)
        print(f"N={N} f_(1,1) を Laurent にする m_(2) = ±r^α s^β（|α| ≤ 12、|β| ≤ 50）: {cands}", flush=True)
        # 4. 決めた m_(2) での f
        f = free_energy_2(N, trefoil_W(N, M2[N]))
        for R in ((1,), (2,), (1, 1)):
            print(f"N={N} f_{R}: {stats(f[R])}", flush=True)
        print(f"N={N} m_(2) = 1 での f_(1,1) の Laurent 性: {laurent(free_energy_2(N, trefoil_W(N, sp.Integer(1)))[(1, 1)])}", flush=True)
        # 5. 1 probe の計算との一致
        old = pickle.load(open("data/refined_kauffman/trefoil_fixedN_f.pkl", "rb"))[N]
        for a, R in ((1, (1,)), (2, (2,))):
            print(f"N={N} f_{R} == 1 probe の f_{a}:", sp.cancel(sp.sympify(old[a], locals={"r": r, "s": s}) - f[R]) == 0, flush=True)
        # 6. hatted：M^{-1} f（参考）
        v = sp.Matrix([f[(2,)], f[(1, 1)]])
        for tt, lab in ((1 / r, "M(q^{-1})"), (s, "M(t)")):
            w = M_matrix(2, tt).inv() * v
            print(f"N={N} {lab}^{{-1}} f の分母:", [sp.factor(sp.fraction(sp.cancel(w[i]))[1]) for i in range(2)], flush=True)
        out[N] = {"m2": str(M2[N]), "f": {str(R): str(v) for R, v in f.items()}}
        print(f"  N={N}: {time.time() - t0:.0f}s", flush=True)
    pickle.dump(out, open("data/refined_kauffman/trefoil_fixedN_2box.pkl", "wb"))
