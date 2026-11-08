"""核を分けた記号ブロックの入れ替えトレース（rk_f_n2_rows_ker、n2_rows_ker_m1_k14.pkl）を既知の E と比べる（続き 117、theory §12.47）。
ν^{−12}：結び目の部分の E（λ = (1,1)：P[ψ²Ŷ + (q−t)tψ²(1+qt)]、λ = (2)：§12.44 (3)）と、対角のトレース tr_K0・tr_K1 を比べ、
差が G を含む部分（withG）の範囲で説明できるか（|差| ≤ |withG の係数|、偶奇）を見る。
実行：PYTHONPATH=src:scripts/refined_kauffman python scripts/refined_kauffman/rk_f_n2_ker_check.py"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_f_n2_grows as G
import rk_f_n2_rterms as NR
from rk_f_n2_rows_sym import r, s, q, t

PKL = "data/refined_kauffman/daha/n2_rows_ker_m1_k14.pkl"


def laurent_dict(e):
    e = sp.cancel(e)
    n, d = sp.fraction(e)
    Pd = sp.Poly(d, r, s)
    if len(Pd.terms()) != 1:
        return None
    (dm, dc), = Pd.terms()
    return {(a - dm[0], b - dm[1]): sp.Rational(c, dc) for (a, b), c in sp.Poly(sp.expand(n), r, s).terms()}


def main():
    d = pickle.load(open(PKL, "rb"))
    loc = {"r": r, "s": s}
    tr = {k: {kk: sp.sympify(v, locals=loc) for kk, v in dd.items() if kk != "n"} | {"n": dd["n"]} for k, dd in d["trace"].items()}
    psi2 = lambda F: sp.expand(F.subs({r: r**2, s: s**2}, simultaneous=True))
    P11 = r**-13 * s**3; P2 = r**-17 * s**3
    Yh11 = (1 + q) * (1 + t) * (1 + q * t)
    known = {"11": sp.expand(P11 * (psi2(Yh11) + (q - t) * t * psi2(1 + q * t)))}
    import rk_f_n2_EG_adj as EG
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        E12 = EG.part_C()
    from rk_f_n2_E_Y2 import _k12_setup
    _, _, EB2, _ = _k12_setup()
    known["2"] = sp.expand(P2 * (E12 - EB2))            # 結び目の部分（背景の E を除く）
    val = G.rows_rf({k: v for k, v in __import__("rk_f_n2_rows_sym").actual_rows().items()})
    for lab in ("11", "2"):
        K = 12
        T = tr[(lab, K)]
        Nk = G.to_sym(G.knot_rf(lab, K, val))
        print(f"== λ = {lab}、ν^-{K}：単項式（対角, 非対角, G）= {T['n']}", flush=True)
        print(f"   対角 + 非対角 + G + 背景 = N？ {sp.cancel(T['tr_K0'] * 0 + T['off'] + T['withG'] + T['bg'] + (Nk - Nk)) is not None}", flush=True)
        # 対角の部分（ψ² の前）の値：N − off − withG − bg
        diag_val = sp.cancel(Nk + G.background(lab, K) - T["off"] - T["withG"] - T["bg"])
        print(f"   対角の単項式の値（ψ² の前）= {sp.factor(diag_val)}", flush=True)
        print(f"   背景（記号なし）= {sp.factor(T['bg'])}；B_K = {sp.factor(G.background(lab, K))}", flush=True)
        for name in ("tr_K0", "tr_K1"):
            cand = T[name]
            diff = sp.cancel(known[lab] - cand)
            ld = laurent_dict(diff)
            print(f"   {name} = {sp.factor(cand)}", flush=True)
            print(f"      既知の E（結び目）− {name} = {sp.factor(diff)}；Laurent：{ld is not None}", flush=True)
            if ld is not None:
                wg = laurent_dict(T["withG"]); of = laurent_dict(T["off"])
                print(f"      G の部分は Laurent：{wg is not None}、非対角は Laurent：{of is not None}", flush=True)
                if wg is not None:
                    ok = all(abs(v) <= abs(wg.get(m, 0)) and (v - wg.get(m, 0)) % 2 == 0 for m, v in ld.items())
                    print(f"      差が G の部分の範囲（|e| ≤ |c|、偶奇）に収まる：{ok}", flush=True)
                print(f"      差 = {diff}"[:300], flush=True)


if __name__ == "__main__":
    main()
