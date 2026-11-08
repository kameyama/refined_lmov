"""ν^{−14} の E（続き 116、theory §12.46）：ν^{−8}..ν^{−12} の E（上端 §12.26、λ = (1,1) §12.27、λ = (2) の全部入れ替え §12.44）を固定し、
ν^{−14} の両 λ の e を自由（|e| ≤ |c|、偶奇）にして、条件 (i)・(i)+(ii) の解の有無、一意な単項式、入れ替わる対の最小を整数計画で求める。
一意な単項式を、全部入れ替え ＋ 符号 (−1)^{k/2+1} の予想（§12.45）と照合する。
実行：PYTHONPATH=src:scripts/refined_kauffman python scripts/refined_kauffman/rk_f_n2_E14.py"""
import sys
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import numpy as np
import sympy as sp
from scipy.optimize import milp, LinearConstraint, Bounds
import rk_f_n2_rterms as NR
from rk_f_n2_rterms import r_, s_, _e_of, _top_E, _E11_12, _n2_rows, _ranges, _poly_of

q, t = r_**2, s_**2


def E12_2():
    """λ = (2) の ν^{−12} の全部入れ替えの E（§12.44 (3)・§12.45）：{key: e}。"""
    import rk_f_n2_EG_adj as EG
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        E12 = EG.part_C()
    return _e_of(r_**-17 * s_**3 * E12, "2", -12)          # E12 は r_, s_ の式（q = r_², t = s_²）


def main():
    fixed = {**_top_E(), **_E11_12(), **E12_2()}
    for with_ii in (False, True):
        mult, keys, A, b, idx, sel, rows_, use = _n2_rows(with_ii, -14)
        # 固定した鍵が台に無い／|e| > |c| なら警告
        bad = [k for k, v in fixed.items() if (k not in idx and v != 0) or (k in idx and abs(v) > abs(mult[k]))]
        print(f"== 条件 {use}：行 {len(rows_)} 本、ν^-8..ν^-14 の単項式 {len(sel)} 個、固定 {len(fixed)}（台の外・違反 {len(bad)}）", flush=True)
        unknown = [k for k in sel if k[4] == -14]
        print(f"   ν^-14 の未知数 {len(unknown)}（λ = (2)：{sum(1 for k in unknown if k[0] == '2')}、(1,1)：{sum(1 for k in unknown if k[0] == '11')}）", flush=True)
        # 固定部分の整合（ν^-14 を含まない行）
        rows_fixed = [i for i in rows_ if all(A[i, idx[k]] == 0 for k in unknown)]
        res = b[rows_fixed] - sum(A[rows_fixed][:, idx[k]] * v for k, v in fixed.items() if k in idx)
        print(f"   ν^-14 を含まない行 {len(rows_fixed)} 本の残差の最大 {np.max(np.abs(res)) if len(rows_fixed) else 0}", flush=True)
        rg = _ranges(mult, keys, A, b, rows_, idx, unknown, {k: v for k, v in fixed.items() if k in idx})
        if rg is None:
            print("   解なし", flush=True)
            continue
        for lab in ("2", "11"):
            un = [k for k in unknown if k[0] == lab]
            uniq = {k: rg[k][0] for k in un if rg[k][0] == rg[k][1]}
            print(f"   λ = {lab}：一意 {len(uniq)}/{len(un)}", flush=True)
            Eu = sum(v * r_**k[2] * s_**k[3] for k, v in uniq.items())
            Nu = sum(mult[k] * r_**k[2] * s_**k[3] for k in uniq)
            print(f"      一意な部分の E = {sp.factor(sp.expand(Eu))}", flush=True)
            print(f"      同じ単項式の N = {sp.factor(sp.expand(Nu))}", flush=True)
            amb = {k: (mult[k], rg[k]) for k in un if rg[k][0] != rg[k][1]}
            print(f"      不定：{ {(k[2], k[3]): v for k, v in amb.items()} }", flush=True)
        # 入れ替わる対の最小
        c = np.array([float(mult[k]) for k in unknown]); sg = np.sign(c); ca = np.abs(c)
        Ar = A[rows_][:, [idx[k] for k in unknown]]
        rhs = b[rows_] - sum(A[rows_][:, idx[k]] * v for k, v in fixed.items() if k in idx)
        M = Ar * (-2 * sg); beq = rhs - Ar @ c
        cons = LinearConstraint(M, beq - 1e-9, beq + 1e-9); bnd = Bounds(np.zeros(len(unknown)), ca)
        pm = milp(np.ones(len(unknown)), constraints=cons, bounds=bnd, integrality=np.ones(len(unknown)))
        print(f"   入れ替わる対の最小 {round(pm.fun)}（Σ|c| = {int(ca.sum())}）", flush=True)


if __name__ == "__main__":
    main()
