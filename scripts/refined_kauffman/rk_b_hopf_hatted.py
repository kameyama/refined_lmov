r"""B：Hopf 絡み目の reformulated（hatted）クロスキャップ不変量 ĝ_{R,[1]} = Σ_S M^{−1}_{RS} g_{S,1}（成分 1 は Schur 基底の多 probe、成分 2 は [1]）。
M は KN の M(x)（src/rcsknots/bps.py の M_matrix、x^{1/2} を渡す）で、t-brane 型の核には M(q^{−1})（KN の f^{\bar t} 側、Mariño の M_{RS} の refined 版）。
判定：Laurent・整数性、t = q で Mariño の ĝ（R ↔ R^T、全体の符号 −1）、report.md §8 の分類、t が ν~ = ν (q/t)^{1/2} だけを通して入るか。
入力：data/refined_kauffman/hopf_g_multi.pkl（scripts/rk_k4_adams.py、|R| ≤ 4、b ≤ 2）。"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
from rcsknots.bps import M_matrix
from rcsknots.partitions import partitions
from rk_b_classify import classify, r, s
from rk_k4_crosscap import nu

nut = sp.Symbol("nut", positive=True)


def hatted(n, tt=None):
    loc = {"r": r, "s": s, "nu": nu}
    g = {k: sp.sympify(v, locals=loc) for k, v in pickle.load(open("data/refined_kauffman/hopf_g_multi.pkl", "rb")).items()}
    parts = partitions(n)
    Minv = M_matrix(n, 1 / r if tt is None else tt).inv()
    hat = Minv * sp.Matrix([g[(R, 1)] for R in parts])
    return {R: sp.cancel(v) for R, v in zip(parts, hat)}


if __name__ == "__main__":
    for n in (2, 3, 4):
        for R, v in hatted(n).items():
            lau = sp.Poly(sp.fraction(v)[1], r, s, nu).is_monomial
            w = sp.cancel(v.subs(nu, nut * s / r))
            print(f"ĝ_{{{R},[1]}}: Laurent={lau}、ν~ を固定して t に依らない={not w.has(s)}、分類={classify(v) if lau else '—'}", flush=True)
            print(f"      = {sp.factor(w)}（ν~ 表示）", flush=True)
