"""G_adj の行の Ĩ（rk_f_n2_EG_adj、theory §12.44、続き 114）の機能テスト：n = 1 の trefoil の保存した E を §12.35 の形に分けると
Ĝ₁ = −t²（= −ψ²(ρ₁)：g₁ = (2ρ₀ρ₁ − ρ₁²)/q の交差項が入れ替えで消え、ρ₁⊗ρ₁ が残る）、Ĝ₂ = −t²(q−1)(t+1)、D₂ = t³（§12.35 の D_{δ+1} = t^{3δ}）。"""
import sys
import sympy as sp
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import rk_f_n2_EG_adj as EG


def test_trefoil_n1_G_hat():
    q, t = EG.q, EG.t
    Gh, g = EG.part_A(1)
    assert sp.expand(Gh[1] + t**2) == 0
    assert sp.expand(Gh[2] + t**2 * (q - 1) * (t + 1)) == 0
    assert sp.expand(g[1] + t * (2 * q * t + t + 2) / q) == 0
    rule2 = sp.expand(g[2].subs(q, -q) * q**2)
    assert sp.expand((Gh[2] - rule2) / 2 - t**3) == 0
