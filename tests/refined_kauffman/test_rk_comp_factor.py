"""composite の大きな N の因子化（rk_f_comp_factor、theory §12.36、続き 106）の機能テスト：
trefoil の [2;1] と [2;2] で、R_{[α;β]} の最上行が R_α R_β の最上行に一致し、差が全行で (q−1)(t−1) で割れる。"""
import sys
import sympy as sp
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import rk_f_comp_factor as CF

r_, s_, nu = CF.r_, CF.s_, CF.nu


def _R():
    w = s_ / r_
    return {c: sp.expand(w**k * CF.R_data(c, "T2_3")) for c, k in (("1", 1), ("2", 2), ("11", 2))}


def test_composite_factorization_trefoil():
    R = _R()
    for name, path, unit, a, b in (("[2;1]", f"{CF.OUT}/Rcomp_A_m1_2.pkl", 1, "2", "1"), ("[2;2]", f"{CF.OUT}/Rcomp22_A_m1_2_2.pkl", 2, "2", "2")):
        Rab = CF.R_comp(path, unit)
        G, ra, rb = CF.G_rows(Rab, R[a], R[b])          # 最上行が比例しなければ assert、割れなければ None
        assert all(g is not None for g in G.values())
        assert max(G) == max(ra) + max(rb) - 2               # G は積の最上行の一つ下から始まる（g_0 = 0）
