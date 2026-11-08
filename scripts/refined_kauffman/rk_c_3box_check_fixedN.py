r"""大きな N の 3 箱の g を U = t², t³ で固定 N の g（T2_{2m+1}_fixedN_3box_g.pkl、rk_k5_torus_fixedN_hatted.py）と全係数で照合する
（progress.md 続き 139。組み立ての後に固定 N の g を作った場合に使う）。
保存した式の文字列は m = 4 で sympy が読めない大きさになるので、組み立て（rk_c_3box_largeN、数十秒）をその場でやり直し、
flint のまま u = s^{2N} を代入して比べる（指数の重なる項は足す）。
使い方：python rk_c_3box_check_fixedN.py m"""
import sys, pickle
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
from fractions import Fraction
from rcsknots.characters import chi
from rcsknots.partitions import partitions
from rcsknots.refined_kauffman import flint_rf as FR
import rk_c_3box_largeN as T3


def assemble(m):
    T3.MK = m
    F = T3.plog(T3.so_Z(T3.so_W(), T3.K_matrix()))
    Wo, pair = T3.cover_parts()
    Fc = T3.plog(T3.cover_Z(Wo, pair))
    return {k: F.get(k, T3.ZERO) - T3.conv(Fraction(1, 2)) * Fc.get(k, T3.ZERO) for k in set(F) | set(Fc)}, F


def subst(x, N):
    def cv(p):
        out = {}
        for e, c in p.to_dict().items():
            k = (e[0], e[1] + 2 * N * e[2])
            out[k] = out.get(k, 0) + int(c)
        return FR.CTX.from_dict({k: c for k, c in out.items() if c})
    return FR.RF(cv(x.n), cv(x.d))


if __name__ == "__main__":
    m = int(sys.argv[1])
    G, _ = assemble(m)
    fx = pickle.load(open(f"data/refined_kauffman/T2_{2 * m + 1}_fixedN_3box_g.pkl", "rb"))
    for N in (2, 3):
        ok = {}
        for R_ in partitions(3):
            gR = sum((G[rho] * T3.conv(chi(R_, rho)) for rho in partitions(3)), T3.ZERO)
            ok[R_] = (subst(gR, N) - FR.from_sympy(fx[N]["g_schur"][str(R_)])).is_zero()
        print(f"T(2,{2 * m + 1})、U = t^{N}：g_R（|R| = 3）が固定 N の g と一致 {ok}", flush=True)
