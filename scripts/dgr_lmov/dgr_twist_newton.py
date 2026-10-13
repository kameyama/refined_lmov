"""段階 D4（docs/dgr_lmov/research_plan.md）：ツイスト結び目 K_p の商 R_ρ = \\hat f_ρ / D の、全ての p での spin 表示の正値性。

\\hat f_ρ(K_p) は p の符号ごとに 5 つの単項式 m_j の指数和（scripts/framed_knot/v4_closed.py）で、D は p によらないので R_ρ も同じ根の指数和。
ニュートン形 R(|p| = n) = Σ_k c_k h_{n−1−k}(m_{π0..πk}) の係数 c_k が全て同じ符号なら（h は正係数）、|p| ≥ 1 の全てで R_ρ は一様な符号【導出】。
根の並び π は 5! 通りを試す。使い方：PYTHONPATH=src python scripts/dgr_lmov/dgr_twist_newton.py
"""
import sys
from itertools import permutations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "framed_knot"))
from v4_closed import M_ROOTS, RHOS, hat_direct, load_inputs, newton_coeffs, newton_eval  # noqa: E402
from rcsknots.dgr_lmov.dgr import quotient  # noqa: E402


def uniform(cs):
    vals = [v for c in cs for v in c.values()]
    return 1 if all(v > 0 for v in vals) else (-1 if all(v < 0 for v in vals) else 0)


def main():
    cache = load_inputs()
    for sign in (1, -1):
        R = {n: {rho: quotient(b, 2) for rho, b in hat_direct(sign * n, cache).items()} for n in range(1, 11)}
        for rho in RHOS:
            seq = [R[n][rho] for n in range(1, 6)]
            found = None
            for perm in permutations(M_ROOTS):
                order = [M_ROOTS[k] if sign == 1 else tuple(-e for e in M_ROOTS[k]) for k in perm]
                cs = newton_coeffs(seq, order)
                if uniform(cs):
                    found = (perm, uniform(cs), cs, order)
                    break
            if found is None:
                print(f"p の符号 {sign:+d} {rho}: 一様な符号のニュートン形が見つからない", flush=True)
                continue
            perm, sg, cs, order = found
            ok = all(newton_eval(cs, order, n - 1) == R[n][rho] for n in range(6, 11))
            print(f"p の符号 {sign:+d} {rho}: 並び {perm}、c_k の符号 {sg:+d}、項数 {[len(c) for c in cs]}、|p| = 6..10 の直接計算と一致 {ok}", flush=True)


if __name__ == "__main__":
    main()
