"""framed knot の refined LMOV 不変量：framing x を掛けた rCS から \\hat f_ρ を作り、G 展開の前後の正整数性を判定する。

framing の因子（ユーザーの以前の計算 Dropbox/Mathematica/marino-vafa の "framed refined Oogri-Vafa ver2" の framing[λ, q, t]）：
  framing[λ] = (−1)^{|λ|} Π_{□∈λ} q^{a(□)+1/2} t^{−l(□)−1/2} = (−(q/t)^{1/2})^{|λ|} q^{n(λ')} t^{−n(λ)}，
  framed の rCS_λ = framing[λ]^x rCS_λ。Adams 作用の写し rCS_λ(a^k, q^k, t^k) には framing[λ](q^k, t^k)^x が付く
  （本研究の reduced の多項式に単項式を掛ければ自動的にそうなる）。KN の (A, Q, T) = (a, q, t)^{1/2} では
  framing[λ]^x = (−1)^{|λ|x} Q^{x(|λ| + 2n(λ'))} T^{−x(|λ| + 2n(λ))}。
判定は第 3 期と同じ（太字の変数、a 次数のずらしは nontorus 規約、G 展開の前と後、全体の符号は良い方。scripts/u4_run.py の judge）。
使い方：python scripts/framing_check.py 結び目 XMIN XMAX [NMAX]
  結び目：unknot、data/reference/notebook_torus.txt の名前（T2_3, 4_1 など）、data/torus/・data/cables/ の名前（T3_4, K1 など）、
          twist:p（ツイスト結び目 K_p の 2 箱まで、GNSSS の予想式。NMAX = 2）
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

from rcsknots.bps import g_expand, to_bold
from rcsknots.bpsnum import Reformulated
from rcsknots.exact import reconstruct_points
from rcsknots.extract import A_OFFSET, A_SHIFT, _arr
from rcsknots.field import V, make_conv
from rcsknots.laurent import load_table
from rcsknots.partitions import conjugate, n_stat, partitions

COLORS = {"1": (1,), "2": (2,), "11": (1, 1), "3": (3,), "111": (1, 1, 1), "21": (2, 1)}


def load_reduced(name, nmax=3):
    """{λ: {(e_A, e_Q, e_T): c}}（reduced、framing 0、KN の (A, Q, T)）。"""
    lams = [lam for k in range(1, nmax + 1) for lam in partitions(k)]
    if name.startswith("twist"):                   # twist:p = ツイスト結び目 K_p の [1], [2], [1,1]（GNSSS の予想式、第 2 期 T3 と同じ入力）
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from t3_pretzel import tc_graded
        from t4_colored import bold_to_small
        p = int(name.split(":")[1])
        P1, P2 = tc_graded(p, 1, 1), tc_graded(p, 1, 2)
        P11 = {(ea, -eq, et - eq): c for (ea, eq, et), c in P2.items()}
        assert nmax <= 2
        return {(1,): bold_to_small(P1), (2,): bold_to_small(P2), (1, 1): bold_to_small(P11)}
    if name == "unknot":
        return {lam: {(0, 0, 0): 1} for lam in lams}
    nb = load_table("data/reference/notebook_torus.txt")
    out = {}
    for c, lam in COLORS.items():
        if sum(lam) > nmax:
            continue
        if (name, c) in nb:
            out[lam] = nb[(name, c)]
            continue
        for d in ("data/torus", "data/cables"):
            p = Path(f"{d}/{name}_{c}.json")
            if p.exists():
                out[lam] = {tuple(map(int, k.split(","))): v for k, v in json.loads(p.read_text())["poly"].items()}
    missing = [lam for lam in lams if lam not in out]
    assert not missing, f"{name}: 入力が無い {missing}"
    return out


def framing_monomial(lam, x):
    """framing[λ]^x の (符号, (e_A, e_Q, e_T))。"""
    k = sum(lam)
    return (-1) ** (k * x), (0, x * (k + 2 * n_stat(conjugate(lam))), -x * (k + 2 * n_stat(lam)))


def framed(reduced, x):
    out = {}
    for lam, poly in reduced.items():
        s, (da, dq, dt) = framing_monomial(lam, x)
        out[lam] = {(ea + da, eq + dq, et + dt): s * c for (ea, eq, et), c in poly.items()}
    return out


def zero_rhos(reduced, n, trials=3):
    """乱数の点（素数 2^31 − 1 を法とする）で \\hat f_ρ が 0 になる ρ。恒等的に 0 の ρ は補間（hat_f_exact）が失敗するので先に除く。"""
    p, rng = 2**31 - 1, np.random.default_rng(11)
    zero = set(partitions(n))
    for _ in range(trials):
        xs, ys = rng.integers(2, 10**6, 4), rng.integers(2, 10**6, 4)
        hat, _ = Reformulated(reduced, V(xs, p), V(ys, p), make_conv(p, 4)).hat_f(n)
        zero = {r for r in zero if not any(np.any(_arr(c, 4, p) != 0) for c in hat[r].terms.values())}
    return zero


def hat_f_nonzero(reduced, n, rhos):
    """rcsknots.extract.hat_f_exact と同じ補間を、指定した ρ だけについて行う（f^t の照合も同じ）。"""
    out = {}
    for rho in rhos:
        def ev(xs, ys, p, rho=rho):
            k = len(xs)
            hat, resid = Reformulated(reduced, V(xs, p), V(ys, p), make_conv(p, k)).hat_f(n)
            for r in resid.values():
                if any(np.any(_arr(c, k, p) != 0) for c in r.terms.values()):
                    raise ArithmeticError("f^t ≠ M(q^{-1}) \\hat f at some point")
            arr = np.zeros((k, 2 * A_OFFSET + 1), dtype=np.int64)
            for e, c in hat[rho].terms.items():
                arr[:, e + A_OFFSET] = (arr[:, e + A_OFFSET] + _arr(c, k, p)) % p
            return arr
        try:
            raw = reconstruct_points(ev, seed=hash(rho) % 1000, log=lambda s: None)
        except RuntimeError:
            out[rho] = None
            continue
        out[rho] = {(d - A_OFFSET, eq, et): c for (d, eq, et), c in raw.items()}
    return out


def judge(hat, n):
    """(G 展開の前の項数, 前の負, 後の項数, 後の負) と G 展開の結果（全体の符号は良い方。閉じなければ後は None）。"""
    if hat is None:
        return "Laurent でない", None
    bold = to_bold(hat, A_SHIFT["nontorus"][n])
    pos, neg = sum(v > 0 for v in bold.values()), sum(v < 0 for v in bold.values())
    sign = 1 if neg <= pos else -1
    try:
        N = g_expand({k: sign * v for k, v in bold.items()})
    except ValueError:
        return f"前 {min(pos, neg)}/{len(bold)} 負、G 展開が閉じない", None
    nneg = sum(v < 0 for v in N.values())
    return f"前 {min(pos, neg)}/{len(bold)} 負、後 {nneg}/{len(N)} 負", N


def run(name, xs, nmax=3, log=print):
    reduced = load_reduced(name, nmax)
    results = {}
    for x in xs:
        R = framed(reduced, x)
        for n in range(1, nmax + 1):
            t0 = time.time()
            zeros = zero_rhos(R, n)
            hats = hat_f_nonzero(R, n, [r for r in partitions(n) if r not in zeros])
            for rho in zeros:
                results[(x, rho)] = ("\\hat f = 0", {}, {})
                log(f"{name} framing {x:+d} ρ={list(rho)}：\\hat f = 0（乱数の 3 点で 0）")
            for rho, hat in hats.items():
                msg, N = judge(hat, n)
                results[(x, rho)] = (msg, N, hat)
                log(f"{name} framing {x:+d} ρ={list(rho)}：{msg}（{time.time() - t0:.1f}s）")
    return results


if __name__ == "__main__":
    name, x0, x1 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    nmax = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    run(name, range(x0, x1 + 1), nmax, log=lambda s: print(s, flush=True))
