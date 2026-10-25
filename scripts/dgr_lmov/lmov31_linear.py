"""4_1 の [3,1]（4 箱）の線形解析：scripts/framed_knot/f41_box4.py の一次の系（t = q、転置、両手型、EGP、U(2)）と、4 箱の \\hat f_ρ の Laurent 性を
**因子ごとの代入**で一次式にしたもの（補助変数なし。未知数は rCS_{[3,1]} の係数 2367 個だけ）に、正値性を課さず零空間の次元を法 p で数える。
さらに hook 性質（a = t/q で rCS_{[3,1]} = c³·A(𝐪⁴𝐭², 𝐭)）を加え、一意になれば厳密解を求めて検算する。

Laurent 性：D·\\hat f_ρ = Σ_λ (Dα unk)_{λρ}·x_λ + D·K_ρ =: N_ρ(x)（λ = [3,1], [2,1,1]、x_{[2,1,1]} は x_{[3,1]} の転置）が D の各既約因子で割り切れること。
D = (Q∓1)²(Q²+1)(T∓1)²(T²+1)(QT∓1)²(QT³∓1)(Q²T²+1)(Q³T∓1)（f41_box4.alpha の D）。重複度 2 の因子は横断方向の微分も課す。
使い方：PYTHONPATH=src:scripts/framed_knot:scripts/dgr_lmov python scripts/dgr_lmov/lmov31_linear.py。出力：data/dgr_lmov/lmovstruct/lmov31_linear.json、logs/dgr_lmov/lmov31_linear.log。
"""
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "framed_knot"))
import f41_box4 as F  # noqa: E402
from lmov21_linear import rank_nullity, solve_unique  # noqa: E402
from rcsknots.bps import to_bold  # noqa: E402
from rcsknots.dgr_lmov.dgr import specialize  # noqa: E402

I4 = {0: (1, 0), 1: (0, 1), 2: (-1, 0), 3: (0, -1)}            # i^k の (実, 虚) 成分


def images4(e):
    """小文字の単項式 (ea, eq, et) の、D の各因子での像 {(因子名, key): 係数}。"""
    ea, eq, et = e
    out = {}
    # (Q ∓ 1)²、(T ∓ 1)²：値と微分
    out[("Q=1", (ea, et))] = 1
    out[("dQ=1", (ea, et))] = eq
    out[("Q=-1", (ea, et))] = (-1) ** (eq % 2)
    out[("dQ=-1", (ea, et))] = eq * (-1) ** ((eq - 1) % 2)
    out[("T=1", (ea, eq))] = 1
    out[("dT=1", (ea, eq))] = et
    out[("T=-1", (ea, eq))] = (-1) ** (et % 2)
    out[("dT=-1", (ea, eq))] = et * (-1) ** ((et - 1) % 2)
    # Q² + 1、T² + 1：Q = i
    for comp, v in enumerate(I4[eq % 4]):
        if v:
            out[("Q^2+1", (ea, et, comp))] = v
    for comp, v in enumerate(I4[et % 4]):
        if v:
            out[("T^2+1", (ea, eq, comp))] = v
    # (QT ∓ 1)²：T = ±1/Q の値と ∂_T
    out[("QT-1", (ea, eq - et))] = 1
    out[("dQT-1", (ea, eq - et + 1))] = et
    out[("QT+1", (ea, eq - et))] = (-1) ** (et % 2)
    out[("dQT+1", (ea, eq - et + 1))] = et * (-1) ** ((et - 1) % 2)
    # QT³ ∓ 1：Q = ±T^{−3}
    out[("QT3-1", (ea, et - 3 * eq))] = 1
    out[("QT3+1", (ea, et - 3 * eq))] = (-1) ** (eq % 2)
    # Q³T ∓ 1：T = ±Q^{−3}
    out[("Q3T-1", (ea, eq - 3 * et))] = 1
    out[("Q3T+1", (ea, eq - 3 * et))] = (-1) ** (et % 2)
    # Q²T² + 1：T = i/Q
    for comp, v in enumerate(I4[et % 4]):
        if v:
            out[("Q2T2+1", (ea, eq - et, comp))] = v
    return out


def laurent_rows4(xs, al, dk):
    """{ρ ごと} の一次式 [(row {列: 係数}, rhs)]。N_ρ(x) の各因子での像 = 0。"""
    rows = []
    for rho in F.partitions(F.N):
        eqs = defaultdict(dict)
        rhs = defaultdict(int)
        for n, e in enumerate(xs):
            for lam, ee in ((F.UNK[0], e), (F.UNK[1], (e[0], -e[2], -e[1]))):
                for k, c in al[(lam, rho)].items():
                    kk = (ee[0] + k[0], ee[1] + k[1], ee[2] + k[2])
                    for key, v in images4(kk).items():
                        eqs[key][n] = eqs[key].get(n, 0) + c * v
        for k, c in dk[rho].items():
            for key, v in images4(k).items():
                rhs[key] -= c * v
                eqs.setdefault(key, {})
        for key, row in eqs.items():
            row = {c: v for c, v in row.items() if v}
            if row or rhs[key]:
                rows.append((row, rhs[key]))
    return rows


def eq_rows(red, h31, al, Dp, dk, nx, log):
    """f41_box4.build の (1) t = q、(3) 両手型、(4)(5) 特殊化の行だけ（h の列を含む行を捨てる）。"""
    xs, Amat, b, lb, ub = F.build(red, h31, al, Dp, dk, sigma=1, positive=False, margin=0, log=log, with_G=False)
    A = Amat.tocsr()
    rows = []
    for i in range(A.shape[0]):
        s, e = A.indptr[i], A.indptr[i + 1]
        cols = A.indices[s:e]
        if len(cols) and cols.max() >= nx:
            continue
        row = {int(c): int(round(v)) for c, v in zip(cols, A.data[s:e]) if round(v)}
        rows.append((row, int(round(b[i]))))
    return xs, rows


def hook_rows(xs, r1_small):
    """a = t/q（A → T/Q）での rCS_{[3,1]} = c³·A(𝐪⁴𝐭², 𝐭)。太字で作って小文字 (Q, T) に直す：𝐪 = T、𝐭 = −Q/T。"""
    r1 = to_bold(r1_small, 1 if any(k[0] % 2 for k in r1_small) else 0)
    (cq, ct), cs = next(iter(specialize(r1, 1, 0).items()))
    assert cq == 0
    Ab = specialize(r1, 1, 1)
    target_bold = defaultdict(int)
    for (eq, et), v in Ab.items():
        target_bold[(4 * eq, et + 2 * eq + 3 * ct)] += cs ** 3 * v
    target = defaultdict(int)
    for (eq, et), v in target_bold.items():                # 𝐪^{eq}𝐭^{et} → (−1)^{et} Q^{et} T^{eq − et}
        target[(et, eq - et)] += v * (-1) ** (et % 2)
    img = defaultdict(dict)
    for n, (ea, eq, et) in enumerate(xs):
        img[(eq - ea, et + ea)][n] = img[(eq - ea, et + ea)].get(n, 0) + 1
    return [(dict(img.get(key, {})), target.get(key, 0)) for key in set(img) | set(target)]


def main():
    Path("data/dgr_lmov/lmovstruct").mkdir(parents=True, exist_ok=True)
    lg = open("logs/dgr_lmov/lmov31_linear.log", "a")

    def log(s):
        print(s, flush=True)
        lg.write(s + "\n")
        lg.flush()

    log(f"=== lmov31_linear（因子ごとの代入）{time.strftime('%Y-%m-%d %H:%M:%S')} ===")
    t0 = time.time()
    red, h31 = F.inputs(log)
    al, Dp, dstr = F.alpha()
    dk = F.dK(red, Dp, log)
    xs0 = F.window([F.mul(red[(1, 1)], red[(1,)], red[(1,)]), F.mul(red[(3,)], red[(1,)]), F.mul(red[(2,)], red[(1, 1)])], h31, 0)
    xs, rows_eq = eq_rows(red, h31, al, Dp, dk, len(xs0), log)
    assert xs == xs0
    nx = len(xs)
    log(f"  未知数 {nx}、特殊化などの式 {len(rows_eq)}（{time.time() - t0:.0f}s）")
    res = {"x": nx}
    rows_L = laurent_rows4(xs, al, dk)
    log(f"  Laurent 性の式 {len(rows_L)}（{time.time() - t0:.0f}s）")
    rows_H = hook_rows(xs, red[(1,)])
    steps = [("eq", rows_eq), ("eq+L", rows_eq + rows_L), ("eq+H", rows_eq + rows_H), ("eq+L+H", rows_eq + rows_L + rows_H)]
    for name, rows in steps:
        rank, feas = rank_nullity(rows, nx)
        res[name] = {"rows": len(rows), "rank": rank, "nullity": nx - rank, "feasible": feas}
        log(f"  {name:8s} 式 {len(rows):6d} 階数 {rank:5d} 零空間 {nx - rank:5d} 実行可能 {feas}（{time.time() - t0:.0f}s）")
        if nx - rank == 0 and feas:
            x = solve_unique(rows, nx)
            if x is not None:
                sol = {",".join(map(str, xs[c])): v for c, v in x.items() if v}
                res[name + "_solution"] = sol
                log(f"  一意解 rCS_[3,1]：{len(sol)} 項、負 {sum(1 for v in sol.values() if v < 0)}")
                Path("data/dgr_lmov/lmovstruct/lmov31_41_solution.json").write_text(json.dumps(sol))
            else:
                log("  一意解の整数での検算に失敗")
    Path("data/dgr_lmov/lmovstruct/lmov31_linear.json").write_text(json.dumps(res, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
