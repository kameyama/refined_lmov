"""DGR / GGS の微分 d_{n|m} の類似を refined LMOV 不変量（spin 表示・genus 表示）で調べる観察の一括実行（docs/dgr_lmov.md）。

出力：logs/dgr_lmov.log（人が読む記録）、data/dgr_lmov/summary.json（判定の表）。
使い方：PYTHONPATH=src python scripts/dgr_lmov_survey.py
"""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from framing_check import load_reduced  # noqa: E402
from framing_twist import hat_framed  # noqa: E402
from v4_closed import inputs, load_inputs, mir  # noqa: E402
from rcsknots.bps import g_expand, to_bold  # noqa: E402
from rcsknots.dgr import (divides, d11_top_genus, hook_rule, hook_shift, judge_quotient, quotient,  # noqa: E402
                          sigma, specialize, u)
from rcsknots.extract import A_SHIFT  # noqa: E402
from rcsknots.knot_table import load_homology, superpoly  # noqa: E402
from rcsknots.lmov import RHOS, load_lmov  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
LMOV_KNOTS = ["T2_3", "T2_5", "T2_11", "T3_4", "T3_5", "T3_7", "T4_5", "4_1",
              "K1", "K3", "K6a", "K7", "K8", "K9", "K2", "K4", "K5a", "K5b", "K5c", "K6b"]
PAIRS = [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2)]
UNKNOT = {(0, 0, 0): 1}


def name(rho):
    return "".join(map(str, rho))


def analyse(hat, log, label):
    """{ρ: 太字の \\hat f_ρ} を判定。"""
    out = {}
    for rho, b in hat.items():
        n = sum(rho)
        r = {"divides": {f"{a}|{c}": divides(b, a, c) for a, c in PAIRS}}
        sn, ln, gn, lg, same = judge_quotient(b, n)
        r.update(R_spin_neg=sn, R_spin_terms=ln, R_genus_neg=gn, R_genus_terms=lg, same_sign=same)
        R = quotient(b, n)
        r["R_genus_max"] = max(int(k[0]) for k in g_expand(R))
        r["R_d11_top_genus"] = d11_top_genus(R)
        out[name(rho)] = r
        divs = "".join("1" if v else "0" for v in r["divides"].values())
        log(f"{label:>8} {name(rho):>3}: 因子 {divs}  R spin 負 {sn}/{ln}  genus 負 {gn}/{lg}  符号一致 {same}  "
            f"genus 最大 {r['R_genus_max']}  d11 後 {r['R_d11_top_genus']}")
    return out


def survey_lmov(log):
    res = {}
    for K in LMOV_KNOTS:
        t0 = time.time()
        F = load_lmov(K)
        res[K] = analyse({rho: to_bold(F[("hat", rho)], A_SHIFT["nontorus"][sum(rho)]) for rho in RHOS}, log, K)
        log(f"         ({time.time() - t0:.0f} 秒)")
    return res


def survey_twist(log, cache):
    res = {}
    for p in [v for v in range(-10, 11) if v]:
        res[str(p)] = analyse(hat_framed(*inputs(p, cache), 0), log, f"K_{p}")
    return res


def survey_thick(log):
    H = load_homology()
    res = {}
    for K in ("10_124", "10_139", "10_152"):
        d = json.loads((ROOT / f"data/b4/{K}_s2_tc.json").read_text())
        P2 = {tuple(map(int, k.split(","))): v for k, v in d["[2]"].items()}
        res[K] = analyse(hat_framed(superpoly(H[K]), P2, 0), log, K)
    return res


def survey_framed(log, cache):
    """\\hat f^{(x)}(K) ≡ c^{|ρ|} \\hat f^{(x)}(○) mod (1 + u_{±1})（c は [1] の d_{1|0} での値）と、差を D で割った商の正値性。"""
    res = {}
    for p in [v for v in range(-8, 9) if v]:
        P1, P2 = inputs(p, cache)
        c = specialize(P1, 1, 0)
        ((cq, ct), cs), = c.items()
        assert cq == 0 and specialize(P1, 0, 1) == c
        for x in range(-6, 7):
            H, U = hat_framed(P1, P2, x), hat_framed(UNKNOT, UNKNOT, x)
            row = {}
            for rho in H:
                E = {(k[0], k[1], k[2] + 2 * ct): cs ** 2 * v for k, v in U[rho].items()}
                diff = dict(H[rho])
                for k, v in E.items():
                    diff[k] = diff.get(k, 0) - v
                diff = {k: v for k, v in diff.items() if v}
                if not diff:
                    row[name(rho)] = "0"
                    continue
                sn, _, gn, _, _ = judge_quotient(diff, 2)
                row[name(rho)] = [sn, gn]
            res[f"{p},{x}"] = row
        log(f"framed K_{p}: c = {cs}𝐭^{ct}  " + " ".join(f"x={x}:{res[f'{p},{x}']}" for x in range(-6, 7)))
    return res


def survey_hook(log, cache):
    """hook λ で rCS_λ|_{d11} = c · A_K(𝐪^{|λ|}𝐭^{arm}, 𝐭)、A_K = rCS_{[1]}|_{d11}。"""
    res = {}
    for K in ("T2_3", "T2_5", "T3_4", "T3_5", "4_1", "K1", "K2", "K4"):
        R = load_reduced(K, 3)
        bold = {lam: to_bold(P, 0 if all(k[0] % 2 == 0 for k in P) else 1) for lam, P in R.items()}
        A1 = specialize(bold[(1,)], 1, 1)
        res[K] = {name(lam): hook_rule(A1, specialize(b, 1, 1), lam) for lam, b in bold.items() if hook_shift(lam)}
        log(f"hook {K}: {res[K]}")
    for p in [v for v in range(-6, 7) if v]:
        P1, P2 = inputs(p, cache)
        A1 = specialize(P1, 1, 1)
        res[f"K_{p}"] = {nm: hook_rule(A1, specialize(P, 1, 1), lam) for nm, P, lam in (("2", P2, (2,)), ("11", mir(P2), (1, 1)))}
        log(f"hook K_{p}: {res[f'K_{p}']}")
    return {k: {kk: (str(v) if v else None) for kk, v in d.items()} for k, d in res.items()}


def survey_orbits(log):
    """σ 軌道 {d_{n|m}, d_{m|n}} の因子 (1 + u)(1 + σu) の genus 表示。"""
    res = {}
    for n in range(4):
        for m in range(n, 4):
            e = u(n, m)
            s = next(iter(sigma({e: 1})))
            prod = {(0, 0, 0): 1, e: 1} if e == s else {(0, 0, 0): 1, e: 1, s: 1, tuple(e[i] + s[i] for i in range(3)): 1}
            N = g_expand(prod)
            res[f"{n}|{m}"] = {"genus_max": max(k[0] for k in N), "negative": sum(1 for v in N.values() if v < 0)}
            log(f"σ 軌道 d_{{{n}|{m}}}: genus 表示 {sorted(N.items())}")
    return res


def main():
    (ROOT / "logs").mkdir(exist_ok=True)
    (ROOT / "data/dgr_lmov").mkdir(parents=True, exist_ok=True)
    logf = open(ROOT / "logs/dgr_lmov.log", "w")

    def log(s):
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    t0 = time.time()
    cache = load_inputs()
    summary = {"orbits": survey_orbits(log), "hook": survey_hook(log, cache), "twist": survey_twist(log, cache),
               "thick": survey_thick(log), "framed": survey_framed(log, cache), "lmov": survey_lmov(log)}
    (ROOT / "data/dgr_lmov/summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=1))
    log(f"完了（{time.time() - t0:.0f} 秒）")


if __name__ == "__main__":
    main()
