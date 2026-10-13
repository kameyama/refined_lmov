"""lmov21_linear.py の一意解（rCS_{[2,1]}、太字）から 3 箱の \\hat f_ρ を作り、genus 表示・spin 表示・既約表示の正値性と特殊化を判定する（予言の検査）。

使い方：PYTHONPATH=src:scripts/framed_knot:scripts/dgr_lmov python scripts/dgr_lmov/lmov21_predict_check.py [data/dgr_lmov/lmovstruct/lmov21_linear.json]
出力：data/dgr_lmov/lmovstruct/lmov21_predictions.json、logs/dgr_lmov/lmov21_predict_check.log。
"""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "framed_knot"))
from dgr_box3 import hat3  # noqa: E402
from lmov21_linear import known_inputs, twist_inputs  # noqa: E402
from t4_colored import bold_to_small  # noqa: E402
from rcsknots.bps import g_expand, to_bold  # noqa: E402
from rcsknots.extract import A_SHIFT  # noqa: E402
from rcsknots.dgr_lmov.lmovstruct import delta_values, is_W_invariant, specialize_lower, spin_count  # noqa: E402


def judge(small_all, log, name):
    hats = hat3(small_all)
    out = {}
    for rho, hf in hats.items():
        hf = {k: v for k, v in hf.items() if v}
        sh = A_SHIFT["nontorus"][3]
        bold = to_bold(hf, sh)
        neg = sum(1 for v in bold.values() if v < 0)
        pos = len(bold) - neg
        sign = 1 if neg <= pos else -1
        bold = {k: sign * v for k, v in bold.items()}
        rec = {"terms": len(bold), "before_neg": min(neg, pos), "W": is_W_invariant(bold)}
        if rec["W"]:
            N = g_expand(bold)
            rec["after_terms"] = len(N)
            rec["after_neg"] = sum(1 for v in N.values() if v < 0)
            rec["irr_terms"], rec["irr_neg"] = spin_count(bold)
        rec["a=t"] = len(specialize_lower(hf, 1))
        rec["a=1"] = len(specialize_lower(hf, 0))
        rec["delta"] = delta_values(bold)[:1] + delta_values(bold)[-1:]
        out[str(rho)] = rec
        log(f"  {name:6s} ρ={str(rho):10s} 前 {rec['terms']:6d} 項（負 {rec['before_neg']}）、W {rec['W']}、後 {rec.get('after_terms')} 項（負 {rec.get('after_neg')}）、既約 {rec.get('irr_terms')} 項（負 {rec.get('irr_neg')}）、a=t:{rec['a=t']} a=1:{rec['a=1']}")
    return out


def main():
    path = Path(sys.argv[1] if len(sys.argv) > 1 else "data/dgr_lmov/lmovstruct/lmov21_linear.json")
    data = json.loads(path.read_text())
    logf = open("logs/dgr_lmov/lmov21_predict_check.log", "a")

    def log(s):
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    log(f"=== lmov21_predict_check {time.strftime('%Y-%m-%d %H:%M:%S')} {path} ===")
    res = {}
    for name, d in data.items():
        sol = d.get("solution_bold")
        if not sol:
            log(f"  {name}: 一意解なし（零空間 {d.get('kernel_dim')}）")
            continue
        bold21 = {tuple(map(int, k.split(","))): v for k, v in sol.items()}
        if name.startswith("K_"):
            small, _, _, _ = twist_inputs(int(name[2:]))
        else:
            small, _, _, _ = known_inputs(name)
        small_all = {**small, (2, 1): bold_to_small(bold21)}
        res[name] = {"rcs21_terms": len(bold21), "rcs21_neg": sum(1 for v in bold21.values() if v < 0), "hat": judge(small_all, log, name)}
    Path("data/dgr_lmov/lmovstruct/lmov21_predictions.json").write_text(json.dumps(res, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
