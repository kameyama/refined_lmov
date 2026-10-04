"""正規化 (N2') を data/cables/*.json に適用する（docs/conventions.md §5）。

(N2'): 彩色 λ の単項式は M_λ = (M_□)^{|λ|} · T_λ^{F}、T_λ = q^{n(λ')} t^{−n(λ)}、F = p₂q₂（一番外側のケーブル、Seifert framing）。
M_□ は □ について (N1)+(N2) で決めたもの（"monomial_applied"）。
旧データ（各 λ に (N1)+(N2) を課したもの）との差は (q/t)^δ でなければならない（(N1) は両者で同じ）。差があれば書き換える。
"""
import json
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage5_cables import KNOTS, COLORS           # noqa: E402
from rcsknots.partitions import n_stat, conjugate   # noqa: E402


def predicted(m1, lam, F):
    k = sum(lam)
    return {"sign": m1["sign"] ** k, "a": Fraction(m1["a"]) * k,
            "q": Fraction(m1["q"]) * k + F * n_stat(conjugate(lam)), "t": Fraction(m1["t"]) * k - F * n_stat(lam)}


def main(names):
    for name in names:
        p1, q1, p2, q2 = KNOTS[name]
        F = p2 * q2
        m1 = json.loads(Path(f"data/cables/{name}_1.json").read_text())["monomial_applied"]
        for cname, lam in COLORS.items():
            path = Path(f"data/cables/{name}_{cname}.json")
            if not path.exists():
                continue
            d = json.loads(path.read_text())
            old = d.get("monomial_N1N2", d["monomial_applied"])
            new = predicted(m1, lam, F)
            dq, dt = new["q"] - Fraction(old["q"]), new["t"] - Fraction(old["t"])
            assert new["sign"] == old["sign"] and new["a"] == Fraction(old["a"]) and dq + dt == 0, (name, cname, old, new)
            shift = 2 * dq                      # in units of Q (and −T)
            assert shift.denominator == 1
            shift = int(shift)
            cur = Fraction(d["monomial_applied"]["q"]) - Fraction(old["q"])     # already applied correction
            extra = shift - int(2 * cur)
            if extra:
                d["poly"] = {",".join(map(str, (ea, eq + extra, et - extra))): v
                             for (ea, eq, et), v in ((tuple(map(int, k.split(","))), v) for k, v in d["poly"].items())}
                d["poly"] = dict(sorted(d["poly"].items()))
            d["monomial_N1N2"] = old
            d["monomial_applied"] = {"sign": new["sign"], "a": str(new["a"]), "q": str(new["q"]), "t": str(new["t"])}
            d["framing"] = "0; normalization (N1) + (N2') of docs/conventions.md §5 (hypothesis for cables)"
            d["N2_differs_by_q_over_t_power"] = str(dq)
            path.write_text(json.dumps(d, indent=1))
            print(name, cname, "change (q/t)^", dq)


if __name__ == "__main__":
    main(sys.argv[1:] or list(KNOTS))
