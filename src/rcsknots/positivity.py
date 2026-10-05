"""Positivity of the refined BPS (refined LMOV) invariants before and after the G-expansion.

before G : the coefficients of \\hat f_ρ (after the a-grading shift) in the bold variables,
           i.e. KN's naive substitution a = −𝐚²𝐭, q = 𝐪²𝐭², t = 𝐪² (with √q = −𝐪𝐭).
after G  : \\hat N_{ρ,g,β,J_r} from the G-expansion.
The overall sign is a convention (a^{±1/2} = (−𝐚²𝐭)^{±1/2}); we report the minimal number of
negative coefficients over the overall sign, and whether one common sign works for both forms.
"""
import json
from pathlib import Path

from .bps import to_bold


def counts(values):
    vals = list(values)
    neg = sum(1 for v in vals if v < 0)
    return neg, len(vals) - neg, len(vals)


def judge(path):
    d = json.loads(Path(path).read_text())
    hf = {tuple(map(int, x.split(","))): v for x, v in d["hat_f"].items()}
    bold = to_bold(hf, d["a_shift_A_units"])
    N = {tuple(map(int, x.split(","))): v * d["overall_sign_applied"] for x, v in d["N"].items()}   # raw sign
    bn, bp, bt = counts(bold.values())
    an, ap, at = counts(N.values())
    common = (bn == 0 and an == 0) or (bp == 0 and ap == 0)
    return {"before_G": {"terms": bt, "min_negative": min(bn, bp)},
            "after_G": {"terms": at, "min_negative": min(an, ap)},
            "both_positive_same_sign": common,
            "integral": all(isinstance(v, int) for v in list(bold.values()) + list(N.values()))}
