"""Positivity and integrality of the refined LMOV invariants (1703.05408 §4), judged from lowercase data.

All quantities f^q_μ, f^{\\bar t}_μ, \\hat f_ρ are computed in KN's lowercase variables (a,q,t) and are
changed to the bold variables only here, immediately before reading off the invariants:
  1. a-grading shift for odd |μ| (A = a^{1/2}; nontorus convention: |μ|=1 → a^{−1/2}, |μ|=3 → a^{+1/2});
  2. KN (cov): a = −𝐚²𝐭, q^{1/2} = −𝐪𝐭, t^{1/2} = 𝐪.
Invariants:
  N^q_μ, N^{\\bar t}_μ : coefficients of f^q_μ, f^{\\bar t}_μ after 1–2   (KN: non-negative for μ = [r] up to ±𝐪^•)
  \\hat N_{μ;i,j,k}    : coefficients of \\hat f_ρ after 1–2 (before the G-expansion)
  \\hat N_{ρ,g,β,J_r}  : after the G-expansion.
The overall sign is a convention; we report the minimal number of negative coefficients over it.
"""
import json
from pathlib import Path

from .bps import to_bold, g_expand
from .extract import A_SHIFT

RHOS = [(1,), (2,), (1, 1), (3,), (2, 1), (1, 1, 1)]


def load_lmov(name):
    out = {}
    for n in (1, 2, 3):
        d = json.loads(Path(f"data/lmov/{name}_{n}.json").read_text())
        for key, poly in d["f"].items():
            kind, rho = key.split("_")
            out[(kind, tuple(int(c) for c in rho))] = {tuple(map(int, k.split(","))): v for k, v in poly.items()}
    return out


def sign_count(values):
    vals = list(values)
    neg = sum(1 for v in vals if v < 0)
    return {"terms": len(vals), "min_negative": min(neg, len(vals) - neg),
            "integral": all(isinstance(v, int) for v in vals), "sign": 1 if neg <= len(vals) - neg else -1}


def judge(name, convention="nontorus"):
    F = load_lmov(name)
    res = {}
    for rho in RHOS:
        sh = A_SHIFT[convention][sum(rho)]
        r = {}
        for kind in ("fq", "ft", "hat"):
            r[kind] = sign_count(to_bold(F[(kind, rho)], sh).values())
        N = g_expand(to_bold(F[("hat", rho)], sh))
        r["hatN"] = sign_count(N.values())
        parity = {}
        for (g, b, Fj) in N:
            parity.setdefault((g, b), set()).add(Fj % 2)
        r["hatN"]["parity_aligned"] = all(len(s) == 1 for s in parity.values())
        r["hat_and_hatN_same_sign"] = (r["hat"]["min_negative"] == 0 and r["hatN"]["min_negative"] == 0
                                      and r["hat"]["sign"] == r["hatN"]["sign"])
        res[rho] = r
    return res
