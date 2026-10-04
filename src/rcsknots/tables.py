"""Save and format refined BPS tables (\\hat N_{ρ,g,β,J_r}) in the layout of 1703.05408 App. C."""
import json
from pathlib import Path


def normalize_sign(raw):
    """Overall sign σ such that the lowest term (g=0, minimal β, then minimal F) is positive.

    The overall sign of \\hat f after the a-grading shift a^{±1/2} = (−𝐚²𝐭)^{±1/2} is a convention.
    Returns (σ, {(g, β, F): σ·c}).
    """
    lowest = min(raw, key=lambda k: (k[0], k[1], k[2]))
    sigma = 1 if raw[lowest] > 0 else -1
    return sigma, {k: sigma * v for k, v in raw.items()}


def positivity(N):
    """Summary of signs and of the 2J_r parity per (g, β)."""
    neg = sorted((k, v) for k, v in N.items() if v < 0)
    mixed = []
    by_gb = {}
    for (g, b, F) in N:
        by_gb.setdefault((g, b), set()).add(F % 2)
    for gb, parities in sorted(by_gb.items()):
        if len(parities) > 1:
            mixed.append(gb)
    return {"n_terms": len(N), "n_negative": len(neg), "negative": [[*k, v] for k, v in neg],
            "all_nonnegative": not neg, "mixed_parity_g_beta": [list(x) for x in mixed]}


def format_table(N, title):
    """Text table: rows g, columns grouped by β, sub-columns 2J_r (as in KN App. C)."""
    betas = sorted({b for _, b, _ in N})
    gs = range(0, max(g for g, _, _ in N) + 1)
    cols = []
    for b in betas:
        Fs = [F for _, bb, F in N if bb == b]
        cols += [(b, F) for F in range(min(Fs), max(Fs) + 1)]
    w = max(3, max(len(str(abs(v))) + 1 for v in N.values()) + 1)
    lines = [title, "β    " + "".join(f"{b:>{w}}" for b, _ in cols), "2J_r " + "".join(f"{F:>{w}}" for _, F in cols)]
    for g in gs:
        lines.append(f"g={g:<3}" + "".join(f"{N.get((g, b, F), 0):>{w}}" for b, F in cols))
    return "\n".join(lines)


def save(path_stem, rho, hat_f, N, meta):
    path = Path(path_stem)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = dict(meta)
    data.update({"rho": list(rho),
                 "hat_f": {",".join(map(str, k)): v for k, v in sorted(hat_f.items())},
                 "N": {",".join(map(str, k)): v for k, v in sorted(N.items())},
                 "positivity": positivity(N)})
    path.with_suffix(".json").write_text(json.dumps(data, indent=1))
    path.with_suffix(".txt").write_text(format_table(N, f"{meta.get('knot')}  rho={list(rho)}  \\hat N_(rho,g,beta,J_r)") + "\n")
