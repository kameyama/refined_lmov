"""[2,1] の構成：一次の条件（U, T, E, U2, [H], L と 𝐭 の窓の規則）の解空間 x = x₀ + K y（K は整数の零空間の基底）の上で、
\\hat f_ρ（ρ ⊢ 3）の正値性（spin 表示 = G 展開前の係数、genus 表示、既約表示）を y の一次不等式にし、多面体が 1 点かを LP で調べる。
1 点なら、それが構成の答え（既知の結び目では正解と照合、ツイスト結び目では予言）。

使い方：PYTHONPATH=src:scripts python scripts/lmov21_polytope.py [--margin=0] [--nohook] T2_5 tw:2 ...
出力：data/lmovstruct/lmov21_polytope.json、logs/lmov21_polytope.log。
"""
import json
import sys
import time
from fractions import Fraction
from math import gcd
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lmov21_linear as LL  # noqa: E402
from dgr_box3 import KEYS, build21, formula3, inputs_R  # noqa: E402
from t4_colored import bold_to_small  # noqa: E402
from v4_closed import L  # noqa: E402
from rcsknots.bps import g_expand, to_bold  # noqa: E402
from rcsknots.lmovstruct import rows as w_rows  # noqa: E402

PRIMES = [2 ** 31 - 1, 2 ** 31 - 19]


def ratrec(a, m):
    """有理数の復元（|分子|, 分母 < sqrt(m/2)）。失敗なら None。"""
    a %= m
    bound = int((m // 2) ** 0.5)
    r0, r1 = m, a
    s0, s1 = 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if abs(s1) > bound or s1 == 0:
        return None
    return Fraction(r1, s1) if s1 > 0 else Fraction(-r1, -s1)


def rref_mod(rows, ncols, p):
    pivots = []
    for row, b in rows:
        r = {c: v % p for c, v in row.items() if v % p}
        rb = b % p
        for pc, prow, pb in pivots:
            f = r.get(pc)
            if f:
                for c, v in prow.items():
                    nv = (r.get(c, 0) - f * v) % p
                    if nv:
                        r[c] = nv
                    else:
                        r.pop(c, None)
                rb = (rb - f * pb) % p
        if not r:
            assert rb == 0, "実行不能"
            continue
        pc = min(r)
        inv = pow(r[pc], p - 2, p)
        pivots.append((pc, {c: v * inv % p for c, v in r.items()}, rb * inv % p))
    for i in range(len(pivots) - 1, -1, -1):
        pc, prow, pb = pivots[i]
        for j in range(i):
            qc, qrow, qb = pivots[j]
            f = qrow.get(pc)
            if f:
                for c, v in prow.items():
                    nv = (qrow.get(c, 0) - f * v) % p
                    if nv:
                        qrow[c] = nv
                    else:
                        qrow.pop(c, None)
                pivots[j] = (qc, qrow, (qb - f * pb) % p)
    piv = {pc: (prow, pb) for pc, prow, pb in pivots}
    free = [c for c in range(ncols) if c not in piv]
    x0 = {pc: pb for pc, (prow, pb) in piv.items()}
    basis = []
    for fc in free:
        vec = {fc: 1}
        for pc, (prow, pb) in piv.items():
            v = prow.get(fc)
            if v:
                vec[pc] = (-v) % p
        basis.append(vec)
    return x0, basis


def exact_affine(rows, ncols):
    """特解 x₀（有理数）と零空間の整数基底。2 つの素数で復元し、整数で検算する。"""
    x0s, bases = zip(*(rref_mod(rows, ncols, p) for p in PRIMES))
    assert len(bases[0]) == len(bases[1]) and set(x0s[0]) == set(x0s[1])
    m = PRIMES[0] * PRIMES[1]

    def crt(a, b):
        p0, p1 = PRIMES
        inv = pow(p0, -1, p1)
        return (a + p0 * ((b - a) * inv % p1)) % m

    def rec(d0, d1):
        out = {}
        for c in set(d0) | set(d1):
            v = ratrec(crt(d0.get(c, 0), d1.get(c, 0)), m)
            if v is None:
                raise ValueError("有理数の復元に失敗")
            if v:
                out[c] = v
        return out
    x0 = rec(x0s[0], x0s[1])
    basis = []
    for v0, v1 in zip(bases[0], bases[1]):
        v = rec(v0, v1)
        den = 1
        for f in v.values():
            den = den * f.denominator // gcd(den, f.denominator)
        vi = {c: int(f * den) for c, f in v.items()}
        g = 0
        for z in vi.values():
            g = gcd(g, abs(z))
        basis.append({c: z // g for c, z in vi.items()})
    # 検算
    for row, b in rows:
        assert sum(Fraction(v) * x0.get(c, 0) for c, v in row.items()) == b, "特解の検算に失敗"
        for vec in basis:
            assert sum(v * vec.get(c, 0) for c, v in row.items()) == 0, "零空間の検算に失敗"
    return x0, basis


def hat_pieces(small, mons, lin_col, x0, basis):
    """\\hat f_ρ = h0_ρ + Σ_j y_j hj_ρ（小文字、Laurent）。h0 は有理数係数（分母 d）、hj は整数。"""
    inv = {v: k for k, v in lin_col.items()}
    d = 1
    for f in x0.values():
        d = d * f.denominator // gcd(d, f.denominator)
    x0i = {inv[c][1]: int(f * d) for c, f in x0.items()}
    R = inputs_R({**small, (2, 1): {(0, 0, 0): 0}})
    out = {}
    for rho, (ns, den) in formula3().items():
        const = None
        for k, n in ns.items():
            if KEYS[k] == "r21":
                continue
            term = n * R[KEYS[k]]
            const = term if const is None else const + term
        num0 = const * L.from_dict({(0, 0, 0): d}) + ns[1] * L.from_dict(bold_to_small(x0i))
        h0 = {k: Fraction(v, d) for k, v in (num0 / den).to_dict().items()}
        hs = []
        for vec in basis:
            vb = {inv[c][1]: z for c, z in vec.items()}
            hs.append((ns[1] * L.from_dict(bold_to_small(vb)) / den).to_dict())
        out[rho] = (h0, hs)
    return out


def bold_frac(P, shift):
    """有理数係数の小文字 → 太字（to_bold は整数前提なので分子で呼ぶ）。"""
    d = 1
    for f in P.values():
        d = d * Fraction(f).denominator // gcd(d, Fraction(f).denominator)
    b = to_bold({k: int(Fraction(v) * d) for k, v in P.items()}, shift)
    return {k: Fraction(v, d) for k, v in b.items()}


def irr_coeffs(bold):
    """既約表示の係数（行が対称でなければ None）。"""
    out = {}
    for (b, f), r in w_rows(bold).items():
        if any(r.get(-k) != c for k, c in r.items()):
            return None
        K = max(abs(k) for k in r)
        for k in range(0, K + 1):
            v = r.get(k, 0) - r.get(k + 2, 0)
            if v:
                out[(b, k, f)] = v
    return out


def genus_frac(bold):
    d = 1
    for f in bold.values():
        d = d * f.denominator // gcd(d, f.denominator)
    try:
        N = g_expand({k: int(v * d) for k, v in bold.items()})
    except ValueError:
        return None
    return {k: Fraction(v, d) for k, v in N.items()}


def polytope(pieces, sigmas, kinds):
    """不等式 A y ≤ b（σ_ρ (h0 + Σ y_j hj) ≥ 0 を各表示で）。"""
    A, b = [], []
    kdim = None
    for rho, (h0, hs) in pieces.items():
        sig = sigmas[rho]
        shift = 1 if any(k[0] % 2 for k in list(h0) + [k for h in hs for k in h]) else 0
        b0 = bold_frac(h0, shift)
        bj = [bold_frac(h, shift) for h in hs]
        kdim = len(bj)
        reps = {"spin": (b0, bj)}
        if "genus" in kinds:
            g0 = genus_frac(b0)
            gj = [genus_frac(h) for h in bj]
            if g0 is None or any(g is None for g in gj):
                raise ValueError("G 展開が閉じない")
            reps["genus"] = (g0, gj)
        if "irr" in kinds:
            i0 = irr_coeffs(b0)
            ij = [irr_coeffs(h) for h in bj]
            if i0 is None or any(g is None for g in ij):
                raise ValueError("既約表示の行が対称でない")
            reps["irr"] = (i0, ij)
        for kind, (c0, cj) in reps.items():
            keys = set(c0)
            for c in cj:
                keys |= set(c)
            for e in keys:
                # σ (c0_e + Σ y_j cj_e) ≥ 0  ⇔  −σ Σ y_j cj_e ≤ σ c0_e
                A.append([-sig * float(c.get(e, 0)) for c in cj])
                b.append(sig * float(c0.get(e, 0)))
    return np.array(A), np.array(b), kdim


def extent(A, b, kdim, trials, rng):
    """多面体の実行可能性と、乱数の方向での幅の最大値、1 点なら y*。"""
    if kdim == 0:
        return (bool(np.all(b >= -1e-9)), 0.0, np.zeros(0)) if np.all(b >= -1e-9) else (None, None, None)
    res = linprog(np.zeros(kdim), A_ub=A, b_ub=b, bounds=[(None, None)] * kdim, method="highs")
    if res.status != 0:
        return None, None, None
    width = 0.0
    ystar = res.x
    for _ in range(trials):
        c = rng.normal(size=kdim)
        lo = linprog(c, A_ub=A, b_ub=b, bounds=[(None, None)] * kdim, method="highs")
        hi = linprog(-c, A_ub=A, b_ub=b, bounds=[(None, None)] * kdim, method="highs")
        if lo.status == 3 or hi.status == 3:
            return True, float("inf"), None
        if lo.status == 0 and hi.status == 0:
            width = max(width, float(-hi.fun - lo.fun))
            ystar = lo.x
    return True, width, ystar


def main():
    margin, hook, trials = 0, True, 12
    names = []
    for arg in sys.argv[1:]:
        if arg.startswith("--margin="):
            margin = int(arg.split("=")[1])
        elif arg == "--nohook":
            hook = False
        elif arg.startswith("--trials="):
            trials = int(arg.split("=")[1])
        else:
            names.append(arg)
    LL.MARGIN = margin
    lg = open("logs/lmov21_polytope.log", "a")

    def log(s):
        print(s, flush=True)
        lg.write(s + "\n")
        lg.flush()

    log(f"=== lmov21_polytope {time.strftime('%Y-%m-%d %H:%M:%S')} margin={margin} hook={hook} ===")
    out_path = Path("data/lmovstruct/lmov21_polytope.json")
    out = json.loads(out_path.read_text()) if out_path.exists() else {}
    rng = np.random.default_rng(0)
    for name in names:
        t0 = time.time()
        if name.startswith("tw:"):
            small, H21u, window, true21 = LL.twist_inputs(int(name[3:]))
            name = f"K_{name[3:]}"
        else:
            small, H21u, window, true21 = LL.known_inputs(name)
        conds = ("U", "T", "E", "U2", "H") if hook else ("U", "T", "E", "U2")
        lin, mons = build21(small, H21u, window, (1, 1, 1), conds=conds)
        rows = [(dict(r), b) for r, b in zip(lin.rows, lin.rhs)] + LL.laurent_rows(small, mons, lin.col)
        x0, basis = exact_affine(rows, len(mons))
        log(f"  {name}: 未知数 {len(mons)}、零空間 {len(basis)}（{time.time() - t0:.0f}s）")
        rec = {"window_t": window, "hook": hook, "n_unknowns": len(mons), "kernel_dim": len(basis)}
        pieces = hat_pieces(small, mons, lin.col, x0, basis)
        inv = {v: k for k, v in lin.col.items()}
        import itertools
        for kinds in (("spin",), ("spin", "genus"), ("spin", "genus", "irr")):
            best = None
            for sg in itertools.product((1, -1), repeat=3):
                sigmas = dict(zip(pieces.keys(), sg))
                try:
                    A, b, kdim = polytope(pieces, sigmas, kinds)
                except ValueError as e:
                    log(f"    {'+'.join(kinds)} σ={sg}: {e}")
                    break
                feas, width, ystar = extent(A, b, kdim, trials, rng)
                if feas:
                    log(f"    {'+'.join(kinds)} σ={sg}: 実行可能、幅 {width:.3g}")
                    if best is None or width < best[1]:
                        best = (sg, width, ystar)
            if best is None:
                rec["+".join(kinds)] = {"feasible": False}
                log(f"    {'+'.join(kinds)}: 実行不能")
                continue
            sg, width, ystar = best
            entry = {"feasible": True, "sigma": sg, "width": width}
            if width < 1e-6 and ystar is not None:
                y = [Fraction(float(v)).limit_denominator(10 ** 6) for v in ystar]
                x = {}
                for c, f in x0.items():
                    x[c] = x.get(c, 0) + f
                for yj, vec in zip(y, basis):
                    for c, z in vec.items():
                        x[c] = x.get(c, 0) + yj * z
                xi = {}
                integral = True
                for c, f in x.items():
                    r = round(f)
                    if abs(f - r) > 1e-6:
                        integral = False
                    if r:
                        xi[inv[c][1]] = int(r)
                ok_rows = all(sum(v * xi.get(inv[c][1], 0) for c, v in row.items()) == bb for row, bb in rows)
                entry.update({"unique": True, "integral": integral, "rows_ok": ok_rows, "terms": len(xi),
                              "matches_known": (true21 is not None) and (xi == true21),
                              "solution_bold": {",".join(map(str, k)): v for k, v in xi.items()}})
                log(f"    → 1 点：{len(xi)} 項、整数 {integral}、式の検算 {ok_rows}、既知と一致 {entry['matches_known']}")
            else:
                entry["unique"] = False
            rec["+".join(kinds)] = entry
        out[name] = rec
        out_path.write_text(json.dumps(out, ensure_ascii=False, indent=1))
        log(f"  {name}: 完了（{time.time() - t0:.0f}s）")


if __name__ == "__main__":
    main()
