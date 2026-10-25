r"""sympy の有理関数体（sympy.polys.fields）の元を係数とする小さい連立一次方程式を、python-flint の fmpz_mpoly で
分母を払ってから分数を使わない Bareiss の消去で解く（sympy の体での Gauss の消去は約分が重い）。"""
from math import lcm

import flint
from concurrent.futures import ProcessPoolExecutor


def _ctx(field):
    names = tuple(str(g) for g in field.gens)
    return flint.fmpz_mpoly_ctx.get(names, "lex")


def _poly_to_flint(p, ctx):
    """sympy の PolyElement（QQ 係数）→（fmpz_mpoly、整数の分母）。"""
    terms = p.terms()
    den = 1
    for _, c in terms:
        den = lcm(den, int(c.denominator))
    d = {tuple(m): int(c.numerator) * (den // int(c.denominator)) for m, c in terms}
    return ctx.from_dict(d), den


def _flint_to_poly(f, ring):
    return ring.from_dict({m: c for m, c in f.to_dict().items()}) if not f.is_zero() else ring.zero


def solve(M, rhs, field):
    """M x = rhs（field の元）。解を field の元のリストで返す。"""
    ctx = _ctx(field)
    n = len(M)
    rows = []
    for i in range(n):
        ents = list(M[i]) + [rhs[i]]
        nums, dens = [], []
        for e in ents:
            pn, dn = _poly_to_flint(e.numer, ctx)
            pd, dd = _poly_to_flint(e.denom, ctx)
            nums.append((pn * dd, pd * dn))          # e = (pn/dn)/(pd/dd) = pn·dd / (pd·dn)
        L = nums[0][1]
        for _, d in nums[1:]:
            g = L.gcd(d)
            L = L * (d / g)
        rows.append([a * (L / b) for a, b in nums])
    one = ctx.from_dict({(0,) * len(ctx.gens()): 1})

    def det(Mx):
        A = [r[:] for r in Mx]
        m = len(A)
        prev, sign = one, 1
        for k in range(m):
            if A[k][k].is_zero():
                for i2 in range(k + 1, m):
                    if not A[i2][k].is_zero():
                        A[k], A[i2] = A[i2], A[k]; sign = -sign
                        break
                else:
                    return ctx.from_dict({})
            piv = A[k][k]
            for i2 in range(k + 1, m):
                for j in range(k + 1, m):
                    A[i2][j] = (A[i2][j] * piv - A[i2][k] * A[k][j]) / prev
            prev = piv
        return A[m - 1][m - 1] * sign

    coef = [r[:n] for r in rows]
    b = [r[n] for r in rows]
    D = det(coef)
    assert not D.is_zero(), "特異"
    ring = field.ring
    sol = []
    for i in range(n):
        Mi = [row[:i] + [b[k]] + row[i + 1:] for k, row in enumerate(coef)]
        Di = det(Mi)
        g = Di.gcd(D)
        sol.append(field(_flint_to_poly(Di / g, ring)) / field(_flint_to_poly(D / g, ring)))
    return sol



def _det_worker(args):
    """(変数名, 行列 [[{単項式: 整数}]]) → 行列式 {単項式: 整数}（別プロセスで呼ぶ）。"""
    names, rows = args
    ctx = flint.fmpz_mpoly_ctx.get(tuple(names), "lex")
    A = [[ctx.from_dict(e) for e in row] for row in rows]
    m = len(A)
    one = ctx.from_dict({(0,) * len(names): 1})
    prev, sign = one, 1
    for k in range(m):
        if A[k][k].is_zero():
            for i2 in range(k + 1, m):
                if not A[i2][k].is_zero():
                    A[k], A[i2] = A[i2], A[k]; sign = -sign
                    break
            else:
                return {}
        piv = A[k][k]
        for i2 in range(k + 1, m):
            for j in range(k + 1, m):
                A[i2][j] = (A[i2][j] * piv - A[i2][k] * A[k][j]) / prev
        prev = piv
    return (A[m - 1][m - 1] * sign).to_dict()


def solve_parallel(M, rhs, field, workers=8):
    """solve と同じだが、Cramer の公式の n + 1 個の行列式を別プロセスで並列に計算する。"""
    ctx = _ctx(field)
    names = [str(g) for g in field.gens]
    n = len(M)
    rows = []
    for i in range(n):
        nums = []
        for e in list(M[i]) + [rhs[i]]:
            pn, dn = _poly_to_flint(e.numer, ctx)
            pd, dd = _poly_to_flint(e.denom, ctx)
            nums.append((pn * dd, pd * dn))
        L = nums[0][1]
        for _, d in nums[1:]:
            L = L * (d / L.gcd(d))
        rows.append([(a * (L / b)).to_dict() for a, b in nums])
    coef = [r[:n] for r in rows]
    b = [r[n] for r in rows]
    jobs = [(names, coef)] + [(names, [row[:i] + [b[k]] + row[i + 1:] for k, row in enumerate(coef)]) for i in range(n)]
    with ProcessPoolExecutor(max_workers=min(workers, len(jobs))) as ex:
        dets = list(ex.map(_det_worker, jobs))
    D = ctx.from_dict(dets[0])
    assert not D.is_zero(), "特異"
    ring = field.ring
    sol = []
    for i in range(n):
        Di = ctx.from_dict(dets[i + 1])
        g = Di.gcd(D)
        sol.append(field(_flint_to_poly(Di / g, ring)) / field(_flint_to_poly(D / g, ring)))
    return sol
