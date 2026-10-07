"""open_problems D1 (I) の組合せ論の 1 段：⟨Hik_{m,n}, e_n⟩ = Σ_D q^{δ−|D|} t^{h₊(D)} を BGLX の定義どおりに確かめる。

BGLX（arXiv:1404.4616, Final-4-17.tex:418-451）の定義をそのまま実装する：
  (m,n)-Dyck 路：m×n の格子で (0,0)→(m,n)、N/E の歩み、lattice diagonal の上（各格子点で m·y ≥ n·x）。
  area(π)：路と lattice diagonal の間のセルの数。λ(π)：路の上（北西）のセルの Ferrers 図形。
  dinv(π)（式 1.8）：λ(π) のセル c で arm/(leg+1) < m/n < (arm+1)/leg を満たすものの数（GN の h₊ と同じ式）。
  駐車関数：N の歩みの東隣のセルに車 1..n を列ごとに上へ増えるように置く。rank：北西のセルを r0 とし、東へ −n、北へ +m。
  tdinv（式 1.9）、maxtdinv（1.10）、dinv(PF) = dinv(π) + tdinv − maxtdinv（1.11）。
  σ(PF)：車を rank の大きい順に読んだ語。ides = σ^{-1} の降下集合。Hik = Σ t^{area} q^{dinv} F_{ides}（defHikita）。
Gessel の基本準対称関数の展開から ⟨Hik, e_n⟩ = Σ_{ides = [n−1]}、⟨Hik, h_n⟩ = Σ_{ides = ∅}。
確かめること：
  K1 ⟨Hik,e_n⟩ に寄与する駐車関数は各 Dyck 路にちょうど 1 つ（車の番号が rank の順に増えるもの、σ = n…1）で、その dinv = dinv(π)。
  K2 area(π) = δ − |λ(π)|（定義から）。
  K3 ⟨Hik,e_n⟩ = Σ_π t^{area} q^{dinv(π)} が GN の C_{m,n}（scripts/d1_catalan.py の catalan）の (q,t) または (t,q) に一致。
使い方：PYTHONPATH=src:scripts python scripts/d1_hik.py（ログ logs/d1_hik.log）。
"""
from itertools import permutations
from math import gcd, factorial

import sympy as sp

from d1_catalan import catalan, q, t


def dyck_paths(m, n):
    """N の歩みの x 座標の列 (x_0, …, x_{n−1})（行 y の N の歩みは (x_y, y)→(x_y, y+1)）。"""
    def rec(y, prev):
        if y == n:
            yield ()
            return
        for x in range(prev, (m * y) // n + 1):
            for rest in rec(y + 1, x):
                yield (x,) + rest
    return list(rec(0, 0))


def path_dinv(xs, m, n):
    cnt = 0
    for y, r in enumerate(xs):
        for x in range(r):
            arm = r - x - 1
            leg = sum(1 for y2 in range(y) if xs[y2] > x)
            if arm * n < m * (leg + 1) and (leg == 0 or m * leg < n * (arm + 1)):
                cnt += 1
    return cnt


def main():
    log = open("logs/d1_hik.log", "w")

    def out(s):
        print(s, flush=True)
        log.write(s + "\n")

    for m in range(2, 8):
        for n in range(2, 7):
            if gcd(m, n) != 1 or m ** (n - 1) > 20000:
                continue
            delta = (m - 1) * (n - 1) // 2
            r0 = (m - 1) * (n - 1) - 1
            paths = dyck_paths(m, n)
            assert len(paths) == factorial(m + n - 1) // (factorial(m) * factorial(n))
            E = H = 0
            k1 = k2 = True
            npf = 0
            for xs in paths:
                lam = sum(xs)
                area = sum((m * y) // n - x for y, x in enumerate(xs))
                k2 &= area == delta - lam
                dpath = path_dinv(xs, m, n)
                rank = [r0 - n * x + m * (y - (n - 1)) for y, x in enumerate(xs)]
                pfs = []
                for perm in permutations(range(1, n + 1)):        # perm[y] = 行 y の車
                    if any(xs[y] == xs[y + 1] and perm[y] > perm[y + 1] for y in range(n - 1)):
                        continue
                    rk = {perm[y]: rank[y] for y in range(n)}
                    td = sum(1 for i in range(1, n + 1) for j in range(i + 1, n + 1) if rk[i] < rk[j] < rk[i] + m)
                    sigma = sorted(range(1, n + 1), key=lambda c: -rk[c])
                    pos = {c: k for k, c in enumerate(sigma)}
                    ides = frozenset(i for i in range(1, n) if pos[i + 1] < pos[i])
                    pfs.append((perm, td, ides, rk))
                npf += len(pfs)
                mx = max(p[1] for p in pfs)
                full = [p for p in pfs if p[2] == frozenset(range(1, n))]
                empty = [p for p in pfs if not p[2]]
                k1 &= len(full) == 1 and full[0][1] == mx and \
                    all(full[0][3][i] < full[0][3][i + 1] for i in range(1, n))
                for p in full:
                    E += t ** area * q ** (dpath + p[1] - mx)
                for p in empty:
                    H += t ** area * q ** (dpath + p[1] - mx)
            assert npf == m ** (n - 1)
            C = catalan(m, n)
            E, H = sp.expand(E), sp.expand(H)
            sw = sp.expand(C.subs({q: t, t: q}, simultaneous=True))
            out(f"(m,n)=({m},{n})：Dyck 路 {len(paths)}、駐車関数 {npf}；K1 {k1}、K2 {k2}、"
                f"⟨Hik,e_n⟩ = C(q,t) {sp.expand(E - C) == 0}、= C(t,q) {sp.expand(E - sw) == 0}、"
                f"⟨Hik,h_n⟩ = C {sp.expand(H - C) == 0}")


if __name__ == "__main__":
    main()
