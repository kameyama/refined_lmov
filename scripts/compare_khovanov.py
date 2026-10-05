"""[1] の refined 不変量（DAHA/EHA の rCS_□）を sl(2) に特殊化し、実際の reduced Khovanov ホモロジーと比べる。

太字の変数で P(𝐚,𝐪,𝐭) とし、𝐚 = 𝐪^2 とおいた P_2(q,t) を Kh(q,t)（data/homology）と比べる。
全体の単項式 ±q^x t^y は、Euler 標数（t = −1）が一致するように x と符号を、y は走査して選ぶ。
P が HOMFLY-PT ホモロジーの Poincaré 多項式なら、HOMFLY → sl(2) のスペクトル系列（微分は q 次数を保ち t 次数を 1 変える）により
  P_2 − Kh = (1 + t) Q,  Q ≥ 0
が成り立つ（必要条件）。各 y についてこれを判定する。
使い方：python scripts/compare_khovanov.py T2_3 K1 K2 ...
"""
import json
import sys
from pathlib import Path

from rcsknots.bps import to_bold
from rcsknots.laurent import load_table


def rcs_box(name):
    if name.startswith("K"):
        d = json.loads(Path(f"data/cables/{name}_1.json").read_text())["poly"]
        return {tuple(map(int, k.split(","))): v for k, v in d.items()}
    return load_table("data/reference/notebook_torus.txt")[(name, "1")]


def specialize(P):
    """bold (e_𝐚, e_𝐪, e_𝐭) → (q, t) with 𝐚 = 𝐪^2."""
    out = {}
    for (ea, eq, et), c in to_bold(P, next(iter(P))[0] % 2).items():
        k = (2 * ea + eq, et)
        out[k] = out.get(k, 0) + c
    return {k: v for k, v in out.items() if v}


def at_minus_one(P):
    out = {}
    for (q, t), c in P.items():
        out[q] = out.get(q, 0) + c * (-1) ** t
    return {k: v for k, v in out.items() if v}


def divide_one_plus_t(R):
    """R = (1 + t) Q in Z[q^±, t^±]; returns Q or None."""
    Q, R = {}, dict(R)
    while R:
        q, t = min(R, key=lambda k: (k[0], k[1]))
        qs = [k for k in R if k[0] == q]
        t = min(k[1] for k in qs)
        c = R[(q, t)]
        Q[(q, t)] = c
        for dt in (0, 1):
            R[(q, t + dt)] = R.get((q, t + dt), 0) - c
            if R[(q, t + dt)] == 0:
                del R[(q, t + dt)]
        if len(Q) > 10000:
            return None
        if any(k[0] == q and k[1] > 200 for k in R):
            return None
    return Q


def main(names):
    for name in names:
        kh = {tuple(map(int, k.split(","))): v for k, v in json.loads(Path(f"data/homology/{name}.json").read_text())["kh"].items()}
        P2 = specialize(rcs_box(name))
        e_kh, e_p = at_minus_one(kh), at_minus_one(P2)
        # match Euler characteristics: e_kh(q) = s · e_p(q − x)
        found = None
        for s in (1, -1):
            for x in range(-80, 81, 2):
                if {q + x: s * c for q, c in e_p.items()} == e_kh:
                    found = (s, x)
        neg = sum(v < 0 for v in P2.values())
        line = [f"{name}: P_2 の項 {len(P2)}（負 {min(neg, len(P2) - neg)}）、|P_2| の総和 {sum(map(abs, P2.values()))}、Kh の総階数 {sum(kh.values())}"]
        if not found:
            print(line[0], "Euler 標数が一致しない")
            continue
        s, x = found
        ok = []
        for y in range(-30, 31):
            if (y % 2) != 0 and s == s:  # parity of y fixed by the sign of t = −1 match
                pass
            P = {(q + x, t + y): s * c * (-1) ** y for (q, t), c in P2.items()}
            if at_minus_one(P) != e_kh:
                continue
            R = {k: P.get(k, 0) - kh.get(k, 0) for k in set(P) | set(kh)}
            R = {k: v for k, v in R.items() if v}
            Q = divide_one_plus_t(R)
            if Q is not None and all(v >= 0 for v in Q.values()):
                ok.append((y, sum(Q.values())))
        line.append(f"(1+t)Q, Q ≥ 0 となる t のずれ y：{ok if ok else 'なし'}")
        print("  ".join(line), flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
