"""第 4 期 V2-2：特殊化 A = s·Q^i·T^j で、rCS_[2,1] = α·X + β·Y（α, β は結び目によらない (Q, T) の関数）となる関係を数値で探す。

ランダムな (Q, T) の点で mod p の値を求め、T(2,3), T(2,5) から α, β を解き（2×2）、T(2,7), T(3,4), 4_1, T(3,5) で同じ点で成り立つかを確かめる。
3 点の全てで成り立った組を候補として出力する。X = 1 または Y = 1 の組も含む（α X + β の形）。結果は logs/framed_knot/v2_linear_relations.log。
"""
import sys
from itertools import combinations, product
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from v2_specializations import load  # noqa: E402

P = 2**31 - 1
SOLVE, CHECK = ["T2_3", "T2_5"], ["T2_7", "T3_4", "4_1", "T3_5"]
NAMES = ["r1·r11", "r1·r2", "r3", "r111", "r1^3", "r1", "r2", "r11", "1"]


def ev(poly, A, Q, T):
    s = 0
    for (a, q, t), c in poly.items():
        s = (s + c * pow(A, a % (P - 1), P) * pow(Q, q % (P - 1), P) * pow(T, t % (P - 1), P)) % P
    return s


def values(r, A, Q, T):
    v = {c: ev(p, A, Q, T) for c, p in r.items()}
    return {"21": v["21"], "r1·r11": v["1"] * v["11"] % P, "r1·r2": v["1"] * v["2"] % P, "r3": v["3"], "r111": v["111"],
            "r1^3": pow(v["1"], 3, P), "r1": v["1"], "r2": v["2"], "r11": v["11"], "1": 1}


def main():
    data = load()
    import random
    rnd = random.Random(1)
    pts = [(rnd.randrange(2, P - 1), rnd.randrange(2, P - 1)) for _ in range(3)]
    for s, i, j in product((1, -1), range(-3, 4), range(-3, 4)):
        if (i, j) == (0, 0):
            continue
        vals = []
        for Q, T in pts:
            A = s * pow(Q, i % (P - 1), P) * pow(T, j % (P - 1), P) % P
            vals.append({k: values(data[k], A, Q, T) for k in SOLVE + CHECK})
        for X, Y in combinations(NAMES, 2):
            ok = True
            for v in vals:
                a1, b1, c1 = v[SOLVE[0]][X], v[SOLVE[0]][Y], v[SOLVE[0]]["21"]
                a2, b2, c2 = v[SOLVE[1]][X], v[SOLVE[1]][Y], v[SOLVE[1]]["21"]
                det = (a1 * b2 - a2 * b1) % P
                if det == 0:
                    ok = False
                    break
                inv = pow(det, P - 2, P)
                al, be = (c1 * b2 - c2 * b1) * inv % P, (a1 * c2 - a2 * c1) * inv % P
                if any((v[k]["21"] - al * v[k][X] - be * v[k][Y]) % P for k in CHECK):
                    ok = False
                    break
            if ok:
                print(f"A = {'+' if s > 0 else '-'}Q^{i} T^{j}: r21 = α·{X} + β·{Y}（α, β は結び目によらない）", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
