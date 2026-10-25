"""四重次数（GGS）の S² 側の構造を、正解のデータから探す（ユーザーの指示「DGR や quadruply の計算をいくつかの例でフォローし…」）。

正解：S² の t_c 次数の Poincaré 多項式 P₂（太字）。ツイスト結び目は GNSSS の予想式（scripts/framed_knot/v4_closed.inputs）、
10_124, 10_139, 10_152 は quiver と KRSS の予想式（data/framed_knot/b4）。[1] は P₁（superpolynomial）。
各 d_{n|m}（0 ≤ n, m ≤ 3）の特殊化 𝐚² = −𝐪^{2(n−m)}𝐭^{−2m−1} について、
  (i)  P₂ の像が単項式か（打ち消す微分で 1 次元のホモロジー）
  (ii) P₂ の像 = 単項式 × P₁ の d_{n'|m'} の像を 𝐪 ↦ 𝐪^α 𝐭^β で付け替えたもの、となる (n', m', α, β) があるか（色付き微分）
を調べる。使い方：PYTHONPATH=src python scripts/dgr_lmov/dgr_s2_discover.py
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "framed_knot"))
from v4_closed import inputs  # noqa: E402
from rcsknots.dgr_lmov.dgr import specialize  # noqa: E402
from rcsknots.knot_table import load_homology, superpoly  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


def regrade(S, al, be):
    return {(al * eq, et + be * eq): c for (eq, et), c in S.items()}


def ratio_mono(A, B):
    if not A or not B or len(A) != len(B):
        return None
    ka, kb = min(A), min(B)
    s = (ka[0] - kb[0], ka[1] - kb[1])
    if A[ka] % B[kb]:
        return None
    c = A[ka] // B[kb]
    return (s, c) if all(A.get((k[0] + s[0], k[1] + s[1])) == c * v for k, v in B.items()) else None


def analyse(name, P1, P2, log):
    P2 = {k: int(v) for k, v in P2.items() if v}
    P1 = {k: int(v) for k, v in P1.items() if v}
    res = {}
    for n in range(4):
        for m in range(4):
            S2 = specialize(P2, n, m)
            if len(S2) == 1:
                res[f"{n}|{m}"] = f"単項式 {S2}"
                continue
            found = []
            for n2 in range(4):
                for m2 in range(4):
                    S1 = specialize(P1, n2, m2)
                    for al in (1, 2):
                        for be in range(-3, 4):
                            r = ratio_mono(S2, regrade(S1, al, be))
                            if r:
                                found.append((n2, m2, al, be, r))
            res[f"{n}|{m}"] = found[:3] if found else "なし"
    log(f"{name}: " + "  ".join(f"d_{{{k}}}: {v}" for k, v in res.items() if v != "なし"))
    return res


def main():
    logf = open(ROOT / "logs/dgr_lmov/dgr_s2_discover.log", "w")

    def log(s):
        print(s, flush=True)
        logf.write(s + "\n")
    out = {}
    for p in (1, -1, 2, -2, 3, -3, 4):
        P1, P2 = inputs(p)
        out[f"K_{p}"] = analyse(f"K_{p}", P1, P2, log)
    H = load_homology()
    for K in ("10_124", "10_139", "10_152"):
        d = json.loads((ROOT / f"data/framed_knot/b4/{K}_s2_tc.json").read_text())
        P2 = {tuple(map(int, k.split(","))): v for k, v in d["[2]"].items()}
        out[K] = analyse(K, superpoly(H[K]), P2, log)
    (ROOT / "data/dgr_lmov/s2_discover.json").write_text(json.dumps({k: {kk: str(vv) for kk, vv in v.items()} for k, v in out.items()}, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
