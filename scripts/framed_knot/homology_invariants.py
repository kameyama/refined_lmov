"""ケーブルの reduced Khovanov ホモロジー（khoca、係数 F_101）と knot Floer ホモロジー（snappy）を計算する。

ブレイドは data/framed_knot/reference/braid_homfly.txt（Stage 4 の HOMFLY-PT の照合に使ったもの）。
khoca の t は通常の規約の逆（三葉結び目 σ1^3 で t^-3 q^8 + t^-2 q^6 + q^2）なので、t → −t で直して保存する。
出力：data/framed_knot/homology/{K}.json（"kh": {"q,t": rank}, "hfk": {"alexander,maslov": rank} と付随情報）
使い方：python scripts/framed_knot/homology_invariants.py K2 K4 ...
"""
import json
import re
import sys
import time
from pathlib import Path

OUT = Path("data/framed_knot/homology")


def braid(name):
    for line in Path("data/framed_knot/reference/braid_homfly.txt").read_text().splitlines():
        cols = line.split("\t")
        if cols[0] == name:
            return [int(x) for x in cols[2].split()]
    raise KeyError(name)


def khovanov(word):
    from khoca import InteractiveCalculator
    KH = InteractiveCalculator(101, (0, 0), 0)
    _, mess = KH(tuple(word), print_messages=True)
    line = mess[mess.index("Reduced Homology:") + 3]
    out = {}
    for c, t, q in re.findall(r"(\d*)t\^(-?\d+)q\^(-?\d+)", line):
        out[(int(q), -int(t))] = int(c) if c else 1
    return out


def knot_floer(word):
    import snappy
    L = snappy.Link(braid_closure=word)
    H = L.knot_floer_homology()
    return H


def main(names):
    OUT.mkdir(parents=True, exist_ok=True)
    for name in names:
        word = braid(name)
        t0 = time.time()
        kh = khovanov(word)
        t1 = time.time()
        H = knot_floer(word)
        t2 = time.time()
        deltas = sorted({q - 2 * t for (q, t) in kh})
        hfk = {f"{a},{m}": r for (a, m), r in H["ranks"].items()}
        hfk_deltas = sorted({a - m for (a, m) in H["ranks"]})
        data = {"knot": name, "crossings": len(word),
                "kh": {f"{q},{t}": r for (q, t), r in sorted(kh.items())},
                "kh_total_rank": sum(kh.values()), "kh_deltas": deltas,
                "hfk": hfk, "hfk_total_rank": H["total_rank"], "hfk_deltas": hfk_deltas,
                "seifert_genus": H["seifert_genus"], "tau": H["tau"], "L_space_knot": H["L_space_knot"],
                "fibered": H["fibered"], "seconds": {"khoca": round(t1 - t0, 1), "hfk": round(t2 - t1, 1)}}
        (OUT / f"{name}.json").write_text(json.dumps(data, indent=1))
        print(name, f"Kh rank {data['kh_total_rank']} δ {deltas}", f"HFK rank {H['total_rank']} δ {hfk_deltas}",
              f"g {H['seifert_genus']} τ {H['tau']} L-space {H['L_space_knot']}", data["seconds"], flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
