"""保存済みの \\hat f（data/lmov/{K}_{n}.json）から、G 展開の後の \\hat N を data/bps/{K}_{ρ}.json と .txt に書き出す。

使い方：python scripts/export_bps.py [--check] K ...
  形式は scripts/stage3_bps.py と同じ（a 次数シフトは nontorus 規約、全体の符号は rcsknots.tables.normalize_sign）。
  --check：書き出さずに、既存の data/bps/ と \\hat f・\\hat N を全係数で比べる。
"""
import json
import lzma
import sys
import time
from pathlib import Path

from rcsknots.extract import A_SHIFT, bps_numbers
from rcsknots.tables import normalize_sign, save

CONVENTION = "nontorus"


def load_hat(name, n):
    path = Path(f"data/lmov/{name}_{n}.json")
    text = path.read_text() if path.exists() else lzma.open(f"{path}.xz", "rt").read()
    f = json.loads(text)["f"]
    return {tuple(int(c) for c in k.split("_")[1]): {tuple(map(int, e.split(","))): v for e, v in d.items()}
            for k, d in f.items() if k.startswith("hat_")}


def main(names, check):
    for name in names:
        for n in (1, 2, 3):
            for rho, hf in sorted(load_hat(name, n).items()):
                t0 = time.time()
                sigma, N = normalize_sign(bps_numbers(hf, rho, CONVENTION))
                stem = f"data/bps/{name}_{''.join(map(str, rho))}"
                neg = sum(v < 0 for v in N.values())
                if check:
                    old = json.loads(Path(stem + ".json").read_text())
                    same_f = old["hat_f"] == {",".join(map(str, k)): v for k, v in sorted(hf.items())}
                    same_N = old["N"] == {",".join(map(str, k)): v for k, v in sorted(N.items())}
                    print(f"{name} ρ={list(rho)}: hat_f 一致 {same_f}, N 一致 {same_N}", flush=True)
                    continue
                meta = {"knot": name, "a_shift_convention": CONVENTION, "a_shift_A_units": A_SHIFT[CONVENTION][n],
                        "overall_sign_applied": sigma, "ft_consistency_checked": True,
                        "variables_hat_f": "exponents of (A=a^(1/2), Q=q^(1/2), T=t^(1/2)) of 1703.05408, before the a-shift",
                        "N_keys": "g,beta,2J_r", "source": f"data/lmov/{name}_{n}.json (scripts/export_bps.py)",
                        "seconds_g_expand": round(time.time() - t0, 1)}
                save(stem, rho, hf, N, meta)
                print(f"{name} ρ={list(rho)}: 前 {len(hf)} 項、後 {len(N)} 項、σ = {sigma}、負 {neg}、"
                      f"{time.time() - t0:.1f} s", flush=True)


if __name__ == "__main__":
    args = sys.argv[1:]
    check = args[:1] == ["--check"]
    main(args[1:] if check else args, check)
