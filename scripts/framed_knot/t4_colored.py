"""ステージ T4：GS（arXiv:1112.0030）の S² 彩色ホモロジー（t_r 次数）から t_c 次数の [2], [1,1] を作り、2 箱の refined LMOV を判定する。

t_c 次数への変換（三葉結び目と 8_19 = T(3,4) で第 1 期の rCS と照合して決めた。docs/framed_knot/conventions.md §11）：
  GS の old 規約の生成元 x = (i, j, k) と δ'(x) から
    [2]   ：(e_𝐚, e_𝐪, e_𝐭) = (2i,  2j, 4i − k + 2j + 4δ'(x))
    [1,1] ：(e_𝐚, e_𝐪, e_𝐭) = (2i, −2j, 4i − k + 4δ'(x))
  δ' は GS の「δ 次数の一般化」で、thin な結び目では −S（三葉結び目で −1）。thick な結び目では GS に規則が無い。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "data/framed_knot/reference/gs_1112.0030"))

import s2_thick as gs  # noqa: E402
from rcsknots.extract import bps_numbers, hat_f_exact  # noqa: E402
from rcsknots.knot_table import load_homology, superpoly  # noqa: E402


def bold_to_small(bold):
    """太字 (e_𝐚, e_𝐪, e_𝐭) → KN の (A, Q, T) = (a, q, t)^{1/2} の指数（rcsknots.bps.to_bold の逆、a 次数シフトなし）。"""
    out = {}
    for (ea, eq, et), c in bold.items():
        m = ea // 2
        key = (2 * m, et - m, eq - et + m)
        out[key] = out.get(key, 0) + c * (-1) ** (et % 2)
    return out


def tc_colored(gens, delta):
    """GS の生成元 {(i, j, k): 重複度} と δ'（関数 (i, j, k) → 値、4δ' は整数）から t_c 次数の [2], [1,1]（太字）。"""
    two, one_one = {}, {}
    for (i, j, k), c in gens.items():
        d4 = round(4 * delta((i, j, k)))
        for D, key in ((two, (2 * i, 2 * j, 4 * i - k + 2 * j + d4)), (one_one, (2 * i, -2 * j, 4 * i - k + d4))):
            D[key] = D.get(key, 0) + c
    return two, one_one


def judge2(P1, P2, P11, log=lambda s: None):
    """[1], [2], [1,1]（太字）から \\hat f_{[2]}, \\hat f_{[1,1]} を作り、G 展開の前後の負の項を数える（全体の符号は良い方）。"""
    reduced = {(1,): bold_to_small(P1), (2,): bold_to_small(P2), (1, 1): bold_to_small(P11)}
    out = {}
    for rho, hat in hat_f_exact(reduced, 2, log=log).items():
        if hat is None:
            out[rho] = "Laurent でない"
            continue
        from rcsknots.bps import to_bold
        from rcsknots.extract import A_SHIFT
        before = to_bold(hat, A_SHIFT["nontorus"][2])
        nb = min(sum(v < 0 for v in before.values()), sum(v > 0 for v in before.values()))
        try:
            N = bps_numbers(hat, rho, "nontorus")
        except ValueError:
            out[rho] = f"前 {nb}/{len(before)} 負、G 展開が閉じない"
            continue
        na = min(sum(v < 0 for v in N.values()), sum(v > 0 for v in N.values()))
        out[rho] = f"前 {nb}/{len(before)} 負、後 {na}/{len(N)} 負"
    return out


def main():
    H = load_homology()
    line2 = gs.terms((1 + gs.q) * (1 + gs.a**-1 * gs.q * gs.t**-1) * (1 + gs.a**-1 * gs.q**-2 * gs.t**-3)
                     * (1 + gs.a**-1 * gs.q**-3 * gs.t**-3) * (gs.a**2 * gs.q**2 * gs.t**6 + gs.a * gs.q**3 * gs.t**5))
    printed = gs.terms(gs.P2_942)
    fixed = {k: printed[k] - line2.get(k, 0) for k in printed if printed[k] - line2.get(k, 0)}
    cases = [
        ("3_1", "3_1", gs.terms(gs.P2_31), -1),
        ("8_19（δ' = −S = −3 一定）", "8_19", gs.terms(gs.P2_819), -3),
        ("9_42 印刷どおり（δ' = 0 一定）", "9_42", printed, 0),
        ("9_42 2 行目を除く（δ' = 0 一定）", "9_42", fixed, 0),
    ]
    for label, name, gens, d in cases:
        P2, P11 = tc_colored(gens, lambda x, d=d: d)
        print(label, judge2(superpoly(H[name]), P2, P11), flush=True)


if __name__ == "__main__":
    main()
