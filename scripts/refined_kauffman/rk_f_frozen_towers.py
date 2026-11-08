r"""凍結した状態の σ 非対称部分に含まれる t の有限の塔（theory §12.5、progress.md 続き 64）。
  ブロック X_β = A⁻_β + D_β、D_β = Σ_{m>0} d_m q_ℓ^m（rk_f_frozen_basis.ell_table の記法）。
  検査：(1) D_β = (1 − q/t) × e_β で、e_β を q の冪ごとの列（t の多項式）に分け、各列を t の指標 [L]_t = 1 + t + … + t^{L−1} の
      シフト t^k [L]_t の和に貪欲に分解する（同符号なら可能。最小の L の個数を記録）。
      (2) 最低の β のブロックの閉じた形：D = −(1 − q/t) w^{−(8m+1)} t ([2n]_t + [2n−1]_t)（n = 2m + 1）の検算。
      (3) 別の読み：D (1 − t)/(t − q) が疎になるか（「t モードの塔 − q モードの塔」の形）。
  対象：T(2,3), T(2,5) の 2 箱 p₁² hatted の c = 1（ν 奇）の σ 非対称なブロック。"""
import sys
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import sympy as sp
import rk_f_frozen_basis as F

r, s, nu = F.r, F.s, F.nu
w = sp.Symbol("w"); t = sp.Symbol("t")


def asym_rs(blk):
    """D = Σ_{m>0} (c_m − c_{−m}) q_ℓ^m を (r, s) で（q_ℓ^m w^j = r^{m+j} s^{m−j}）。"""
    d = F.ell_table(blk)
    D = sp.Integer(0)
    for (m, j) in {(abs(m), j) for (m, j) in d}:
        if m > 0:
            D += (d.get((m, j), 0) - d.get((-m, j), 0)) * r**(m + j) * s**(m - j)
    return sp.expand(D)


def columns(e):
    """e（r, s の Laurent 多項式）を r の冪ごとに {a: {b: c}}（r^a s^b）に。"""
    out = {}
    for tm in sp.Add.make_args(sp.expand(e)):
        c, mon = tm.as_coeff_Mul()
        a = int(mon.as_coeff_exponent(r)[1]) if mon != 1 else 0
        b = int(mon.as_coeff_exponent(s)[1]) if mon != 1 else 0
        out.setdefault(a, {})[b] = int(c)
    return out


def sign_changes(col):
    """列 {b: c} を t の次数の順に見たときの符号の変わる回数（0 は飛ばす）。"""
    sg = [1 if c > 0 else -1 for b, c in sorted(col.items()) if c != 0]
    return sum(1 for a, b in zip(sg, sg[1:]) if a != b)


def truncated_tower_g(L, eps, kind):
    """長さ L の z₂ の塔をもつ状態（電荷 1）の g = f − ½f^{cov} への寄与（theory §12.5）。
    kind：'pair'（Ĩ で入れ替わる対）、'J2'（Ĩ 不変、塔の l 番目に (−1)^l ε）、'SR'（Ĩ 不変、塔に ε が一様にかかる）、'frozen'。"""
    tw = sum(t**l for l in range(L))
    if kind == "pair":
        return sp.expand(tw - sp.Rational(1, 2) * 2 * tw)
    if kind == "J2":
        return sp.expand(sum(t**l for l in range(L) if (-1)**l * eps == 1) - sp.Rational(1, 2) * tw)
    if kind == "SR":
        return sp.expand((tw if eps == 1 else 0) - sp.Rational(1, 2) * tw)
    if kind == "frozen":
        return tw
    raise ValueError(kind)


def greedy_characters(col):
    """列 {b: c}（s^b、b は同じ偶奇）を t^k [L]_t（連続する b に同じ係数 1）の和に貪欲に分解する。
    係数がすべて同符号なら必ず分解でき、結果は (符号, [(k, L, 重複度)])。混在なら None。"""
    sgn = F.sign_of(col.values())
    if sgn in ("±", "0"):
        return None
    rem = {b: abs(c) for b, c in col.items()}
    parts = []
    while any(rem.values()):
        bs = sorted(b for b, c in rem.items() if c > 0)
        # 最長の連続区間（b が 2 ずつ）を取る
        start = bs[0]; L = 1
        while start + 2 * L in rem and rem[start + 2 * L] > 0:
            L += 1
        mult = min(rem[start + 2 * i] for i in range(L))
        for i in range(L):
            rem[start + 2 * i] -= mult
        parts.append((start, L, mult))
    return sgn, parts


if __name__ == "__main__":
    print("== 長さ L の塔の g への寄与（[L]_{−t} = (1 − (−t)^L)/(1 + t)）", flush=True)
    for L in (1, 2, 3, 4):
        for eps in (1, -1):
            j2 = truncated_tower_g(L, eps, "J2")
            ok = sp.expand(j2 - sp.Rational(eps, 2) * sp.cancel((1 - (-t)**L) / (1 + t))) == 0
            print(f"  L={L} ε={eps:+d}：対 {truncated_tower_g(L, eps, 'pair')}、J₂ の塔 {j2}（= ½ε[L]_{{−t}} {ok}）、"
                  f"S_R の塔 {truncated_tower_g(L, eps, 'SR')}、凍結 {truncated_tower_g(L, eps, 'frozen')}", flush=True)
    for m_ in (1, 2):
        n = 2 * m_ + 1
        X = F.knot_two_box(m_)
        blocks = F.nu_blocks(F.parity_split(X, False))
        print(f"== T(2,{n}) 2 箱 p₁² hatted の c=1（ν 奇）の σ 非対称部分", flush=True)
        for b, blk in sorted(blocks.items()):
            D = asym_rs(blk)
            if D == 0:
                continue
            e = sp.cancel(D / (1 - r**2 / s**2))
            assert sp.Poly(sp.fraction(e)[1], r, s).is_monomial
            cols = columns(e)
            print(f"  ν^{b}: e = D/(1 − q/t) の q の冪ごとの列（r^a、t の指標 t^k[L]_t への分解）:")
            for a in sorted(cols):
                res = greedy_characters(cols[a])
                if res is None:
                    print(f"     r^{a}: 係数の符号が混在（t の次数の順の符号の変化 {sign_changes(cols[a])} 回） {dict(sorted(cols[a].items()))}")
                else:
                    sg, parts = res
                    print(f"     r^{a}: 符号 {sg}、" + " + ".join(f"{mult}·s^{k}[{L}]" if mult > 1 else f"s^{k}[{L}]" for k, L, mult in parts))
            # (3) D(1 − t)/(t − q)
            alt = sp.expand(sp.cancel(D * (1 - s**2) / (s**2 - r**2)))
            print(f"     項数：e {len(sp.Add.make_args(sp.expand(e)))}、D(1−t)/(t−q) {len(sp.Add.make_args(alt))}")
        # (2) 最低のブロックの閉じた形
        bmin = min(b for b in blocks if asym_rs(blocks[b]) != 0)
        D = asym_rs(blocks[bmin])
        tower = sum(t**i for i in range(2 * n)) + sum(t**i for i in range(2 * n - 1))
        pred = -(1 - w**2) * w**(-(8 * m_ + 1)) * t * tower
        pred_rs = sp.expand(pred.subs({w: r / s, t: s**2}))
        print(f"  最低のブロック ν^{bmin}: D = −(1 − q/t) w^{{−{8*m_+1}}} t ([{2*n}]_t + [{2*n-1}]_t) と一致: {sp.expand(D - pred_rs) == 0}", flush=True)
