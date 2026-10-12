"""4 箱の \\hat f の Laurent 性が rCS_{[2,1]} の自由度を削るか（D11 の「追加の仮定の候補」、D12 段階 3）。

[2,1] の自由度：3 箱の一次の条件（U, T, E, U2, H, L、𝐭 の窓は r₂r₁ の範囲）の零空間の基底 v_j（scripts/lmov21_linear）。
4 箱の \\hat f_ρ（ρ ⊢ 4）に [2,1] は C_{21}·C_1 の項で一次に入る：Δ\\hat f_ρ = Σ_j y_j n_{ρ,(21)(1)}·r₁·v_j / D_ρ。
Laurent 性 ⇔ 分子が D_ρ で割り切れる ⇔ Σ_j y_j rem_{D_ρ}(P_j) = 0（単一の多項式による剰余は線形）。
  known31：[3,1], [2,1,1], [2,2], [4], [1⁴] を既知（GN から補間）として、y の零空間の次元を数える。
  joint  ：[3,1]（と転置の [2,1,1]）も未知数（箱の係数）にして、連立の零空間を [2,1] の座標に射影した次元を数える。
使い方：PYTHONPATH=src:scripts python scripts/lmov4_joint.py known31 T2_5
記録：logs/lmov4_joint.log、data/lmovstruct/lmov4_joint.json。
"""
import json
import sys
import time
from pathlib import Path

import flint
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lmov21_linear as LL  # noqa: E402
from dgr_box3 import build21  # noqa: E402
from framing_check import load_reduced  # noqa: E402
from lmov4_formula import formula_n, inputs_n  # noqa: E402
from t4_colored import bold_to_small  # noqa: E402
from v4_closed import L  # noqa: E402
from rcsknots.partitions import partitions  # noqa: E402

LOG = Path("logs/lmov4_joint.log")
OUT = Path("data/lmovstruct/lmov4_joint.json")
TORUS = {"T2_3": (2, 3), "T2_5": (2, 5), "T2_7": (2, 7), "T3_4": (3, 4)}


def log(s):
    print(s, flush=True)
    with LOG.open("a") as f:
        f.write(s + "\n")


def reduced_4box(K):
    """|λ| ≤ 4 の reduced の rCS（小文字）。|λ| ≤ 3 は保存値、|λ| = 4 は GN から補間（キャッシュ data/lmovstruct/rcs4_<K>.json）。"""
    path = Path(f"data/lmovstruct/rcs4_{K}.json")
    small = load_reduced(K, 3)
    if path.exists():
        for k, v in json.loads(path.read_text()).items():
            small[tuple(map(int, k.split(",")))] = {tuple(map(int, kk.split(","))): c for kk, c in v.items()}
        return small
    from rcsknots.exact import reconstruct_points
    from rcsknots.extract import A_OFFSET, _arr
    from rcsknots.field import V
    from rcsknots.bpsgen import torus_rcs_factory
    m, n = TORUS[K]
    fac = torus_rcs_factory(m, n)
    new = {}
    for lam in partitions(4):
        def ev(xs, ys, p, lam=lam):
            k = len(xs)
            Qv, Tv = V(xs, p), V(ys, p)
            lo, cs = fac(lam, Qv, Tv)
            arr = np.zeros((k, 2 * A_OFFSET + 1), dtype=np.int64)
            for i, c in enumerate(cs):
                arr[:, lo + i + A_OFFSET] = _arr(c, k, p) % p
            return arr
        raw = reconstruct_points(ev, seed=7, log=lambda s: None)
        new[lam] = {(d - A_OFFSET, eq, et): c for (d, eq, et), c in raw.items()}
    path.write_text(json.dumps({",".join(map(str, k)): {",".join(map(str, kk)): c for kk, c in v.items()} for k, v in new.items()}))
    small.update(new)
    return small


def nullspace21(K):
    """[2,1] の一次の条件（U, T, E, U2, H, L、窓の規則）の零空間の整数基底（太字 {m: c}）。lmov21_polytope.exact_affine を使う。"""
    from lmov21_polytope import exact_affine
    LL.MARGIN = 0
    small, H21u, window, true21 = LL.known_inputs(K)
    lin, mons = build21(small, H21u, window, (1, 1, 1), conds=("U", "T", "E", "U2", "H"))
    rows = [(dict(r), b) for r, b in zip(lin.rows, lin.rhs)] + LL.laurent_rows(small, mons, lin.col)
    x0, basis = exact_affine(rows, len(mons))
    inv = {v: k for k, v in lin.col.items()}
    vecs = [{inv[c][1]: z for c, z in vec.items()} for vec in basis]
    return vecs, true21


def remainder_rows(polys_small, D):
    """polys_small：[{(e_A, e_Q, e_T): c}]（Laurent）。共通の単項式でずらして非負にし、D（L）で割った剰余を {単項式: {列: 係数}} にする。"""
    ctx = flint.fmpz_mpoly_ctx.get(("A", "Q", "T"), "lex")
    dd = D.to_dict()
    lo = [min(min(k[i] for k in P) for P in polys_small if P) for i in range(3)]
    lo = [min(lo[i], min(k[i] for k in dd)) for i in range(3)]
    Dp = ctx.from_dict({(k[0] - lo[0], k[1] - lo[1], k[2] - lo[2]): c for k, c in dd.items()})
    rows = {}
    for j, P in enumerate(polys_small):
        if not P:
            continue
        Pp = ctx.from_dict({(k[0] - lo[0], k[1] - lo[1], k[2] - lo[2]): c for k, c in P.items()})
        _, r = divmod(Pp, Dp)
        for e, c in r.to_dict().items():
            rows.setdefault(tuple(int(v) for v in e), {})[j] = int(c)
    return rows


def dense_nullspace_modp(rows, ncols, p=1048573, batch=2000):
    """法 p（< 2²⁰、float64 で厳密）の密な消去で零空間の基底を求める。行は {列: 係数}（整数）。
    numpy の行列積で、まとまりごとに現在のピボット行（RREF）で簡約してから、まとまりの中で消去する。"""
    import numpy as np
    piv_cols = []                                            # ピボット列（挿入順）
    B = np.zeros((0, ncols))                                  # RREF の行（ピボット列は単位ベクトル）
    for start in range(0, len(rows), batch):
        chunk = rows[start:start + batch]
        M = np.zeros((len(chunk), ncols))
        for i, (row, _) in enumerate(chunk):
            for c, v in row.items():
                M[i, c] = v % p
        if len(piv_cols):
            M = (M - (M[:, piv_cols] @ B) % p) % p             # 既存のピボットで簡約（M[:, piv] は簡約後 0）
        # まとまりの中の消去
        for i in range(len(chunk)):
            nz = np.nonzero(M[i])[0]
            if len(nz) == 0:
                continue
            c = int(nz[0])
            inv = pow(int(M[i, c]), p - 2, p)
            M[i] = (M[i] * inv) % p
            # 他の行（まとまり内と B）から列 c を消す
            col = M[:, c].copy(); col[i] = 0
            M = (M - np.outer(col, M[i])) % p
            if len(piv_cols):
                B = (B - np.outer(B[:, c], M[i])) % p
            B = np.vstack([B, M[i:i + 1]])
            piv_cols.append(c)
    pivset = set(piv_cols)
    free = [c for c in range(ncols) if c not in pivset]
    basis = []
    for fc in free:
        vec = {fc: 1}
        for r, pc in enumerate(piv_cols):
            v = int(B[r, fc])
            if v:
                vec[pc] = (-v) % p
        basis.append(vec)
    return basis


def known31(K):
    t0 = time.time()
    small = reduced_4box(K)
    vecs, true21 = nullspace21(K)
    log(f"{K}: [2,1] の零空間 {len(vecs)} 次元 ({time.time() - t0:.0f}s)")
    fo, keys = formula_n(4)
    R = inputs_n(small, 4)
    idx21 = keys.index(((2, 1), 1)); idx1 = keys.index(((1,), 1))
    allrows = []
    ncols = len(vecs)
    for rho, (ns, D) in fo.items():
        expo = tuple(1 if i in (idx21, idx1) else 0 for i in range(len(keys)))
        n = ns.get(expo)
        if n is None:
            log(f"  ρ={rho}: C_21·C_1 の項が無い")
            continue
        base = n * R[((1,), 1)]
        polys = [(base * L.from_dict(bold_to_small(v))).to_dict() for v in vecs]
        rows = remainder_rows(polys, D)
        nz = sum(1 for r in rows.values() if r)
        allrows += [(r, 0) for r in rows.values() if r]
        log(f"  ρ={rho}: 剰余の式 {nz} ({time.time() - t0:.0f}s)")
        # 正解の検算：\\hat f_ρ（正解）が Laurent であること
    rank = ncols - len(dense_nullspace_modp(allrows, ncols)) if ncols > 100 else LL.rank_nullity(allrows, ncols)[0]
    log(f"{K}: 4 箱の Laurent 性（[3,1] 既知）が [2,1] の零空間 {ncols} のうち殺す次元 {rank}、残り {ncols - rank} ({time.time() - t0:.0f}s)")
    return {"nullity21": ncols, "rank_4box": rank, "remaining": ncols - rank}


if __name__ == "__main__":
    cmd, K = sys.argv[1], sys.argv[2]
    log(f"=== lmov4_joint {time.strftime('%Y-%m-%d %H:%M:%S')} {' '.join(sys.argv[1:])} ===")
    if cmd == "known31":
        res = known31(K)
        prev = json.loads(OUT.read_text()) if OUT.exists() else {}
        prev[f"{K} known31"] = res
        OUT.write_text(json.dumps(prev, ensure_ascii=False, indent=1))


# ---- joint：[3,1]（と転置の [2,1,1]）も未知数 --------------------------------------------------------

def bold_of(P):
    from rcsknots.bps import to_bold
    return {k: v for k, v in to_bold(P, 1 if any(k[0] % 2 for k in P) else 0).items() if v}


def nullspace31(K, small, true31_bold):
    """[3,1] の一次の条件（t = q、EGP（q = 1：r₁₁r₁r₁、t = 1：r₃r₁）、U(2)（a = t²：r₂r₁₁）、hook（a = t/q：c³·A(𝐪⁴𝐭², 𝐭)））の
    零空間の整数基底（太字）。窓：(𝐚, 𝐪) は非 refined の台の箱、𝐭 は正解の範囲（検査用）。"""
    from dgr_d7_construct import Lin
    from dgr_box3 import _mul
    from lmov21_polytope import exact_affine
    from rcsknots.dgr import specialize
    b = {l: bold_of(P) for l, P in small.items() if sum(l) <= 3}
    H31u = {}
    for (a, q, t), c in true31_bold.items():
        H31u[(a, q)] = H31u.get((a, q), 0) + c * (-1) ** (t % 2)
    H31u = {k: v for k, v in H31u.items() if v}
    As = sorted({a for a, _ in H31u}); Qs = sorted({q for _, q in H31u})
    ts = sorted({t for _, _, t in true31_bold})
    mons = [(a, q, t) for a in range(As[0], As[-1] + 1, 2) for q in range(Qs[0], Qs[-1] + 1, 2) for t in range(ts[0], ts[-1] + 1)]
    lin = Lin()
    for m in mons:
        lin.var(("x", m))
    # U：𝐭 = −1 で非 refined
    for a in range(As[0], As[-1] + 1, 2):
        for q in range(Qs[0], Qs[-1] + 1, 2):
            lin.eq({("x", (a, q, t)): (-1) ** (t % 2) for t in range(ts[0], ts[-1] + 1)}, H31u.get((a, q), 0))

    def special(n, m, target):
        img = {}
        for mm in mons:
            for k, v in specialize({mm: 1}, n, m).items():
                img.setdefault(k, {})[("x", mm)] = v
        for k in set(img) | set(target):
            lin.eq(img.get(k, {}), target.get(k, 0))
    # E（q = 1）：build21 と同じ（a, t − q）の和 = r₁₁·r₁·r₁
    prod = _mul(_mul(b[(1, 1)], b[(1,)]), b[(1,)])
    target = {}
    for (a, q, t), v in prod.items():
        target[(a, t - q)] = target.get((a, t - q), 0) + v
    img = {}
    for m in mons:
        img.setdefault((m[0], m[2] - m[1]), {})[("x", m)] = 1
    for k in set(img) | set(target):
        lin.eq(img.get(k, {}), target.get(k, 0))
    # E'（t = 1 ⇔ 𝐪 = ±1）：Σ_q x (±1)^q = r₃·r₁ の同じ特殊化
    prod = _mul(b[(3,)], b[(1,)])
    for sgn in (1, -1):
        target = {}
        for (a, q, t), v in prod.items():
            target[(a, t)] = target.get((a, t), 0) + v * sgn ** (q % 2)
        img = {}
        for m in mons:
            img.setdefault((m[0], m[2]), {})[("x", m)] = sgn ** (m[1] % 2)
        for k in set(img) | set(target):
            lin.eq(img.get(k, {}), target.get(k, 0))
    # U2（a = t²）：r₂·r₁₁
    special(2, 0, _mul(specialize(b[(2,)], 2, 0), specialize(b[(1, 1)], 2, 0)))
    # hook（a = t/q）：c³·A(𝐪⁴𝐭², 𝐭)
    (cq, ct), cs = next(iter(specialize(b[(1,)], 1, 0).items()))
    A = specialize(b[(1,)], 1, 1)
    special(1, 1, {(4 * eq, et + 2 * eq + 3 * ct): cs ** 3 * v for (eq, et), v in A.items()})
    rows = [(dict(r), bb) for r, bb in zip(lin.rows, lin.rhs)]
    # 正解の検算
    xv = {lin.col[("x", m)]: true31_bold.get(m, 0) for m in mons}
    ok = all(sum(v * xv.get(c, 0) for c, v in row.items()) == bb for row, bb in rows)
    x0, basis = exact_affine(rows, len(mons))
    inv = {v: k for k, v in lin.col.items()}
    vecs = [{inv[c][1]: z for c, z in vec.items()} for vec in basis]
    return vecs, ok, len(mons)


def transpose_bold(P):
    return {(a, -q, t - q): c for (a, q, t), c in P.items()}


def joint(K):
    t0 = time.time()
    small = reduced_4box(K)
    vecs21, true21 = nullspace21(K)
    true31 = bold_of(small[(3, 1)])
    vecs31, ok31, n31 = nullspace31(K, small, true31)
    log(f"{K}: [2,1] の零空間 {len(vecs21)}、[3,1] の未知数 {n31}・零空間 {len(vecs31)}（正解が条件を満たす {ok31}）({time.time() - t0:.0f}s)")
    fo, keys = formula_n(4)
    R = inputs_n(small, 4)
    idx = {k: i for i, k in enumerate(keys)}
    e21 = tuple(1 if i in (idx[((2, 1), 1)], idx[((1,), 1)]) else 0 for i in range(len(keys)))
    e31 = tuple(1 if i == idx[((3, 1), 1)] else 0 for i in range(len(keys)))
    e211 = tuple(1 if i == idx[((2, 1, 1), 1)] else 0 for i in range(len(keys)))
    ncols = len(vecs21) + len(vecs31)
    import pickle
    cache = Path(f"data/lmovstruct/joint_rows_{K}.pkl")           # 再開用：ρ ごとの剰余の行を保存する
    cached = pickle.loads(cache.read_bytes()) if cache.exists() else {}
    allrows = []
    for rho, (ns, D) in fo.items():
        if str(rho) in cached:
            rows_list = cached[str(rho)]
        else:
            polys = []
            base21 = ns[e21] * R[((1,), 1)]
            for v in vecs21:
                polys.append((base21 * L.from_dict(bold_to_small(v))).to_dict())
            for v in vecs31:
                P = ns[e31] * L.from_dict(bold_to_small(v)) + ns[e211] * L.from_dict(bold_to_small(transpose_bold(v)))
                polys.append(P.to_dict())
            rows = remainder_rows(polys, D)
            rows_list = [r for r in rows.values() if r]
            cached[str(rho)] = rows_list
            cache.write_bytes(pickle.dumps(cached))
        allrows += [(r, 0) for r in rows_list]
        log(f"  ρ={rho}: 剰余の式 {len(rows_list)} ({time.time() - t0:.0f}s)")
    log(f"  消去：式 {len(allrows)}、列 {ncols}")
    basis = dense_nullspace_modp(allrows, ncols)
    # [2,1] の座標への射影の階数
    proj = [({c: v for c, v in vec.items() if c < len(vecs21)}, 0) for vec in basis]
    rank21, _ = LL.rank_nullity(proj, len(vecs21))
    log(f"{K}: 連立の零空間 {len(basis)}（[2,1] 20 ＋ [3,1] {len(vecs31)} 次元から）、[2,1] に残る自由度 {rank21} ({time.time() - t0:.0f}s)")
    return {"nullity21": len(vecs21), "n31": n31, "nullity31": len(vecs31), "joint_nullity": len(basis), "remaining21": rank21}


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "joint":
    K = sys.argv[2]
    res = joint(K)
    prev = json.loads(OUT.read_text()) if OUT.exists() else {}
    prev[f"{K} joint"] = res
    OUT.write_text(json.dumps(prev, ensure_ascii=False, indent=1))
