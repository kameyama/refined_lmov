r"""固定セクターの因数分解に Ĩ がテンソル積で働くとした E の整数計画（theory §12.10、progress.md 続き 74）。
  観察（trefoil の n = 1）：P = s²(t − ν²)(ν² − 1)(q − 1)(ν²q − 1)·N、§12.7 の一意な E = s²(t − ν²)(ν² + 1)(q + 1)(ν²q + 1)·Ñ
  （Ñ は N の係数の符号を付け替えたもの、全 8 状態が不変）。
  仮説：固定セクターの状態の空間は 2 状態の因子（境界のモード q、ν²、ν²q、t/ν²）と残りの部分 N の積で、
    状態の数は N の係数の絶対値（打ち消しは N の中では起きない）。Ĩ は因子ごとに働き（(t − ν²) は不変、(ν² − 1)、(q − 1)、
    (ν²q − 1) は和に替わる）、N の部分では e_m（|e_m| ≤ |c_m|、e_m ≡ c_m mod 2）。
  n = 1：E = s²(t − ν²)(ν² + 1)(q + 1)(ν²q + 1) Σ e_m m、条件 (i) 2(1 + t)g_(2)（t = −1）、(ii) −(1 + q)(ψ²g_(1) + 1)（t = q）。
    因子の捻りの 16 通り（各因子の ±）を全部試す。
  n = 2：P̂_λ = s(t − ν²)(ν² − 1)(ν²q − 1)·N_λ（境界のモードは hatted で除いてある、rk_f_twisted22_lam）、
    E_λ = s(t − ν²)(ν² + 1)(ν²q + 1) Σ e^λ_m m、条件は rk_f_twisted22_lam と同じ（i22, i4, ii）。因子が ν のブロックをまたぐので全体で解く。
使い方：python rk_f_twisted_tensor.py n1 | n2 [--all-inv] | n2series（N_λ の系列の同定）（--all-inv：N の全状態が不変、e = ±c）。"""
import sys, os, io, contextlib, pickle, itertools, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import numpy as np
import sympy as sp
from scipy.optimize import milp, LinearConstraint, Bounds
import rk_b_twisted_multiplet as TM

r, s, nu = TM.r, TM.s, TM.nu
t, q, I = s**2, r**2, sp.I
OUT = "data/refined_kauffman/daha"


def solve(mult, basis_fn, conds, all_inv=False):
    """TM.solve と同じ整数計画。all_inv なら e = ±c（全状態が不変）。一意でない変数も返す。"""
    if not all_inv:
        return TM.solve(mult, basis_fn, conds, uniq=True)
    keys = sorted(mult)
    basis = [basis_fn(k) for k in keys]
    A, b = [], []
    for target, sub, *w in conds:
        a_, b_ = TM.rows_for([x * w[0](k) if w else x for k, x in zip(keys, basis)], target, sub)
        A += a_; b += b_
    A, b = np.array(A), np.array(b)
    c = np.array([float(mult[k]) for k in keys]); n = len(keys)
    cons = LinearConstraint(A * (-2 * c), b - A @ c - 1e-9, b - A @ c + 1e-9)       # e = c(1 − 2z)、z ∈ {0, 1}
    opt = lambda cc: milp(cc, constraints=cons, bounds=Bounds(0, 1), integrality=np.ones(n))
    res = opt(np.zeros(n))
    if res.status != 0:
        return dict(feasible=False, n=n)
    free = []
    for i in range(n):
        cc = np.zeros(n); cc[i] = 1
        if round(opt(cc).fun) != round(-opt(-cc).fun):
            free.append(keys[i])
    z = np.round(res.x)
    return dict(feasible=True, n=n, free=free, e={k: int(c[i] * (1 - 2 * z[i])) for i, k in enumerate(keys)})


def n1(all_inv=False):
    import rk_b_knot_fixed_sector as K
    import rk_b_cover_split as C
    u = sp.Symbol("u", positive=True)
    out = {}
    for m in (1, 2):
        with contextlib.redirect_stdout(io.StringIO()):
            sec = C.cover_sectors(m)
        P = sp.cancel((1 - t) * sp.cancel(sec["x"].to_sympy().subs(u, nu * s)))
        N = sp.cancel(P / (s**2 * (t - nu**2) * (nu**2 - 1) * (q - 1) * (nu**2 * q - 1)))
        mult = TM.monomials(N)
        g = K.load_g(m)
        R = sp.cancel(2 * (1 + t) * g[(2,)])
        Eq = sp.cancel(-(1 + q) * (g[(1,)].subs({r: r**2, s: s**2, nu: nu**2}, simultaneous=True) + 1))
        print(f"T(2,{2 * m + 1})：N の単項式 {len(mult)}、状態の数 Σ|c| = {sum(abs(v) for v in mult.values())}、"
              f"指数の単項式（P/(s z_q)）{len(TM.monomials(sp.cancel(P / (s * (r - 1 / r)))))}", flush=True)
        for sg in itertools.product((1, -1), repeat=4):
            pre = s**2 * (t - sg[0] * nu**2) * (nu**2 + sg[1]) * (q + sg[2]) * (nu**2 * q + sg[3])
            res = solve(mult, lambda k: pre * TM.twisted(k, "mono"), [(R, {s: I}), (Eq, {s: r})], all_inv)
            if res["feasible"]:
                inv = sum(1 for k in mult if abs(res["e"][k]) == abs(mult[k]))
                print(f"  捻り (t {'−' if sg[0] > 0 else '+'} ν²)(ν² {'+' if sg[1] > 0 else '−'} 1)(q {'+' if sg[2] > 0 else '−'} 1)"
                      f"(ν²q {'+' if sg[3] > 0 else '−'} 1)：解あり、一意でない {len(res['free'])}/{res['n']}、|e| = |c| が {inv}", flush=True)
                out[(m, sg)] = res
        print("  （上に無い捻りは解なし）", flush=True)
    return out


def n2_setup():
    import rk_f_twisted22_lam as L
    bl = L.blocks()
    full = lambda name: sum((v * nu**k for k, v in bl[name].items()), sp.Integer(0))
    d = pickle.load(open(L.COVER, "rb"))
    mult = {}
    for lab in L.LAMS:
        X = sp.sympify(d["hat" + lab], locals={"r": r, "s": s, "nu": nu})
        for key, c in TM.monomials(sp.cancel(X / (s * (t - nu**2) * (nu**2 - 1) * (nu**2 * q - 1)))).items():
            mult[(lab,) + key] = c
    conds = {"i22": (full("R22"), {s: I}, lambda key: L.W22),
             "i4": (full("R4"), {s: I}, lambda key: L.CHI[key[0]][(2,)] * L.W4),
             "ii": (full("Eq"), {s: r}, lambda key: L.W22)}
    return mult, conds


def n2_scan():
    """n = 2：因子の捻り 8 通り × 条件の組で、解の有無（整数計画、uniq なし）と一次方程式の解の有無。"""
    import rk_f_twisted22_lam as L
    mult, conds = n2_setup()
    for sg in itertools.product((1, -1), repeat=3):
        pre = s * (t - sg[0] * nu**2) * (nu**2 + sg[1]) * (nu**2 * q + sg[2])
        basis = lambda key: pre * TM.twisted(key[1:], "mono")
        row = []
        for use in (("i22", "i4"), ("ii",), ("i22", "ii"), ("i22", "i4", "ii")):
            cs = [conds[u_] for u_ in use]
            res = TM.solve(mult, basis, cs, uniq=False)
            lin = L.linear_ok(mult, basis, cs)
            row.append(f"{'+'.join(use)}: 整数計画 {'○' if res['feasible'] else '×'} 一次 {'○' if lin else '×'}")
        print(f"捻り (t {'−' if sg[0] > 0 else '+'} ν²)(ν² {'+' if sg[1] > 0 else '−'} 1)(ν²q {'+' if sg[2] > 0 else '−'} 1)："
              + "、".join(row), flush=True)


def rows(mult, basis, conds, beta_odd_only=()):
    """条件の行。beta_odd_only に入る条件は、ν^{2β} の β が奇数（c = 1）の単項式の行だけを使う。"""
    keys = sorted(mult)
    bs = [basis(k) for k in keys]
    A, b = [], []
    for name, (target, sub, w) in conds.items():
        cols = [TM.lp(sp.expand((x * w(k)).subs(sub)), (r, s, nu)) for k, x in zip(keys, bs)]
        T_ = TM.lp(sp.expand(target.subs(sub)), (r, s, nu)) if target != 0 else {}
        for mo in set(T_).union(*cols):
            if name in beta_odd_only and (mo[2] // 2) % 2 == 0:
                continue
            for f in (sp.re, sp.im):
                row = [float(f(c.get(mo, 0))) for c in cols]
                v = float(f(T_.get(mo, 0)))
                if any(row) or v:
                    A.append(row); b.append(v)
    return keys, np.array(A), np.array(b)


def freedom(mult, A, b, uniq=True):
    """|e| ≤ |c|・偶奇の整数計画の解の有無、一次方程式の零空間の次元、（uniq なら）一意でない変数の数。"""
    keys = sorted(mult)
    c = np.array([float(mult[k]) for k in keys]); n = len(keys)
    cons = LinearConstraint(-2 * A, b - A @ c - 1e-9, b - A @ c + 1e-9)
    bnd = Bounds(np.minimum(0, c), np.maximum(0, c))
    opt = lambda cc: milp(cc, constraints=cons, bounds=bnd, integrality=np.ones(n))
    res = opt(np.zeros(n))
    nullity = n - np.linalg.matrix_rank(A) if len(A) else n
    if res.status != 0:
        return dict(feasible=False, n=n, nullity=nullity)
    free = None
    if uniq:
        free = []
        for i in range(n):
            cc = np.zeros(n); cc[i] = 1
            if round(opt(cc).fun) != round(-opt(-cc).fun):
                free.append(keys[i])
    y = np.round(res.x)
    return dict(feasible=True, n=n, nullity=nullity, free=free, e={k: int(c[i] - 2 * y[i]) for i, k in enumerate(keys)})


def n2_partial(uniq=True):
    """n = 2、捻り (t − ν²)(ν² + 1)(ν²q + 1)：(a) t = −1 だけ、(b) t = −1 ＋ (ii) の c = 1（β 奇）の部分。"""
    mult, conds = n2_setup()
    pre = s * (t - nu**2) * (nu**2 + 1) * (nu**2 * q + 1)
    basis = lambda key: pre * TM.twisted(key[1:], "mono")
    for lab, use, odd in (("t = −1 だけ（i22 + i4）", ("i22", "i4"), ()),
                          ("t = −1 ＋ (ii) の β 奇の部分", ("i22", "i4", "ii"), ("ii",))):
        t0 = time.time()
        keys, A, b = rows(mult, basis, {k: conds[k] for k in use}, odd)
        res = freedom(mult, A, b, uniq)
        msg = (f"解あり、一意でない {len(res['free'])}/{res['n']}" if res["feasible"] and res["free"] is not None
               else ("解あり" if res["feasible"] else "解なし"))
        print(f"n = 2 {lab}：{msg}、一次方程式の零空間の次元 {res['nullity']}（{time.time() - t0:.0f} 秒）", flush=True)


def qrule(mult, A, b, free_blocks=(), lab_of=lambda k: ""):
    """e_m = σ_{β} (−1)^{a/2} c_m（a は r の指数、σ は ν のブロック（と λ）ごとの符号）。free_blocks の ν のブロックは一般の e。
    σ と free_blocks の e を整数計画で解き、解の有無と σ の一意性を返す。"""
    keys = sorted(mult)
    grp = sorted({(lab_of(k), k[-1]) for k in keys})
    gi = {g: i for i, g in enumerate(grp)}
    rule = [k for k in keys if k[-1] not in free_blocks]
    fr = [k for k in keys if k[-1] in free_blocks]
    assert all(k[-3] % 2 == 0 for k in rule), "r の指数が奇数"
    col = {k: i for i, k in enumerate(keys)}
    # 変数：z_g ∈ {0,1}（σ_g = 1 − 2z_g）、自由な単項式の y_m（e = c − 2y）
    ng, nf = len(grp), len(fr)
    M = np.zeros((len(b), ng + nf)); rhs = b.copy()
    for k in rule:
        v = float(mult[k]) * (-1)**(k[-3] // 2)
        rhs = rhs - A[:, col[k]] * v
        M[:, gi[(lab_of(k), k[-1])]] += -2 * v * A[:, col[k]]
    for j, k in enumerate(fr):
        c = float(mult[k]); rhs = rhs - A[:, col[k]] * c
        M[:, ng + j] = -2 * A[:, col[k]]
    cf = np.array([float(mult[k]) for k in fr])
    lo = np.concatenate([np.zeros(ng), np.minimum(0, cf)]); hi = np.concatenate([np.ones(ng), np.maximum(0, cf)])
    cons = LinearConstraint(M, rhs - 1e-9, rhs + 1e-9)
    opt = lambda cc: milp(cc, constraints=cons, bounds=Bounds(lo, hi), integrality=np.ones(ng + nf))
    res = opt(np.zeros(ng + nf))
    if res.status != 0:
        return None
    sig = {}
    for g, i in gi.items():
        cc = np.zeros(ng + nf); cc[i] = 1
        mn, mx = round(opt(cc).fun), round(-opt(-cc).fun)
        sig[g] = "±" if mn != mx else ("+" if mn == 0 else "−")
    return sig


def lower_series(m=2):
    """T(2,2m+1) の N の下のブロック ν^{−(6+2k)}（k = 0..2m−1）は Sym^k V / r^{4+2k}、V = 1 + t + qt（基底 1, t, qt）。
    V の対角な対合 diag(ε₁, ε₂, ε₃) の Sym^k への誘導作用 × σ_k（ブロックごと）を下のブロックに課し、上のブロックは一般の e として
    (i)+(ii) を解く。σ_k と上のブロックの変数の一意性も調べる。"""
    import rk_b_knot_fixed_sector as K
    import rk_b_cover_split as C
    import rk_f_frozen_basis as FB
    u = sp.Symbol("u", positive=True)
    with contextlib.redirect_stdout(io.StringIO()):
        sec = C.cover_sectors(m)
    P = sp.cancel((1 - t) * sp.cancel(sec["x"].to_sympy().subs(u, nu * s)))
    N = sp.cancel(P / (s**2 * (t - nu**2) * (nu**2 - 1) * (q - 1) * (nu**2 * q - 1)))
    low = {-(6 + 2 * k): k for k in range(2 * m)}
    sym = lambda k: sp.expand(sum(sp.Integer(1) * t**(b + c) * q**c for a in range(k + 1) for b in range(k + 1 - a) for c in [k - a - b]))
    blocks_ = FB.nu_blocks(N)
    print(f"T(2,{2 * m + 1})：下のブロックが Sym^k V / r^(4+2k) に一致："
          + " ".join(f"ν^{b_}:{sp.expand(blocks_[b_] * r**(4 + 2 * k) - sym(k)) == 0}" for b_, k in low.items()), flush=True)
    mult = TM.monomials(N)
    g = K.load_g(m)
    conds = {"i": (sp.cancel(2 * (1 + t) * g[(2,)]), {s: I}, lambda k: 1),
             "ii": (sp.cancel(-(1 + q) * (g[(1,)].subs({r: r**2, s: s**2, nu: nu**2}, simultaneous=True) + 1)), {s: r}, lambda k: 1)}
    pre = s**2 * (t - nu**2) * (nu**2 + 1) * (q + 1) * (nu**2 * q + 1)
    keys, A, b = rows(mult, lambda k: pre * TM.twisted(k, "mono"), conds)
    col = {k: i for i, k in enumerate(keys)}
    for eps in ((1, 1, -1), (1, -1, 1), (1, -1, -1), (1, 1, 1)):
        rule = [(k, low[k[-1]]) for k in keys if k[-1] in low]
        free = [k for k in keys if k[-1] not in low]
        ng, nf = len(low), len(free)
        Mx = np.zeros((len(b), ng + nf)); rhs = b.copy()
        for k, kk in rule:
            qe = (k[1] + 4 + 2 * kk) // 2; te = k[2] // 2
            v = float(mult[k]) * eps[0]**(kk - te) * eps[1]**(te - qe) * eps[2]**qe
            rhs -= A[:, col[k]] * v; Mx[:, kk] += -2 * v * A[:, col[k]]
        for j, k in enumerate(free):
            rhs -= A[:, col[k]] * float(mult[k]); Mx[:, ng + j] = -2 * A[:, col[k]]
        cf = np.array([float(mult[k]) for k in free])
        bnd = Bounds(np.concatenate([np.zeros(ng), np.minimum(0, cf)]), np.concatenate([np.ones(ng), np.maximum(0, cf)]))
        cons = LinearConstraint(Mx, rhs - 1e-9, rhs + 1e-9)
        opt = lambda cc: milp(cc, constraints=cons, bounds=bnd, integrality=np.ones(ng + nf))
        if opt(np.zeros(ng + nf)).status != 0:
            print(f"  diag{eps}：解なし", flush=True); continue
        nfree = []
        for i in range(ng + nf):
            cc = np.zeros(ng + nf); cc[i] = 1
            if round(opt(cc).fun) != round(-opt(-cc).fun):
                nfree.append(i)
        print(f"  diag{eps}：解あり、σ_k の一意でない数 {sum(1 for i in nfree if i < ng)}/{ng}、"
              f"上のブロックの一意でない変数 {sum(1 for i in nfree if i >= ng)}/{nf}", flush=True)


def knot_rule(m=2, use=("i", "ii"), save=False):
    """T(2,2m+1) の n = 1：下の系列（ν^{−(6+2k)}、k < 2m、Sym^k V）と、k = 2m のブロックの正の部分（T(2,5) では t³V）に
    V の対合 diag(1, 1, −1) の誘導作用 e = σ_g (−1)^{q の指数} c を課し（σ_g は群ごと）、残りは一般の e で解く。
    一意なら E と、全単項式で「e = σ (−1)^q (c − 2p)」（p は Ĩ で入れ替わる対の数、σ はブロックと正負の部分ごと）の形かを表示する。"""
    import rk_b_knot_fixed_sector as K
    import rk_b_cover_split as C
    import rk_f_frozen_basis as FB
    u = sp.Symbol("u", positive=True)
    with contextlib.redirect_stdout(io.StringIO()):
        sec = C.cover_sectors(m)
    P = sp.cancel((1 - t) * sp.cancel(sec["x"].to_sympy().subs(u, nu * s)))
    N = sp.cancel(P / (s**2 * (t - nu**2) * (nu**2 - 1) * (q - 1) * (nu**2 * q - 1)))
    mult = TM.monomials(N)
    g = K.load_g(m)
    allc = {"i": (sp.cancel(2 * (1 + t) * g[(2,)]), {s: I}, lambda k: 1),
            "ii": (sp.cancel(-(1 + q) * (g[(1,)].subs({r: r**2, s: s**2, nu: nu**2}, simultaneous=True) + 1)), {s: r}, lambda k: 1)}
    pre = s**2 * (t - nu**2) * (nu**2 + 1) * (q + 1) * (nu**2 * q + 1)
    keys, A, b = rows(mult, lambda k: pre * TM.twisted(k, "mono"), {k_: allc[k_] for k_ in use})
    col = {k: i for i, k in enumerate(keys)}
    kk = lambda k: (-k[-1] - 6) // 2                                      # ブロックの番号 k
    qe = lambda k: (k[1] + 4 + 2 * kk(k)) // 2                             # q の指数（r^{−(4+2k)} で正規化）
    grp = lambda k: (kk(k), "+" if mult[k] > 0 else "−")                   # 符号の群：ブロックと正負の部分
    fixed = [k for k in keys if kk(k) < 2 * m or (kk(k) == 2 * m and mult[k] > 0)]
    free = [k for k in keys if k not in fixed]
    gl = sorted({grp(k) for k in fixed}); gi = {g_: i for i, g_ in enumerate(gl)}
    ng, nf = len(gl), len(free)
    Mx = np.zeros((len(b), ng + nf)); rhs = b.copy()
    for k in fixed:
        v = float(mult[k]) * (-1)**qe(k); rhs -= A[:, col[k]] * v; Mx[:, gi[grp(k)]] += -2 * v * A[:, col[k]]
    for j, k in enumerate(free):
        rhs -= A[:, col[k]] * float(mult[k]); Mx[:, ng + j] = -2 * A[:, col[k]]
    cf = np.array([float(mult[k]) for k in free])
    bnd = Bounds(np.concatenate([np.zeros(ng), np.minimum(0, cf)]), np.concatenate([np.ones(ng), np.maximum(0, cf)]))
    cons = LinearConstraint(Mx, rhs - 1e-9, rhs + 1e-9)
    opt = lambda cc: milp(cc, constraints=cons, bounds=bnd, integrality=np.ones(ng + nf))
    res = opt(np.zeros(ng + nf))
    tag = f"T(2,{2 * m + 1}) 条件 {'+'.join(use)}"
    if res.status != 0:
        print(f"{tag}：解なし", flush=True); return None
    nfree = []
    for i in range(ng + nf):
        cc = np.zeros(ng + nf); cc[i] = 1
        if round(opt(cc).fun) != round(-opt(-cc).fun): nfree.append(i)
    print(f"{tag}：解あり、一意でない 符号の群 {sum(1 for i in nfree if i < ng)}/{ng}、一般の変数 {sum(1 for i in nfree if i >= ng)}/{nf}", flush=True)
    if nfree:
        return None
    x = np.round(res.x)
    e = {}
    for k in fixed:
        e[k] = int(mult[k] * (-1)**qe(k) * (1 - 2 * x[gi[grp(k)]]))
    for j, k in enumerate(free):
        e[k] = int(mult[k] - 2 * x[ng + j])
    # 符号の規則の確認：群ごとに σ = e (−1)^q / c の符号がそろうか（|e| < |c| は入れ替わる対を含む）
    sig = {}
    for k in keys:
        if e[k] != 0:
            sig.setdefault(grp(k), set()).add(int(np.sign(e[k] * (-1)**qe(k) / mult[k])))
    pairs = sum((abs(int(mult[k])) - abs(e[k])) // 2 for k in keys)
    print(f"   全単項式で σ_g = sign(e (−1)^q / c) が群ごとにそろう：{all(len(v) == 1 for v in sig.values())}"
          f"（σ：{', '.join(f'k={g_[0]}{g_[1]}:{list(v)[0]:+d}' for g_, v in sorted(sig.items()) if len(v) == 1)}）、"
          f"入れ替わる対の数 Σ(|c| − |e|)/2 = {pairs}（その単項式：{[(kk(k), qe(k), k[2] // 2, int(mult[k]), e[k]) for k in keys if abs(e[k]) < abs(mult[k])]}）", flush=True)
    for name, ex in (("(−1)^q", qe), ("(−1)^t", lambda k: k[2] // 2)):
        sg = {}
        for k in keys:
            if e[k] != 0:
                sg.setdefault(grp(k), set()).add(int(np.sign(e[k] * (-1)**ex(k) / mult[k])))
        print(f"   {name} で σ がそろう群：{[f'k={g_[0]}{g_[1]}' for g_, v in sorted(sg.items()) if len(v) == 1]}、"
              f"そろわない群：{[f'k={g_[0]}{g_[1]}' for g_, v in sorted(sg.items()) if len(v) > 1]}", flush=True)
    E = sp.expand(pre * sum(v * TM.twisted(k, "mono") for k, v in e.items()))
    if save:
        pickle.dump({"E": str(E), "e": {str(k): v for k, v in e.items()}}, open(f"{OUT}/twisted_E_n1_m{m}.pkl", "wb"))
        print(f"   保存：{OUT}/twisted_E_n1_m{m}.pkl", flush=True)
    return E


def full_rule(m=2, use=("i",)):
    """全ブロックに符号の規則を課す：e_m = σ_g (−1)^{d_m} sign(c_m) w_m、0 ≤ w_m ≤ |c_m|、w_m ≡ c_m (mod 2)
    （|c| − w は Ĩ で入れ替わる対）。d_m は k ≤ 2m で q の指数、k > 2m で t の指数。σ_g は群（ブロックと c の符号）ごと。
    σ_g w_m の積は z_g ∈ {0,1} と v_m = z_g w_m（McCormick の不等式）で線形にする。"""
    import rk_b_knot_fixed_sector as K
    import rk_b_cover_split as C
    u = sp.Symbol("u", positive=True)
    with contextlib.redirect_stdout(io.StringIO()):
        sec = C.cover_sectors(m)
    P = sp.cancel((1 - t) * sp.cancel(sec["x"].to_sympy().subs(u, nu * s)))
    N = sp.cancel(P / (s**2 * (t - nu**2) * (nu**2 - 1) * (q - 1) * (nu**2 * q - 1)))
    mult = TM.monomials(N)
    g = K.load_g(m)
    allc = {"i": (sp.cancel(2 * (1 + t) * g[(2,)]), {s: I}, lambda k: 1),
            "ii": (sp.cancel(-(1 + q) * (g[(1,)].subs({r: r**2, s: s**2, nu: nu**2}, simultaneous=True) + 1)), {s: r}, lambda k: 1)}
    pre = s**2 * (t - nu**2) * (nu**2 + 1) * (q + 1) * (nu**2 * q + 1)
    keys, A, b = rows(mult, lambda k: pre * TM.twisted(k, "mono"), {k_: allc[k_] for k_ in use})
    kk = lambda k: (-k[-1] - 6) // 2
    d = lambda k: (k[1] + 4 + 2 * kk(k)) // 2 if kk(k) <= 2 * m else k[2] // 2
    grp = lambda k: (kk(k), int(np.sign(float(mult[k]))))
    res = rule_milp(keys, mult, A, b, grp, d)
    tag = f"T(2,{2 * m + 1}) 全ブロックに規則、条件 {'+'.join(use)}"
    if res is None:
        print(f"{tag}：解なし", flush=True); return
    print(f"{tag}：解あり、e が一意でない単項式 {res['nfree']}/{len(keys)}、入れ替わる対の数 {res['pairs']}（最小 {res['min_pairs']}）", flush=True)


def rule_milp(keys, mult, A, b, grp, d, pairs=True, uniq=True, minp=True, max_pairs=None, ranges=None, samples=None, flip=None):
    """符号の規則 e_m = σ_g (−1)^{d_m} sign(c_m) w_m、0 ≤ w_m ≤ |c_m|、w_m ≡ c_m (mod 2) の下で A e = b を解く
    （|c| − w は Ĩ で入れ替わる対、pairs=False なら w = |c|）。σ_g は群 grp(k) ごと。
    σ_g w_m の積は z_g ∈ {0,1} と v_m = z_g w_m（McCormick の不等式）で線形にする。
    解が無ければ None、あれば e、（uniq なら）e が一意でない単項式、入れ替わる対の数、（minp なら）対の数の最小。
    max_pairs：Σ y ≤ max_pairs を加える。ranges（単項式の列）：その各単項式の e の最小と最大 {単項式: (最小, 最大)} だけを返す。
    samples = (個数, 乱数の種)：Σ r_m e_m（r_m は ±1 の乱数）を最小にする解を個数だけ求め、e の列を返す（解の集合の標本）。
    flip = {単項式: ±1}：Σ flip_m e_m を最小にする解（flip_m = sign(e0_m) なら、なるべく多くの符号を反転させる）の (最小値, e) を返す。"""
    n = len(keys)
    gl = sorted({grp(k) for k in keys}, key=str); gi = {g_: i for i, g_ in enumerate(gl)}; ng = len(gl)
    sgn = np.array([float(np.sign(float(mult[k]))) * (-1)**d(k) for k in keys])
    cab = np.array([abs(float(mult[k])) for k in keys])
    # 変数 [y (n), v (n), z (ng)]：w = |c| − 2y、e = sgn (w − 2v)
    nv = 2 * n + ng
    # 等式：A e = b → A diag(sgn) (|c| − 2y − 2v) = b
    AS = A * sgn
    Aeq = np.hstack([-2 * AS, -2 * AS, np.zeros((len(b), ng))]); beq = b - AS @ cab
    cons = [LinearConstraint(Aeq, beq - 1e-9, beq + 1e-9)]
    # McCormick：v ≤ w = |c| − 2y、v ≤ |c| z、v ≥ w − |c|(1 − z)
    Ms = []
    for i, k in enumerate(keys):
        gj = 2 * n + gi[grp(k)]
        r1 = np.zeros(nv); r1[n + i] = 1; r1[i] = 2; Ms.append((r1, -np.inf, cab[i]))            # v + 2y ≤ |c|
        r2 = np.zeros(nv); r2[n + i] = 1; r2[gj] = -cab[i]; Ms.append((r2, -np.inf, 0))          # v − |c| z ≤ 0
        r3 = np.zeros(nv); r3[n + i] = 1; r3[i] = 2; r3[gj] = -cab[i]; Ms.append((r3, 0, np.inf))  # v + 2y − |c| z ≥ 0
    if max_pairs is not None:
        r4 = np.zeros(nv); r4[:n] = 1; Ms.append((r4, -np.inf, max_pairs))
    cons.append(LinearConstraint(np.array([x[0] for x in Ms]), [x[1] for x in Ms], [x[2] for x in Ms]))
    bnd = Bounds(np.zeros(nv), np.concatenate([np.floor(cab / 2) if pairs else np.zeros(n), cab, np.ones(ng)]))
    opt = lambda cc: milp(cc, constraints=cons, bounds=bnd, integrality=np.ones(nv))
    if ranges is not None:
        col = {k: i for i, k in enumerate(keys)}; out = {}
        for k in ranges:
            i = col[k]; cc = np.zeros(nv); cc[i] = -2 * sgn[i]; cc[n + i] = -2 * sgn[i]     # e_i − sgn|c| の係数
            lo_, hi_ = opt(cc), opt(-cc)
            if lo_.status != 0 or hi_.status != 0:
                return None
            out[k] = (round(sgn[i] * cab[i] + lo_.fun), round(sgn[i] * cab[i] - hi_.fun))
        return out
    if flip is not None:
        col = {k: i for i, k in enumerate(keys)}; cc = np.zeros(nv); const = 0.0
        for k, f in flip.items():
            i = col[k]; cc[i] = -2 * sgn[i] * f; cc[n + i] = -2 * sgn[i] * f; const += f * sgn[i] * cab[i]
        res = opt(cc)
        if res.status != 0:
            return None
        x = res.x
        return round(const + res.fun), {k: int(round(sgn[i] * (cab[i] - 2 * x[i] - 2 * x[n + i]))) for i, k in enumerate(keys)}
    if samples is not None:
        rng = np.random.default_rng(samples[1]); out = []
        for _ in range(samples[0]):
            rr = rng.choice([-1.0, 1.0], n); cc = np.zeros(nv); cc[:n] = -2 * sgn * rr; cc[n:2 * n] = -2 * sgn * rr
            x = opt(cc).x
            out.append({k: int(round(sgn[i] * (cab[i] - 2 * x[i] - 2 * x[n + i]))) for i, k in enumerate(keys)})
        return out
    res = opt(np.zeros(nv))
    if res.status != 0:
        return None
    e0 = np.round(sgn * (cab - 2 * res.x[:n] - 2 * res.x[n:2 * n]))
    free = None
    if uniq:
        free = []
        for i in range(n):           # e_i の最小と最大（e は y, v の線形関数）
            cc = np.zeros(nv); cc[i] = -2 * sgn[i]; cc[n + i] = -2 * sgn[i]
            if round(opt(cc).fun) != round(-opt(-cc).fun): free.append(keys[i])
    cc = np.zeros(nv); cc[:n] = 1                                     # 入れ替わる対の数 Σ y の最小
    return dict(e={k: int(e0[i]) for i, k in enumerate(keys)}, free=free, nfree=None if free is None else len(free),
                pairs=int(sum(cab - np.abs(e0)) // 2), min_pairs=round(opt(cc).fun) if minp else None)


def n2_rule(uniq=True, opts=("q", "t"), uniform=False):
    """n = 2（trefoil）の N_λ に n = 1 の符号の規則を課し、t = −1 の条件（i22 + i4）だけで解く（progress.md 続き 76）。
    部分（λ、ν のブロック、c の符号）ごとに e = σ (−1)^{q の指数} sign(c) w か σ (−1)^{t の指数} sign(c) w を選ぶ
    （w ≤ |c| は入れ替わる対を許す）。(1) 部分ごとに一つだけ規則を課して（他の単項式は一般の e）q・t の可否を調べ、
    (2) 一方だけ可の部分はそれに固定し、両方可の部分は組を全部試し、(3) 全部分に規則を課した解の一意性と対の数の最小を調べる。
    opts は部分ごとの規則の候補：q（(−1)^a）、t（(−1)^b）、qt（(−1)^{a+b}）、0（すべて同符号）、a, b は q, t の指数。
    塔のセル (q, t, qt) の対角な対合 ε_qt = ε_q ε_t の誘導作用はこの 4 通りで、鏡映 (q, t) → (t^{−1}, q^{−1}) で対称なのは 0 と qt
    （theory §12.14、progress.md 続き 78）。uniform なら全部分に同じ規則を課した結果も表示する。"""
    mult, conds = n2_setup()
    pre = s * (t - nu**2) * (nu**2 + 1) * (nu**2 * q + 1)
    keys, A, b = rows(mult, lambda key: pre * TM.twisted(key[1:], "mono"), {k: conds[k] for k in ("i22", "i4")})
    part = lambda k: (k[0], k[-1], int(np.sign(float(mult[k]))))
    parts = sorted({part(k) for k in keys}, key=lambda x: (x[0], -x[1], -x[2]))
    exps = {"q": lambda k: k[2] // 2, "t": lambda k: k[3] // 2, "qt": lambda k: k[2] // 2 + k[3] // 2, "0": lambda k: 0}

    def run(assign, uniq_=False, minp=False):
        grp = lambda k: part(k) if part(k) in assign else k               # 規則の無い単項式は自分だけの群（一般の e と同じ）
        d = lambda k: exps[assign[part(k)]](k) if part(k) in assign else 0
        return rule_milp(keys, mult, A, b, grp, d, True, uniq_, minp)

    base = run({}, minp=True)
    print(f"規則なし（t = −1 だけ）：入れ替わる対の数の最小 {base['min_pairs']}（状態の数 {sum(abs(int(v)) for v in mult.values())}）", flush=True)
    if uniform:
        for x in opts:
            res = run({p_: x for p_ in parts}, minp=True)
            print(f"全部分に規則 {x}：" + ("解なし" if res is None else f"解あり、入れ替わる対の数の最小 {res['min_pairs']}"), flush=True)
    ok = {}
    for p_ in parts:
        n_ = sum(1 for k in keys if part(k) == p_)
        ok[p_] = [opts[0]] if n_ == 1 else [x for x in opts if run({p_: x}) is not None]   # 単項式 1 個では規則は σ に吸収される
        print(f"λ = {p_[0]}、ν^{p_[1]}、c の符号 {p_[2]:+d}（単項式 {n_}）：規則 {'・'.join(ok[p_]) or 'なし'}", flush=True)
    if any(not v for v in ok.values()):
        print("規則を課せない部分がある", flush=True); return None
    both = [p_ for p_ in parts if len(ok[p_]) >= 2]
    ncomb = int(np.prod([len(ok[p_]) for p_ in both]))
    if ncomb > 4096:
        print(f"複数可の部分 {len(both)} 個の組 {ncomb} 通り：多すぎるので組は調べない", flush=True); return None
    feas = []
    for ch in itertools.product(*[ok[p_] for p_ in both]):
        assign = {p_: ok[p_][0] for p_ in parts if len(ok[p_]) == 1} | dict(zip(both, ch))
        if run(assign) is not None:
            feas.append((ch, assign))
    print(f"複数可の部分 {len(both)} 個の組 {ncomb} 通りのうち解ありは {len(feas)} 通り", flush=True)
    for i, p_ in enumerate(both):
        print(f"   λ = {p_[0]}、ν^{p_[1]}、{p_[2]:+d}：解ありの組での選択 {sorted({ch[i] for ch, _ in feas})}", flush=True)
    out = []
    for ch, assign in feas:
        res = run(assign, False, True)
        out.append((res["min_pairs"], ch, assign))
    for mp in sorted({x[0] for x in out}):
        print(f"   入れ替わる対の数の最小 {mp}：{sum(1 for x in out if x[0] == mp)} 通り（例 {'/'.join(next(x[1] for x in out if x[0] == mp))}）", flush=True)
    if uniq and out:
        mp, ch, assign = min(out, key=lambda x: x[0])
        t0 = time.time()
        res = run(assign, True)
        print(f"対の数の最小が最も小さい組 {'/'.join(ch)} の一意性：e が一意でない単項式 {res['nfree']}/{len(keys)}（{time.time() - t0:.0f} 秒）", flush=True)
        for bk in sorted({(k[0], k[-1]) for k in keys}, key=lambda x: (x[0], -x[1])):
            n_ = sum(1 for k in keys if (k[0], k[-1]) == bk)
            print(f"   λ = {bk[0]}、ν^{bk[1]}：一意でない {sum(1 for k in res['free'] if (k[0], k[-1]) == bk)}/{n_}", flush=True)
    return out


def choice_milp(keys, mult, A, b, part, pats, minimize_pairs=True, max_pairs=None, ranges=None, log=None):
    """部分ごとに符号の型 π（pats[部分] の ±1 の辞書の列）を一つ選び、e_m = π_m sign(c_m) w_m、w_m = |c_m| − 2y_m で A e = b を解く。
    v_{m,π} = u_{p,π} w_m（McCormick）で線形にする。入れ替わる対の数 Σ y を最小にし、(最小値, 選んだ型の番号, e) を返す。解なしは None。
    max_pairs：Σ y ≤ max_pairs を加える。ranges（単項式の列）：その各単項式について e の最小と最大を求め、{単項式: (最小, 最大)} を返す。"""
    n = len(keys)
    cab = np.array([abs(float(mult[k])) for k in keys]); sc = np.array([float(np.sign(float(mult[k]))) for k in keys])
    plist = sorted(pats, key=str)
    uidx, i0 = {}, n
    for p_ in plist:
        for j in range(len(pats[p_])):
            uidx[(p_, j)] = i0; i0 += 1
    vidx = {}
    for i, k in enumerate(keys):
        for j in range(len(pats[part(k)])):
            vidx[(i, j)] = i0; i0 += 1
    nv = i0
    Aeq = np.zeros((len(b), nv))
    for (i, j), vi in vidx.items():
        Aeq[:, vi] = A[:, i] * sc[i] * pats[part(keys[i])][j][keys[i]]
    rows_, lo, hi = [], [], []
    def row(coef, l, h):
        x = np.zeros(nv)
        for idx, c in coef: x[idx] += c
        rows_.append(x); lo.append(l); hi.append(h)
    for p_ in plist:
        row([(uidx[(p_, j)], 1) for j in range(len(pats[p_]))], 1, 1)
    for (i, j), vi in vidx.items():
        u = uidx[(part(keys[i]), j)]
        row([(vi, 1), (u, -cab[i])], -np.inf, 0)                      # v ≤ |c| u
        row([(vi, 1), (i, 2)], -np.inf, cab[i])                        # v ≤ w = |c| − 2y
        row([(vi, 1), (i, 2), (u, -cab[i])], 0, np.inf)                # v ≥ w − |c|(1 − u)
    if max_pairs is not None:
        row([(i, 1) for i in range(n)], -np.inf, max_pairs)
    cons = [LinearConstraint(Aeq, b - 1e-9, b + 1e-9), LinearConstraint(np.array(rows_), lo, hi)]
    ub = np.zeros(nv); ub[:n] = np.floor(cab / 2)
    for k_, vi in uidx.items(): ub[vi] = 1
    for (i, j), vi in vidx.items(): ub[vi] = cab[i]
    opt = lambda cc: milp(cc, constraints=cons, bounds=Bounds(np.zeros(nv), ub), integrality=np.ones(nv))
    if ranges is not None:
        col = {k: i for i, k in enumerate(keys)}
        out = {}
        for c_, k in enumerate(ranges):
            i = col[k]; cc = np.zeros(nv)
            for j in range(len(pats[part(k)])):
                cc[vidx[(i, j)]] = sc[i] * pats[part(k)][j][k]               # e_k = Σ_π sign(c) π v
            lo_, hi_ = opt(cc), opt(-cc)
            if lo_.status != 0 or hi_.status != 0:
                return None
            out[k] = (round(lo_.fun), round(-hi_.fun))
            if log and (c_ + 1) % 100 == 0:
                log(f"   e の範囲：{c_ + 1}/{len(ranges)}")
        return out
    cc = np.zeros(nv)
    if minimize_pairs: cc[:n] = 1
    res = opt(cc)
    if res.status != 0:
        return None
    x = np.round(res.x)
    chosen = {p_: next(j for j in range(len(pats[p_])) if x[uidx[(p_, j)]] == 1) for p_ in plist}
    e = {k: int(sc[i] * sum(pats[part(k)][j][k] * x[vidx[(i, j)]] for j in range(len(pats[part(k)])))) for i, k in enumerate(keys)}
    return int(round(sum(x[:n]))), chosen, e


def n2_choice(opts=("0", "q", "t", "qt"), by="sign"):
    """n = 2（trefoil）、t = −1 だけ：部分（λ、ν のブロック、c の符号）ごとに規則 opts のどれかと σ を選ぶ整数計画で、
    入れ替わる対の数の最小と、そのときの選択を求める（theory §12.14、progress.md 続き 78）。
    by = "sign"：部分は（λ、ブロック、c の符号）。by = "stripe"：（λ、ブロック、t の指数 − q の指数）で、qt 方向の弦（縞）ごと。"""
    mult, conds = n2_setup()
    pre = s * (t - nu**2) * (nu**2 + 1) * (nu**2 * q + 1)
    keys, A, b = rows(mult, lambda key: pre * TM.twisted(key[1:], "mono"), {k: conds[k] for k in ("i22", "i4")})
    if by == "sign":
        part = lambda k: (k[0], k[-1], int(np.sign(float(mult[k]))))
    else:
        part = lambda k: (k[0], k[-1], (k[3] - k[2]) // 2)
    exps = {"q": lambda k: k[2] // 2, "t": lambda k: k[3] // 2, "qt": lambda k: k[2] // 2 + k[3] // 2, "0": lambda k: 0}
    pats, names = {}, {}
    for p_ in {part(k) for k in keys}:
        ks = [k for k in keys if part(k) == p_]
        seen = {}
        for o in opts:
            for sg in (1, -1):
                pat = tuple(sg * (-1)**exps[o](k) for k in ks)
                seen.setdefault(pat, []).append(("+" if sg > 0 else "−") + o)
        pats[p_] = [dict(zip(ks, pat)) for pat in seen]
        names[p_] = ["/".join(v) for v in seen.values()]
    t0 = time.time()
    res = choice_milp(keys, mult, A, b, part, pats)
    if res is None:
        print(f"規則 {','.join(opts)}、部分 {by}：解なし（{time.time() - t0:.0f} 秒）", flush=True); return None
    mp, chosen, e = res
    print(f"規則 {','.join(opts)}、部分 {by}（{len(pats)} 個）：入れ替わる対の数の最小 {mp}（{time.time() - t0:.0f} 秒）", flush=True)
    if by != "sign":
        return res
    for p_ in sorted(chosen, key=lambda x: (x[0], -x[1], -x[2])):
        ks = [k for k in keys if part(k) == p_]
        pr = sum(abs(int(mult[k])) - abs(e[k]) for k in ks) // 2
        print(f"   λ = {p_[0]}、ν^{p_[1]}、{p_[2]:+d}（単項式 {len(ks)}）：選択 {names[p_][chosen[p_]]}、対 {pr}", flush=True)
    return res


def _n2_model():
    mult, conds = n2_setup()
    pre = s * (t - nu**2) * (nu**2 + 1) * (nu**2 * q + 1)
    keys, A, b = rows(mult, lambda key: pre * TM.twisted(key[1:], "mono"), {k: conds[k] for k in ("i22", "i4")})
    return mult, keys, A, b


def _n2_sample(arg):
    mp, k_, seed = arg
    mult, keys, A, b = _n2_model()
    return rule_milp(keys, mult, A, b, lambda k: k, lambda k: 0, max_pairs=mp, samples=(k_, seed))


def _n2_range(arg):
    mp, key = arg
    mult, keys, A, b = _n2_model()
    return rule_milp(keys, mult, A, b, lambda k: k, lambda k: 0, max_pairs=mp, ranges=[key])


def n2_pattern():
    """n = 2（trefoil）、t = −1 だけ、規則なし：入れ替わる対の数が最小の解の集合で、各単項式の e の範囲を求め（data に保存して再開できる）、
    符号の型 e/c を (q, t) の格子で表示し、一意に決まる単項式で規則 0・q・t・q+t が部分ごとに成り立つかを調べる（theory §12.15、progress.md 続き 79）。"""
    path = f"{OUT}/twisted_n2_pattern.pkl"
    mult, keys, A, b = _n2_model()
    pats = {k: [{k: 1}, {k: -1}] for k in keys}
    st = pickle.load(open(path, "rb")) if os.path.exists(path) else {}
    if "mp" not in st:
        t0 = time.time()
        st["mp"], _, st["e"] = choice_milp(keys, mult, A, b, lambda k: k, pats)
        pickle.dump(st, open(path, "wb"))
        print(f"規則なしの対の数の最小 {st['mp']}（{time.time() - t0:.0f} 秒）", flush=True)
    if "samples" not in st:                                              # (1) 標本：解の間で e が違えば一意でない
        import multiprocessing as mp_                                    # ProcessPoolExecutor はサンドボックスで使えない
        t0 = time.time()
        with mp_.get_context("fork").Pool(8) as ex:
            got = ex.map(_n2_sample, [(st["mp"], 2, seed) for seed in range(8)])
        st["samples"] = [st["e"]] + [x for g_ in got for x in g_]
        pickle.dump(st, open(path, "wb"))
        print(f"標本の解 {len(st['samples'])} 個（{time.time() - t0:.0f} 秒）", flush=True)
    var = {k for k in keys if len({e_[k] for e_ in st["samples"]}) > 1}
    # (2) 行空間：単位ベクトル e_m が A の行空間にあれば、A e = b だけで e_m が決まる
    _, sv, Vt = np.linalg.svd(A, full_matrices=True)
    rank = int((sv > 1e-8 * sv[0]).sum()); Nul = Vt[rank:].T
    lin = {k for i, k in enumerate(keys) if np.abs(Nul[i]).max() < 1e-8}
    print(f"A の階数 {rank}（単項式 {len(keys)}）、標本で e が変わる単項式 {len(var)}、一次方程式だけで決まる単項式 {len(lin)}", flush=True)
    st.setdefault("rg", {})
    for k in lin:
        st["rg"][k] = (st["e"][k], st["e"][k])
    for k in var:
        vs = [e_[k] for e_ in st["samples"]]
        st["rg"].setdefault(k, (min(vs), max(vs)))                       # 一意でないことだけが確か（範囲は下限）
    # (3) まとめて反転：未判定の集合 U で Σ_{m∈U} sign(e0_m) e_m を最小にする。反転した単項式は一意でない。
    #     反転が無く、U の全単項式で |e0| = |c| かつ最小値 = Σ|c| なら、各項が上限にあるので U は全部一意（証明）。
    e0 = st["e"]
    U = [k for k in keys if k not in st["rg"] and abs(e0[k]) == abs(int(mult[k]))]
    t0 = time.time()
    while U:
        val, e1 = rule_milp(keys, mult, A, b, lambda k: k, lambda k: 0, max_pairs=st["mp"], flip={k: int(np.sign(e0[k])) for k in U})
        moved = [k for k in U if e1[k] != e0[k]]
        print(f"   まとめて反転：U {len(U)} 個、反転 {len(moved)}、最小値 {val}、Σ|c| {sum(abs(int(mult[k])) for k in U)}（{time.time() - t0:.0f} 秒）", flush=True)
        for k in moved:
            st["rg"][k] = (min(e0[k], e1[k]), max(e0[k], e1[k]))
        pickle.dump(st, open(path, "wb"))
        if not moved:
            assert val == sum(abs(int(mult[k])) for k in U)
            for k in U:
                st["rg"][k] = (e0[k], e0[k])
            break
        U = [k for k in U if k not in moved]
    todo = [k for k in keys if k not in st["rg"]]                        # |e0| < |c| の単項式：1 個ずつ正確に
    print(f"1 個ずつの判定が要る単項式 {len(todo)}", flush=True)
    for k in todo:
        st["rg"].update(rule_milp(keys, mult, A, b, lambda k: k, lambda k: 0, max_pairs=st["mp"], ranges=[k]))
        pickle.dump(st, open(path, "wb"))
    mp, e, rg = st["mp"], st["e"], st["rg"]
    fixed = [k for k in keys if rg[k][0] == rg[k][1]]
    print(f"規則なし、対の数 {mp} の解：e が一意な単項式 {len(fixed)}/{len(keys)}、"
          f"そのうち e = c が {sum(1 for k in fixed if rg[k][0] == mult[k])}、e = −c が {sum(1 for k in fixed if rg[k][0] == -mult[k])}", flush=True)
    # 格子：行 = q の指数、列 = t の指数。記号は e/c の符号（+：e = c、−：e = −c、小文字 p/m：|e| < |c| で同符号/逆符号、0：e = 0）。
    # 一意でない単項式は記号の後に ?（最小の解の e で表示）
    sym = lambda k, v: "0" if v == 0 else ("+" if v == mult[k] else "−" if v == -mult[k] else ("p" if v * mult[k] > 0 else "m"))
    part = lambda k: (k[0], k[-1], int(np.sign(float(mult[k]))))
    for bk in sorted({(k[0], k[-1]) for k in keys}, key=lambda x: (x[0], -x[1])):
        ks = [k for k in keys if (k[0], k[-1]) == bk]
        a0 = min(k[2] for k in ks); b0 = min(k[3] for k in ks)
        g = {((k[2] - a0) // 2, (k[3] - b0) // 2): k for k in ks}
        A_ = max(a for a, _ in g); B_ = max(b_ for _, b_ in g)
        nf = sum(1 for k in ks if k in fixed)
        print(f"λ = {bk[0]}、ν^{bk[1]}（r^{a0} s^{b0}、一意 {nf}/{len(ks)}）", flush=True)
        for a in range(A_ + 1):
            print("   " + " ".join((sym(g[(a, b_)], e[g[(a, b_)]]) + ("" if g[(a, b_)] in fixed else "?")).ljust(2)
                                   if (a, b_) in g else ". " for b_ in range(B_ + 1)), flush=True)
    # 一意な単項式で、部分ごとに e/c の符号が σ(−1)^d の形か
    exps = {"0": lambda k: 0, "q": lambda k: k[2] // 2, "t": lambda k: k[3] // 2, "qt": lambda k: k[2] // 2 + k[3] // 2}
    print("一意で e ≠ 0 の単項式だけで、部分（λ、ブロック、c の符号）ごとに e/c の符号 = σ(−1)^d が成り立つ規則 d", flush=True)
    for p_ in sorted({part(k) for k in keys}, key=lambda x: (x[0], -x[1], -x[2])):
        ks = [k for k in fixed if part(k) == p_ and rg[k][0] != 0]
        ok = [o for o, d in exps.items() if len({int(np.sign(rg[k][0] / mult[k])) * (-1)**d(k) for k in ks}) <= 1]
        print(f"   λ = {p_[0]}、ν^{p_[1]}、{p_[2]:+d}：一意 {len(ks)}/{sum(1 for k in keys if part(k) == p_)}、成り立つ規則 {'・'.join(ok) if ks else '（一意な単項式なし）'}", flush=True)
    # λ = (2) の r^a s^b ν^c と λ = (1,1) の r^{a+4} s^b ν^c（§12.13 で N が一致する低い q の領域）で e も一致するか
    pr = [(k, ("11",) + k[1:2] + (k[2] + 4,) + k[3:]) for k in keys if k[0] == "2"]
    pr = [(k, k2) for k, k2 in pr if k2 in mult and mult[k] == mult[k2] and k in fixed and k2 in fixed]
    print(f"λ = (2) と q^(−2)·λ = (1,1) で c が等しく両方一意な単項式 {len(pr)}、そのうち e も等しい {sum(1 for k, k2 in pr if rg[k] == rg[k2])}", flush=True)
    for lab in ("2", "11"):
        cnt = {"下": [0, 0], "上": [0, 0]}
        for c in {k[-1] for k in keys if k[0] == lab}:
            ks = [k for k in keys if k[0] == lab and k[-1] == c]
            mid = (min(k[2] for k in ks) + max(k[2] for k in ks)) / 2
            for k in ks:
                h = "下" if k[2] <= mid else "上"; cnt[h][0] += k in fixed; cnt[h][1] += 1
        print(f"λ = {lab}：一意な単項式は q の指数がブロックの下半分で {cnt['下'][0]}/{cnt['下'][1]}、上半分で {cnt['上'][0]}/{cnt['上'][1]}", flush=True)
    must = [k for k in keys if max(abs(rg[k][0]), abs(rg[k][1])) < abs(int(mult[k]))]
    print("どの解でも入れ替わる対を含む単項式（λ、ν、r、s の指数、c、e）：" + "、".join(f"{k[0]} ν^{k[-1]} r^{k[2]} s^{k[3]} c={int(mult[k])} e={rg[k][0]}" for k in must), flush=True)
    return st


def n1_partial():
    """n = 1、捻り (t − ν²)(ν² + 1)(q + 1)(ν²q + 1)：(i) だけ、(i) ＋ (ii) の β 奇の部分、(i) ＋ (ii)。"""
    import rk_b_knot_fixed_sector as K
    import rk_b_cover_split as C
    u = sp.Symbol("u", positive=True)
    for m in (1, 2):
        with contextlib.redirect_stdout(io.StringIO()):
            sec = C.cover_sectors(m)
        P = sp.cancel((1 - t) * sp.cancel(sec["x"].to_sympy().subs(u, nu * s)))
        mult = TM.monomials(sp.cancel(P / (s**2 * (t - nu**2) * (nu**2 - 1) * (q - 1) * (nu**2 * q - 1))))
        g = K.load_g(m)
        conds = {"i": (sp.cancel(2 * (1 + t) * g[(2,)]), {s: I}, lambda k: 1),
                 "ii": (sp.cancel(-(1 + q) * (g[(1,)].subs({r: r**2, s: s**2, nu: nu**2}, simultaneous=True) + 1)), {s: r}, lambda k: 1)}
        pre = s**2 * (t - nu**2) * (nu**2 + 1) * (q + 1) * (nu**2 * q + 1)
        basis = lambda k: pre * TM.twisted(k, "mono")
        for lab, use, odd in (("(i) だけ", ("i",), ()), ("(i) ＋ (ii) の β 奇", ("i", "ii"), ("ii",)), ("(i) ＋ (ii)", ("i", "ii"), ())):
            keys, A, b = rows(mult, basis, {k: conds[k] for k in use}, odd)
            res = freedom(mult, A, b)
            msg = f"解あり、一意でない {len(res['free'])}/{res['n']}" if res["feasible"] else "解なし"
            print(f"T(2,{2 * m + 1}) {lab}：{msg}、零空間の次元 {res['nullity']}", flush=True)
            if use == ("i", "ii") and not odd:
                blocks_ = sorted({k[-1] for k in mult})
                for fb in [()] + [(k,) for k in blocks_]:
                    sig = qrule(mult, A, b, fb)
                    if sig is not None or not fb:
                        print(f"   規則 e = σ_β (−1)^(q の指数) c（自由なブロック {fb}）："
                              + ("解なし" if sig is None else "解あり、σ = " + " ".join(f"ν^{g[1]}:{v}" for g, v in sig.items())), flush=True)


def times(P, a, b, c, sg=-1):
    """P·(1 + sg r^a s^b ν^c)。P は {(r, s, ν の指数): 係数}。"""
    out = dict(P)
    for (x, y, z), v in P.items():
        k = (x + a, y + b, z + c); out[k] = out.get(k, 0) + sg * v
        if out[k] == 0: del out[k]
    return out


def greedy_cells(P, cand, cost, steps=6):
    """(1 + sg m)（cand の (a, b, c, sg)）を cost が減る限り貪欲に掛ける。選んだ因子と、各段の cost を返す。"""
    seq = []
    for _ in range(steps):
        best = min(cand, key=lambda x: cost(times(P, *x)))
        P2 = times(P, *best)
        if cost(P2) >= cost(P):
            break
        seq.append((best, cost(P), cost(P2))); P = P2
    return seq, P


def mirror_sign(P):
    """ν のブロックごとの (q, t) → (t^{−1}, q^{−1}) の対称性：単項式のずらしの後に +P、−P、どちらでもない（×）。"""
    out = {}
    for c in sorted({k[2] for k in P}, reverse=True):
        g = {k[:2]: v for k, v in P.items() if k[2] == c}
        m = {(-b, -a): v for (a, b), v in g.items()}
        sh = lambda h: {(a - min(x for x, _ in h), b - min(y for _, y in h)): v for (a, b), v in h.items()}
        g, m = sh(g), sh(m)
        out[c] = "+" if m == g else ("−" if m == {k: -v for k, v in g.items()} else "×")
    return out


def n2_series():
    """N_λ（n = 2、trefoil）の系列の同定（theory §12.13、progress.md 続き 77）。
    (a) ν 方向の系列：端から 4 ブロックの状態数を減らす (1 − ν^{∓2} r^a s^b) を貪欲に探す（n = 1 の T(2,5) を対照にする）。
    (b) ブロックの中の塔：(1 ± m)、m = q^a t^b（a, b ≤ 2）を貪欲に探す。
    (c) 鏡映 (q, t) → (t^{−1}, q^{−1})。(d) 端のブロックの閉じた形。(e) λ = (2) と q^{−2}·λ = (1,1) の差。"""
    import rk_b_cover_split as C
    u = sp.Symbol("u", positive=True)
    st = lambda P, lim=lambda c: True: sum(abs(v) for k, v in P.items() if lim(k[2]))
    nu_cand = lambda dc: [(a, b, dc, -1) for a in range(-14, 15, 2) for b in range(-14, 15, 2)]
    nm = lambda a, b: (f"q{'²' if a == 4 else ''}" if a else "") + (f"t{'²' if b == 4 else ''}" if b else "")

    def nu_series(P, tag):
        top, bot = max(k[2] for k in P), min(k[2] for k in P)
        for dc, lim, end in ((-2, lambda c: c > top - 8, "上"), (2, lambda c: c < bot + 8, "下")):
            seq, _ = greedy_cells(P, nu_cand(dc), lambda X: st(X, lim))
            print(f"  {tag}、{end}の端から：" + (" ".join(f"(1 − ν^{dc} r^{a} s^{b})[{c0}→{c1}]" for (a, b, _, _), c0, c1 in seq) or "なし"), flush=True)

    print("(a) ν 方向の系列（端から 4 ブロックの状態数）", flush=True)
    with contextlib.redirect_stdout(io.StringIO()):
        sec = C.cover_sectors(2)
    P1 = sp.cancel((1 - t) * sp.cancel(sec["x"].to_sympy().subs(u, nu * s)))
    N1 = sp.cancel(P1 / (s**2 * (t - nu**2) * (nu**2 - 1) * (q - 1) * (nu**2 * q - 1)))
    nu_series({k[1:]: int(v) for k, v in TM.monomials(N1).items()}, "対照 n = 1 の T(2,5)")
    mult, _ = n2_setup()
    N = {lab: {k[2:]: int(v) for k, v in mult.items() if k[0] == lab} for lab in ("2", "11")}
    tower = lambda P: times(times(times(P, 2, 0, 0), 0, 2, 0), 2, 2, 0)          # ×(1 − q)(1 − t)(1 − qt)
    for lab, P in N.items():
        nu_series(P, f"N_{lab}")
        nu_series(tower(P), f"N_{lab}(1 − q)(1 − t)(1 − qt)")
    print("(b) ブロックの中の塔 (1 ± q^a t^b) と (c) 鏡映", flush=True)
    cand = [(a, b, 0, sg) for a in (0, 2, 4) for b in (0, 2, 4) if (a, b) != (0, 0) for sg in (1, -1)]
    for lab, P in N.items():
        mir = mirror_sign(P)
        for c in sorted({k[2] for k in P}, reverse=True):
            g = {k: v for k, v in P.items() if k[2] == c}
            seq, _ = greedy_cells(g, cand, st)
            print(f"  λ = {lab}、ν^{c}：鏡映 {mir[c]}、{st(g)} → {seq[-1][2] if seq else st(g)}  "
                  + "".join(f"(1{'+' if x[3] > 0 else '−'}{nm(x[0], x[1])})" for x, _, _ in seq), flush=True)
        print(f"  N_{lab} の状態数 {st(P)} → ×(1 − q)(1 − t)(1 − qt) で {st(tower(P))}", flush=True)
    print("(d) 端のブロック", flush=True)
    blk = lambda lab, c: sp.expand(sum(v * r**k[0] * s**k[1] for k, v in N[lab].items() if k[2] == c))
    V, Vp, W = 1 + t + q * t, 1 + q + q * t, 1 + q**2 + q**2 * t**2
    h = lambda k, x: sp.expand(sum(x[0]**i * x[1]**j * x[2]**(k - i - j) for i in range(k + 1) for j in range(k + 1 - i)))
    print(f"  N_2 の ν^-8 = r^-7 s^3 (t − q)：{sp.expand(blk('2', -8) - r**-7 * s**3 * (t - q)) == 0}、"
          f"ν^-10 = r^-11 s^3 (t − q) V V′：{sp.expand(blk('2', -10) - r**-11 * s**3 * (t - q) * V * Vp) == 0}", flush=True)
    for lab, K, a0 in (("2", 5, -41), ("11", 4, -37)):
        ok = sp.expand(blk(lab, -24) + r**a0 * s**17 * (1 + q * t) * h(K, (1, q**2, q**2 * t**2))) == 0
        print(f"  N_{lab} の ν^-24 = −r^{a0} s^17 (1 + qt) h_{K}(W) = −r^{a0} s^17 (1 + qt) ψ²h_{K}(V′)：{ok}", flush=True)
    print("(e) λ = (2) と q^{−2}·λ = (1,1) の差（ブロックごとの状態数、差の q の指数の範囲は N_11 の範囲と比べる）", flush=True)
    for c in sorted({k[2] for k in N["2"]}, reverse=True):
        g2 = {k[:2]: v for k, v in N["2"].items() if k[2] == c}
        g11 = {k[:2]: v for k, v in N["11"].items() if k[2] == c}
        D = {(a + 4, b): v for (a, b), v in g2.items()}
        for k, v in g11.items():
            D[k] = D.get(k, 0) - v
        D = {k: v for k, v in D.items() if v}
        rng = lambda g: (min(a for a, _ in g) // 2, max(a for a, _ in g) // 2) if g else None
        print(f"  ν^{c}：差 {sum(abs(v) for v in D.values())}（N_2 {sum(abs(v) for v in g2.values())}、N_11 {sum(abs(v) for v in g11.values())}）、"
              f"r の指数/2 の範囲：差 {rng(D)}、N_11 {rng(g11)}", flush=True)


def n2(all_inv=False):
    import rk_f_twisted22_lam as L
    bl = L.blocks()
    full = lambda name: sum((v * nu**k for k, v in bl[name].items()), sp.Integer(0))
    d = pickle.load(open(L.COVER, "rb"))
    mult = {}
    for lab in L.LAMS:
        X = sp.sympify(d["hat" + lab], locals={"r": r, "s": s, "nu": nu})
        N = sp.cancel(X / (s * (t - nu**2) * (nu**2 - 1) * (nu**2 * q - 1)))
        ml = TM.monomials(N)
        print(f"N_({','.join(lab)})：単項式 {len(ml)}、状態の数 Σ|c| = {sum(abs(v) for v in ml.values())}、"
              f"指数の単項式（P̂/s）{len(TM.monomials(sp.cancel(X / s)))}", flush=True)
        for key, c in ml.items():
            mult[(lab,) + key] = c
    pre = s * (t - nu**2) * (nu**2 + 1) * (nu**2 * q + 1)
    basis = lambda key: pre * TM.twisted(key[1:], "mono")
    conds = [(full("R22"), {s: I}, lambda key: L.W22),
             (full("R4"), {s: I}, lambda key: L.CHI[key[0]][(2,)] * L.W4),
             (full("Eq"), {s: r}, lambda key: L.W22)]
    t0 = time.time()
    res = solve(mult, basis, conds, all_inv)
    if res["feasible"]:
        fl = {lab: sum(1 for k in res["free"] if k[0] == lab) for lab in L.LAMS}
        nl = {lab: sum(1 for k in mult if k[0] == lab) for lab in L.LAMS}
        inv = sum(1 for k in mult if abs(res["e"][k]) == abs(mult[k]))
        print(f"n = 2：解あり、一意でない (2): {fl['2']}/{nl['2']}、(1,1): {fl['11']}/{nl['11']}、|e| = |c| が {inv}/{len(mult)}"
              f"（{time.time() - t0:.0f} 秒）", flush=True)
    else:
        print(f"n = 2：解なし（変数 {res['n']}、{time.time() - t0:.0f} 秒）", flush=True)
    return res


if __name__ == "__main__":
    ai = "--all-inv" in sys.argv
    if sys.argv[1] == "full":
        for m_ in (1, 2):
            for use in (("i",), ("i", "ii")):
                full_rule(m_, use)
    elif sys.argv[1] == "n2series":
        n2_series()
    elif sys.argv[1] == "n2rule":
        n2_rule("--nouniq" not in sys.argv)
    elif sys.argv[1] == "n2pattern":        # 規則なしで対が最小の解の符号の型と一意性
        n2_pattern()
    elif sys.argv[1] == "n2choice":         # 部分ごとの規則の選択も変数にして対の数を最小にする
        n2_choice(tuple(sys.argv[2].split(",")), sys.argv[3] if len(sys.argv) > 3 else "sign")
    elif sys.argv[1] == "n2tower":           # 塔のセルの規則（0, q, t, qt）と、鏡映で対称な規則（0, qt）だけの場合
        n2_rule("--nouniq" not in sys.argv, tuple(sys.argv[2].split(",")), uniform=True)
    elif sys.argv[1] == "rule":
        knot_rule(1, save=False); knot_rule(2, ("i",)); knot_rule(2, save=True)
    elif sys.argv[1] == "lower":
        lower_series(1); lower_series(2)
    elif sys.argv[1] == "n1partial":
        n1_partial()
    elif sys.argv[1] == "n2partial":
        n2_partial("--nouniq" not in sys.argv)
    elif sys.argv[1] == "n2scan":
        n2_scan()
    elif sys.argv[1] == "n1":
        out = n1(ai)
        pickle.dump(out, open(f"{OUT}/twisted_tensor_n1{'_allinv' if ai else ''}.pkl", "wb"))
    else:
        res = n2(ai)
        pickle.dump(res, open(f"{OUT}/twisted_tensor_n2{'_allinv' if ai else ''}_m1.pkl", "wb"))
