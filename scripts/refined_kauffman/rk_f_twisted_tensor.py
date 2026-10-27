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
使い方：python rk_f_twisted_tensor.py n1 | n2 [--all-inv]（--all-inv：N の全状態が不変、e = ±c）。"""
import sys, io, contextlib, pickle, itertools, time
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
    n = len(keys)
    kk = lambda k: (-k[-1] - 6) // 2
    d = lambda k: (k[1] + 4 + 2 * kk(k)) // 2 if kk(k) <= 2 * m else k[2] // 2
    grp = lambda k: (kk(k), int(np.sign(float(mult[k]))))
    gl = sorted({grp(k) for k in keys}); gi = {g_: i for i, g_ in enumerate(gl)}; ng = len(gl)
    sgn = np.array([float(np.sign(float(mult[k]))) * (-1)**d(k) for k in keys])
    cab = np.array([abs(float(mult[k])) for k in keys])
    # 変数 [y (n), v (n), z (ng)]：w = |c| − 2y、e = sgn (w − 2v)
    nv = 2 * n + ng
    rowsA, lo_, hi_ = [], [], []
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
    cons.append(LinearConstraint(np.array([x[0] for x in Ms]), [x[1] for x in Ms], [x[2] for x in Ms]))
    bnd = Bounds(np.zeros(nv), np.concatenate([np.floor(cab / 2), cab, np.ones(ng)]))
    opt = lambda cc: milp(cc, constraints=cons, bounds=bnd, integrality=np.ones(nv))
    res = opt(np.zeros(nv))
    tag = f"T(2,{2 * m + 1}) 全ブロックに規則、条件 {'+'.join(use)}"
    if res.status != 0:
        print(f"{tag}：解なし", flush=True); return
    def e_of(x):
        return np.round(sgn * (cab - 2 * x[:n] - 2 * x[n:2 * n]))
    e0 = e_of(res.x)
    nfree = 0
    for i in range(n):           # e_i の最小と最大（e は y, v の線形関数）
        cc = np.zeros(nv); cc[i] = -2 * sgn[i]; cc[n + i] = -2 * sgn[i]
        if round(opt(cc).fun) != round(-opt(-cc).fun): nfree += 1
    print(f"{tag}：解あり、e が一意でない単項式 {nfree}/{n}、入れ替わる対の数 {int(sum(cab - np.abs(e0)) // 2)}", flush=True)


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
