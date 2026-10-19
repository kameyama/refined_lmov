r"""閉弦の refined orientifold（SO(2N) の S³）の自由エネルギーを、凍結した状態と不変な状態に分ける（progress.md 続き 38）。
(A) D_N の積表示 Z = Π_{m<β} Π_{α>0}(q^{−m/2}t^{−(α,ρ)/2} − q^{m/2}t^{(α,ρ)/2})（AgSc 1202.4456 §6、t = q^β）から
    F_d = (1/d)(1 − t^d)/(1 − q^d) Σ_{α>0} t^{(α,ρ)d}（−log Z の log 部分の d 番目）を厳密に作り、
    AgSc §6 の 3 項（½ 向き付けられた項、RP² の項、(t^{1/2}+t^{−1/2}) の項）を引いた残りが N によらないかを見る。
(B) Avetisyan–Mkrtchyan 2205.12832 (3.4)–(3.7) の D 型の c = 1（奇数の n）と Klein 瓶（c = 2、全 n）の和が、
    全 n の凍結した状態の項と、(1 + t^n) を分母にもつ不変な状態の項の和に等しいこと（n ごとの有理式の恒等式）。
    AM の変数 (q, t) は本研究の (q^{1/2}, t^{1/2})。"""
import sympy as sp

q, t, n_ = sp.symbols("q t n", positive=True)


def F_product(N):
    """d = 1 の F_d（q, t をそれぞれ q^d, t^d に置き換えれば一般の d）。"""
    rho = [N - 1 - i for i in range(N)]
    hts = [rho[i] - rho[j] for i in range(N) for j in range(i + 1, N)] + [rho[i] + rho[j] for i in range(N) for j in range(i + 1, N)]
    return (1 - t) / (1 - q) * sum(t**h for h in hts)


def F_AgSc(N):
    rq, rt = sp.sqrt(q), sp.sqrt(t)
    zq, zt = rq - 1 / rq, rt - 1 / rt
    pre = rt / rq
    return (sp.Rational(1, 2) * t**(2 * N - 1) * pre / (zq * zt) + t**(N - sp.Rational(1, 2)) * pre / zq
            - sp.Rational(1, 2) * t**(2 * N - 1) * pre / (zq * (rt + 1 / rt)))


def check_A(Ns=range(2, 8)):
    """残りが N の 1 次式（t^N を含まない古典的な項と定数写像）であることを確かめる。"""
    rests = [sp.simplify(F_product(N) - F_AgSc(N)) for N in Ns]
    lin = [sp.simplify(rests[i + 2] - 2 * rests[i + 1] + rests[i]) == 0 for i in range(len(rests) - 2)]
    slope = sp.factor(sp.simplify(rests[1] - rests[0]))
    print(f"(A) 積表示 − AgSc §6 の 3 項 は N = {list(Ns)} で N の 1 次式：{all(lin)}；N の係数 = {slope}、"
          f"N = 0 での値 = {sp.factor(sp.simplify(rests[0] - Ns[0] * slope))}", flush=True)
    return all(lin)


def check_B():
    """AM の変数 Q, T（= 本研究の q^{1/2}, t^{1/2}）。n 番目の項。X = e^{τ1/2}、Y = e^{τ2} = X²/Q（τ1 = g_s(1+ay−y), τ2 = (1+2ay−2y)g_s/2）。"""
    Q, T, X = sp.symbols("Q T X", positive=True)
    c1_all = X / (Q - 1 / Q)                                    # c = 1 の n 番目を全 n に
    c1_psi2 = X**2 / (2 * (Q**2 - Q**-2))                       # 奇数の n に直すときに引く ½ψ²
    Y = X**2 / Q
    kb = Y * (Q / T - T / Q) / (2 * (Q**2 - Q**-2) * (T + 1 / T))   # Klein 瓶
    Q0 = X**2 * T / Q                                           # e^{τ0}：τ0 = (1+2ay−y)g_s/2
    inv = -Q0 / (2 * (Q - 1 / Q) * (T + 1 / T))                  # 不変な状態（ε = +1、電荷 e^{τ0} の (1+t) の項）
    ok = sp.simplify(-c1_psi2 + kb - inv) == 0
    print(f"(B) AM の Σ_奇数 c=1 + Klein 瓶 = Σ_全 凍結 X/(Q−Q⁻¹) − ½ e^(nτ0)/((Q^n−Q^−n)(T^n+T^−n))：{ok}", flush=True)
    even = sp.simplify(Q0 / (2 * (Q - 1 / Q) * (T - 1 / T)) + inv - Q0 / (T * (Q - 1 / Q) * (T**2 - T**-2)))
    print(f"    ½F_A の n 番目と合わせると z₂ の偶数のモードだけの e^(nτ0) T^(−n)/((Q^n−Q^−n)(T^2n−T^−2n))：{even == 0}；"
          f"凍結した状態の電荷 X と e^(τ0) の比 X²/e^(τ0) = {sp.simplify(X**2 / Q0)}", flush=True)
    return ok


if __name__ == "__main__":
    check_A()
    check_B()
