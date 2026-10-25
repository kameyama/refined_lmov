r"""被覆の随伴（composite (1;1)）の reduced な値 R_adj（data/refined_kauffman/daha/Radj_A_m{m}.pkl、GL の DAHA、大きな N）が、
ホモロジーの Poincaré 多項式として読めるか（係数が正になる符号の規則があるか）を調べる（progress.md 続き 40）。
(1) KN の太字変換を鏡像で当てる：ν² = −1/(𝐚²𝐭)、r = −1/(𝐪𝐭)、s = 1/𝐪。
(2) 単項式 r^a s^b ν^c に i^{α a/2 + β b + γ c/2}（α, β, γ = 0..3）を掛ける偶奇の規則で、係数が一斉にそろうか。"""
import itertools, pickle
import sympy as sp

r, s = sp.symbols("r s", positive=True)
nu = sp.Symbol("nu", positive=True)


def load(m):
    d = pickle.load(open(f"data/refined_kauffman/daha/Radj_A_m{m}.pkl", "rb"))
    return sp.expand(sum(sp.sympify(v, locals={"r": r, "s": s}) * nu**j for j, v in d.items()))


def bold_signs(R):
    n2, Q, T, A = sp.symbols("n2 Q T A")
    e = sp.expand(sp.cancel(sp.expand(R.subs(nu, sp.sqrt(n2))).subs({n2: -1 / (A**2 * T), r: -1 / (Q * T), s: 1 / Q}, simultaneous=True)))
    cs = sp.Poly(sp.expand(e * (Q * T * A)**200), Q, T, A).coeffs()
    return sum(1 for c in cs if c > 0), sum(1 for c in cs if c < 0)


def parity_rules(R):
    terms = sp.Poly(sp.expand(sp.cancel(R * (nu * r * s)**40)), r, s, nu).terms()
    hits = []
    for al, be, ga in itertools.product(range(4), repeat=3):
        sg = [c * sp.I**(al * (e[0] // 2) + be * e[1] + ga * (e[2] // 2)) for e, c in terms]
        if all(x.is_real for x in sg) and (all(x > 0 for x in sg) or all(x < 0 for x in sg)):
            hits.append((al, be, ga))
    return hits, len(terms)


if __name__ == "__main__":
    for m in (1, 2):
        R = load(m)
        pos, neg = bold_signs(R)
        hits, n = parity_rules(R)
        print(f"T(2,{2*m+1})：R_adj の項 {n}、太字（鏡像の KN 規約）で正 {pos}・負 {neg}、係数がそろう偶奇の規則 {hits}", flush=True)
