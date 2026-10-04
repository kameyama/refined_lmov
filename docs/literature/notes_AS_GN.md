# Notes: refined CS / DAHA / refined Rosso–Jones conventions (AS, GN, DMMSS, MMSS)

Sources (LaTeX only):

- **AS-long** = `1105.5117/src/Final_revised7arxivb.tex` (Aganagic–Shakirov, harvmac; equation labels `\eqn\label`, auto-numbered, so labels are cited instead of numbers)
- **AS-short** = `1202.2489/src/short_final.tex` (ISO-8859 encoded; use `grep -a`)
- **GN** = `1304.3328/src/main.tex` (Gorsky–Negut). Sec. 2.6 = "Refined Chern-Simons theory" (line 473)
- **DMMSS** = `1106.4305/src/main.tex` (Dunin-Barkowski, Mironov, Morozov, Sleptsov, Smirnov, "Superpolynomials for torus knots from evolution induced by cut-and-join operators")
- **MMSS** = `1201.3339/src/main.tex` (Mironov, Morozov, Shakirov, Sleptsov, "Interplay between MacDonald and Hall-Littlewood expansions of extended torus superpolynomials")

Remarks marked **[remark]** are my own checks or comments, not quotations.

---

## 1. Refined S and T matrices; the torus knot operator

### 1.1 AS-long, Sec. 4 (first form, with unspecified normalizations)

AS-long:1533–1544
```latex
\eqn\ssone{
S_{{\bar i}j}/S_{00} = M_{R_i}(t^{\rho}) M_{R_j}(q^{\lambda_{R_i}} t^{\rho}),
}
...
\eqn\sstwo{
S_{00} =c_S \prod_{m=0}^{\beta -1} \prod_{i=1}^{N-1}(q^{-m/2} t^{-i/2}-q^{m/2}t^{i/2})^{N-i}.
}
...
\eqn\ttone{
{T_{{\bar i} j}} = T_i \delta^i_ j =c_T  q^{{1\over 2}(\lambda_{R_i} + \beta \rho)^2} \delta^i_j.
}
```
Here $\lambda_R$ is the $SU(N)$ highest weight and $\rho$ the Weyl vector. AS-long:2733 (App. A) gives the concrete conventions: "Weyl vector $\rho$ ... (namely, $\rho_j = (N+1)/2 - j$) ... $\lambda_R$ ... (namely, $(\lambda_R)_j = R_j - |R|/N$)".

### 1.2 AS-long, Sec. 5.2 (normalizations fixed); same as AS-short Sec. 2.2

AS-long:1805–1846 (AS-short:771–804 is identical, except that it writes $R_i$ for $\lambda_{R_i}$ and $T^i_j$ for $T_{\bar i j}$):
```latex
\eqn\metd{
g_i = \prod_{m=0}^{\beta -1}
\prod_{\alpha>0}{q^{-{1\over 2}( \alpha, \lambda_R)}t^{ -{1\over 2} (\alpha,\rho)} q^{-{m\over 2}}-
q^{{1\over 2} (\alpha, \lambda_R)} t^{ {1\over 2} (\alpha,\rho)} q^{m\over 2}\over
q^{-{1\over 2}( \alpha, \lambda_R)}t^{ -{1\over 2} (\alpha,\rho)} q^{m\over 2}-
q^{{1\over 2} (\alpha, \lambda_R)} t^{ {1\over 2} (\alpha,\rho)} q^{-{m\over 2}}}
}
...
\eqn\rtu{
q = e^{2\pi i \over k+\beta N}, \qquad t = e^{2 \pi i \beta \over k+\beta N}.
}
...
\eqn\nsm{
S_{{\bar i}j} =d_{R_i}(t^{\rho})M_{R_j}(t^{\rho} q^{ \lambda_{R_i}}),
}
\eqn\ntm{T_{{\bar i} j} = {\delta^i}_j \;q^{{1\over 2}(\lambda_{R_i}, \lambda_{R_i})}\,t^{(\lambda_{R_i}, \rho)}t^{{\beta-1\over 2}(\rho, \rho)}q^{-{k\over 2 N}(\rho,\rho)},
}
d_{R_i}(t^{\rho}) = {{i}^{N(N-1)/2}\over N^{1\over 2}(k+\beta N)^{N-1\over 2}}
\; \prod_{m=0}^{\beta-1}\prod_{\alpha>0}(q^{-m/2} t^{-(\alpha,  \rho)/2} - q^{m/2} t^{(\alpha,  \rho)/2}) \;M_{R_i}(t^{\rho})
```
Index raising: AS-long:1871 "$S_{{\bar i} {j}} = \sum_{k} g_{{\bar i}k} {S^{k}}_j = g_{i} {S^{i}}_j$". The relations are $S^4=1$, $(ST)^3=S^2$ (AS-long:1880–1882).

Quantum dimension (AS-long:1897–1901):
```latex
S_{{\bar i}0}/S_{00} =
\prod_{m=0}^{\beta -1} \prod_{1\leq i<j\leq N}
{q^{{R_j - R_i\over 2}} t^{{i-j\over 2}} q^{-{m\over 2}}-q^{R_i - R_j\over 2} t^{j-i\over 2} q^{m\over 2}\over
 t^{{i-j\over 2}} q^{-{m\over 2}}-t^{j-i\over 2} q^{m\over 2}}
=   M_{R_i}(t^{\rho}).
```

**Macdonald normalization.** AS uses monic $M_R = P_R(q,t)$; it does not use the orthonormal form (AS-long:1798 footnote; AS-short:761 footnote). The commented-out examples at AS-long:2760–2761 confirm this:
```latex
%M_1 = p_1
%M_2 = \frac{(1-t)(1+q)}{(1-tq)} \frac{p_1^2}{2} + \frac{(1+t)(1-q)}{(1-tq)} \frac{p_2}{2}, \ \ \ M_{11} = \frac{p_1^2}{2} - \frac{p_2}{2}
```
**[remark]** The $m_2$ coefficient of this $M_2$ is 1, so it is the monic $P_2$. The Macdonald norm is AS-long:2722
```latex
g_R = g_{\emptyset} \prod\limits_{(i,j) \in R} \frac{1 - t^{R^T_j-i} q^{R_i-j+1} }{1 - t^{R^T_j-i+1} q^{R_i-j} } \ \frac{1 - (t/q) t^{N-i} q^{j} }{1 - t^{N-i} q^{j}}
```
The Cauchy kernel is AS-long:2780–2784 ($m_R$ = large-$N$ limit of $g_R^{-1}$).

**Gamma factors:** AS has none. The only extra structure is the non-trivial metric $g_i$ and the prefactor $S_{00}$ / $d_R$ given above.

**[remark] T-eigenvalue in partition language.** With $\rho_j=(N+1)/2-j$ and $(\lambda_R)_j=R_j-|R|/N$:
- $(\lambda,\rho)=|R|(N-1)/2-\sum_i (i-1)R_i$
- $\tfrac12(\lambda,\lambda)=\sum_i R_i(R_i-1)/2+|R|/2-|R|^2/(2N)$

So
$$T_R \propto q^{n(R')}\,t^{-n(R)}\times q^{|R|/2-|R|^2/2N}\,t^{|R|(N-1)/2},$$
with $n(R)=\sum_i (i-1)R_i$ and $n(R')=\sum_i R_i(R_i-1)/2$. Note the **minus** sign on the $t$-exponent (compare GN below).

### 1.3 Torus knot operator (AS-long Sec. 5.4)

AS-long:2001–2050:
```latex
{\cal O}^{(n,m)}_{R_i} =K{\cal O}^{(0,1)}_{R_i} K^{-1}
...
K = \pmatrix{ a & n \cr
              b &  m \cr} \;\;  \in \;\; SL(2,{\bf Z}).
...
S = \pmatrix{ 0 & -1 \cr
              1 &  0 \cr} , \;\; T = \pmatrix{1& 1\cr 0&1} \;\; \in \;\; SL(2,{\bf Z}).
...
\eqn\tn{
K{\cal O}^{(0,1)}_{R_i}K^{-1} |0\rangle = \sum_{j,k,\ell} {{(K)^{\ell}}_k}{N^{k}}_{ij}{(K^{-1})^{j}}_0|R_{\ell}\rangle.
}
...
\langle 0 |{\cal O}^{(n,m)}_{R_i} S |0\rangle =
\langle 0 | K {\cal O}^{(0,1)}_{R_i}{K}^{-1} S|0\rangle =
\sum_{j,k,\ell} K_{0k}{N^{k}}_{ij}{(K^{-1})^{j}}_{\ell} {S^{\ell}}_0,
```
Here $K$ "takes the $(1,0)$ cycle to the cycle $(a,b)$ cycle dual to $(n,m)$", with $am-nb=1$. The Wilson line is multiplication by the Macdonald polynomial, with fusion coefficients ${\cal O}^{(0,1)}_{R_i}|R_j\rangle=\sum_k N^k_{ij}|R_k\rangle$; App. A defines $M_PM_Q=\sum_Y N^Y_{PQ}M_Y$ (AS-long:2793).

AS-short:904–925 (eq. `\knots`) gives the same formula, quoted verbatim:
```latex
Z(S^3, K_{n,m}, R_i) = \sum_{j,k,\ell} K_{0k}{N^{k}}_{ij}{(K^{-1})^{j}}_{\ell} {S^{\ell}}_p,
```
Here "$K$ is an element of $SL(2,\IZ)$ that takes the $(0,1)$ cycle to $(n,m)$ cycle ... with the only ambiguity being the choice of framing $K \mapsto K T^f$ with integer $f$".

Trefoil word (AS-long:2468–2477):
```latex
K =              S^{-1} T^{-2} S^{-1} T^{-2}   = \pmatrix{ -1 & 2 \cr -2 & 3 \cr}
\langle 0 |{\cal O}^{(2,3)}_{\tableau{1}} S |0\rangle =   \sum_{j,k} {{(S^{-1} T^{-2} S^{-1} T^{-2} )^{0}}_k} \ {N^{k}}_{\tableau{1}, j} \ {(T^2 S^{-1} T^2)^{j}}_0.
```
Other words used:
- $(2,5)$: $K=S^{-1}T^{-3}S^{-1}T^{-2}$
- $(2,2m-1)$: $K=S^{-1}T^{-m}S^{-1}T^{-2}$
- $(3,4)$: $K=S^{-1}T^{-1}ST^{3}$

(AS-long:2523, 2560, 2594.)

### 1.4 GN Sec. 2.6 (Kirillov / Etingof–Kirillov form; this is the cleanest statement)

GN:479–505:
```latex
W_{n,m}^{\lambda}|P_{\mu}\rangle=v_{n,m,\mu}^{\lambda}.
...
\begin{equation}
\label{conj}
W_{n,m}^{\lambda}=K^{-1}W_{1,0}^{\lambda}K,
\end{equation}
where $W_{1,0}^{\lambda}$ is the operator of multiplication by $P_{\lambda}$ and $K$ is any element of $\slz$ taking $(1,0)$ to $(n,m)$, seen as an endomorphism of $Z(T^2)$. The action of $\slz$ was introduced by A. Kirillov Jr. in \cite{K} and is given by the two matrices  (\cite{AS,EK,K}) written in the Macdonald polynomial basis \footnote{The operators in \cite{AS} differ from these by overall scalar factors, which are not important for us}:
\begin{equation}
\label{S and T}
S_{\lambda}^{\mu}=P_{\lambda}(t^{\rho}q^{\mu})P_{\mu}(t^{\rho}),\qquad T_{\lambda}^{\mu}=\delta_{\lambda}^{\mu}q^{\frac{1}{2}\sum_{i}\lambda_i(\lambda_i-1)}t^{\sum_{i}\lambda_i(i-1)}.
\end{equation}
which correspond to the matrices:
\sigma = \left( \begin{array}{cc} 0 & 1 \\ -1 & 0 \end{array} \right), \qquad \tau = \left( \begin{array}{cc} 1 & 1 \\ 0 & 1 \end{array} \right) \in \slz.
 The refined Chern-Simons knot invariant is then defined as
\begin{equation}
\label{CS evaluation}
\CP^{CS}_{n,m,\lambda}(q,t):=\langle 1|SW_{n,m}^{\lambda}|1 \rangle_{q,t}.
\end{equation}
```
Further GN conventions:
- $P_\lambda$ is monic: GN:233 "$P_\lambda = m_\lambda + \sum_{\mu < \lambda} d^\mu_\lambda m_\mu$".
- Inner product (GN:228–230): $\langle p_{\lambda},p_{\mu}\rangle_{q,t}= \delta_{\lambda}^{\mu} z_{\lambda}\prod_{i}\frac{1-q^{\lambda_{i}}}{1-t^{\lambda_{i}}}$.
- $\rho_N=(N-1,N-2,\ldots,1,0)$ (GN:312).

Evaluation (GN:316–332):
```latex
\e_N: V_N \longrightarrow \BK, \qquad \e_N(f) = f(t^{\rho}) = f(t^{N-1},t^{N-2},\ldots,1)
\e_N(p_{k}) = \frac{1-t^{kN}}{1-t^k}=\frac{1-u^{k}}{1-t^{k}},
\e_N(P_{\lambda}) = \frac 1{h_\lambda} \prod_{\square\in \lambda} (t^{l'(\square)}-uq^{a'(\square)})
 \qquad \Longrightarrow \qquad \e_N(J_{\lambda}) = \prod_{\square\in \lambda}\left(t^{l'(\square)}-uq^{a'(\square)}\right)
```
with $h_{\lambda} = \prod_{\square}(1-q^{a(\square)}t^{l(\square)+1})$ and $h'_{\lambda}=\prod_{\square}(1-q^{a(\square)+1}t^{l(\square)})$ (GN:242). The covector $\langle 1|S|\cdot\rangle_{q,t}$ equals $\varepsilon_N$ (GN:547), so GN's invariant is unnormalized: the unknot gives $\varepsilon_N(P_\lambda)$, the principal specialization.

DAHA form (GN:419–449):
```latex
P_{n,m}^{\lambda,N} := \gamma_{n,m}(P^{\lambda,N}) \in \sh_N
P^N_{kn,km} := \gamma_{n,m}(p_k^N) \in \sh_N
\mathcal{P}_{n,m}^{\lambda,N}(q,t) =\varepsilon_N(P_{n,m}^{\lambda,N}\cdot 1)
```
The Corollary at GN:537–550 states $\CP^{CS}_{n,m,\lambda}(q,t)=\CP^{\lambda,N}_{n,m}(q,t)$ (AS = Cherednik).

**[remark, possible sign issue]** GN writes $T=q^{n(\lambda')}t^{+n(\lambda)}$, but AS's \ntm gives $q^{n(\lambda')}t^{-n(\lambda)}$ up to $|\lambda|$-dependent factors (see 1.2). Kirillov's formula $q^{(\lambda,\lambda+2k\rho)/2}$ with $t=q^k$ and decreasing $\rho$ also gives $t^{-n(\lambda)}$. Before using GN's $T$ literally, check $(ST)^3\propto S^2$ numerically.

DMMSS uses yet another normalization, quoting AS (DMMSS:1504–1506; their $q,t$ are the squares of the standard ones):
```latex
S_{R,Q}=S_{[],[]}\frac{M_{R}(t^{2 \rho}) M_{Q}(t^{2\rho} q^{2R})}{ q^{2 |R|
|Q|/N}}, \ \ \  T_{R,Q}=\frac{q^{2 ||R||} t^{N|R|}}{t^{ ||R^{\prime}||} q^{2
|R|^2/N}} \delta_{Q,R}
```
with "$||R||=\sum_{i} R_{i}^2$". DMMSS:1508–1514 derives the unknot ${\cal P}^{\rm Unknot}_R=M_R(t^{2\rho})=M_R^*$ from it.

---

## 2. "Refined Rosso–Jones": how the coefficients $c$ are defined

**Bottom line.** None of these papers states that $c_{\lambda Q}$ comes from expanding $P_\lambda[p_k\to p_{mk}]$ in $P_Q$ with $q,t$ unchanged. DMMSS and MMSS say explicitly that the naive Adams/plethysm rule fails in the refined case. Extra "$\gamma$-factors" are needed, and the coefficients depend on the residue $r=n \bmod m$, not only on $m$ and $\lambda$. Neither paper suggests replacing $q,t$ by $q^m,t^m$. GN and AS do not use a Rosso–Jones form at all: they use $K^{-1}W_{1,0}K$, or the DAHA/elliptic Hall algebra.

### 2.1 HOMFLY (unrefined) Rosso–Jones, which DMMSS takes as the starting point

DMMSS:669–694:
```latex
\chi_R\{p^{(m)}\}\ = \sum_{Q\vdash m|R|} c_R^Q \chi_Q\{p\}
\label{cfromshift}
p^{(m)}_k = p_{mk}
\label{pmk}
\boxed{
H_R^{T[m,n]}\{p\} =
\sum_{Q\vdash m|R|} q^{-2\frac{n}{m}\varkappa_Q} c^Q_R s_Q\{p\}
= q^{-2\frac{n}{m}\hat W_{[2]}} \sum_{Q\vdash m|R|} c^Q_R s_Q\{p\}
= q^{-2\frac{n}{m}\hat W_{[2]}} s_R\{p^{(m)}\}
}
\label{Wrep}
p_k^* = \frac{a^k-a^{-k}}{q^k-q^{-k}}
```
with (DMMSS:381–392)
```latex
\varkappa_Q = \nu_{Q'} - \nu_{Q}, \ \ \ \ \
\hat W_{[2]} s_Q = \varkappa_Q s_Q
\lambda_Q = q^{-2\varkappa_Q/m},
\ \ \ \ \
\nu_Q = \sum_i (i-1)Q_i,
\ \ \ \
\varkappa_Q = {1\over 2}\sum_i Q_i(Q_i-2i+1)=\sum_{(i,j)\in Q}(i-j)
```

### 2.2 Refined version (DMMSS)

DMMSS:403–416 (introduction):
```latex
In order to lift HOMFLY to the superpolynomials one needs to deform
from the Schur to MacDonald characters
$\chi_Q = M_Q(q^2|t^2)$ separating the two constituents $\nu_Q$ and $\nu_{Q'}$
of $\varkappa_Q$:
\be
\lambda_Q = t^{2\nu_Q/m} q^{-2\nu_{Q'}/m}
\label{evms}
\ee
...
p_k =p^*_k = \frac{A^k-A^{-k}}{t^k-t^{-k}}
\label{pstar}
```
DMMSS:360–362: $P = \sum_{Q}c_{R_1\ldots R_m}^Q\lambda_Q^n \chi^*_Q$ (`\Ksum`), with the $c$ "specified by the initial condition" `\boco`: $T[m,p]=T[p,m]$ for $p<m$ (DMMSS:320–323).

DMMSS:898–915 (Sec. 3.1):
```latex
\be\label{McDd}
{\cal P}_R^{[m,n]}(A|q,t)\sim\sum c_R^Q\ t^{2\nu_Q/m} q^{-2\nu_{Q'}/m} M_Q^*
\ee
where the coefficients $c^Q_R$ now depend on not only on the number of strands
$m$ and the representation $R$, but also on
$p$. ...
The biggest loophole is the lack of the simple rule (\ref{pmk}): this makes
evaluation of
$c_R^Q$ for $|Q|\ge 4$ a piece of art.
```
DMMSS:425–437: "The rule involves a $p$-dependent exchange of $q$ and $t$ at some places".

For $T[m,km+1]$, fundamental representation (DMMSS:1240–1268; repeated in App. 6253–6280):
```latex
\boxed{ P^{T[m,k
m+1]}_{[1]}=\sum\limits_{|Q|=m} c^{Q}_{[1]} M_{Q}^{\ast} q^{-2
\nu(Q^{\prime})k} t^{2 \nu(Q)k} }
c^{Q}_{[1]}=\bar{c}^{Q}_{[1]} \gamma^{Q}_{[1]}
p_{m}=\sum_{|Q|=m}
\bar{c}^{Q}_{[1]} M_{Q}(p_{k}),
\gamma^{Q}_{[1]}=\dfrac{q^{\alpha_{Q}} }{[m]_{q}}
\sum\limits_{(i,j)\in Q} t^{2 (i-1)} q^{2Q_{1}-2Q_{i}+2 j-2}
\gamma^{[2]}=\dfrac{1+q^2}{1+q^2}=1, \ \ \ \
\gamma^{[1,1]}=\dfrac{1+t^2}{1+q^2}
```
So $\bar c$ is the naive plethysm coefficient (in Macdonald $M_Q(q^2,t^2)$, with $q,t$ unchanged), and $\gamma$ is the correction. Note also that the evolution here is in **integer** powers $k$ (relative to the $T[m,1]$ initial point), not in $n/m$.

$T[2,n]$, fundamental (DMMSS:1150–1156 and App. 2802–2830):
```latex
\boxed{
P_{[1]}^{T[2,n]}(A|q|t)
= q^{-n}\left(M_2^* - (qt)^n\frac{t^2-t^{-2}}{qt - (qt)^{-1}}M_{11}^*\right)
}
\boxed{
c^{[2]}_{[1]}=1,\ \ \ c^{[1,1]}_{[1]}=-\frac{1-t^4}{1-q^2 t^2}
\Big(\frac{q}{t}\Big)
=-{\{t^2\}\over \{qt\}}}
```
($P^{T[2,n]}_{[1]}=c^{[2]}_{[1]} M_{[2]}^{\ast} q^{-n}+ c^{[1,1]}_{[1]}M_{[1,1]}^{\ast} t^{n}$.)

**[remark]** In standard variables ($q_{\rm std}=q^2$, $t_{\rm std}=t^2$), the naive plethysm $p_2=P_2-\frac{(1+q)(1-t)}{1-qt}P_{11}$ gives $\bar c^{[1,1]}=-\frac{(1+q)(1-t)}{1-qt}$. The refined coefficient is $c^{[1,1]}=-(q/t)^{1/2}\frac{1-t^2}{1-qt}$, so the ratio is $\gamma=(q/t)^{1/2}(1+t)/(1+q)$. This is not a monomial, so it cannot be absorbed into a framing or eigenvalue redefinition.

$T[2,2k+1]$, $R=[2]$ (DMMSS:6400–6416, boxed):
```latex
P^{T[2,n]}=c^{[4]}_{[2]} M_{[4]}^{\ast} q^{-6n} + \Big(\frac{q}{t}\Big)
c^{[3,1]}_{[2]} M_{[3,1]}^{\ast} q^{-3n}t^{n} + \Big(\frac{q}{t}\Big)^2
c^{[2,2]}_{[2]} M_{[2,2]}^{\ast} q^{-2n}t^{2n}
c^{[4]}_{[2]}=1,\ \ \ c^{[3,1]}_{[2]}=-\frac{(1-t^2)(1+q^2)(1+t^2q^2)}{1-q^6
t^2} \ \ \ c^{[2,2]}_{[2]}=\frac{(1-t^4)(1-t^4q^2)}{(1-q^4 t^2)(1-q^2t^2)}
```
**[remark]** The same table in Sec. 4.1 (DMMSS:1390–1392) prints $c^{[2,2]}_{[2]}=\frac{(1-t^4)(1+t^4q^2)}{(1-q^4 t^2)(1-q^2t^2)}$, with $+$ instead of $-$. One of the two is a typo; check against AS's $[2]$-trefoil (Sec. 4).

### 2.3 MMSS (1201.3339)

MMSS:113–124:
```latex
\boxed{\ \
{\cal P}_{m,mk+r}\{p\}
 = \sum_{Q\vdash m} c_Q^{(m,r)} q^{-k\nu(Q^T)}t^{k\nu(Q)} M_Q\{p\}
= e^{k\hat W} {\cal P}_{m,r}\{p\}\ \
}
```
These are standard (asymmetric) variables: "$[x]_t = \frac{1-t^x}{1-t}$ ... differing in particular by the change $t\rightarrow t^2$, $q\rightarrow q^2$" from DMMSS (MMSS:139–147).

MMSS:157–163:
```latex
For the HOMFLY polynomials these were essentially
$q$-independent integers,
read off from the Adams decomposition at the "initial" point $n=0$ \cite{chi}.
After $\beta$-deformation they become non-trivial rational
functions of $q$ and $t$, and are {\it not} given by just a
naive deformation of the Adams rule \cite{DMMSS}.
```
MMSS:183–188 gives the Hall–Littlewood re-expansion at $k=0$:
${\cal P}_{m,r}\{p\} = \sum_{Q\vdash m} c_Q^{(m,r)}M_Q\{p\} = \sum_{Q\vdash m,\ l(Q)\leq r} h_Q^{(m,r)} L_Q\{p\}$, with $L_Q=M_Q|_{q=0}$ and $h_Q = q^{\nu(Q)} (1-t)^{l(Q)-1} \hat h_Q$. For example $\hat h^{(m,1)}_Q=\delta_{l(Q),1}$ and $\hat h^{(m,2)}_Q=1$ for $l(Q)\le2$ (MMSS:250–256, 296–303). MMSS:304–312 restates DMMSS's split $c=C\gamma$ into an Adams coefficient times $\gamma$. Example (MMSS:313–321):
```latex
C_{[4,3]}^{(7,2)} = \dfrac{(1-t)(q-t)(q^2-t)}{(1-q^2 t)(1-q^3 t)(1-q^4 t)}
```

### 2.4 GN: colored torus knots via plethysm in the elliptic Hall algebra (not Rosso–Jones)

The colored operator is a plethystic substitution in the **generators** $P_{kn,km}$, not in the power sums $p_k$ of a symmetric function. GN:419–443:
```latex
P^N_{kn,km} := \gamma_{n,m}(p_k^N) \in \sh_N
...Since the Macdonald polynomials $P^{\lambda,N}$ are polynomials in the power-sum functions $p_k^N$, then the $P_{n,m}^{\lambda,N}$ are polynomials in $P^N_{kn,km}$.
```
GN:916–926 (modified version):
```latex
\TP_{n,m}^{\la}:= \frac{t^{-n(\lambda)}}{h_{\lambda}(q,t^{-1})}\Hl\left[p_k\to \TP_{kn,km}\right].
```
GN:955–967 (colored case): $\TCP^\lambda_{n,m} = \sum_{\mu \vdash n|\lambda|} ( I_\mu | \TP^\lambda_{n,m} | 1 ) \prod_{\square \in \mu} (1 - u\chi(\square))/g_\mu$. It is computed by iterating eq. (eqn:fix3); GN gives "no nice formula".

**[remark]** This is the correct "refined Rosso–Jones" mechanism. $P_\lambda[p_k\to P_{km,kn}]$ with $q,t$ unchanged, where $P_{km,kn}=\gamma_{m,n}(p_k)$ is the $SL_2(\mathbb Z)$ image of $p_k$. Only when $(m,n)=(1,0)$ does $P_{k,0}$ act as multiplication by $p_k$. The operators $P_{km,kn}$ are not $T^{n/m}$-conjugates of multiplication by $p_{km}$ with fixed $Q$-independent $c$. This is consistent with the $\gamma$-factors of DMMSS and MMSS.

Uncolored GN theorem (GN:121–135):
```latex
\CP_{n,m}(u,q,t) =\sum_{\mu \vdash n} \frac {\widetilde{\gamma}^n}{\widetilde{g}_\mu} \sum^{\syt}_{\text{of shape }\mu}  \frac {\prod_{i=1}^{n} \chi_i^{S_{m/n}(i)} (1-u\chi_i)(q \chi_i - t)}{\left(1-\frac {q\chi_{2}}{t\chi_{1}}\right) \ldots  \left(1-\frac {q\chi_{n}}{t\chi_{n-1}}\right)}\prod_{1\leq i < j\leq n} \frac {(\chi_j-q\chi_i)(t\chi_j-\chi_i)}{(\chi_j-\chi_i)(t\chi_j-q\chi_i)}
S_{m/n}(i)=\left\lfloor \frac {im}n\right\rfloor - \left\lfloor \frac {(i-1)m}n \right\rfloor, \qquad \qquad \widetilde{\gamma}=\frac {(t-1)(q-1)}{(q-t)}
\widetilde{g}_\mu =\prod_{\square \in \lambda}(1-q^{a(\square)}t^{l(\square)+1}) \prod_{\square \in \lambda}(1-q^{-a(\square)-1}t^{-l(\square)})
```
Here "$\chi_i$ denotes the $q,t^{-1}$-weight of the box labeled by $i$". The modified form is (GN:933–940):
```latex
\TCP^\lambda_{n,m}(u,q,t) =\left ( \Lambda(u) |\TP^{\lambda}_{n,m}|1 \right) ,\ \text{where}\ \Lambda(u) = \sum_{\mu \vdash n} \frac {I_\mu}{g_\mu}  \prod_{\square\in \mu}\left(1 - u q^{a'(\square)}t^{l'(\square)}\right)
\TCP^\lambda_{n,m}(u,q,t) = (-q)^{-n|\lambda|}\CP^\lambda_{n,m}(u,q,t^{-1})
```
GN:886 and 760 give $\gamma = \frac {(q-1)(t-1)}{qt(qt-1)}$ and $g_\lambda = \prod(1-q^{a}t^{-l-1}) \prod(1-q^{-a-1}t^{l})$. Box weight $\chi(\square)=q^{i-1}t^{j-1}$ (GN:869).

---

## 3. T-eigenvalue, framing correction, normalization

- **AS** ($R$-dependent part): $q^{\frac12(\lambda_R,\lambda_R)}t^{(\lambda_R,\rho)}$, i.e. $\propto q^{n(R')}t^{-n(R)}$ (see 1.2). The full \ntm also has $R$-independent factors $t^{\frac{\beta-1}{2}(\rho,\rho)}q^{-\frac{k}{2N}(\rho,\rho)}$.
- **GN**: $q^{\frac12\sum\lambda_i(\lambda_i-1)}t^{\sum\lambda_i(i-1)}$, the same $q$-part as AS (the $t$-sign is flagged in 1.4). GN applies no framing correction; it simply defines $\langle 1|S W_{n,m}|1\rangle$.
- **DMMSS**: $\lambda_Q=t^{2\nu_Q/m}q^{-2\nu_{Q'}/m}$ (squared variables), i.e. the inverse of AS's $T^{1/m}$ in $R$-dependent content.
- **MMSS**: $q^{-k\nu(Q^T)}t^{k\nu(Q)}$, standard variables.

**AS framing correction (fundamental):** AS-long:2489–2515:
```latex
\frac{Z(K_{2,3})}{Z(\bigcirc)} \Big|_{SU(N)} = t^{2-2N} q^{-3 + 3/N} + t^{3-2N} q^{-2+3/N} - t^{2-N} q^{-2+3/N}
...
{\bf a}^2 = t^N \sqrt{t/q}, \ \ \ {\bf q} = \sqrt{t}, \ \ {\bf t} = -\sqrt{q/t}
...
\frac{Z(K_{2,3})}{Z(\bigcirc)} \Big|_{SU(N)} = \left(\frac{t}{q}\right)^{5/2} \frac{q^{3/N}}{t^{3N}} \; {\cal P}(K_{2,3})
...
T_{\tableau{1}} = \sqrt{\frac{t^{N+1}}{q^{1/N+1}}}
so the prefactor corresponds to -6 units of framing
\left(\frac{t}{q}\right)^{5/2} \frac{q^{3/N}}{t^{3N}} = \sqrt{\frac{q}{t}} \ T_{\tableau{1}}^{-6}
modulo a $\sqrt{q/t}$ factor, which vanishes in the unrefined case.
```
Further prefactors (AS-long:2546, 2577, 2618, 2641):
- $(2,5)$: $\frac{q}{t}\,T^{-10}_{\tableau{1}}$
- $(2,2m-1)$: $\left(\sqrt{\frac{q}{t}}\right)^{m-1} T^{-2(2m-1)}_{\tableau{1}}$
- $(3,4)$: $\left(\frac{q}{t}\right)^{3/2} T^{-12}_{\tableau{1}}$
- $[2]$-colored trefoil: $\frac{q}{t}\,T_{\tableau{2}}^{-6}$

**[remarks]**
1. **Pattern.** For $R=\square$, every quoted case fits $Z/Z(\bigcirc)=(q/t)^{(m-1)(n-1)/4}\,T_\square^{-mn}\,{\cal P}$. The $[2]$-trefoil prefactor $(q/t)^1$ fits $(q/t)^{|R|(m-1)(n-1)/4}T_R^{-mn}$. So the framing-0 normalization is: multiply by $T_R^{mn}$, the standard $-mn$ framing of $T(m,n)$, together with a $(t/q)^{\ldots}$ monomial.
2. **Typo in the printed $T_\square$.** The printed $T_{\tableau{1}}=\sqrt{t^{N+1}/q^{1/N+1}}$ does **not** reproduce AS's own prefactor identity. From \ntm one gets $T_\square/T_\emptyset = q^{\frac12-\frac1{2N}}t^{\frac{N-1}2}=\sqrt{q^{1-1/N}t^{N-1}}$, and with this value $\sqrt{q/t}\,T_\square^{-6}=(t/q)^{5/2}q^{3/N}t^{-3N}$ holds exactly. So the printed $T_{\tableau{1}}$ is a typo; use $T_R/T_\emptyset$ from \ntm.
3. **Normalization.** AS normalizes by the unknot in the same representation, $Z_R(K)/Z_R(\bigcirc)$, with $Z_R(\bigcirc)=M_R(t^\rho)$ (reduced). GN and DMMSS keep the unreduced invariant, whose unknot value is the principal specialization $\varepsilon_N(P_\lambda)$ resp. $M_R^*$. GN Conj. 2 (GN:462–469) defines $\CP^{\lambda,red}_{n,m}=\CP^\lambda_{n,m}/\CP^\lambda_{1,0}$.

AS $[2]$-unknot (AS-long:2637):
```latex
Z_{\tableau{2}}(\bigcirc) = \frac{(t^{N/2}-t^{-N/2})(q t^{N/2}-t^{-N/2})}{(t^{1/2}-t^{-1/2})(q t^{1/2}-t^{-1/2})}
```
DMMSS unknots (DMMSS:2644–2660): $M_1^*=\frac{A-A^{-1}}{t-t^{-1}}$ and $M_2^*=\frac{(A-A^{-1})(Aq-A^{-1}q^{-1})}{(t-t^{-1})(qt-q^{-1}t^{-1})}$. DMMSS:1352–1365 discusses that for $R\neq[1^r]$ the unknot is a choice; DMMSS and AS use $M_R^*$.

---

## 4. Stable limit $a=t^N$ and DGR variables; test polynomials

### 4.1 Variable changes

**AS-long** (steps 1–3, AS-long:2416–2440):
```latex
\eqn\mapagain{\lambda = t^{N} t^{1/2}q^{-1/2}.}
\eqn\changeofv{{\bf q} = \sqrt{t}, \qquad {\bf t} = -\sqrt{q/t}, \qquad {\bf a} = \sqrt{\lambda },}
```
In the trefoil example this reads ${\bf a}^2 = t^N\sqrt{t/q}$. "For ${\bf t} = -1$, the ${\cal P}(K)$ manifestly reduces to Jones polynomial by setting ${\bf a} = {\bf q}^2$".

**AS-short** (AS-short:931–933 and App. 1566–1568):
```latex
{\bf t } = \sqrt{t} , \ \ \ {\bf q} = -\sqrt{q/t}, \ \ \ {\bf a}^2 = t^{N} \sqrt{t/q}.
```
**[remark]** AS-short has ${\bf q}$ and ${\bf t}$ swapped relative to AS-long, and its polynomials use halved ${\bf a},{\bf q}$ exponents: trefoil ${\bf a}{\bf q}^{-1}+\ldots$ vs AS-long ${\bf a}^2{\bf q}^{-2}+\ldots$. The appendix normalizes ${\cal P}=1+\ldots$. Treat AS-short's variable line as a typo, and use AS-long / DGR.

**DMMSS** (squared variables; DMMSS:1143–1148):
```latex
\boxed{
{\bf q} = t, \ \ \ \ \ \ {\bf t} = -q/t,
\ \ \ \ \ \ A \ \stackrel{(\ref{Avsatq})}{=}\ {\bf a}\sqrt{-{\bf t}}
}
```
**MMSS** (standard variables, MMSS:149–152): ${\bf t} = -\sqrt{q/t},\ {\bf q} = \sqrt{t},\ {\bf a}^2 = A\sqrt{t/q}$, with $p^*_k = \frac{1-A^k}{1-t^k}$ and $A=t^N$.

**GN**: $u=t^N$ (GN:324), and Sec. 3.4 (proof of Conj. 1, stabilization) gives $\CP^\lambda_{n,m}=\langle \Vv(u)|P^\lambda_{n,m}|1\rangle_{q,t}$ (GN:671–676), where (GN:563)
```latex
\Vv(u) := \sum_{\lambda} \frac {J_\lambda }{h_\la h'_\la}  \prod_{\square\in \lambda}\left(t^{l'(\square)}-u q^{a'(\square)}\right)
```
GN gives no explicit map to DGR $({\bf a},{\bf q},{\bf t})$. For the trefoil, the AS-long map applied to $u=t^N$ gives ${\bf a}^2=u\sqrt{t/q}$.

### 4.2 Explicit test polynomials (verbatim)

AS-long trefoil (AS-long:2504):
```latex
{\cal P}(K_{2,3})({\bf q}, {\bf a}, {\bf t}) = {\bf a}^2 {\bf q}^{-2} + {\bf a}^2 {\bf q}^2 {\bf t}^2 + {\bf a}^4 {\bf t}^3
```
$(2,5)$ (AS-long:2546–2550):
```latex
\frac{Z(K_{2,5})}{Z(\bigcirc)} \Big|_{SU(N)} = \frac{q}{t} \ T^{-10}_{\tableau{1}} \ \big( t^{2N-1} q^{-1} + t^{2N} + t^{2N+1} q - t^{3N-1} - t^{3N} q \big)
{\cal P}(K_{2,5})\left({\bf q}, {\bf a}, {\bf t} \right) = {\bf a}^4 {\bf q}^{-4} + {\bf a}^4 {\bf t}^2 + {\bf a}^4 {\bf q}^4 {\bf t}^4 + {\bf a}^6 {\bf q}^{-2} {\bf t}^3 + {\bf a}^6 {\bf q}^2 {\bf t}^5
```
$(2,2m-1)$ (AS-long:2572, 2584):
```latex
\frac{Z(K_{2,2m-1})}{Z(\bigcirc)} \Big|_{SU(N)} = T^{-2(2m-1)}_{\tableau{1}} \ t^{Nm-m} \ \left[ t^{N+1} \sum\limits_{i=0}^{m-1} (tq)^i - \sum\limits_{i=1}^{m-1} (tq)^i \right].
{\cal P}(K_{2,2m-1})\left({\bf q}, {\bf a}, {\bf t} \right) = {a}^{2m-2} \sum\limits_{i=0}^{m-1} {\bf q}^{4m-2i+2} {\bf t}^{2i} + {\bf a}^{2m} \sum\limits_{i=1}^{m-1} {\bf q}^{4i-2m} {\bf t}^{2i+1}.
```
**[remark]** The first sum is misprinted. To agree with the $(2,3)$ and $(2,5)$ results it must be ${\bf q}^{4i-2m+2}$ instead of ${\bf q}^{4m-2i+2}$.

$(3,4)$ (AS-long:2622–2625):
```latex
{\cal P}(K_{3,4})\left({\bf q}, {\bf a}, {\bf t} \right) = {\bf a}^6 {\bf t}^4+{\bf a}^6 {\bf q}^6 {\bf t}^6+{\bf a}^6 {\bf q}^2 {\bf t}^4+{\bf a}^6 {\bf q}^{-2} {\bf t}^2+{\bf a}^6 {\bf q}^{-6}+
+ {\bf a}^8 {\bf q}^4 {\bf t}^7+{\bf a}^8 {\bf q}^2 {\bf t}^7+{\bf a}^8 {\bf t}^5+{\bf a}^8 {\bf q}^{-2} {\bf t}^5+{\bf a}^8 {\bf q}^{-4} {\bf t}^3+{\bf a}^{10} {\bf t}^8
```
**$[2]$-colored trefoil** (AS-long:2641–2651), normalized by $Z_{\tableau{2}}(\bigcirc)$:
```latex
\frac{Z_{\tableau{2}}(K_{2,3})}{Z_{\tableau{2}}(\bigcirc)} \Big|_{SU(N)} = \frac{q}{t} \ T_{\tableau{2}}^{-6} {\cal P}_{\tableau{2}}(K_{2,3})
{\cal P}_{\tableau{2}}(K_{2,3}) = q t^{2 N-3}+q^3 t^{2 N-2}-q^3 t^{3 N-3}-q^4 t^{3 N-3}+
+ q^4 t^{2 N-2}-q^5 t^{3 N-2}+q^5 t^{2 N-1}-q^6 t^{3 N-2}+q^6 t^{4 N-3} =
= {\bf q}^4 {\bf t}^6 {\bf a}^4 + {\bf q}^{-4} {\bf a}^4+{\bf q}^2 {\bf t}^4 {\bf a}^4+{\bf t}^{8} {\bf q}^8 {\bf a}^4+{\bf q}^8 {\bf t}^{11} {\bf a}^6+{\bf q}^6 {\bf t}^{9} {\bf a}^6+{\bf t}^5 {\bf a}^6+{\bf q}^2 {\bf t}^7 {\bf a}^6+{\bf q}^6 {\bf t}^{12} {\bf a}^8
```
AS-short trefoil and $(4,5)$ (AS-short:936–948, AS-short conventions):
```latex
\frac{Z(S^3, K_{2,3}, \tableau{1})}{Z(S^3, \bigcirc, \tableau{1})} = {\bf a} {\bf q}^{-1} + {\bf a} {\bf q} {\bf t}^2 + {\bf a}^2 {\bf t}^3.
={\bf a}^3 {\bf t}^4 + {\bf a}^3 {\bf q}^3 {\bf t}^6 + {\bf a}^3 {\bf q} {\bf t}^4 + {\bf a}^3 {\bf q}^{-1} {\bf t}^2 + {\bf a}^3 {\bf q}^{-3}
+ {\bf a}^4 {\bf q}^2 {\bf t}^7 +
+ {\bf a}^4 {\bf q} {\bf t}^7 + {\bf a}^4 {\bf t}^5 + {\bf a}^4 {\bf q}^{-1} {\bf t}^5 + {\bf a}^4 {\bf q}^{-2} {\bf t}^3 + {\bf a}^5 {\bf t}^4.
```
The AS-short appendix (AS-short:1563–1600) gives $(3,7)$, $(4,7)$, $(5,8)$, $(7,10)$, … normalized to $1+O(\cdot)$ (long; read them in the file).

DMMSS $T[2,n]$, $[1]$ (DMMSS:2828–2832), in DMMSS bold variables:
```latex
P^{T[2,3]}_{[1]}=\dfrac{\{A\} t}{\{t\} A q^3}
\Big({q}^{2}(-A^{2})+{q}^{2}{t}^{2}+1\Big)=\bf{{\dfrac {{a}^{2}t+1}{ \left(
-1+{q}^{2} \right) {a}^{2}{t}^{4}q}} \Big(
{q}^{4}{t}^{2}+{q}^{2}{t}^{3}{a}^{2}+1 \Big)}
```
DMMSS $T[2,3]$, $[2]$ (DMMSS:6424–6432), unreduced, with unknot $M_{[2]}^*$:
```latex
P^{T[2,3]}_{[2]}=\dfrac{\{A\}\{qA\} t^2}{\{t\}\{tq\} A^2 q^{18}}
\Big(q^{10}A^4+(q^{10}t^2+q^8t^2+q^6+q^4)(-A^2)+q^8t^4+q^6t^2+q^4t^2+1\Big)
=\\= \bf{\frac { \left( {a}^{2}t+1 \right)  \left( {q}^{2}{t}^{3}{a}^{2}+1
\right)}{{q}^{17}
\mbox{} \{q\}  \{{q}^{2}t\} {a}^{4}{t}^{21}}\left( {q}^{10}{t}^{12}{a}^{4}+
\left( {t}^{5}{q}^{4}+{t}^{11}{q}^{12}+{t}^{9}{q}^{10}+{t}^{7}{q}^{6} \right)
{a}^{2}
\mbox{}+{t}^{6}{q}^{8}+1+{q}^{6}{t}^{4}+{q}^{12}{t}^{8} \right)}
```
DMMSS:6465 states: "Our answer for $P^{T[2,3]}_{[2]}$ is in full agreement with the known results \cite{ASh}". The bracket in DMMSS bold variables is ${\bf a}^4{\bf q}^{10}{\bf t}^{12} + {\bf a}^2({\bf q}^4{\bf t}^5+{\bf q}^6{\bf t}^7+{\bf q}^{10}{\bf t}^9+{\bf q}^{12}{\bf t}^{11}) + 1+{\bf q}^6{\bf t}^4+{\bf q}^8{\bf t}^6+{\bf q}^{12}{\bf t}^8$. Up to a monomial, this equals AS's ${\cal P}_{\tableau{2}}$.

GN $t=q^{-1}$ check (GN:1536–1541): $\TCP_{n,m}(0,q,q^{-1}) = \frac {[m+n-1]!}{[m]!  [n]!}$.

DMMSS large-$n$ limit (DMMSS:1236–1238):
```latex
q^{(m-1)n}{\cal P}^{[m,n]}_{[1]}\ \stackrel{n\to\infty}{=}\ {M_{[m]}^*\over
M_{[1]}^*}
```

---

## 5. Cables / iterated torus knots

None of these five papers discusses cables, iterated torus knots or satellites. A grep for "cabl", "iterated" and "satellite" finds nothing relevant; the only hits are bibliography/footnote text in DMMSS:2501 and MMSS:693, about Shakirov's $\gamma$-ansatz. Related statements:

- AS-long:2050–2054: links of torus knots via $\langle 0 | {\cal O}^{(n,m)}_{R_i} S{\cal O}^{(0,1)}_{R_p}|0\rangle$ (torus knot linked with an unknot/core). This is the operator building block for a satellite/cable approach.
- AS-long and AS-short restrict the refined theory to knots preserved by the $U(1)$ (Seifert) action, i.e. torus knots/links in $S^3$ (AS-long:1965). Cables are not Seifert-fibered in this sense, so AS gives no prescription for them.
- GN: $W^\lambda_{n,m}$ acting on $|P_\mu\rangle$ gives the vector of the torus knot linked with a $\mu$-colored meridian, $v^\lambda_{n,m,\mu}$ (GN:479–481). It is defined, but GN does not use it for cables.
- DMMSS Sec. 4.4 (DMMSS:1613ff) "Non-torus knots": the evolution method for $5_2$ / $10_{139}$ works, but for $4_1$ "the same procedure leads to a superpolynomial with some negative coefficients" (DMMSS:494–499).
