# DAHA superpolynomials of iterated torus knots — extracted formulas and test data

Sources (LaTeX, line numbers refer to the files below, CRLF stripped):

| key | arXiv | file | content |
|---|---|---|---|
| **CD14** | 1408.4348 (Cherednik–Danilenko, "DAHA and iterated torus knots", = [ChD1]) | `1408.4348/src/cable.tex` | definitions for knots, all knot examples |
| **CD15** | 1509.08351 (Cherednik–Danilenko, "DAHA approach to iterated torus links", = [ChD2]) | `1509.08351/src/links.tex` | links; hat-normalization; Connection Conjecture for links |
| **C24** | 2410.14703 (Cherednik 2024, superpolynomials / motivic) | `2410.14703/src/superpolynoms.tex` | trefoil computation, cable dictionary, a few cable examples |
| **CP16** | 1605.00978 (Cherednik–Philipp, flagged Jacobian factors) | `1605.00978/src/jacobian-fin-full.tex` | positivity remarks, more cable superpolynomials |

Everything inside ```latex blocks is copied verbatim (by script, from the line ranges given).
Text outside the blocks is my summary. "[checked]" means I re-verified it numerically
(scratch parser: LaTeX -> Laurent polynomial in q,t,a; see Section 6).

---------------------------------------------------------------------------------------------------

## 0. Executive summary

* **Definition (knots).** For Newton-pair data $\vec r=(r_1,\dots,r_\ell)$, $\vec s=(s_1,\dots,s_\ell)$,
  $$JD_{\vec r,\vec s}(b;q,t)=\Big\{\hat\gamma_1\Big(\cdots\hat\gamma_{\ell-1}\Big(\big(\hat\gamma_\ell(P_b)/P_b(q^{-\rho_k})\big)\!\Downarrow\Big)\!\Downarrow\cdots\Big)\Big\}_{ev},$$
  $\gamma_i\in SL_2(\mathbb Z)$ any matrix with **first column $(r_i,s_i)^T$**, lifted to a word in $\tau_\pm$
  ($\tau_+\mapsto\begin{pmatrix}1&1\\0&1\end{pmatrix}$, $\tau_-\mapsto\begin{pmatrix}1&0\\1&1\end{pmatrix}$), $H\Downarrow:=H(1)$ (action on the polynomial representation).
  **Order:** the *last* pair $(r_\ell,s_\ell)$ = the *outermost* cable acts first on $P_b$; the *first* pair
  $(r_1,s_1)$ = the companion torus knot is applied last, just before the coinvariant $\{\cdot\}_{ev}$.
  Normalization: divide by $P_b(q^{-\rho_k})$ (= $P_b(t^{-\rho})$) once, at the innermost step (reduced invariant);
  then tilde/hat-normalize (divide by the lowest monomial, so constant term 1 for algebraic knots).
  Framing factors are dropped: everything is "up to $q^\bullet t^\bullet$".
* **Superpolynomial.** $\mathcal H_{\vec r,\vec s}(b;q,t,a=-t^{N+1})=\widetilde{JD}^{A_N}_{\vec r,\vec s}(b;q,t)$ for all $N\ge n-1$ (algebraic case; up to $\pm q^\bullet t^\bullet$ in general).
* **Topology.** $a_1=s_1,\ a_i=a_{i-1}r_{i-1}r_i+s_i$; knot $=Cab(a_\ell,r_\ell)\cdots Cab(a_2,r_2)\,T(r_1,s_1)$,
  where $Cab(a,r)(K)$ has **winding number $r$** (along $K$) and **$a$ meridional turns, measured in the zero (Seifert) framing of $K$**,
  i.e. $Cab(a,r)K$ = the standard $(r,a)$-cable $K_{r,a}$. So for the $(p_2,q_2)$-cable of $T(p_1,q_1)$:
  $\vec r=(p_1,p_2)$, $\vec s=(q_1,\ q_2-p_1q_1p_2)$. Algebraic iff all $r_i,s_i>0$.
  Trefoil cables: $(2,13)\to[(3,2),(2,1)]$, $(2,15)\to[(3,2),(2,3)]$, $(2,3)\to[(3,2),(2,-9)]\equiv[(3,-2),(2,9)]$,
  $(2,1)\to[(3,2),(2,-11)]\equiv[(3,-2),(2,11)]$ (CD14 uses the $r_2=-2$ versions). [checked via Alexander polynomials]
* **Variables.** DAHA $(q,t,a)$ vs. standard/DGR $(q_{st},t_{st},a_{st})$: $t=q_{st}^2,\ q=q_{st}^2t_{st}^2,\ a=a_{st}^2t_{st}$.
  Trefoil: $\mathcal H_{3,2}=1+qt+aq$ $\Rightarrow$ $1+q_{st}^4t_{st}^2+a_{st}^2q_{st}^2t_{st}^3$ $=(a^2q^{-2})^{-1}\cdot$DGR$(a^2q^{-2}+a^2q^2t^2+a^4t^3)$. [checked]
* **HOMFLY-PT:** $\mathcal H(q,t\mapsto q,a\mapsto -a)=\widetilde{HOM}$ in the variables $(q,a)=(q_{st}^2,a_{st}^2)$. Alexander: $\mathcal H(q,q,-1)$.

---------------------------------------------------------------------------------------------------

## 1. Definition of the DAHA-Jones polynomial and of the DAHA superpolynomial

### 1.1 Automorphisms $\tau_\pm$, matrices (CD14 §1.3)

`1408.4348/src/cable.tex:513-535`
```latex
\subsection{\bf The automorphisms}\label{sect:Aut}
The following maps can be (uniquely) extended to
an automorphism of $\HH\,$, where 
$q^{1/(2\mm)}$ must be added to $\Z_{q,t}$
(see \cite{C101}, (3.2.10)-(3.2.15)):
\begin{align}\label{tauplus}
& \tau_+:\  X_b \mapsto X_b, \ T_i\mapsto T_i\, (i>0),\
\ Y_r \mapsto X_rY_r q^{-\frac{(\om_r,\om_r)}{2}}\,,
\\
& \tau_+:\ T_0\mapsto  q^{-1}\,X_\vth T_0^{-1},\
\pi_r \mapsto q^{-\frac{(\om_r,\om_r)}{2}}X_r\pi_r\
(r\in O'),\notag\\
& \label{taumin}
\tau_-:\ Y_b \mapsto \,Y_b, \ T_i\mapsto T_i\, (i\ge 0),\
\ X_r \mapsto Y_r X_r q^\frac{(\om_r,\om_r)}{ 2},\\
&\tau_-(X_{\vth})= 
q T_0 X_\vth^{-1} T_{s_{\vth}}^{-1};\ \
\si\equal \tau_+\tau_-^{-1}\tau_+\, =\,
\tau_-^{-1}\tau_+\tau_-^{-1},\notag\\
&\si(X_b)=Y_b^{-1},\   \si(Y_b)=
T_{w_0}^{-1}X_{b^\iota}^{-1}T_{w_0},\ \si(T_i)=T_i (i>0).
\label{taux}
\end{align}
```

`1408.4348/src/cable.tex:566-601`
```latex
Let us list the matrices corresponding to the automorphisms and
anti-automorphisms above upon the natural projection 
onto $GL_2(\Z)$, corresponding to  
$t^{\frac{1}{2\mm}}_\nu=1=q^{\frac{1}{2\mm}}$. 
The matrix {\tiny 
$\begin{pmatrix} \al & \be \\ \ga & \de\\ \end{pmatrix}$}
will represent the map $X_b\mapsto X_b^\al Y_b^\ga,
Y_b\mapsto X_b^\be Y_b^\de$ for $b\in P$. One has:
\smallskip

$\!\!\!\tau_+\mapsto$ 
{\tiny 
$\begin{pmatrix}1 & 1 \\0 & 1 \\ \end{pmatrix}$},\ 
$\tau_-\mapsto$ 
{\tiny 
$\begin{pmatrix}1 & 0 \\1 & 1 \\ \end{pmatrix}$},\
$\si\mapsto$ 
{\tiny 
$\begin{pmatrix}0 & 1 \\-1 & 0 \\ \end{pmatrix}$},\
$\vph\mapsto$ 
{\tiny 
$\begin{pmatrix}0 & -1 \\-1 & 0 \\ \end{pmatrix}$},\
$\eta\mapsto$ 
{\tiny 
$\begin{pmatrix}-1 & 0 \\0 & 1 \\ \end{pmatrix}$}.\
\smallskip

{\sf Projective $GL_2$.}
We define the {\em projective\,} $GL_{\,2}^\wedge(\Z)\,$ as
the group generated by $\tau_{\pm}, \eta$ subject to the
relations $\tau_+\tau_-^{-1}\tau_+\, =\,
\tau_-^{-1}\tau_+\tau_-^{-1},\ \eta^2=1$ and 
$\eta\tau_{\pm}\eta=\tau_{\pm}^{-1}.$ The span 
of $\tau_{\pm}$ is the projective $PSL_{\,2}^\wedge(\Z)\,$
(due to Steinberg), which is isomorphic to the braid 
group $B_3$.
```

Same formulas in C24 (`2410.14703/src/superpolynoms.tex:1021-1103`); there the group is called $PSL_2^\wedge(\mathbb Z)$ and only $\tau_\pm$ are used (no $\eta$).

**A_1 version** (CD14 §3.2, also C24 trefoil example):

`1408.4348/src/cable.tex:3686-3707`
```latex
The {\em double affine Hecke algebra
$\HH$\,} is generated by invertible elements
$Y=Y_{\om}, T=T_1, X=X_{\om}$
subject to the quadratic relation $(T-t^{1/2})(T+t^{-1/2})=0$
and the cross-relations:
\begin{align}\label{dahaone}
&TXT=X^{-1},\ T^{-1}YT^{-1}=Y^{-1},\ Y^{-1}X^{-1}YXT^2q^{1/2}=1.
\end{align}
%Setting $\pi\equal YT^{-1}$, the second relation 
%becomes $\pi^2=1$.
The field of definition will be $\Q(q^{1/4},t^{1/2})$.
% though the ring $\Z[q^{1/4},t^{1/2}]$ is sufficient 
%for DAHA-Jones polynomials.
Here $q^{\pm 1/4}$ is needed for the automorphisms $\tau_{\pm}$:
\begin{align}\label{taudef}
&\tau_+(X)=X,\ \ \tau_+(T)=T,\ \ \tau_+(Y)=q^{-1/4}XY,\ \\
%\tau_+(\pi)=q^{-1/4}X\pi,\\
%\label{tau-def}
&\tau_-(Y)=Y,\ \ \,\tau_-(T)=T,\ \ \,\tau_-(X)=q^{1/4}YX.\ 
%\tau_-(\pi)=\pi.
\notag
\end{align}
```

### 1.2 Coinvariant (evaluation) $\{\cdot\}_{ev}$, polynomial representation, $P_b$, evaluation formula (CD14 §1.4)

`1408.4348/src/cable.tex:605-625`
```latex
\subsection{\bf Macdonald polynomials}
Following \cite{C101}, 
we use the PBW Theorem to express any $H\in \HH$ in the form 
\,$\sum_{a,w,b} c_{a,w,b}\, X_a T_{w} Y_b$\, for $w\in W$,
$a,b\in P$ (this presentation is unique). Then we substitute:
\begin{align}\label{evfunct}
\{\,\}_{ev}:\ X_a \ \mapsto\  q^{-(\rho_k,a)},\ 
Y_b \ \mapsto\  q^{(\rho_k,b)},\ 
T_i \ \mapsto\  t_i^{1/2}. 
\end{align}

{\sf Polynomial representation.}
The functional $\,\HH\ni H\mapsto \{H\}_{ev}$, 
called {\em coinvariant\,}, acts via the projection 
$H\mapsto H\!\Downarrow \equal H(1)$ of $\HH\,$
onto the {\em polynomial representation \,}$\v$, which is
the $\HH$\~module induced from the one-dimensional
character $T_i(1)=t_i^{-1/2}=Y_i(1)$ for $1\le i\le n$ and
$T_0(1)=t_0^{-1/2}$. Recall  that $t_0=t_{\sht}$; 
see \cite{C101,CJ}.
\smallskip
```

`1408.4348/src/cable.tex:681-697`
```latex
\begin{align}\label{macdsym}
&P_b(X^{-1})\,=\,P_{b^\iota}(X)\,=\,
P_b^\star(X),\ \, P_{b}(q^{-\rho_k})=
P_{b}(q^{\rho_k})\\
\label{macdeval}
=&(P_{b}(q^{-\rho_k}))^\star=
q^{-(\rho_k,b)}
\prod_{\al>0}\prod_{j=0}^{(\al^{\!\vee},b)-1}
\Bigl(
\frac{
1- q_\al^{j}t_\al X_\al(q^{\rho_k})
 }{
1- q_\al^{j}X_\al(q^{\rho_k})
}
\Bigr).
\end{align}
Recall that $\iota(b)=b^\iota=-w_0(b)$ for $b\in P$.
```

`1408.4348/src/cable.tex:726-731`
```latex
The product $\Pi_R^b$ becomes the {\em Poincar\'e polynomial\,} 
$\Pi_R$ for generic $b\in P_+$. See formula (3.3.23) from 
\cite{C101}.
We set  $P_b^\circ\equal P_b/P_b(q^{-\rho_k})$ for $b\in P_+$,
this is the so-called {\em spherical normalization\,}.
\medskip
```

C24 also writes $\{\mathcal P\}_{ev}=\mathcal P(t^{-\rho})=\mathcal P(t^{\rho})$ for symmetric $\mathcal P$ (`C24:1908-1909`), and $J$-polynomials $J_\lambda=h_\lambda P_\lambda$ (needed only for links):

`2410.14703/src/superpolynoms.tex:1289-1296`
```latex
The {\em $J$\~polynomials\,}
are as follows:
\begin{align}\label{P-arms-legs}
J_\la\equal h_\la P_b \for \la=\la(b),\,
h_\la=\prod_{\Box\in\la}
(1-q^{arm(\Box)}t^{leg(\Box)+1});
\end{align}
they are $q,t$\~integral.
```

`2410.14703/src/superpolynoms.tex:1324-1333`
```latex
The following the simplest instance of the 
{\em $a$-stabilization\,}; $J_\la(t^{-\rho})$ for $A_n$ is the value
of the formula in (\ref{stabevalx}) for  $a\!=\!-t^{n+1}$:\,
\begin{align}\label{stabevalx}
J_\la(t^{-\rho})\!=\!
(a^2)^{-\frac{|\la|}{4}}\,t^
{\sum_{i=0}^{n-1}(i+1/2)\la_{i+1}}\,
&\prod_{p=1}^n \prod_{j=0}^{\la_{p\,}-1}
\Bigl( 1\!+\!q^{j}\, a\, t^{-p+1}\Bigr).
\end{align}
```

### 1.3 Tilde-normalization and the DAHA-Jones polynomial of an iterated torus knot (CD14 §2.3)

`1408.4348/src/cable.tex:1112-1174`
```latex
Torus knots $T(\rr,\ss)$ are naturally represented   
by $\ga_{\rr,\ss}\in GL_{\,2}(\Z)$ with the 
first column $(\rr,\ss)^{tr}$ ($tr$ is the transposition)
for $\,\rr,\ss\in \Z$ 
assuming that \,gcd$(\rr,\ss)=1$. 
Then let $\hat{\ga}_{\rr,\ss}\in GL_{\,2}^{\wedge}(\Z)$ be
any pullback of $\ga_{\rr,\ss}$.

Note that $\,(\rr,\ss)\,$ can be obviously lifted
to $\,\ga\,$ of determinant $1$ and, accordingly, 
to $\hat{\ga}\,$ from the subgroup $PSL^\wedge_{\,2}(\Z)$
generated by $\{\tau_{\pm}\}$, i.e. without using
$\eta$. This is actually sufficient for the construction 
of the DAHA-Jones polynomials below.
% using $GL_2(\Z)$ 
%instead of $PSL_2(\Z)\subset GL_2(\Z)$ 
%is actually equivalent to the conjugation relation 
%from (\ref{jones-sym-iter}).
\smallskip

For a polynomial $F$ in terms of 
fractional powers of $q$ and $t_\nu$, 
the {\em tilde-normalization}
$\tilde{F}$ will be the result of the division of $F$
by the lowest $q,t_\nu$\~monomial, assuming that it 
is well defined. We put $q^\bullet t^\bullet$ for
a monomial factor (possibly fractional)
in terms of $q,t_\nu$.
%thus $\tilde{F}$ contains only non-negative powers (possibly
%rational) of $q,t_\nu$. 
%We will constantly
%use the tilde-normalization below, though not all properties of 
%the DAHA-Jones polynomials really require such a 
%normalization (see \cite{CJJ}). 


\begin{theorem}\label{JONITER}
Recall that $H\mapsto H\!\!\Downarrow\, \equal H(1)$, where the
action of $H\in \HH$ in $\v$ is used.
Given two sequences $\vec\rr,\vec\ss$ of length $\ell$ as 
in (\ref{iterrss}), a reduced root system $R$ and a weight 
$b\in P_+\,$, we lift $(\rr_i,\ss_i)^{tr}$ to
$\ga_i$ and then to $\hat{\ga_i}\in GL_{\,2}^{\wedge}(\Z)$ as above.
The {\sf DAHA-Jones polynomial} is
\begin{align}\label{jones-dit}
&J\!D_{\,\vec\rr,\,\vec\ss\,}^{R}(b\,;\,q,t)\ =\ 
J\!D_{\,\vec\rr,\,\vec\ss\,}(b\,;\,q,t)\, \equal\\ 
\Bigl\{\hat{\ga_{1}}&\Bigr(
\cdots\Bigl(\hat{\ga}_{\ell-1}
\Bigl(\bigl(\hat{\ga}_\ell(P_b)/P_b(q^{-\rho_k})\bigr)\!\Downarrow
\Bigr)\!\Downarrow\Bigr) \cdots\Bigr)\Bigr\}_{ev}.\notag
\end{align}
It does not depend
on the particular choice of the lifts \,$\ga_i\in GL_2(\Z)$
and $\hat{\ga}_i\in GL_{\,2}^{\wedge}(\Z)$\, for $1\le  i
\le\ell.$ Assuming that $\rr_i,\ss_i>0$ for all $i$,
i.e. $\t(\vec\rr,\vec\ss)$ is algebraic,
the tilde-normalization 
$\tilde{J\!D}_{\vec\rr,\vec\ss}\,(b\,;\,q,t)$ is well defined
and is a polynomial in terms of $\,q,t_\nu$ with
constant term $1$. 
\sq
\end{theorem}
```

The same formula is (2.12) of CD14, repeated verbatim in CD15 and C24:

`1509.08351/src/links.tex:4984-4999`
```latex
{\sf Iterated knots.\,}
In the case of torus iterated {\em knots\,} (when there is 
only one path) and in the absence of $'\!\l$,
we arrive at formula (2.12) from \cite{ChD}:
\begin{align}\label{jones-ditx}
& J\!D_{\vec\rr,\vec\ss}
(b;q,t)\! =\!  
\Bigl\{\hat{\ga_{1}}\Bigr(
\cdots\Bigl(\hat{\ga}_{\ell-1}
\Bigl(\bigl(\hat{\ga}_\ell(P_b)/
P_b(q^{-\rho_k})\bigr)\!\Downarrow
\Bigr)\!\Downarrow\Bigr) \cdots\Bigr)\Bigr\}_{ev}.
\end{align}
It includes only one $b\in P_+$ and therefore does not depend
on the choice of the integral form $\hat{P}_b$ of $P_b$.

```

`2410.14703/src/superpolynoms.tex:1911-1925`
```latex
In the case of iterated torus {\em knots\,},
when $\l$ is one path with one $b$ (and there is no
$'\!\l$), 
we arrive at formula (2.12) from \cite{ChD1}:
\begin{align}\label{jones-ditx}
& J\!D_{\vec\rr,\vec\ss}
(b;q,t)\! =\!
\Bigl\{\hat{\ga_{1}}\Bigr(
\cdots\Bigl(\hat{\ga}_{\ell-1}
\Bigl(\bigl(\hat{\ga}_\ell(P_b)/
P_b(t^{-\rho})\bigr)\!\Downarrow
\Bigr)\!\Downarrow\Bigr) \cdots\Bigr)\Bigr\}_{ev}.
\end{align}
Only $P_b^\circ=P_b/P_b(t^{-\rho})$ are needed in this case
(no $J$-polynomials). 
```

C24's choice of lift (only $PSL_2$):

`2410.14703/src/superpolynoms.tex:1825-1832`
```latex
We represent torus knots $T(\rr,\ss)$
by the matrices $\ga[\rr,\ss]=\ga_{\rr,\ss}\in PSL_{\,2}(\Z)$
with the
first column $(\rr,\ss)^{tr}$ ($tr$ is the transposition)
for $\,\rr,\ss\in \N$ such that \,gcd$(\rr,\ss)=1$.
Let $\hat{\ga}_{\rr,\ss}\in PSL_{\,2}^{\wedge}(\Z)$ be
{\em any\,} pullback of $\ga_{\rr,\ss}$ to 
$PSL_{\,2}(\Z)$; see Section \ref{sect:Aut}.
```

**Iteration order (my reading, unambiguous from the formula):** start with $P_b$ (or $E_b$, see C24 below);
apply $\hat\gamma_\ell$ (last Newton pair = outermost cable), divide by $P_b(q^{-\rho_k})$, project $\Downarrow$ (act on $1$) to get a
Laurent polynomial in $X$; apply $\hat\gamma_{\ell-1}$ to *that* polynomial viewed as an element of $\mathcal H$
(multiplication operator), project $\Downarrow$; ... ; apply $\hat\gamma_1$ (companion torus knot $T(r_1,s_1)$) and finally take
$\{\cdot\}_{ev}$. Topologically the first cabling is the companion (see 2.1: "the first iteration ... is for $\{a_1,r_1\}$ (not for the last pair!)").

Explicit lifts used in CP16 for $\vec r=\{3,3\},\vec s=\{2,4\}$:

`1605.00978/src/jacobian-fin-full.tex:5204-5216`
```latex
\subsection{\bf Superpolynomial}\label{sec:sup6-9-22}
Let us provide the corresponding DAHA superpolynomial,
which is for $\vec\rr=\{3,3\},\,
\vec\ss=\{2,4\}$. Recall that
the DAHA construction requires the matrices from
$PSL_2(\Z)$ with the first columns $(\rr_i,\ss_i)^{tr}$,
namely:
{\small $
\begin{pmatrix}3,*\\2,*\end{pmatrix}\,=\,
\tau_+\tau_-^2\,,\ \ 
 \begin{pmatrix}3,*\\4,*\end{pmatrix}\,=\,
\tau_-\tau_+^2\tau_-\,.
$} Then:
```

(So $\gamma_{3,2}=\tau_+\tau_-^2=\begin{pmatrix}3&1\\2&1\end{pmatrix}$, $\gamma_{3,4}=\tau_-\tau_+^2\tau_-=\begin{pmatrix}3&2\\4&3\end{pmatrix}$.)

Non-symmetric simplification and worked trefoil computation (C24 §2.4):

`2410.14703/src/superpolynoms.tex:2257-2346`
```latex
\subsection{\bf Some examples}\label{sec:exa}
An important simplification of the construction of the $J\!D$-polynomials
and $\hat{\h}$-polynomials for iterated torus {\em knots}, is due to
the usage of nonsymmetric Macdonald polynomials $E_\la$ for dominant $\la$
instead of $P_\la$. One simply replaces $P_\la^\circ
=P_\la/P_\la(t^{-\rho})$ by  $E_\la^\circ
=E_\la/E_\la(t^{-\rho})$; no other changes are
necessary. Recall, that $E_\la$ are eigenvectors of $Y$-operators.
Also, one can do the calculations for $GL_{n+1}$, which
is somewhat simpler than those for $SL_{n+1}$ (for $A_n$).  
\vskip 0.2cm


{\sf The case of trefoil}.
Let us calculate $\hat{\h}_{3,2}$ for uncolored trefoil.
We begin with the case of $A_1$. 

In this case,
$\HH$ is generated by $X^{\pm1},Y^{\pm1},T$ subject to  
group relations $TXTX=1=TY^{-1}TY^{-1}, Y^{-1}X^{-1}YXT^2=q^{-1/2}$
and the quadratic one 
$(T-t^{1/2})(T+t^{-1/2})=1$. The action of $\tau_{\pm}$ is:
\begin{align*}
&\tau_+\!:\!
Y\!\mapsto\! q^{-\frac{1}{4}}XY,\ X\!\mapsto\!X,\ T\!\mapsto\!T,\ 
\tau_-\!:\! X\!\mapsto\! q^{\frac{1}{4}}YX,\  
Y\!\mapsto\!Y,\  T\!\mapsto\!T.
\end{align*}
The polynomial representation is given by the formulas
$T\mapsto t^{1/2}s\!+\! 
\frac{t^{1/2}\!- t^{-1/2}}{X^2-1}(s\!-\!1)$, 
$X\mapsto X,\ \, Y\mapsto spT$, where $ s(X)=X^{-1},\,
\, p(X)=q^{1/2}X.$
\vskip 0.2cm

By $\,\sim$\,, 
we mean ``\,up to $q^{\bullet}t^{\bullet}$\,"
in the following calculation: 
\begin{align*}
&J\!D_{3,2}\!=\!
\{\tau_+\tau_-^2(X)\}_{ev}\!\sim\!\{(XY)(XY)X(1)\}_{ev}\!\sim\!
\{Y(X^2)\}_{ev}\\
&=t^{-\frac{1}{2}}q^{-1}X^2-
t^{\frac{1}{2}}+t^{-\frac{1}{2}}|_{X^2\mapsto t^{-1}}\sim
1+qt-qt^2,
\end{align*}
where we replaced  $P_1$ 
by the {\em nonsymmetric} Macdonald polynomial $E_1$.
One has:  $Y(X)\!=\!(qt)^{-\frac{1}{2}}X$ and $E_1=X$.
Recall that 
$\{H\}_{ev}\!\equal\! H(1)(X\!\mapsto\! t^{-\rho})$.
The {\em Jones hat-polynomial} is
$\hat{J\!D}_{3,2}(q\mapsto t)= 1\!+\!t^2\!-\!t^3.$


\vskip 0.2cm

For $GL_{n+1}$, the corresponding
$\HH$ is generated by $X_i^{\pm 1}, Y_j^{\pm 1}, T_k$,  where 
$1\le i,j\le n+1$ and
$1\le k \le n$ for pairwise commutative $\{X_i\}$ and 
$\{Y_j\}$. One has:   
 $\tau_+(Y_1)\!=\!q^{-1/2} 
X_1 Y_1$,\  $\tau_-(X_1)\!=\!q^{+1/2}
Y_1 X_1$ and so on. The action of $Y_1$ in the polynomial
representation is via
the formula 
$Y_1\!=\!\pi T_{n}\cdots T_1$, where 
$\pi: X_1\!\mapsto\! X_2, X_2\!\mapsto\! X_3,
\ldots, X_n\!\mapsto\! q^{-1} X_1$. Here $T_k$ are
those for $A_n$, and $X_{\al}=X_i X_j^{-1}$ for
$\al=\ep_i-\ep_j$. 
These formulas are almost equally simple 
for any $Y_i$, but we need only $X_1, Y_1$.

The calculation of $J\!D$
remains practically the same as it was for $A_1$.
The polynomial $E_1=X_1$ can be used instead of $P_1$.
We obtain that  $\hat{J\!D}^{A_n}_{3,2}
=1+qt-qt^{n+1}$ and 
$\hat{\h}_{3,2}=1+qt+aq$. Indeed, the values of the latter 
at $a=-t^{n+1}$ are $\hat{J\!D}^{A_n}_{3,2}$.

Note that the relations
$\hat{\h}_{3,2}(a\!\mapsto\! -t)=1$ and    
$\hat{\h}_{3,2}(a\!\mapsto\! -t^{2})=1+qt-qt^2$ (only two)
are sufficient to fix $\hat{\h}_{3,2}$ uniquely using
that $deg_a=1$. Generally, 
$deg_a\hat{\h}_{\rr,\ss}(\la)\!=\!
|\la|\bigl(\text{Min}(\rr,\ss)\!-\!1\bigr)$. 
```

The DAHA-Jones polynomial of a single path with hat-normalization (C24, Theorem THM-integr-Jones):

`2410.14703/src/superpolynoms.tex:1940-1962`
```latex
\begin{theorem}\label{THM-integr-Jones}
The invariants  $J\!D_{\l,\,'\!\l}\,$ 
are polynomials in terms of $q,t$
up to a factor $q^\bullet t^\bullet$,
where the powers $\bullet$ can be rational.
Modulo such factors,
it does not depend
on the particular choice of the lifts \,$\ga\in PSL_2(\Z)$
and $\hat{\ga}\in PSL_{\,2}^{\wedge}(\Z)$\ for the vertices
of $\Up,'\Up$. 

Up to the $q^\bullet t^\bullet$\~equivalence
and a possible factor from $\Q$, the following
{\em\sf hat-normalization} $\hat{J\!D}_{\l,\,'\!\l}\,$ exists.
It is a $q,t$\~polynomial not divisible by
$\,q\,$ and by $\,t\,$ with integral
coefficients of its $q,t$\~monomials such that 
their total $GCD$ equals $1$, and  
the coefficient of the minimal pure
power of $t$ must be positive. For algebraic links colored
by rows $m\om_1$, this normalization becomes
$\hat{J\!D}_{\l,\,'\!\l}\,(q\!=\!0,t\!=\!0)=1$. \sq
\end{theorem}
```

### 1.4 Symmetries (CD14 Thm 2.2)

`1408.4348/src/cable.tex:1180-1234`
```latex
\begin{theorem} \label{MAINTHM} 
(i) The polynomials
$J\!D_{\vec\rr,\vec\ss}(b\,;\,q,t)$  
coincide up to $q^\bullet t^\bullet$
for the pairs from (\ref{iter-triv}-\ref{iter-s1}), i.e.
for $\{\vec{\rr},\vec{\ss}\}$ 
with isotopic $\t(\vec\rr,\vec\ss)$. Also,
(\ref{iter-sym}) corresponds to
\begin{align}\label{iter-sym-j}
&J\!D_{\,\vec{\rr},\,\vec{\ss}}\,(b\,;\,q,t)=
J\!D_{\{\ss_1,\rr_{2},\cdots,\rr_{\ell}\},
\{\rr_{1},\ss_2\cdots,\ss_{\ell}\}}\,(b\,;\,q,t),
\end{align}
and the following DAHA counterparts of
(\ref{iterorien}) and (\ref{iter2min}) hold:
\begin{align}
\label{jones-sym-iter} 
&J\!D_{\,\vec\rr,\,-\vec\ss}\,(b\,;\,q,t) = 
\bigr(J\!D_{\,\vec\rr,\,\vec\ss}\,
(b\,;\,q,t)\bigl)^\star \hbox{\,for\, } \star \hbox{\, 
from (\ref{eta-conj})},
\\
&J\!D_{\,\vec\rr,\,\vec\ss}=
J\!D_{\{\ldots,\rr_{i\!-\!1},
-\rr_i,\rr_{i\!+\!1},\ldots\},
\{\ldots,\ss_{i\!-\!1},-\ss_i,\ss_{i\!+\!1},\ldots\}}
\hbox{\,\ for\,\, } 1\!\le\! i\!\le\! \ell.\notag
\end{align}

(ii) Let us assume that for $b,c\in P_+$ and certain 
\,id\,$\neq w\in W$,  
\begin{align}\label{bcw-rel}
q^{\,(b+\rho_k-w(\rho_k)-w(c)\,,\,\al)\,}=1=
q_{\al}^{\,(b-w(c)\,,\,\al^{\!\vee})}\,
t_{\sht}^{\,(\rho^w_{\sht}\,,\,\al^{\!\vee})}\,
t_{\lng}^{\,(\rho^w_{\lng}\,,\,\al^{\!\vee})}
\end{align}
for any $\al\in R_+$, where we set\,
$\rho_{\nu}^w\equal w(\rho_{\nu})-\rho_{\nu}$. Then 
\begin{align}\label{jones-bc}
J\!D_{\,\vec\rr,\,\vec\ss}\,(b\,;\,q,t)=
J\!D_{\,\vec\rr,\,\vec\ss}\,(c\,;\,q,t)
\hbox{\, for such \, } q,\{t_\nu\} 
\end{align}
and any $\vec\rr,\vec\ss$. Also, for 
$b=\hbox{\small$\sum$}_{i=1}^n b_i \om_i$,
\begin{align}\label{jones-eval}
J\!D_{\,\vec\rr,\,\vec\ss}
\,\bigl(b\,;\,q\!=\!1,t
\bigr)\!=\!
\hbox{\small$\prod$}_{i=1}^n J\!D_{\,\vec\rr,\,\vec\ss}\,
(\om_i\,;\,q\!=\!1,t)^{b_i} \hbox{\, for any \,} \vec\rr,\vec\ss.
\end{align}
\vskip -1.25cm \sq
 \end{theorem}
```

### 1.5 Stabilization: the DAHA superpolynomial (CD14 Thm 2.3)

`1408.4348/src/cable.tex:1293-1381`
```latex
\begin{theorem}\label{STABILIZ}
We switch to $A_n$, setting $t=t_{\sht}=q^k$. 
Let us consider $P_+\ni b=$ $\sum_{i=1}^n b_i \om_i$ as 
a (dominant) weight for any $A_N$ with $N\ge n-1$,
where we set $\om_{n}=0$ upon the restriction to $A_{n-1}$.

(i) {\sf Stabilization.} Given $\t(\vec\rr,\vec\ss)$, there exists 
$\h_{\,\vec\rr,\,\vec\ss}\,(b\,;\,q,t,a)$,
a polynomial from $\Z[q,t^{\pm 1},a]$, satisfying the relations
%such that
%its coefficient of $a^0$ is tilde-normalized 
%(i.e. in the form $\sum_{u,v\ge 0}C_{u,v}q^u t^v$ 
%with $C_{0,0}=1$) and the following holds:
\begin{align}\label{jones-sup}
\h_{\,\vec\rr,\,\vec\ss}\,(b;q,t,a\!=\!-t^{N\!+\!1})\!=\!
\pm q^\bullet t^\bullet\, 
J\!D_{\,\vec\rr,\,\vec\ss\,}^{A_{N}}(b;q,t)
\hbox{\, for any } N\!\ge\! n\!-\!1, 
\end{align} 
and normalized as follows. The polynomial
$\,\h(a\!=\!0)\,$ is assumed minimal, i.e. not divisible by 
$\,q,\,t\,$ and 
any prime number; the coefficient of the minimal power of $t$ in
$\,\h(q\!=\!0,a\!=\!0)\,$ is selected positive.

(ii) {\sf Symmetries.}  The corresponding normalization of
$J\!D$\~polynomials, denoted by $ \tilde{J\!D}$, 
extends the tilde-normalization for algebraic knots
from Theorem \ref{JONITER}. For sufficiently large $N$
(for all $N\ge n-1$ if $\rr_i,\ss_i>0$),  
\begin{align}\label{tilde-sup}
\h_{\,\vec\rr,\,\vec\ss}\,(b;q,t,a\!=\!-t^{N\!+\!1})\!=\! 
\tilde{J\!D}_{\,\vec\rr,\,\vec\ss\,}^{A_{N}}(b;q,t).
\end{align} 
The symmetries from Part (i) of Theorem \ref{MAINTHM}
hold for $\h_{\,\vec\rr,\,\vec\ss}$ under such a normalization,
but (\ref{jones-sym-iter}); the latter holds  
up to $\,a^\bullet q^\bullet t^\bullet\,$
for $a^\star=a^{-1}$.
%For the latter we set 
%$\,a^\star= a^{-1}\,$ and then multiply 
%$\h_{\,\vec\rr,\,\vec\ss}^\star$ by $\,a^{deg}$ and proper
%$q^\bullet t^\bullet$,
%where\, deg=deg$_a$\, is the $a$\~degree of  
%$\h_{\,\vec\rr,\,\vec\ss}$.

(iii) {\sf Color Exchange.} Imposing (\ref{bcw-rel}), 
we consider $\,w\,$ 
%from Part (ii) of Theorem  \ref{MAINTHM} 
as an element of $\S_{N+1}$ for every
$N\ge n$ (the Weyl group for $A_N$) naturally acting in 
the corresponding $P$. Then for any \, $\vec\rr,\vec\ss$,
\begin{align}\label{jones-bca}
\h_{\,\vec\rr,\,\vec\ss}\,(b\,;\,q,t,a) = 
\h_{\,\vec\rr,\,\vec\ss}\,(c\,;\,q,t,a)
\hbox{\, for such \, } q,t.
\end{align}
In particular, let $w=s_i=(i,i+1)$ with $i<n$
and $\,b=c-\bigl(k+(c,\al_i)\bigr)\al_i$ for
a dominant $c$. 
Then the inequalities
$\,c_i/2\!\le\! -k\!\le\! c_i+\min\{c_{i\pm 1>0}\}\,$
for $k\!\in\! -\Z_+$ are sufficient.

(iv) {\sf Specialization and deg${}_a$}.
Making $q=1$, one has:
\begin{align}\label{jones-seval}
\h_{\,\vec\rr,\,\vec\ss}\,(b\,;\,q\!=\!1,t,a)\!=\!
\prod_{i=1}^n\h_{\,\vec\rr,\,\vec\ss}\,
(\om_i\,;\,q\!=\!1,t,a)^{b_i}
\hbox{\, for\, } b\!=\!\sum_{i=1}^n b_i \om_i.
\end{align}
Assuming  that $|\rr_1|>|\ss_1|$, the $a$-degree
deg${}_a \h_{\,\vec\rr,\,\vec\ss}\,(b\,;\,q,t,a)\,$
is no greater than $(|\ss_1\rr_2\cdots\rr_\ell|-1)$ times 
\, ord $(\la_b)$, the 
number of boxes in the Young diagram $\la_{b}$
associated with $b\in P_+$. 

(v) {\sf Super-duality \cite{GS,CJ,GoN}.} 
Up to a power of $q$ and $t$,
\begin{align}\label{iter-duality}
\h_{\,\vec\rr,\,\vec\ss}\,(\la\,;q,t,a)=q^{\bullet}t^{\bullet}
\h_{\,\vec\rr,\,\vec\ss}\,(\la^{tr}\,;t^{-1},q^{-1},a),
\end{align}
where we switch from (arbitrary) dominant
weights $\,b\,$ to the corresponding Young diagrams $\la=\la_b$,
and $\la^{tr}$ is the transposition of $\la$. 
\end{theorem} 
```

`1408.4348/src/cable.tex:1424-1438`
```latex
{\sf The $a$-degree.}
We expect that \,deg${}_a \h_{\,\vec\rr,\,\vec\ss}\,
(b\,;\,q,t,a)\,$
from Part $(iv)$ equals 
$\,(\ss_1\rr_2\cdots\rr_\ell-1)\,\hbox{ord}(\la_b)$
for algebraic knots (there is no such a
coincidence in general). 
Due to (\ref{jones-seval}), it suffices to check the
equality for the fundamental weights only.
This equality is likely to hold
for the HOMFLYPT polynomials of algebraic 
knots (the inequality from Part $(iv)$ holds for all
iterated knots), which would formally result in such
an equality for $\h$ assuming Part $(i)$ of 
Conjecture \ref{CONCONJ} below.
```

Hat-normalization (CD15, for links; for knots it coincides with the above):

`1509.08351/src/links.tex:5513-5530`
```latex
They are normalized as follows ({\sf the hat-normalization}).
 We multiply $\hat{\h}$
by $\pm q^\bullet t^\bullet$ to ensure that $\hat{\h}(a\!=\!0)$
is not divisible by 
$\,q\,$ and $\,t\,$ and that the coefficient of the 
minimal 
power of $t$ in $\hat{\h}(a\!=\!0)$ is positive. If 
$\hat{J\!D}^{A_m}$\~polynomials are considered under the
same hat-normalization, then 
relations (\ref{jones-sup-hat}) will automatically 
hold for sufficiently large $\,m$
without any correction factors $\pm q^\bullet t^\bullet$
(and one sufficiently large
$\,m$\, is actually sufficient to fix $\,\hat{\h}\,$ 
uniquely). The $a$\~stabilization and other claims here
hold for
$\hat{\h}_{\l,\,'\!\l^\vee}$; see (\ref{jones-hat-vee}).

```

C24 version of the stabilization theorem and the a-degree conjecture:

`2410.14703/src/superpolynoms.tex:2142-2172`
```latex
\begin{theorem}\label{STABILIZ}\cite[Theorem 2.3]{ChD2}
Given a pair $\{\l,\, '\!\l\}$
colored by $\mathbf{b}=(b^j),\, '\mathbf{b}=('\!b^j)$,
there exists a
unique polynomial from $\Z[q,t^{\pm 1},a]$
\begin{align}\label{h-polynoms-hat}
&\hat{\h}_{\l,\,'\!\l}\ =\
\hat{\h}_{(\vec\rr^j,\,\vec\ss^j),
('\vec\rr^j,\,'\vec\ss^j)}(\mathbf{b},\,
'\mathbf{b};\,q,t,a)
\end{align}
such that for any  $m\!\ge\! n\!-\!1$ and proper powers
of $q,t$ (possibly rational)\,:
\begin{align}\label{jones-sup-hat}
&\hat{\h}_{\l,\,'\!\l}(q,t,a\!=\!-t^{m+1})\,=\,
\pm\, q^\bullet t^\bullet\,
\hat{J\!D}^{A_m}_{\l,\,'\!\l}(q,t).
\end{align}
Here $\hat{\h}(a\!=\!0)$
is automatically 
hat-normalized. We note that this normalization becomes
$\hat{\h}(q\!=\!0,t\!=\!0,a\!=\!0)=1$ in Theorem
\ref{con:motknot} and further  coincidence
conjectures below.
 The correction factors 
$\pm\, q^\bullet t^\bullet$ are not needed for sufficiently
large $\,m$ for links and for any  $m\!\ge\! n\!-\!1$ in the case of 
iterated torus knots. Also,  the consideration of
one sufficiently large
$\,m$\, is sufficient to fix $\,\hat{\h}\,$
uniquely. \sq
```

`2410.14703/src/superpolynoms.tex:2205-2234`
```latex
We conjectured in \cite{ChD2} for algebraic links that
%\vskip -0.5cm
\begin{align}\label{deg-a-jj}
\hbox{deg}_a\hat{\h}_{\l,\,'\!\l}
=&\,\hbox{\small $\sum$}_{j=1}^\kap\
\min\{\,\rr^j_1,\ss^j_1\}\, \rr^j_2\cdots
\rr^j_{\ell^j}\,|\la^j|\!\\
&+
\hbox{\small $\sum$}_{j=1}^{\,'\!\kap}
\ \min\{\,'\rr^j_1,\,'\ss^j_1\}\,
\,'\rr^j_2\cdots
\,'\rr^j_{\,'\!\ell^j} \,|\,'\!\la^j|-\De\notag\\
\hbox{for\ } \De&=\!|\la^1\!\vee\!
\ldots\!\vee\!\la^\kap\!\vee\,
'\!\la^1\!\vee\!\ldots\!\vee\, '\!\la^{'\!\kap}|
\hbox{\, from \,(\ref{P-arms-legs-union})},
\notag
%\hbox{\, for \,} \hat{\h}^{min},\notag
%\De\!=\!|\la^{j_o}| \hbox{\, for \,} \hat{\h}^{j_o},\notag
\end{align}
where $|\la|$ is the number of boxes in $\la$.
A somewhat weaker statement can be justified. This is
compatible with the specialization $q\!=\!1$, which gives that
the $a$\~degree is no smaller than that in (\ref{deg-a-jj}).

The right-hand side of this
formula is the multiplicity
of the corresponding singularity generalized to the colored
case. For instance, $deg_a=|\la|(\ss-1)$ for torus knots
$T(\rr,\ss)$ colored by $\la$, where $\ss<\rr$. 
```

CP16 version (note: here the relation (jones-sup) is stated with no $\pm q^\bullet t^\bullet$ factor, algebraic case):

`1605.00978/src/jacobian-fin-full.tex:1031-1076`
```latex
\begin{theorem}\label{JONITER}
Let $R$ be a reduced irreducible root system.
Recall that $H\mapsto H\!\!\Downarrow\, \equal H(1)$, where the
action of $H\in \HH$ in $\v$ is used.

(i) Given two strictly
positive sequences $\vec\rr,\vec\ss$ of length $\ell$ as 
in (\ref{iterrss}), we lift $(\rr_i,\ss_i)^{tr}$ to
$\ga_i$ and then to $\hat{\ga_i}$ (acting in $\HH$) as above.
For a weight $b\in P_+\,$, the {\sf DAHA-Jones polynomial} is
\begin{align}\label{jones-dit}
&J\!D_{\,\vec\rr,\,\vec\ss\,}^{R}(b\,;\,q,t)\ =\ 
J\!D_{\,\vec\rr,\,\vec\ss\,}(b\,;\,q,t)\, \equal\\ 
\Bigl\{\hat{\ga_{1}}&\Bigr(
\cdots\Bigl(\hat{\ga}_{\ell-1}
\Bigl(\bigl(\hat{\ga}_\ell(P_b)/P_b(q^{-\rho_k})\bigr)\!\Downarrow
\Bigr)\!\Downarrow\Bigr) \cdots\Bigr)\Bigr\}_{ev}.\notag
\end{align}
It does not depend
on the particular choice of the lifts \,$\ga_i$
and $\hat{\ga}_i\in GL_{\,2}^{\wedge}(\Z)$.
The {\sf tilde-normalization} 
$\tilde{J\!D}_{\vec\rr,\vec\ss}\,(b\,;\,q,t)$ is well defined
and is a polynomial in terms of $\,q,t_\nu$ with the
constant term $1$. 

(ii) Let us switch to the root system $A_n$ for $sl_{n+1}$,
setting $t\!=\!t_{\sht}\!=\!q^k$ and  considering $P_+\ni b=$ 
$\sum_{i=1}^n b_i \om_i$ as 
(dominant) weights for any $A_m$ (for $sl_{m+1}$)
with $m\ge n-1$, where we assume that $\om_{n}=0$ upon the 
restriction to $A_{n-1}$.

Then given $\t(\vec\rr,\vec\ss)$ as above, there exists a 
{\sf DAHA\~superpolynomial\,} 
$\h_{\,\vec\rr,\,\vec\ss}\,(b\,;\,q,t,a)$
 from $\Z[q,t^{\pm 1},a]$ satisfying the relations
\begin{align}\label{jones-sup}
\h_{\,\vec\rr,\,\vec\ss}\,(b;q,t,a\!=\!-t^{m+1})\!=\!
\tilde{J\!D}_{\,\vec\rr,\,\vec\ss\,}^{A_{m}}(b;q,t)
\hbox{\, for any } m\!\ge\! n\!-\!1;
\end{align} 
then its $a$\~constant term %(the coefficient of $a^0$) 
is automatically tilde-normalized.
\sq
\end{theorem}
```

### 1.6 Framing / $q^\bullet t^\bullet$ factors

* The invariants are defined "up to $q^\bullet t^\bullet$" and then tilde/hat-normalized; framing is the standard (zero) one and
  **no framing factor is kept** (CD15:3899-3906, quoted in 2.1 below: "Since the DAHA-invariants are considered in this work up to powers of q,t, this will be sufficient").
* The only place where an explicit framing-type factor appears is the non-reduced $A_1$ comparison with the
  Jones polynomial (CD14 Prop. 3.2), where the factor $q^{-\frac14 a_\ell r_\ell b(b+2)}$ appears and also $q\mapsto q^{-1}$:

`1408.4348/src/cable.tex:3710-3785`
```latex
We will prove the Connection Conjecture 
for the {\em non-reduced\,} version of the
DAHA-Jones polynomials. Namely, we modify
(\ref{jones-dit}) as follows:

\begin{align}\label{nonrjone}
&J\!D^\#_{\,\vec\rr,\,\vec\ss\,}(b\,;\,q,t) \equal 
\Bigl\{\hat{\ga_{1}}\Bigr(
\cdots\Bigl(\hat{\ga}_{\ell-1}
\Bigl(\bigl(\hat{\ga}_\ell(P_b^{(k)})\bigr)\!\Downarrow
\Bigr)\!\Downarrow\Bigr) \cdots\Bigr)\Bigr\}^{(k)}, \\
&\for t \equal q^k\hbox{\,\, and\,\, } 
\ga=\ga_{r,s} \hbox{\,\, lifted to\,\, }
\hat{\ga}\in PSL_{\,2}^{\wedge}(\Z),
\notag
\end{align}
where we show explicitly the dependence of the $P$\~polynomial
on the parameter $k$ and also use ${}^{(k)}$ for the
corresponding coinvariant. 
Generally, the non-reduced DAHA-Jones polynomials are 
inconvenient because of their non-trivial $q,t$\~denominators. 
However they become $q$\~polynomials as $t=q$ and simplify 
the considerations. 

\smallskip

\subsection{\bf Jones polynomials}
We are going to use the
relation for the colored Jones polynomials of the cable 
$C\!ab(\aa,\rr)K$ with that 
of the initial knot $K$. See \cite{Mo,Tr} and references therein.
Note that $t$ in \cite{Tr} is our $q^{1/4}$,  
%we use $(\aa,\rr)$ parameters for the cabling to distinguish them 
%from the $(\rr,\ss)$ parameters used in the construction of the 
%of the DAHA-Jones polynomials.  
and $C\!ab(r,s)$ there is our $C\!ab(\aa,\rr)$  (in this order!). 

Assuming here and below that the $j$-summation has step $1$, one has:
\begin{align}\label{jones-formula}
&\j^{\#}_{C\!ab(\aa,\rr)K}(b\,;\,q) = q^{-\frac{1}{4}\aa\rr\, b(b+2)} 
\sum_{j=-b/2}^{b/2} q^{\aa\, j(j\,\rr +1)} \j^{\#}_{K}(2\rr j\,;\,q).
\end{align}

This formula gives a recursive definition of the Jones 
polynomials of iterated knots; the unknot is normalized as
follows:

\begin{align}
\j^{\#}_{\unknot}(b\,;\,q) = (q^{(b+1)/2}-q^{-(b+1)/2})/
(q^{1/2}-q^{-1/2}).
\end{align}
\smallskip

\begin{proposition}\label{RELATIONTOJONES}
For any (possibly negative)
$\{\rr,\ss\}$ from Theorem \ref{JONITER}, 
\begin{align}
&J\!D^\#_{\{\rr_1,\ldots \rr_\ell\},\,\{\ss_1,\ldots \ss_\ell\}\,}
(b\,;\,q,t\mapsto q)\\
= q^{-\frac{1}{4}\aa_\ell \rr_\ell\, b(b+2)}
\j&^{\#}_{T(\{\rr_1,\ldots 
\rr_\ell\},\,\{\ss_1,\ldots \ss_\ell\})\,}(b\,;\, q\mapsto q^{-1}),
\where 
\notag
%&= q^{-\frac{1}{4}\aa_\ell \rr_\ell b(b+2)} 
%\j^{\#}_{\bigl(C\!ab(\aa_\ell,\rr_\ell)
%\cdots C\!ab(\aa_1,\rr_1)\bigr)
%(\unknot )}(b\,;\, q\mapsto q^{-1}), \notag
\end{align}
$T(\{\rr_1,\ldots \rr_\ell\},\,\{\ss_1,\ldots \ss_\ell\}) = 
C\!ab(\vec\aa,\vec\rr)=\bigl(C\!ab(\aa_\ell,\rr_\ell)\cdots 
C\!ab(\aa_1,\rr_1)\bigr)(\unknot )$,
and $\aa_i$ are recursively defined by relations $\aa_{i+1} = 
\rr_i \rr_{i+1} \aa_i + \ss_{i+1}$ for $1 \leq i \leq \ell - 1$ 
and the initial condition $\aa_1 = \ss_1$. 
\end{proposition}
```

---------------------------------------------------------------------------------------------------

## 2. Newton pairs ↔ topological cables

### 2.1 Definitions (CD14 §2.1; CD15 §3.1; C24 §2.1)

`1408.4348/src/cable.tex:746-802`
```latex
The {\em iterated torus knots\,}, denoted in this paper by 
$\t(\vec{\rr},\vec{\ss})$,
will be associated with two sequences of
integers {\em of any signs\,}: 
\begin{align}\label{iterrss} 
\vec{\rr}=\{\rr_1,\ldots \rr_\ell\}, \ 
\vec{\ss}=\{\ss_1,\ldots \ss_\ell\} \hbox{\, such that\, 
gcd}(\rr_i,\ss_i)=1;
\end{align}
$\ell$ will be called the {\em length\,} of $\,\vec\rr,\vec\ss$.
The pairs $\{\rr_i,\ss_i\}$ will become 
{\em characteristic\,} or {\em Newton pairs\,}
for algebraic knots. We will use this name, but
generally consider arbitrary iterated knots in this paper. 

These are combinatorial data. The (topological) definition
of iterated torus knots requires one more sequence: 
\begin{align}\label{Newtonpair}
\aa_1=\ss_1,\,\aa_{i}=\aa_{i-1}\rr_{i-1}\rr_{i}+\ss_{i}\,\ 
(i=2,\ldots,m). 
\end{align}
See  e.g. \cite{EN} and \cite{Pi}. Then, in terms of
the cabling defined below:
\begin{align}\label{Knotsiter}
\t(\vec{\rr},\!\vec{\ss})\equal
C\!ab(\vec{\aa},\!\vec{\rr})(\unknot)=
\bigl(C\!ab(\aa_\ell,\rr_\ell)\cdots C\!ab(\aa_2,\rr_2)\bigr)
\bigl(T(\aa_1,\!\rr_1)\bigr).
\end{align}
Note that the first iteration (application of 
$C\!ab$) is for $\{\aa_1,\rr_1\}$ (not for the last pair!).
Recall that $\,\unknot$ is the
{\em unknot\,} and knots are considered up to 
{\em ambient isotopy\,} (we use \,``$=$" for it).
\smallskip

{\sf Cabling}.
The {\em cabling\,} $C\!ab(\aa,\bb)(K)$ of any oriented
knot $K$ in (oriented) $S^3$ is defined as follows; 
see e.g. \cite{Mo,EN} and references therein. 
We consider a small $2$\~dimensional torus
around $K$ and put there the torus knot $T(\aa,\bb)$
in the direction of $K$, 
which is $C\!ab(\aa,\bb)(K)$ (up to ambient isotopy).
We will sometimes set $C\!ab(\vec\aa,\vec\rr)= 
C\!ab(\vec\aa,\vec\rr)(\unknot)$.  

This procedure depends on 
the order of $\aa,\bb$ and orientation of
$K$. We choose the latter in the standard way (compatible 
with almost all sources, including  the Mathematica package
``KnotTheory"); the parameter $\aa$ gives the number of
turns around $K$. This construction also depends
on the framing of the cable knots; we take the natural
one, associated with the parallel copy of the torus where 
a given cable knot sits (its parallel copy 
has zero linking number with this knot). 
```

`1408.4348/src/cable.tex:817-850`
```latex
By construction,
$C\!ab(\aa,0)(K)=\unknot$ and  $C\!ab(\aa,1)(K)=K$ for any
knot $K$ and $\aa\neq 0$.
Accordingly, we have the following {\em reduction cases\,}:
\begin{align}\label{iter-triv}
&\hbox{When\, }
\rr_i=0\,\ (i<\ell),\ \, \t(\vec{\ss},\vec{\rr})=\\
\t(\{\rr_{i+1},\cdots,\rr_{\ell}&\},
\{\ss_{i+1},\cdots,\ss_{\ell}\}),\  
\t(\vec{\rr},\vec{\ss})=\unknot
\for i=\ell.\notag\\
&\hbox{When\, } \rr_i=1,\ \ss_i=0,\ 
\t(\vec\rr,\vec\ss)=\label{iter-s0}\\
\t(\{\rr_1,\cdots,\rr_{i-1}&,
\rr_{i+1},\cdots,\rr_\ell\},\{\ss_1,\cdots,\ss_{i-1},
\ss_{i+1},\cdots,\ss_\ell\}).\notag\\
&\hbox{When\, } \rr_i=1,\, \ss_i\neq 0,\ \,
\t(\vec\rr,\vec\ss)=\label{iter-s1}\\
\t(\{\rr_1,\cdots,\rr_{i-1}&,
\rr_{i+1},\cdots,\rr_\ell\},\{\cdots,\ss_{i-1},
\,\ss'_{i+1}\,,\ss_{i+2},\cdots\}),\notag\\
\hbox{where\,\,\, } \ss_{i+1}'\!=\!\ss&_{i+1}\!+\!\ss_{i}\rr_{i+1} 
\hbox{\,\, if\,\, } i <\ell\hbox{\,\, (no\, $s'_{\ell+1}$\, 
for\, } i=\ell). 
\end{align}
\smallskip

Also $T(\rr,\ss)=T(\ss,\rr)$, which results in
the {\em transposition property\,}:
\begin{align}\label{iter-sym}
&\t(\vec{\rr},\vec{\ss})=
\t(\{\ss_1,\rr_{2},\cdots,\rr_{\ell}\},
\{\rr_{1},\ss_2\cdots,\ss_{\ell}\}).
\end{align}
```

`1408.4348/src/cable.tex:898-921`
```latex
Switching from $K$ to its {\em mirror image\,}, 
denoted by $K^\star$, one has:
\begin{align}\label{iterorien}
C\!ab(\aa,\bb)(K^\star)&=
\bigl(C\!ab(-\aa,\bb)(K)\bigr)^\star \hbox{\,\, for any }
\aa,\bb \hbox{ with\,\,   gcd}(\aa,\bb)\!=\!1,\notag\\
\hbox{and\, }&C\!ab(-\vec\aa,\vec\rr)\!=\!
\bigl(C\!ab(\vec\aa,\vec\rr)\bigr)^\star,\ 
\t(\vec\rr, -\vec\ss)\!=\!
\bigl(\t(\vec\rr,\vec\ss)\bigr)^\star.
\end{align}

Changing the {\em orientation\,}, denoted by ``$-$", at the 
$i$\hbox{\tiny th}\, step, we obtain 
that for any $1\le i\le \ell$,
\begin{align}\label{iter2min}
-C\!ab(\vec\aa,\vec\rr)\!=\!
C\!ab&(\{\ldots\!,\aa_{i\!-\!1},\!-\aa_i,\aa_{i\!+\!1},\!\ldots\},
\{\ldots\!,\rr_{i\!-\!1},\!-\rr_i,\rr_{i\!+\!1},\!\ldots\}),
\notag\\
-\t(\vec\rr,\vec\ss)\!=\!\t&(\{\ldots,\rr_{i\!-\!1},
\!-\rr_i,\rr_{i\!+\!1},\ldots\},
\{\ldots\!,\ss_{i\!-\!1},\!-\ss_i,\ss_{i\!+\!1},\!\ldots\} ).
\end{align} 
```

CD15 restates it (note: here the first knot is written $T(r_1,s_1)$, with the remark that $r,s$ are transposed w.r.t. $Cab$):

`1509.08351/src/links.tex:3851-3906`
```latex
This presentation, referred to as the {\em $\tax$-presentation\,}, 
will be exactly the one needed 
in the DAHA approach. However it is not optimal for establishing
the symmetries of our polynomials and the justification 
that our construction depends only on the 
corresponding knot/link. We actually need 
the {\em cable presentation\,} for this, which is actually from
the definition of the corresponding iterated torus knots. 
It requires one more sequence of integers (possibly negative):
\begin{align}\label{Newtonpair}
\aa_1=\ss_1,\,\aa_{i}=\aa_{i-1}\rr_{i-1}\rr_{i}+\ss_{i}\,\ 
(i=2,\ldots,m). 
\end{align}
See  e.g., \cite{EN}. In terms of
the {\em cabling\,} discussed below, the corresponding knots
are as follows. First, $T(\rr,\ss)=C\!ab(\ss,\rr)(\unknot)$
(note that we transpose $\rr,\ss$ here), and then we set:
\begin{align}\label{Knotsiter}
\t(\vec{\rr},\!\vec{\ss})\rightsquigarrow
C\!ab(\vec{\aa},\!\vec{\rr})(\unknot)=
\bigl(C\!ab(\aa_\ell,\rr_\ell)\cdots C\!ab(\aa_2,\rr_2)\bigr)
\bigl(T(\rr_1,\!\ss_1)\bigr).
\end{align}

The first iteration here (application of 
$C\!ab$) is  $C\!ab(\aa_1,\rr_1)$, {\em not that for the last 
pair\,}, and then we proceed inverting the natural ordering.
\smallskip

{\sf Cabling}.
Knots and links will be always considered up to 
{\em ambient isotopy\,}; we use \,``$=$" for it.
The {\em cabling\,} $C\!ab(\aa,\bb)(K)$ of any oriented
knot $K$ in (oriented) $\S^3$ is defined as follows; 
see e.g., \cite{Mo,EN} and references therein. 
We consider a small $2$\~dimensional torus
around $K$ and put there the torus knot $T(\aa,\bb)$
in the direction of $K$, 
which is $C\!ab(\aa,\bb)(K)$ (up to ambient isotopy);
we set $C\!ab(\vec\aa,\vec\rr)\equal
C\!ab(\vec\aa,\vec\rr)(\unknot)$.   

This procedure depends on 
the order of $\aa,\bb$ and orientation of
$K$. We choose the latter in the standard way, matching  
the Mathematica package
Knot Atlas \cite{KA}; the parameter $\aa$ gives the number of
turns around $K$. This construction also depends
on the {\em framing\,} of the cable knots; we take the natural
one, associated with the parallel copy of the 
torus where a given cable knot sits (its parallel copy 
has zero linking number with this knot). It will be the 
same {\em standard framing\,} for iterated torus {\em links\,};
see below. Since the DAHA-invariants
are considered in this work up to powers of $q,t$, this 
will be sufficient. 
```

`1509.08351/src/links.tex:3933-3943`
```latex
Let us comment on the last relation; see
(\ref{Newtonpair}). Since $r_i=1$, one has:
$$
\aa_i\!=\!\aa_{i-1}\rr_{i-1}\!+\!\ss_i,\ \,
\aa_{i+1}\!=\!\aa_i \rr_{i+1}\!+\!\ss_{i+1}\!=\!
\aa_{i-1} \rr_{i-1}\rr_{i+1}\!+ (\ss_{i+1}\!+\ss_i \rr_{i+1}).$$
The pairs $\{\aa_i,\rr_i\}$ are sometimes called 
{\em topological or cable parameters\,}. Indeed,   
the isotopy equivalence of iterated knots generally can be
seen only at the level of these parameters (not at the level 
of Newton pairs).
```

`1509.08351/src/links.tex:3980-3999`
```latex
Furthermore,
changing the {\em orientation\,}, denoted by ``$-$", at the 
$i$\hbox{\tiny th}\, step, we obtain 
that for any $1\le i\le \ell$,
\begin{align}\label{iter2min}
-C\!ab(\vec\aa,\vec\rr)\!=\!
C\!ab&(\{\ldots\!,\aa_{i\!-\!1},\!-\aa_i,\aa_{i\!+\!1},\!\ldots\},
\{\ldots\!,\rr_{i\!-\!1},\!-\rr_i,\rr_{i\!+\!1},\!\ldots\}),
\\
(-)^{\ell-i+1}\,\t(\vec\rr,\vec\ss)\!=
\!\t&(\{\ldots,\rr_{i\!-\!1},
\!-\rr_i,\rr_{i\!+\!1},\ldots\},
\{\ldots\!,\ss_{i\!-\!1},\!-\ss_i,\ss_{i\!+\!1},\!\ldots\} ).
\notag
\end{align} 
The second transformation here results in the following cable: 
$$C\!ab(\{\ldots\!,\aa_{i\!-\!1},\!-\aa_i,-\aa_{i\!
+\!1},\!\ldots\},
\{\ldots\!,\rr_{i\!-\!1},\!-\rr_i,-\rr_{i\!+\!1},\!\ldots\}),$$
which explains the sign.
```

C24 (clearest statement of which number is the winding number):

`2410.14703/src/superpolynoms.tex:1404-1445`
```latex
{\sf Cabling parameters.} 
The {\em cabling\,} $C\!ab(\aa,\bb)(K)$ of 
any oriented
knot $K$ in (oriented) $S^3$ is defined as follows;
see e.g., \cite{EN} and references therein.
We consider a small $2$\~dimensional torus
around $K$ and put there the torus knot $T(\bb,\aa)$
in the direction of $K$,
which is $C\!ab(\aa,\bb)(K)$. The parameter
$\aa$ is the number of
turns around $K$, $\bb$ is that along $K$.
For instance, $C\!ab(a,0)K=\unknot$ and $Cab(a,1)K=K$.

This procedure depends on
the order of numbers $\aa$ and $\bb$, the orientation of $K$,
and its {\em framing}. We take the natural
zero-framing  when  the parallel copy of the knot has
zero linking number with a given knot. 


Using this operation, iterated torus knots (also called
satellite knots)
are as follows.
Given two sequences
$\vec\aa=(\aa_1,\ldots,\aa_\ell),\vec\rr=(\rr_1,\ldots,\rr_{\ell})\}$, 
called a {\em cable presentation}, 
we begin with  $T(\rr_1,\aa_1)=C\!ab(\aa_1,\rr_1)(\unknot)$
(note that we transpose $\aa$ and $\rr$ here), and then set:
\begin{align}\label{Knotsiter} 
\t(\vec\rr,\vec\ss)=C\!ab(\vec{\aa},\!\vec{\rr})(\unknot)=
\Bigl(C\!ab(\aa_\ell,\rr_\ell)\cdots \bigl(C\!ab(\aa_2,\rr_2)
(T(\rr_1,\!\aa_1))\bigr)\Bigr),
\end{align}
where the connection with the $\tax$-presentation is as follows:
\begin{align}\label{Newtonpair}
\aa_1=\ss_1,\,\aa_{i}=\aa_{i-1}\rr_{i-1}\rr_{i}+\ss_{i}\,\
(i=2,\ldots,\ell).
\end{align}
See  e.g., \cite{EN}. The $\tax$-presentation will be the
one used later. However, to prove the topological invariance
of the DAHA construction (a theorem) the passage to the cable presentation
is necessary. 
```

C24 explicit cable/Newton-pair examples and the semigroup:

`2410.14703/src/superpolynoms.tex:2694-2750`
```latex
The simplest ``non-torus" family 
is $\r\!=\!\C[[z^4,z^6\!+\!z^{7+2m}]]$ for $m\!\in\! \Z_+$,
which are 
of multiplicity $4$ and with $\de_m\!=\!8\!+\!m$. The corresponding
cables are $C\!ab(13+2m,2)C\!ab(2,3)$, where 
$C\!ab(2,3)=T(3,2)$. One can replace here 
$T(3,2)$ by $T(2,3)$ because it does not change
the isotopy type of the cable; the cable parameters
are topological invariants. 
\vskip 0.2cm

The simplest ``triple" cable is 
$C\!ab(53,2)C\!ab(13,2)C\!ab(3,2)$.
 The corresponding
ring is $\r\!=\!$
$
\!\C[[x\!\!=\!z^8,y\!=\!z^{12}\!+\!z^{14}\!+\!z^{15}]].
$ 

let us provide the corresponding {\em Newton's pairs.}
Generally, they are
$\{\rr_1,\ss_1\},\{\rr_2,\ss_2\}, \ldots\,$ such that
the Puiseux-type equation of the singularity is
{  $y=x^{\frac{\ss_1}{\rr_1}}\Bigl(1+ c_1 
x^{\frac{\ss_2}{\rr_1 \rr_2}}
\bigl(1+c_2 x^{\frac{\ss_3}{\rr_1 \rr_2 \rr_3}}(\cdots)\bigr)\Bigr)$}
for generic $c_i$. 
The Newton's pairs are $\{(2,3),(2,1),(2,1)\}$ in this example.

The valuation semigroup is $\Gamma=\lan 0,8,12,26,53\ran$, the
linear span of these numbers over  $\Z_+=\{\Z\ni n\ge 0\}$.
The Serre number is $\de\!=\!42$  and $mult=8$.
Generally, $\Ga=\lan 0,\, \rr_1 \rr_2 \rr_3,\,
\aa_1 \rr_2 \rr_3,\, \aa_2 \rr_3,\, \aa_3\ran$ for the cable parameters 
$(\aa_i,\rr_i)$ above (here $\rr_1<\ss_1$ is used). Recall that
$\aa_1\!=\!\ss_1$, $a_2\!=\!\rr_1 \ss_1 \rr_2+\ss_2,$\, 
$\aa_3\!=\!\aa_2 \rr_2 \rr_3+\ss_3\,$, and so on for any number 
of $\{\rr_i,\ss_i\}$. 
\vskip 0.2cm

As another example, let us calculate $\Ga$ and $\de$ for 
the ring 
$\r=\F[[z^{\upsilon \rr}, z^{\upsilon \ss} +z^{\upsilon \ss+p}]]$,
where, gcd$(\rr,\ss)=1$, $\rr<\ss$,
$\upsilon>1$ and gcd$(\upsilon,p)=1$ for $p\ge 1$. The
corresponding cable is $C\!ab(\upsilon \rr \ss\!+\!p,\,\upsilon)\,
C\!ab(\ss,\rr)$. One has:
\begin{align}\label{Ga-cab}
&\Ga=\lan \upsilon \rr,\,
\upsilon \ss,\,\upsilon \rr \ss\!+\!p\ran,\ 
2\de=\upsilon^2 \rr \ss\!-\!\upsilon (\rr\!+\!\ss)+
(\upsilon\!-\!1)p+1.
\end{align}

We will mostly assume below that the first $C\!ab(\rr_1,\ss_1)$ in any 
cables satisfies the inequality $\rr_1\le \ss_1$, as it was in these
examples. 
```

### 2.2 Algebraic knots (CD14 §2.2)

`1408.4348/src/cable.tex:933-979`
```latex
\subsection{\bf Algebraic knots}\label{sec:ALG-KNOTS}
There exists a deep connection of the iterated torus knots
for strictly positive $\rr_i,\ss_i$
with the germs of isolated {\em plane curve singularities\,}. 
Its origin is the Newton's successive approximations 
for $y$ in terms
of $x$ satisfying a polynomial equation $f(x,y)=0$
in a neighborhood of $0=(x=0,y=0)$. The main claim is that
the knot $\t(\vec\rr,\vec\ss)$ for the 
{\em characteristic pairs\,} $\{\rr_i,\ss_i\}$
with $\rr_i,\ss_i>0\,$ is the link 
of the {\em germ of the singularity\,} 
\begin{align}\label{yxcurve}
y = x^{\frac{\ss_1}{\rr_1}}
(c_1 + x^{\frac{\ss_2}{\rr_1\rr_2}}
\bigl(c_2 + \ldots + 
x^{\frac{\ss_\ell}{\rr_1\rr_2\cdots\rr_\ell}}\bigr))
\hbox{\, at\, } 0,
\end{align} 
which is the intersection of the corresponding plane curve 
$f(x,y)=0$ with a small $3$-dimensional sphere in 
$\C^2$ around $0$. We will always
assume that this germ is {\em unibranch}.

The inequality  $\ss_1<\rr_1$ is commonly
imposed here (otherwise $x$ and $y$ can be switched).
Formula (\ref{yxcurve}) is the celebrated 
{\em Newton-Puiseux expansion}, though Puiseux 
performed the multiplication here:
\begin{align*}
y\!=\!b_1\, x^{\frac{\mm_1}{\rr_1}}\!+\!
b_2\, x^{\,\frac{\mm_2}{\rr_1\rr_2}}\!+\!
b_3\, x^{\,\frac{\mm_3}{\rr_1\rr_2\rr_3}}\!+\ldots,\ \,
\mm_1\!=\!\ss_1,\, \mm_i\!=\!\ss_i\!+\!\rr_i\mm_{i-1},\ldots\,,
\end{align*} 
where $\{\rr_i,\mm_i\}$ are called the  
{\em Puiseux pairs\,} (which we will not use).

We will denote the germ of
the curve in (\ref{yxcurve}) by $\,\c_{\vec\rr,\vec\ss}$. 
If it has more than one branch (near $0$), then
its complete topological description requires 
knowing the sequences $\{\vec{\rr},\vec{\ss}\}$ for all components 
together with the pairwise linking numbers. See e.g. \cite{Mi}
and pg. 49 of \cite{EN}. All algebraic knots can be obtained
in such a way; the positivity of $\,\rr_i,\ss_i\,$ is 
necessary for $\,\t(\vec\rr,\vec\ss)\,$ to be an algebraic link.
```

### 2.3 Dictionary to "the $(p_2,q_2)$-cable of $T(p_1,q_1)$, $q_2$ in the Seifert framing of the companion"

From $Cab(a,r)K$ := torus knot $T(r,a)$ on the boundary torus of $K$ "in the direction of $K$", with $r$ = number of turns *along* $K$
(winding number), $a$ = number of turns *around* $K$, zero framing (C24:1404-1421, CD14:782-802):

* $Cab(a_2,r_2)\,T(r_1,s_1)$ = standard $(p_2,q_2)=(r_2,a_2)$-cable of $T(p_1,q_1)=T(r_1,s_1)$, with
  $$\boxed{\,r_2=p_2,\qquad s_2=a_2-a_1r_1r_2=q_2-p_1q_1p_2\,}$$
  ($a_1=s_1$; for $T(p_1,q_1)$ either order may be used for $(r_1,s_1)$ because of (iter-sym); $p_1q_1p_2$ is the
  torus-framing vs. Seifert-framing correction: the companion's torus-framing has linking $p_1q_1$).
* General: $s_i=a_i-a_{i-1}r_{i-1}r_i$.
* **Algebraicity:** all $r_i,s_i>0$ (CD14:941-979; "positivity of $r_i,s_i$ is necessary for $\mathbb T(\vec r,\vec s)$ to be an algebraic link").
  For a 2-step cable of $T(p_1,q_1)$ ($p_i,q_i>0$): algebraic iff $q_2>p_1q_1p_2$. Convention $s_1<r_1$ is "commonly imposed"
  (CD14:957); the a-degree statements assume $|r_1|>|s_1|$, i.e. $\deg_a=(s_1r_2\cdots r_\ell-1)|\lambda|$, e.g. $=3$ for trefoil 2-cables with $(r_1,s_1)=(3,2)$.
* Sign flips: $JD$ is invariant under $(r_i,s_i)\mapsto(-r_i,-s_i)$ for a single $i$ (CD14 Thm MAINTHM, second line of (jones-sym-iter), the DAHA counterpart of (iter2min) = orientation change);
  this changes $a_i\mapsto -a_i$ and $Cab(-a,-r)=\pm Cab(a,r)$. CD14 uses this to write cables with $a_2<12$ of the trefoil with $r_2=-2$:
  $\{3,-2\},\{2,s\}$ gives $a_2=s-12$, i.e. $Cab(12-s,2)T(3,2)$ (up to orientation).
* Mirror: $\mathbb T(\vec r,-\vec s)=\mathbb T(\vec r,\vec s)^\star$, $\mathcal H\mapsto \mathcal H(q^{-1},t^{-1},a^{-1})$ up to monomial (CD14:898-908, C24:2538-2543).

**Worked dictionary (all [checked] by $\mathcal H(q,q,a=-1)=\Delta_{T(p_1,q_1)}(x^{p_2})\Delta_{T(p_2,q_2)}(x)$, see §6):**

| topological knot | DAHA $[\vec r,\vec s]$ (as in the sources) | equivalent | algebraic? | source of $\mathcal H$ |
|---|---|---|---|---|
| (2,13)-cable of $T(2,3)$ = $Cab(13,2)T(3,2)$ | $\{3,2\},\{2,1\}$ | $\{2,2\},\{3,1\}$ (transposition) | yes (simplest non-torus alg.) | CD14 (1), CP16, C24 |
| (2,15)-cable | $\{3,2\},\{2,3\}$ | | yes | CD14 (2), CP16 |
| (2,11)-cable | $\{3,2\},\{2,-1\}$ | | no (pseudo-alg.) | CD14 (7) |
| (2,9)-cable | $\{3,2\},\{2,-3\}$ | | no (pseudo-alg.) | CD14 (6), C24 |
| (2,7)-cable | $\{3,-2\},\{2,5\}$ | $\{3,2\},\{2,-5\}$ | no (pseudo-alg., deg$_a$=2) | CD14 (a) |
| (2,5)-cable | $\{3,-2\},\{2,7\}$ | $\{3,2\},\{2,-7\}$ | no, non-positive | CD14 (b), C24 |
| **(2,3)-cable** | $\{3,-2\},\{2,9\}$ | $\{3,2\},\{2,-9\}$ | no, non-positive | CD14 (c) |
| **(2,1)-cable** | $\{3,-2\},\{2,11\}$ | $\{3,2\},\{2,-11\}$ | no, non-positive | CD14 (d) |
| (2,-1)-cable | $\{3,-2\},\{2,13\}$ | | no, non-positive | CD14 (e) |
| (2,-3)-cable | $\{3,-2\},\{2,15\}$ | | no, non-positive | CD14 (f) |
| (2,-5)-cable | $\{3,-2\},\{2,17\}$ | | no, positive | CD14 (g) |
| (3,19)-cable of $T(2,3)$ | $\{3,3\},\{2,1\}$ | | yes | CD14 (3) |
| (3,20)-cable | $\{3,3\},\{2,2\}$ | | yes | CD14 Betti only (header misprinted as $\{4,2\},\{3,1\}$, CD14:2436 vs 2455) |
| (3,22)-cable | $\{3,3\},\{2,4\}$ | | yes | CP16 |
| (4,25)-cable | $\{3,4\},\{2,1\}$ | | yes | CD14 Betti only |
| (2,53)-cable of (2,13)-cable | $\{3,2,2\},\{2,1,1\}$ | | yes | CD14 Betti only |
| (2,25)-cable of $T(3,4)$ | $\{4,2\},\{3,1\}$ | | yes | CD14 (4) |
| (2,17)-cable of $T(3,4)$ | $\{4,-2\},\{3,7\}$ | $\{4,2\},\{3,-7\}$ | no, positive | CD14 (h) |
| (2,11)-cable of $T(3,4)$ | ($\{4,2\},\{3,-13\}$) | | no, non-positive | C24 (no Newton pairs given) |
| (2,21)-cable of $T(2,5)$ | $\{5,2\},\{2,1\}$ | | yes | CD14 Betti only |
| (3,32)-cable of $T(2,5)$ | $\{5,2\},\{3,2\}$ | | yes | CD14 Betti only |
| (2,31)-cable of $T(2,7)$ | $\{7,2\},\{2,3\}$ | | yes | CP16 |

CD15 example of the cable notation with a negative $s_1$ (here $a_1=-2$, $a_2=-11$; the first knot is the mirror trefoil although written $T(3,2)$):

`1509.08351/src/links.tex:9162-9180`
```latex
For example,
for the tree $\Up=\{\twoone\}$ with 2 vertices and
2 arrowheads, the link
$\l^{\twoonetiny}_{\{3,2\},\{-2,1\}}$ means 
combinatorially and topologically:
$$
\l^{\twoonetiny,\, (\la^1,\la^2)}_{(\bigl\{\{3,2\},\{-2,1\}\bigr\},
\bigl\{\{3,2\}\bigr\})}\, \sim\,
\bigl(C\!ab(-11,2)^{\la^1}\coprod\, C\!ab(0,1)^{\la^2}\bigr)T(3,2),
$$
where we count the paths from the bottom to the top. Recall that 
the labels $[\rr,\ss]$ corresponding here
to $\vec\rr=\{3,2\},\vec\ss=\{-2,1\}$ are $[3,-2]$ and $[2,1]$,
so they can be naturally considered as the columns of 
{\tiny $\begin{pmatrix}
3&\!\!2\\ 
-2&\!\!1
\end{pmatrix}$}; these columns become 
the first columns of the corresponding $\ga$. 
```

---------------------------------------------------------------------------------------------------

## 3. Variables: DAHA $(q,t,a)$ vs. standard / DGR

`1408.4348/src/cable.tex:1488-1498`
```latex
\begin{conjecture}{\sf Connection Conjecture.}\label{CONCONJ}
Let 
$\,\h_{\,\vec\rr,\,\vec\ss}\,(b;q,t,a)_{st}\,$
denote  $\,\h_{\,\vec\rr,\,\vec\ss}\,(b;q,t,a)\,$ expressed
in terms of the standard parameters (in the Khovanov-Rozansky 
theory; see \cite{CJ} and Section 1 in \cite{ORS}):
\begin{align}\label{qtareli}
&t=q^2_{st},\  q=(q_{st}t_{st})^2,\  a=a_{st}^2 t_{st},\notag\\
&q^2_{st}=t,\  t_{st}=\sqrt{q/t},\  a_{st}^2=a\sqrt{t/q}.
\end{align}

```

(identical in CD15:6203-6206 and CP16:1169-1172.)

Khovanov–Rozansky $sl_N$ specialization and the differentials:

`1509.08351/src/links.tex:6066-6082`
```latex
\begin{align}\label{KhR-subs}
a_{st}\mapsto q_{st}^{N}, \hbox{\, equivalently, \, }
a\mapsto t^{N}\sqrt{q/t} \hbox{\ \,in our parameters}.
\end{align} 
This actually depends on the grading used in the theory
and there are variations here in different works.

These relations are applied as such only
for sufficiently large (stable) $N$,
otherwise the theory of differentials is necessary.
The differentials correspond to a different
substitution: $a=-t^N$ in DAHA parameters 
(see the conjecture below). 
Let us mention 
the relation to the 
{\em Heegaard-Floer homology\,}: $N=0$. 
Also, the 
```

`1605.00978/src/jacobian-fin-full.tex:1160-1172`
```latex
The passage to the Khovanov-Rozansky polynomials
for $sl_N$ for sufficiently large $N$   
is the substitution
$a\mapsto\! t^{N}\sqrt{q/t}$. Note the relation 
to the {\em Heegaard-Floer homology\,} for $N=0$.
Equivalently, this passage is
$a_{st}\mapsto q_{st}^{N}$ in the {\em standard topological
parameters} (also used in the {\em ORS Conjecture\,}),
which are related to the DAHA parameters as follows:
\begin{align}\label{qtareli}
&t=q^2_{st},\  q=(q_{st}t_{st})^2,\  a=a_{st}^2 t_{st},\notag\\
&q^2_{st}=t,\  t_{st}=\sqrt{q/t},\  a_{st}^2=a\sqrt{t/q}.
\end{align}
```

Note: CD14:1582 states "$a\mapsto t^{n+1}\sqrt{t/q}$", which is inconsistent with (qtareli); CD15/CP16 have $a\mapsto t^N\sqrt{q/t}$
(this follows from $a=a_{st}^2t_{st}$, $a_{st}=q_{st}^N$). Treat CD14:1582 as a misprint.

`1408.4348/src/cable.tex:1578-1586`
```latex
{\sf Khovanov-Rozansky polynomials.}
Due to $(ii)$, the conjectural relation of the
(uncolored) DAHA- superpolynomials 
to {\em Khovanov-Rozansky polynomials\,} for $sl_{n\!+1}$
is $a\mapsto\! t^{n\!+1}\sqrt{t/q}$ 
in DAHA parameters in the stable case (for sufficiently large
$n$); then one uses the substitution 
from (\ref{qtareli}). 
Also, the relation to the 
```

**Trefoil in DAHA conventions** (CD15:634, C24:2337, C24:2497): $\mathcal H_{3,2}=1+qt+aq$.

`1509.08351/src/links.tex:632-637`
```latex
topology or physics papers. The superpolynomials are
relatively simple  for the simplest torus knots, say 
$1+qt+aq$ for the uncolored trefoil; they are
generally very non-trivial apart from
the family $T(2\mm+1,2)$ even  in the 
absence of colors. %%%!!!
```

`2410.14703/src/superpolynoms.tex:2492-2503`
```latex
For the uncolored trefoil, i.e. for $T(3,2)$ and when 
$\la=\om_1=\yng(1)\, $:\,
$H\!O\!M\!=\!\aa(t+t^{-1}-\aa),\, 
H\!O\!M\,\,\widehat{}\!=\! 1+t^2-t\aa$; 
recall that  
$\h\!=\!1+qt+qa$. 

The Alexander polynomials $Alex(t)$ are, generally,
the following specialization: 
$H\!O\!M(t,\aa\!=\!1)(1-t)^{1-\kappa}$, which is  for links
with $\kappa$ components. In particular,
$Alex\!=\!t^{-1}-1+t,\, Alex\,\widehat{}\!=\!1-t+t^2$ for trefoil.
```

**Conversion to DGR** [checked]: substituting $t=q_{st}^2,\ q=q_{st}^2t_{st}^2,\ a=a_{st}^2t_{st}$:
$1+qt+aq\ \mapsto\ 1+q_{st}^4t_{st}^2+a_{st}^2q_{st}^2t_{st}^3$, and
$\mathcal P^{DGR}_{3_1}=a^2q^{-2}+a^2q^2t^2+a^4t^3=a^2q^{-2}\,(1+q^4t^2+a^2q^2t^3)$ with $(a,q,t)_{DGR}=(a_{st},q_{st},t_{st})$.
So: DGR $=$ (monomial) $\times\ \mathcal H(q\to q^2t^2,\ t\to q^2,\ a\to a^2t)$.
The monomial is not part of the DAHA output. Recipe (mine): choose $\pm a^\alpha q^\beta t^\gamma$ so that $\mathcal P(a,q,t=-1)$ is the reduced
HOMFLY normalized by $\mathcal P(a=q,t=-1)=1$ and $q\leftrightarrow q^{-1}$ symmetric. For algebraic knots this gives $(a/q)^{2\delta}$,
$\delta$ = genus (e.g. $(a/q)^{16}$ for $Cab(13,2)T(3,2)$) [checked for (1),(2)].

HOMFLY-PT reduction (DAHA $(q,a)$ = standard $(q_{st}^2,a_{st}^2)$):

`1408.4348/src/cable.tex:1500-1515`
```latex
(i) We conjecture that for any $\,\vec{\rr},\vec{\ss}\,$
as above (possibly negative), 
\begin{align}\label{conjhomfly}
&\h_{\,\vec\rr,\,\vec\ss}\,(b\,;
\,q,t\!\mapsto\! q,a\!\mapsto\! -a)=
\tilde{
\hbox{\small H\!O\!M}}_{\,\vec\rr,\,\vec\ss}\,(\la_b\,;\,q,a),
\end{align}
where the latter is the tilde-normalization of the (reduced)
HOMFLYPT polynomial; $b$ is an arbitrary dominant weight 
represented by the Young diagram $\la_b$. 
We conjecture that the $a$-degrees of these polynomials 
equal to $(\ss_1 \rr_2\cdots \rr_\ell-1)\hbox{\,ord\,}(\la_b)$
provided the positivity of $\ss_i,\rr_i$
and $\ss_1<\rr_1$.

```

`1408.4348/src/cable.tex:1566-1576`
```latex
The case of $A_1$ was already checked:
\begin{align}\label{homjon2}
\h(m\om_1;\,q,q,-q^2)=
\tilde{\hbox{\small\em H\!O\!M}}(m\om_1;\,q,q^2)
=\tilde{J\!D}^{A_1}(m\om_1;\,q,q),
\end{align}
where the latter is the reduced tilde-normalized
Jones polynomial for any iterated knot 
$\t(\,\vec\rr,\,\vec\ss)$. We use
the explicit formulas from \cite{RJ,Mo}.
\smallskip
```

`1408.4348/src/cable.tex:1773-1778`
```latex
Up to $q^\bullet$, this polynomial is equal to the one considered 
in Section 7 of \cite{ObS} upon its division by $a^{16}$ 
and the substitution $z=q-q^{-1}$. We need to 
change $a^2,q^2$ there by our $a,q$, i.e. switch
from the {\em standard parameters\,} to our ones.

```

`2410.14703/src/superpolynoms.tex:2456-2469`
```latex
Basically, 
$H\!O\!M(t,\aa;\la)\,\,\widehat{}\!=\!\hat{\h}(q=t,t,\aa=-a; \la)$
for {\em iterated torus knots}, 
where $\,\widehat{}\,$ is the {\em hat-normalization} 
of the reduced $H\!O\!M$: must be $1$ for the unknot.
We note that the parameter $\sqrt{\aa}$ is mostly used in
topology.
% Spherical $P_\la^\circ$ can be used for knots;
%moreover, they can be replaced without any other
%changes by nonsymmetric $E_\al^\circ$ (see above). 
The coincidence is due to the author 
for torus knots, Morton-Samuelson (iterated torus knots), and
Cherednik-Danilenko (iterated torus links); see \cite{CJ, MoS,ChD2}.

```

Alexander polynomial: $\hat{\mathcal H}(q,q,a=-1)$ for knots (CD15:6082-6091). Super-duality: $\mathcal H(\lambda;q,t,a)=q^\bullet t^\bullet\mathcal H(\lambda^{tr};t^{-1},q^{-1},a)$.

---------------------------------------------------------------------------------------------------

## 4. Explicit superpolynomials of iterated torus knots (verbatim)

All are **reduced**, tilde/hat-normalized ($\mathcal H(a=0)$ has constant term 1 in the algebraic cases), DAHA variables $(q,t,a)$.
$\square=\omega_1$ (uncolored). Only CD14 (5) is colored ($2\omega_1$). For $\omega_2$ only $q=1,a=0$ is given (super-dual to $2\omega_1$).

### 4.1 (2,13)-cable of trefoil, $\vec r=\{3,2\},\vec s=\{2,1\}$, uncolored — CD14 (1)

`1408.4348/src/cable.tex:1691-1713`
```latex
\subsection{\bf Simplest algebraic iterations}
$$
(\mathbf 1):\ \vec\rr=\{3,2\},\,
\vec\ss=\{2,1\},\,
\t= C\!ab(13,2)T(3,2);\ 
\h_{\,\vec\rr,\,\vec\ss}\,(\square\,;\,\,q,t,a)=
$$
\renewcommand{\baselinestretch}{0.5} 
{\small
\(
1+q t+q^2 t+q^3 t+q^2 t^2+q^3 t^2+2 q^4 t^2+q^3 t^3+q^4 t^3+
2 q^5 t^3+q^4 t^4+q^5 t^4+2 q^6 t^4+q^5 t^5+q^6 t^5+q^7 t^5+
q^6 t^6+q^7 t^6+q^7 t^7+q^8 t^8+a^3 \bigl(q^6+q^7 t+q^8 t^2\bigr)+
a^2 \bigl(q^3+q^4+q^5+q^4 t+2 q^5 t+2 q^6 t+q^5 t^2+2 q^6 t^2
+2 q^7 t^2+q^6 t^3+2 q^7 t^3+q^8 t^3+
q^7 t^4+q^8 t^4+q^8 t^5\bigr)+
a \bigl(q+q^2+q^3+q^2 t+2 q^3 t+3 q^4 t+q^5 t+q^3 t^2+2 q^4 t^2+
4 q^5 t^2+q^6 t^2+q^4 t^3+2 q^5 t^3+4 q^6 t^3+q^7 t^3+q^5 t^4+
2 q^6 t^4+3 q^7 t^4+q^6 t^5+2 q^7 t^5+q^8 t^5+q^7 t^6+q^8 t^6
+q^8 t^7\bigr).
\)
}
\renewcommand{\baselinestretch}{1.2} 
```

DAHA-Betti polynomial and HOMFLY-PT (CD14):

`1408.4348/src/cable.tex:1735-1741`
```latex
For $\vec\rr=\{3,2\},\, \vec\ss=\{2,1\}\,$, 
the DAHA-Bette polynomial equals
\begin{align}\label{bett1}
\h(q\!=\!1,t,a\!=\!0)=
1\!+\!3 t\!+\!4 t^2\!+\!4 t^3\!+\!4 t^4\!+\!3 t^5\!+\!2 t^6
\!+\!t^7\!+\!t^8,
\end{align}
```

`1408.4348/src/cable.tex:1759-1772`
```latex
{\sf HOMFLYPT Polynomial.}
It coincides with the  corresponding DAHA-generated 
tilde-normalized HOMFLYPT polynomial:
\begin{align*}
&\h(q,t\!\mapsto q,a\!\mapsto\!-a)=
1+q^2+q^3+2 q^4+q^5+3 q^6+q^7+3 q^8+q^9\\
&+3 q^{10}+q^{11}+2 q^{12}+q^{13}+q^{14}+q^{16}\,\, -a^3 
\bigl(q^6+q^8+q^{10}\bigr)\\
&+a^2 \bigl(q^3+q^4+2 q^5+2 q^6+3 q^7+2 q^8
+3 q^9+2 q^{10}+2 q^{11}+q^{12}+q^{13}\bigr)\\
&-a \bigl(q+q^2+2 q^3+2 q^4+4 q^5+3 q^6+5 q^7+3 q^8+5 q^9\\
&\ \ \ \ \ \ \ \ +3 q^{10}+4 q^{11}+2 q^{12}+2 q^{13}
+q^{14}+q^{15}\bigr).
\end{align*}
```

Same polynomial, factored by $q$-powers, in CP16 and C24 (C24 calls it $\mathcal H^{mot}$, conjecturally $=\mathcal H^{daha}$) [checked: all three identical]:

`1605.00978/src/jacobian-fin-full.tex:2855-2884`
```latex
The first DAHA superpolynomial is as follows:
$$
\vec\rr=\{3,2\},\,
\vec\ss=\{2,1\},\,
\t= C\!ab(13,2)T(3,2);\ 
\h_{\,\vec\rr,\,\vec\ss}\,(\square\,;\,\,q,t,a)=
$$
\renewcommand{\baselinestretch}{0.5} 
{\small
\(
1 + q t + q^8 t^8 + q^2 \bigl(t + t^2\bigr) 
+ a^3 \bigl(q^6 + q^7 t + q^8 t^2\bigr) + 
 q^3 \bigl(t + t^2 + t^3\bigr) + q^4 \bigl(2 t^2 + t^3 + t^4\bigr) 
+  q^5 \bigl(2 t^3 + t^4 + t^5\bigr) + q^6 \bigl(2 t^4 + t^5 
+ t^6\bigr) + 
 q^7 \bigl(t^5 + t^6 + t^7\bigr) + 
 a^2 \bigl(q^3 + q^4 \bigl(1 + t\bigr) + q^5 \bigl(1 + 2 t 
+ t^2\bigr) + 
    q^6 \bigl(2 t + 2 t^2 + t^3\bigr) + q^7 \bigl(2 t^2 + 2 t^3 
+ t^4\bigr) + 
    q^8 \bigl(t^3 + t^4 + t^5\bigr)\bigr) + 
 a \bigl(q + q^2 \bigl(1 + t\bigr) + q^3 \bigl(1 + 2 t + t^2\bigr) 
+ q^4 \bigl(3 t + 2 t^2 + t^3\bigr) + 
    q^5 \bigl(t + 4 t^2 + 2 t^3 + t^4\bigr) + q^6 \bigl(t^2 
+ 4 t^3 + 2 t^4 + t^5\bigr) + 
    q^7 \bigl(t^3 + 3 t^4 + 2 t^5 + t^6\bigr) + q^8 \bigl(t^5 
+ t^6 + t^7\bigr)\bigr).
\)
}
\renewcommand{\baselinestretch}{1.2} 
```

`2410.14703/src/superpolynoms.tex:3561-3563`
```latex
{\sf The case of $\F[[z^4,z^6+z^7]]$.}
This is $\r$ with the simplest non-torus algebraic 
cable, which is  $K=C\!ab(13,2)C\!ab(2,3)$.
```

`2410.14703/src/superpolynoms.tex:3624-3646`
```latex
The whole (uncolored) superpolynomial is: 
$\h^{mot}(q,t,a)=$
\renewcommand{\baselinestretch}{1.2} 
{\small
\(
1 + q t + q^8 t^8 + q^2 \bigl(t + t^2\bigr) 
 + 
 q^3 \bigl(t + t^2 + t^3\bigr) + q^4 \bigl(2 t^2 + t^3 + t^4\bigr) 
+  q^5 \bigl(2 t^3 + t^4 + t^5\bigr) + q^6 \bigl(2 t^4 + t^5 
+ t^6\bigr) + 
 q^7 \bigl(t^5 + t^6 + t^7\bigr) + 
a \bigl(q + q^2 \bigl(1 + t\bigr) + q^3 \bigl(1 + 2 t + t^2\bigr) 
+ q^4 \bigl(3 t + 2 t^2 + t^3\bigr) + 
    q^5 \bigl(t + 4 t^2 + 2 t^3 + t^4\bigr) + q^6 \bigl(t^2 
+ 4 t^3 + 2 t^4 + t^5\bigr) + 
    q^7 \bigl(t^3 + 3 t^4 + 2 t^5 + t^6\bigr) + q^8 \bigl(t^5 
+ t^6 + t^7\bigr)\bigr)+  
a^2 \bigl(q^3 + q^4 \bigl(1 + t\bigr) + q^5 \bigl(1 + 2 t + t^2\bigr) 
+ q^6 \bigl(2 t + 2 t^2 + t^3\bigr) + q^7 \bigl(2 t^2 + 2 t^3 
+ t^4\bigr) + q^8 \bigl(t^3 + t^4 + t^5\bigr)\bigr) +
a^3 \bigl(q^6 + q^7 t + q^8 t^2\bigr).
\)
}
```

C24 $(1+aq^i)$-decomposition of the same polynomial:

`2410.14703/src/superpolynoms.tex:3723-3736`
```latex
Let us provide the $(1\!+\!a q^i)$-decomposition for this 
superpolynomial: $\h^{mot}=$
\renewcommand{\baselinestretch}{1.2} 
{\small
\(
q^8 t^8\ +\  (1\!+\!a q) (1\!+\!a q^2) (1\!+\!a q^3) 
\bigl(1+q t+q^2 t^2\bigr)\ +\ 
(1\!+\!a q) (1\!+\!a q^2)\bigl(q^2 (1+q) t+q^2 (q+q^2) t^2
+q^2 (q+2 q^2+q^3) t^3
+q^2 (q^2+q^3) t^4+q^5 t^5\bigr)\ +\ (1\!+\!a q) \bigl(q^4 t^2+
(-1+q) q^4 t^3+2 q^6 t^4+q^4 (q^2+q^3) t^5+q^4 (q^2+q^3) t^6
+q^7 t^7\bigr).
\)
}
```

$q=1$ specialization (CD15, as a factor; the second factor is $\mathcal H_{Cab(13,2)T(3,2)}(q=1,t,a)$) [checked]:

`1509.08351/src/links.tex:11208-11226`
```latex
\Yboxdim7pt
This polynomial is self-dual. The positivity of
$\lxi\hat{\h}{}^{min}_{\l,\,'\!\l}(q,t,a)/(1-t)$ holds.
One has: $\lxi\hat{\h}{}^{min}_{\l,\,'\!\l}(q=1,t,a)=$
\renewcommand{\baselinestretch}{0.5} 
\noindent
{\small
\(
(1+a)\bigl(1+3 a+3 a^2+a^3+3 t+7 a t+5 a^2 t+a^3 t
+4 t^2+8 a t^2+5 a^2 t^2+a^3 t^2+4 t^3+8 a t^3+4 a^2 t^3+4 t^4
+6 a t^4+2 a^2 t^4+3 t^5+4 a t^5+a^2 t^5+2 t^6+2 a t^6+t^7
+a t^7+t^8\bigr),
\)
}
\renewcommand{\baselinestretch}{1.0} 
%\smallskip
where the second factor is that for $C\!ab(13,2)T(3,2)$,
which appears when making $\l$ empty.

```

### 4.2 (2,13)-cable of trefoil colored by $2\omega_1$ (two-row $\square\square$) — CD14 (5)

`1408.4348/src/cable.tex:2553-2856`
```latex
\subsection{\bf The case of two boxes}
We will consider now the iterated knot from the first
example with the weight $\,2\om_1$. We calculated
the one for $\om_2$ too; it is super-dual to the one we 
provide and will be omitted. Also, we will
omit the corresponding DAHA- Khovanov 
polynomial (calculated using the algorithm from
\cite{CJ}), since the Khovanov polynomial
are uncolored.

$$
(\mathbf 5):\ \vec\rr=\{3,2\},\,
\vec\ss\!=\!\{2,1\},\,
\t\!=\! C\!ab(13,2)T(3,2);\, 
\h_{\vec\rr,\vec\ss}(\square\square;q,t,a)=
$$
\renewcommand{\baselinestretch}{0.5} 
{\small
\(
1+t^{16} q^{32}+t^{15} q^{31}+t^{14} q^{31}
+t^{13} q^{31}
+t^{15} q^{30}+2 t^{14} q^{30}+2 t^{13} q^{30}
+2 t^{12} q^{30}
+t^{14} q^{29}+3 t^{13} q^{29}+4 t^{12} q^{29}
+3 t^{11} q^{29}
+t^{14} q^{28}+2 t^{13} q^{28}+5 t^{12} q^{28}
+6 t^{11} q^{28}
+4 t^{10} q^{28}+t^{13} q^{27}+3 t^{12} q^{27}
+7 t^{11} q^{27}
+8 t^{10} q^{27}+4 t^9 q^{27}+t^{13} q^{26}
+2 t^{12} q^{26}
+5 t^{11} q^{26}+10 t^{10} q^{26}+9 t^9 q^{26}
+2 t^8 q^{26}
+t^{12} q^{25}+3 t^{11} q^{25}+7 t^{10} q^{25}
+12 t^9 q^{25}
+7 t^8 q^{25}+t^7 q^{25}+t^{12} q^{24}+2 t^{11} q^{24}
+5 t^{10} q^{24}+10 t^9 q^{24}+14 t^8 q^{24}
+4 t^7 q^{24}
+t^{11} q^{23}+3 t^{10} q^{23}+7 t^9 q^{23}
+12 t^8 q^{23}
+11 t^7 q^{23}+t^6 q^{23}+t^{11} q^{22}+2 t^{10} q^{22}
+5 t^9 q^{22}+10 t^8 q^{22}+13 t^7 q^{22}
+6 t^6 q^{22}
+t^{10} q^{21}+3 t^9 q^{21}+7 t^8 q^{21}
+12 t^7 q^{21}
+10 t^6 q^{21}+t^5 q^{21}+t^{10} q^{20}+2 t^9 q^{20}
+5 t^8 q^{20}+10 t^7 q^{20}+13 t^6 q^{20}+3 t^5 q^{20}
+t^9 q^{19}+3 t^8 q^{19}+7 t^7 q^{19}+11 t^6 q^{19}
+9 t^5 q^{19}+t^9 q^{18}+2 t^8 q^{18}+5 t^7 q^{18}
+10 t^6 q^{18}+11 t^5 q^{18}+2 t^4 q^{18}+t^8 q^{17}
+3 t^7 q^{17}+7 t^6 q^{17}+10 t^5 q^{17}+4 t^4 q^{17}
+t^8 q^{16}+2 t^7 q^{16}+5 t^6 q^{16}+9 t^5 q^{16}
+8 t^4 q^{16}+t^7 q^{15}+3 t^6 q^{15}+7 t^5 q^{15}
+8 t^4 q^{15}+2 t^3 q^{15}+t^7 q^{14}+2 t^6 q^{14}
+5 t^5 q^{14}+8 t^4 q^{14}+3 t^3 q^{14}+t^6 q^{13}
+3 t^5 q^{13}+6 t^4 q^{13}+5 t^3 q^{13}+t^6 q^{12}
+2 t^5 q^{12}+5 t^4 q^{12}+6 t^3 q^{12}+t^2 q^{12}
+t^5 q^{11}+3 t^4 q^{11}+5 t^3 q^{11}+t^2 q^{11}
+t^5 q^{10}+2 t^4 q^{10}+4 t^3 q^{10}+3 t^2 q^{10}
+t^4 q^9+3 t^3 q^9+3 t^2 q^9+t^4 q^8+2 t^3 q^8
+3 t^2 q^8+t^3 q^7+2 t^2 q^7+t q^7+t^3 q^6+2 t^2 q^6
+t q^6+t^2 q^5+t q^5+t^2 q^4+t q^4+t q^3+t q^2
\)

\vfil
\noindent
\(
+a^6 \bigl(t^4 q^{35}+t^3 q^{34}+t^3 q^{33}+t^2 q^{33}
+t^2 q^{32}+t^2 q^{31}+t q^{30}+t q^{29}+q^{27}\bigr)
\)

\vfil
\noindent
\(
+a^5 \bigl(t^7 q^{35}+t^6 q^{35}+t^5 q^{35}+t^7 q^{34}
+2 t^6 q^{34}+2 t^5 q^{34}+t^4 q^{34}+2 t^6 q^{33}
+4 t^5 q^{33}+4 t^4 q^{33}+t^3 q^{33}+t^6 q^{32}
+4 t^5 q^{32}+5 t^4 q^{32}+3 t^3 q^{32}+2 t^5 q^{31}
+5 t^4 q^{31}+5 t^3 q^{31}+t^2 q^{31}+t^5 q^{30}
+4 t^4 q^{30}+6 t^3 q^{30}+3 t^2 q^{30}+2 t^4 q^{29}
+6 t^3 q^{29}+5 t^2 q^{29}+t^4 q^{28}+4 t^3 q^{28}
+5 t^2 q^{28}+t q^{28}+2 t^3 q^{27}+5 t^2 q^{27}
+3 t q^{27}+t^3 q^{26}+4 t^2 q^{26}+3 t q^{26}
+2 t^2 q^{25}+3 t q^{25}+q^{25}+t^2 q^{24}+3 t q^{24}
+q^{24}+2 t q^{23}+q^{23}+t q^{22}+q^{22}+q^{21}
+q^{20}\bigr)
\)

\vfil
\noindent
\(
+a^4 \bigl(t^9 q^{35}+t^8 q^{35}
+t^7 q^{35}+t^{10} q^{34}+2 t^9 q^{34}+4 t^8 q^{34}
+3 t^7 q^{34}+2 t^6 q^{34}+2 t^9 q^{33}+5 t^8 q^{33}
+8 t^7 q^{33}+6 t^6 q^{33}+3 t^5 q^{33}+t^9 q^{32}
+5 t^8 q^{32}+11 t^7 q^{32}+13 t^6 q^{32}+8 t^5 q^{32}
+2 t^4 q^{32}+2 t^8 q^{31}+9 t^7 q^{31}+16 t^6 q^{31}
+15 t^5 q^{31}+6 t^4 q^{31}+t^3 q^{31}+t^8 q^{30}
+5 t^7 q^{30}+16 t^6 q^{30}+21 t^5 q^{30}+14 t^4 q^{30}
+3 t^3 q^{30}+2 t^7 q^{29}+9 t^6 q^{29}+22 t^5 q^{29}
+20 t^4 q^{29}+8 t^3 q^{29}+t^7 q^{28}+5 t^6 q^{28}
+16 t^5 q^{28}+25 t^4 q^{28}+14 t^3 q^{28}+2 t^2 q^{28}
+2 t^6 q^{27}+9 t^5 q^{27}+21 t^4 q^{27}+20 t^3 q^{27}
+5 t^2 q^{27}+t^6 q^{26}+5 t^5 q^{26}+16 t^4 q^{26}
+22 t^3 q^{26}+10 t^2 q^{26}+2 t^5 q^{25}+9 t^4 q^{25}
+20 t^3 q^{25}+13 t^2 q^{25}+2 t q^{25}+t^5 q^{24}
+5 t^4 q^{24}+15 t^3 q^{24}+17 t^2 q^{24}+4 t q^{24}
+2 t^4 q^{23}+9 t^3 q^{23}+15 t^2 q^{23}+6 t q^{23}
+t^4 q^{22}+5 t^3 q^{22}+13 t^2 q^{22}+8 t q^{22}
+q^{22}+2 t^3 q^{21}+8 t^2 q^{21}+9 t q^{21}+q^{21}
+t^3 q^{20}+5 t^2 q^{20}+8 t q^{20}+2 q^{20}+2 t^2 q^{19}
+6 t q^{19}+2 q^{19}+t^2 q^{18}+4 t q^{18}+3 q^{18}
+2 t q^{17}+2 q^{17}+t q^{16}+2 q^{16}+q^{15}+q^{14}\bigr)
\)

\vfil
\noindent
\(
+a^3 \bigl(t^{10} q^{35}+t^{12} q^{34}+2 t^{11} q^{34}
+3 t^{10} q^{34}+3 t^9 q^{34}+t^8 q^{34}+t^{12} q^{33}
+3 t^{11} q^{33}+6 t^{10} q^{33}+8 t^9 q^{33}
+6 t^8 q^{33}+2 t^7 q^{33}+2 t^{11} q^{32}
+8 t^{10} q^{32}+14 t^9 q^{32}+16 t^8 q^{32}
+10 t^7 q^{32}+3 t^6 q^{32}+t^{11} q^{31}
+5 t^{10} q^{31}+15 t^9 q^{31}+25 t^8 q^{31}
+22 t^7 q^{31}+11 t^6 q^{31}+2 t^5 q^{31}
+2 t^{10} q^{30}+10 t^9 q^{30}+26 t^8 q^{30}
+36 t^7 q^{30}+26 t^6 q^{30}+9 t^5 q^{30}
+t^4 q^{30}+t^{10} q^{29}+5 t^9 q^{29}
+18 t^8 q^{29}+39 t^7 q^{29}+42 t^6 q^{29}
+21 t^5 q^{29}+4 t^4 q^{29}+2 t^9 q^{28}
+10 t^8 q^{28}+29 t^7 q^{28}+50 t^6 q^{28}
+39 t^5 q^{28}+12 t^4 q^{28}+t^3 q^{28}
+t^9 q^{27}+5 t^8 q^{27}+18 t^7 q^{27}
+43 t^6 q^{27}+52 t^5 q^{27}+25 t^4 q^{27}
+4 t^3 q^{27}+2 t^8 q^{26}+10 t^7 q^{26}
+29 t^6 q^{26}+53 t^5 q^{26}+41 t^4 q^{26}
+10 t^3 q^{26}+t^8 q^{25}+5 t^7 q^{25}
+18 t^6 q^{25}+43 t^5 q^{25}+50 t^4 q^{25}
+20 t^3 q^{25}+t^2 q^{25}+2 t^7 q^{24}
+10 t^6 q^{24}+29 t^5 q^{24}+50 t^4 q^{24}
+32 t^3 q^{24}+5 t^2 q^{24}+t^7 q^{23}
+5 t^6 q^{23}+18 t^5 q^{23}+41 t^4 q^{23}
+39 t^3 q^{23}+10 t^2 q^{23}+2 t^6 q^{22}
+10 t^5 q^{22}+28 t^4 q^{22}+41 t^3 q^{22}
+18 t^2 q^{22}+t q^{22}+t^6 q^{21}+5 t^5 q^{21}
+18 t^4 q^{21}+36 t^3 q^{21}+24 t^2 q^{21}
+3 t q^{21}+2 t^5 q^{20}+10 t^4 q^{20}+26 t^3 q^{20}
+27 t^2 q^{20}+6 t q^{20}+t^5 q^{19}
+5 t^4 q^{19}+17 t^3 q^{19}+26 t^2 q^{19}
+9 t q^{19}+2 t^4 q^{18}+10 t^3 q^{18}
+21 t^2 q^{18}+12 t q^{18}+q^{18}+t^4 q^{17}
+5 t^3 q^{17}+15 t^2 q^{17}+13 t q^{17}
+q^{17}+2 t^3 q^{16}+9 t^2 q^{16}+12 t q^{16}
+2 q^{16}+t^3 q^{15}+5 t^2 q^{15}+10 t q^{15}
+3 q^{15}+2 t^2 q^{14}+7 t q^{14}+3 q^{14}
+t^2 q^{13}+4 t q^{13}+3 q^{13}+2 t q^{12}
+3 q^{12}+t q^{11}+2 q^{11}+q^{10}+q^9\bigr)
\)

\vfil
\noindent
\(
+a^2 \bigl(t^{13} q^{34}+t^{12} q^{34}
+t^{11} q^{34}+t^{14} q^{33}+2 t^{13} q^{33}
+4 t^{12} q^{33}+4 t^{11} q^{33}+3 t^{10} q^{33}
+2 t^{13} q^{32}+5 t^{12} q^{32}+10 t^{11} q^{32}
+10 t^{10} q^{32}+6 t^9 q^{32}+t^8 q^{32}
+t^{13} q^{31}+5 t^{12} q^{31}+12 t^{11} q^{31}
+20 t^{10} q^{31}+18 t^9 q^{31}+8 t^8 q^{31}
+t^7 q^{31}+2 t^{12} q^{30}+9 t^{11} q^{30}
+21 t^{10} q^{30}+32 t^9 q^{30}+24 t^8 q^{30}
+9 t^7 q^{30}+t^6 q^{30}+t^{12} q^{29}
+5 t^{11} q^{29}+17 t^{10} q^{29}+35 t^9 q^{29}
+45 t^8 q^{29}+25 t^7 q^{29}+7 t^6 q^{29}
+2 t^{11} q^{28}+9 t^{10} q^{28}+27 t^9 q^{28}
+49 t^8 q^{28}+49 t^7 q^{28}+20 t^6 q^{28}
+3 t^5 q^{28}+t^{11} q^{27}+5 t^{10} q^{27}
+17 t^9 q^{27}+42 t^8 q^{27}+62 t^7 q^{27}
+44 t^6 q^{27}+11 t^5 q^{27}+t^4 q^{27}
+2 t^{10} q^{26}+9 t^9 q^{26}+27 t^8 q^{26}
+57 t^7 q^{26}+63 t^6 q^{26}+28 t^5 q^{26}
+3 t^4 q^{26}+t^{10} q^{25}+5 t^9 q^{25}
+17 t^8 q^{25}+42 t^7 q^{25}+70 t^6 q^{25}
+48 t^5 q^{25}+11 t^4 q^{25}+2 t^9 q^{24}
+9 t^8 q^{24}+27 t^7 q^{24}+56 t^6 q^{24}
+66 t^5 q^{24}+24 t^4 q^{24}+2 t^3 q^{24}
+t^9 q^{23}+5 t^8 q^{23}+17 t^7 q^{23}
+42 t^6 q^{23}+67 t^5 q^{23}+42 t^4 q^{23}
+6 t^3 q^{23}+2 t^8 q^{22}+9 t^7 q^{22}
+27 t^6 q^{22}+55 t^5 q^{22}+54 t^4 q^{22}
+15 t^3 q^{22}+t^8 q^{21}+5 t^7 q^{21}
+17 t^6 q^{21}+41 t^5 q^{21}+60 t^4 q^{21}
+26 t^3 q^{21}+2 t^2 q^{21}+2 t^7 q^{20}
+9 t^6 q^{20}+27 t^5 q^{20}+50 t^4 q^{20}
+38 t^3 q^{20}+5 t^2 q^{20}+t^7 q^{19}
+5 t^6 q^{19}+17 t^5 q^{19}+39 t^4 q^{19}
+43 t^3 q^{19}+12 t^2 q^{19}+2 t^6 q^{18}
+9 t^5 q^{18}+26 t^4 q^{18}+41 t^3 q^{18}
+18 t^2 q^{18}+t q^{18}+t^6 q^{17}
+5 t^5 q^{17}+17 t^4 q^{17}+34 t^3 q^{17}
+25 t^2 q^{17}+2 t q^{17}+2 t^5 q^{16}
+9 t^4 q^{16}+24 t^3 q^{16}+25 t^2 q^{16}
+5 t q^{16}+t^5 q^{15}+5 t^4 q^{15}
+16 t^3 q^{15}+25 t^2 q^{15}+8 t q^{15}
+2 t^4 q^{14}+9 t^3 q^{14}+19 t^2 q^{14}
+10 t q^{14}+t^4 q^{13}+5 t^3 q^{13}
+14 t^2 q^{13}+11 t q^{13}+q^{13}
+2 t^3 q^{12}+8 t^2 q^{12}+11 t q^{12}
+q^{12}+t^3 q^{11}+5 t^2 q^{11}+9 t q^{11}
+2 q^{11}+2 t^2 q^{10}+6 t q^{10}+2 q^{10}
+t^2 q^9+4 t q^9+3 q^9+2 t q^8+2 q^8+t q^7
+2 q^7+q^6+q^5\bigr)
\)

\vfil
\noindent
\(
+a \bigl(t^{15} q^{33}
+t^{14} q^{33}+t^{13} q^{33}+t^{15} q^{32}
+2 t^{14} q^{32}+3 t^{13} q^{32}+3 t^{12} q^{32}
+t^{11} q^{32}+2 t^{14} q^{31}+5 t^{13} q^{31}
+8 t^{12} q^{31}+7 t^{11} q^{31}+2 t^{10} q^{31}
+t^{14} q^{30}+4 t^{13} q^{30}+9 t^{12} q^{30}
+14 t^{11} q^{30}+11 t^{10} q^{30}+3 t^9 q^{30}
+2 t^{13} q^{29}+7 t^{12} q^{29}+16 t^{11} q^{29}
+22 t^{10} q^{29}+15 t^9 q^{29}+3 t^8 q^{29}
+t^{13} q^{28}+4 t^{12} q^{28}+12 t^{11} q^{28}
+24 t^{10} q^{28}+29 t^9 q^{28}+14 t^8 q^{28}
+2 t^7 q^{28}+2 t^{12} q^{27}+7 t^{11} q^{27}
+19 t^{10} q^{27}+34 t^9 q^{27}+32 t^8 q^{27}
+11 t^7 q^{27}+t^6 q^{27}+t^{12} q^{26}
+4 t^{11} q^{26}+12 t^{10} q^{26}+28 t^9 q^{26}
+42 t^8 q^{26}+27 t^7 q^{26}+5 t^6 q^{26}
+2 t^{11} q^{25}+7 t^{10} q^{25}+19 t^9 q^{25}
+38 t^8 q^{25}+43 t^7 q^{25}+17 t^6 q^{25}
+t^5 q^{25}+t^{11} q^{24}+4 t^{10} q^{24}
+12 t^9 q^{24}+28 t^8 q^{24}+46 t^7 q^{24}
+34 t^6 q^{24}+6 t^5 q^{24}+2 t^{10} q^{23}
+7 t^9 q^{23}+19 t^8 q^{23}+38 t^7 q^{23}
+46 t^6 q^{23}+17 t^5 q^{23}+t^4 q^{23}
+t^{10} q^{22}+4 t^9 q^{22}+12 t^8 q^{22}
+28 t^7 q^{22}+45 t^6 q^{22}+31 t^5 q^{22}
+4 t^4 q^{22}+2 t^9 q^{21}+7 t^8 q^{21}
+19 t^7 q^{21}+37 t^6 q^{21}+42 t^5 q^{21}
+11 t^4 q^{21}+t^9 q^{20}+4 t^8 q^{20}
+12 t^7 q^{20}+28 t^6 q^{20}+42 t^5 q^{20}
+21 t^4 q^{20}+t^3 q^{20}+2 t^8 q^{19}
+7 t^7 q^{19}+19 t^6 q^{19}+35 t^5 q^{19}
+32 t^4 q^{19}+5 t^3 q^{19}+t^8 q^{18}
+4 t^7 q^{18}+12 t^6 q^{18}+27 t^5 q^{18}
+35 t^4 q^{18}+11 t^3 q^{18}+2 t^7 q^{17}
+7 t^6 q^{17}+19 t^5 q^{17}+31 t^4 q^{17}
+18 t^3 q^{17}+t^2 q^{17}+t^7 q^{16}
+4 t^6 q^{16}+12 t^5 q^{16}+25 t^4 q^{16}
+23 t^3 q^{16}+3 t^2 q^{16}+2 t^6 q^{15}
+7 t^5 q^{15}+18 t^4 q^{15}+24 t^3 q^{15}
+7 t^2 q^{15}+t^6 q^{14}+4 t^5 q^{14}
+12 t^4 q^{14}+21 t^3 q^{14}+10 t^2 q^{14}
+2 t^5 q^{13}+7 t^4 q^{13}+16 t^3 q^{13}
+13 t^2 q^{13}+t q^{13}+t^5 q^{12}
+4 t^4 q^{12}+11 t^3 q^{12}+14 t^2 q^{12}
+2 t q^{12}+2 t^4 q^{11}+7 t^3 q^{11}
+12 t^2 q^{11}+3 t q^{11}+t^4 q^{10}
+4 t^3 q^{10}+9 t^2 q^{10}+5 t q^{10}
+2 t^3 q^9+6 t^2 q^9+6 t q^9+t^3 q^8
+4 t^2 q^8+5 t q^8+2 t^2 q^7+4 t q^7+q^7
+t^2 q^6+3 t q^6+q^6+2 t q^5+q^5+t q^4
+q^4+q^3+q^2\bigr).
\)
}
\renewcommand{\baselinestretch}{1.2} 
\smallskip

{\sf Betti numbers} (presumably):
\begin{align} \label{bett5}
&\h(q\!=\!1,t,a\!=\!0)\!=\!
1+6 t+17 t^2+32 t^3+48 t^4+62 t^5+70 t^6\\
+70& t^7\!+\!64 t^8\!+\!54 t^9\!+\!41 t^{10}\!+\!28 t^{11}\!+
\!18 t^{12}
\!+\!10 t^{13}\!+\!5 t^{14}\!+\!2 t^{15}\!+\!t^{16}\notag\\
&=\bigl(1+3 t+4 t^2+4 t^3+4 t^4+3 t^5+2 t^6+t^7+t^8\bigr)^2.
\notag
\end{align}

We put here (and below) ``presumably", since at the moment
it is not clear which Betti numbers this polynomial calculates
(even conjecturally).
The last relation is a particular case of Part $(iv)$ 
of Theorem \ref{STABILIZ} (the specialization at $q=1$). 
The full specialization formula in this case reads
as follows:

\begin{align*}
&\h_{\{3,2\},\{2,1\}}\,(2\om_1;q\!=\!1,t,a)
\!=\!\bigl(1\!+\!3 a\!+\!3 a^2\!+\!a^3\!+\!3 t\!+\!7 a t\!+\!
5 a^2 t\!+\!a^3 t\!+\!4 t^2\\
&+8 a t^2
+5 a^2 t^2+a^3 t^2+4 t^3+8 a t^3+4 a^2 t^3+4 t^4+6 a t^4+2 a^2 t^4
+3 t^5+\\
&4 a t^5+a^2 t^5+2 t^6+2 a t^6+t^7+a t^7+t^8\bigr)^2
=\h_{\{3,2\},\{2,1\}}\,(\om_1\,;\, q\!=\!1,t,a)^2.
\end{align*}
\smallskip
```

### 4.3 (2,13)-cable of trefoil, $\omega_2$ (only $q=1,a=0$) — CD14

`1408.4348/src/cable.tex:2858-2876`
```latex
{\mathversion{bold}$\h_{\{3,2\},\{2,1\}}(\om_2)$.}
The superpolynomial for the same iterated knot and
the weight $\om_2$ is connected with the one for
$2\om_1$ by the Super-Duality, so we will not list it.
However, let us provide its value at $a=0,q=1$ (presumably
a sequence of  Betti numbers):
\begin{align} \label{bett55}
&\h_{\{3,2\},\{2,1\}}\,(\om_2;q\!=\!1,t,a\!=\!0)\!=\!
1\! +\! 3 t\! + 7 t^2\! + 11 t^3\! + 18 t^4\! + 23 t^5\\
+ 29 &t^6+ 31 t^7 + 
36 t^8 + 35 t^9 + 37 t^{10} + 34 t^{11}+ 34 t^{12} 
+31 t^{13} + 31 t^{14}\notag\\
+ 25 &t^{15}
+ 25 t^{16} + 21 t^{17}+ 19 t^{18} + 15 t^{19} + 15 t^{20} + 
10 t^{21} + 10 t^{22}\notag\\ 
+ 7 t&^{23}+ 6 t^{24}+ 4 t^{25}+ 4 t^{26} + 2 t^{27} + 
 2 t^{28} + t^{29} + t^{30} + t^{32}.\notag
\end{align}

```

### 4.4 (2,15)-cable, $\vec r=\{3,2\},\vec s=\{2,3\}$ — CD14 (2) (identical in CP16:3131-3167 [checked])

`1408.4348/src/cable.tex:1843-1877`
```latex
{\mathversion{bold}$\h_{\{3,2\},\{2,3\}}$.}
Let us switch from $C\!ab(13,2)$ to
the ``next" $C\!ab(15,2)$.
$$
(\mathbf 2):\ \vec\rr=\{3,2\},\,
\vec\ss=\{2,3\},\,
\t= C\!ab(15,2)T(3,2);\ 
\h_{\,\vec\rr,\,\vec\ss}\,(\square\,;\,\,q,t,a)=
$$
\renewcommand{\baselinestretch}{0.5} 
{\small
\(
1+q t+q^2 t+q^3 t+q^2 t^2+q^3 t^2+2 q^4 t^2+q^3 t^3+q^4 t^3
+2 q^5 t^3
+q^4 t^4+q^5 t^4+2 q^6 t^4+q^5 t^5+q^6 t^5+2 q^7 t^5+q^6 t^6+q^7 t^6
+q^8 t^6+q^7 t^7+q^8 t^7+q^8 t^8+q^9 t^9+a^3 \bigl(q^6+q^7 t+q^8 t^2
+q^9 t^3\bigr)+a^2 \bigl(q^3+q^4+q^5+q^4 t+2 q^5 t+2 q^6 t+q^5 t^2
+2 q^6 t^2+2 q^7 t^2+q^6 t^3+2 q^7 t^3+2 q^8 t^3+q^7 t^4+2 q^8 t^4
+q^9 t^4+q^8 t^5+q^9 t^5+q^9 t^6\bigr)+a \bigl(q+q^2+q^3+q^2 t
+2 q^3 t+3 q^4 t+q^5 t+q^3 t^2+2 q^4 t^2+4 q^5 t^2+q^6 t^2
+q^4 t^3+2 q^5 t^3+4 q^6 t^3+q^7 t^3+q^5 t^4+2 q^6 t^4
+4 q^7 t^4+q^8 t^4+q^6 t^5+2 q^7 t^5+3 q^8 t^5+q^7 t^6
+2 q^8 t^6+q^9 t^6+q^8 t^7+q^9 t^7+q^9 t^8\bigr).
\)
}
\renewcommand{\baselinestretch}{1.2} 
\smallskip

{\sf Betti Numbers.}
\begin{align} \label{bett2}
\h(q\!=\!1,t,a\!=\!0)\!=\!
1\!+\!3 t\!+\!4 t^2\!+\!4 t^3\!+\!4 t^4\!+\!4 t^5\!+\!3 t^6
\!+\!2 t^7\!+\!t^8\!+\!t^9.
\end{align}
\smallskip
```

### 4.5 (3,19)-cable of trefoil, $\vec r=\{3,3\},\vec s=\{2,1\}$ — CD14 (3)

`1408.4348/src/cable.tex:1906-2126`
```latex
\subsection{\bf Degree 5 and beyond}
The $a$\~degree was $3$ in the previous examples;
let us consider two examples of uncolored DAHA superpolynomials
of degree $5$ (which is $\ss_1\rr_2-1$ in these cases) and
then some further examples of DAHA-Betti polynomials, which
are beyond known or conjectured formulas. We note 
that $D\!A\!H\!A'K\!h\!R$ do not coincide with 
$\tilde{K\!h}$ in these two examples, but the difference
matches our conjecture.


$$
(\mathbf 3):\ \vec\rr=\{3,3\},\,
\vec\ss=\{2,1\},\,
\t= C\!ab(19,3)T(3,2);\ 
\h_{\,\vec\rr,\,\vec\ss}\,(\square\,;\,\,q,t,a)=
$$
\renewcommand{\baselinestretch}{0.5} 
{\small
\(
1+q t+q^2 t+q^3 t+q^4 t+q^5 t+q^2 t^2+q^3 t^2+2 q^4 t^2
+2 q^5 t^2+3 q^6 t^2+2 q^7 t^2+q^8 t^2+q^3 t^3+q^4 t^3
+2 q^5 t^3+3 q^6 t^3+4 q^7 t^3+4 q^8 t^3+4 q^9 t^3
+q^{10} t^3+q^4 t^4+q^5 t^4+2 q^6 t^4+3 q^7 t^4
+5 q^8 t^4+5 q^9 t^4+6 q^{10} t^4+4 q^{11} t^4
+q^{12} t^4+q^5 t^5+q^6 t^5+2 q^7 t^5+3 q^8 t^5
+5 q^9 t^5+6 q^{10} t^5+7 q^{11} t^5+6 q^{12} t^5
+3 q^{13} t^5+q^6 t^6+q^7 t^6+2 q^8 t^6+3 q^9 t^6
+5 q^{10} t^6+6 q^{11} t^6+8 q^{12} t^6+7 q^{13} t^6
+3 q^{14} t^6+q^{15} t^6+q^7 t^7+q^8 t^7+2 q^9 t^7
+3 q^{10} t^7+5 q^{11} t^7+6 q^{12} t^7+8 q^{13} t^7
+8 q^{14} t^7+3 q^{15} t^7+q^8 t^8+q^9 t^8+2 q^{10} t^8
+3 q^{11} t^8+5 q^{12} t^8+6 q^{13} t^8+8 q^{14} t^8
+7 q^{15} t^8+3 q^{16} t^8+q^9 t^9+q^{10} t^9
+2 q^{11} t^9+3 q^{12} t^9+5 q^{13} t^9+6 q^{14} t^9
+8 q^{15} t^9+6 q^{16} t^9+q^{17} t^9+q^{10} t^{10}
+q^{11} t^{10}+2 q^{12} t^{10}+3 q^{13} t^{10}
+5 q^{14} t^{10}+6 q^{15} t^{10}+7 q^{16} t^{10}
+4 q^{17} t^{10}+q^{11} t^{11}+q^{12} t^{11}
+2 q^{13} t^{11}+3 q^{14} t^{11}+5 q^{15} t^{11}
+6 q^{16} t^{11}+6 q^{17} t^{11}+q^{18} t^{11}
+q^{12} t^{12}+q^{13} t^{12}+2 q^{14} t^{12}
+3 q^{15} t^{12}+5 q^{16} t^{12}+5 q^{17} t^{12}
+4 q^{18} t^{12}+q^{13} t^{13}+q^{14} t^{13}
+2 q^{15} t^{13}+3 q^{16} t^{13}+5 q^{17} t^{13}
+4 q^{18} t^{13}+q^{19} t^{13}+q^{14} t^{14}
+q^{15} t^{14}+2 q^{16} t^{14}+3 q^{17} t^{14}
+4 q^{18} t^{14}+2 q^{19} t^{14}+q^{15} t^{15}
+q^{16} t^{15}+2 q^{17} t^{15}+3 q^{18} t^{15}
+3 q^{19} t^{15}+q^{16} t^{16}+q^{17} t^{16}
+2 q^{18} t^{16}+2 q^{19} t^{16}+q^{20} t^{16}
+q^{17} t^{17}+q^{18} t^{17}+2 q^{19} t^{17}
+q^{20} t^{17}+q^{18} t^{18}+q^{19} t^{18}
+q^{20} t^{18}+q^{19} t^{19}+q^{20} t^{19}
+q^{20} t^{20}+q^{21} t^{21}
\)

\vfil
\noindent
\(
+a^5 \bigl(q^{15}
+q^{16} t+q^{17} t+q^{17} t^2+q^{18} t^2
+q^{19} t^2+q^{18} t^3+q^{19} t^3+q^{19} t^4
+q^{20} t^4+q^{20} t^5+q^{21} t^6\bigr)
\)

\vfil
\noindent
\(
+a^4 \bigl(q^{10}+q^{11}+q^{12}+q^{13}+q^{14}
+q^{11} t+2 q^{12} t+3 q^{13} t+3 q^{14} t
+3 q^{15} t+q^{16} t+q^{12} t^2+2 q^{13} t^2
+4 q^{14} t^2+5 q^{15} t^2+5 q^{16} t^2+3 q^{17} t^2
+q^{18} t^2+q^{13} t^3+2 q^{14} t^3+4 q^{15} t^3
+6 q^{16} t^3+6 q^{17} t^3+3 q^{18} t^3+q^{19} t^3
+q^{14} t^4+2 q^{15} t^4+4 q^{16} t^4+6 q^{17} t^4
+6 q^{18} t^4+3 q^{19} t^4+q^{15} t^5+2 q^{16} t^5
+4 q^{17} t^5+6 q^{18} t^5+5 q^{19} t^5+q^{20} t^5
+q^{16} t^6+2 q^{17} t^6+4 q^{18} t^6+5 q^{19} t^6
+3 q^{20} t^6+q^{17} t^7+2 q^{18} t^7+4 q^{19} t^7
+3 q^{20} t^7+q^{21} t^7+q^{18} t^8+2 q^{19} t^8
+3 q^{20} t^8+q^{21} t^8+q^{19} t^9+2 q^{20} t^9
+q^{21} t^9+q^{20} t^{10}+q^{21} t^{10}+q^{21} t^{11}\bigr)
\)

\vfil
%\smallskip
\noindent
\(
+a^3 \bigl(q^6+q^7+2 q^8+2 q^9+2 q^{10}+q^{11}+q^{12}
+q^7 t+2 q^8 t+4 q^9 t+6 q^{10} t+7 q^{11} t+6 q^{12} t
+4 q^{13} t+2 q^{14} t+q^8 t^2+2 q^9 t^2+5 q^{10} t^2
+8 q^{11} t^2+12 q^{12} t^2+12 q^{13} t^2+10 q^{14} t^2
+5 q^{15} t^2+2 q^{16} t^2+q^9 t^3+2 q^{10} t^3
+5 q^{11} t^3+9 q^{12} t^3+14 q^{13} t^3+17 q^{14} t^3
+15 q^{15} t^3+9 q^{16} t^3+3 q^{17} t^3+q^{18} t^3
+q^{10} t^4+2 q^{11} t^4+5 q^{12} t^4+9 q^{13} t^4
+15 q^{14} t^4+19 q^{15} t^4+18 q^{16} t^4+10 q^{17} t^4
+3 q^{18} t^4+q^{11} t^5+2 q^{12} t^5+5 q^{13} t^5
+9 q^{14} t^5+15 q^{15} t^5+20 q^{16} t^5+18 q^{17} t^5
+9 q^{18} t^5+2 q^{19} t^5+q^{12} t^6+2 q^{13} t^6
+5 q^{14} t^6+9 q^{15} t^6+15 q^{16} t^6+19 q^{17} t^6
+15 q^{18} t^6+5 q^{19} t^6+q^{13} t^7+2 q^{14} t^7
+5 q^{15} t^7+9 q^{16} t^7+15 q^{17} t^7+17 q^{18} t^7
+10 q^{19} t^7+2 q^{20} t^7+q^{14} t^8+2 q^{15} t^8
+5 q^{16} t^8+9 q^{17} t^8+14 q^{18} t^8+12 q^{19} t^8
+4 q^{20} t^8+q^{15} t^9+2 q^{16} t^9+5 q^{17} t^9
+9 q^{18} t^9+12 q^{19} t^9+6 q^{20} t^9+q^{21} t^9
+q^{16} t^{10}+2 q^{17} t^{10}+5 q^{18} t^{10}
+8 q^{19} t^{10}+7 q^{20} t^{10}+q^{21} t^{10}
+q^{17} t^{11}+2 q^{18} t^{11}+5 q^{19} t^{11}
+6 q^{20} t^{11}+2 q^{21} t^{11}+q^{18} t^{12}
+2 q^{19} t^{12}+4 q^{20} t^{12}+2 q^{21} t^{12}
+q^{19} t^{13}+2 q^{20} t^{13}+2 q^{21} t^{13}
+q^{20} t^{14}+q^{21} t^{14}+q^{21} t^{15}\bigr)
\)

\vfil
%\smallskip
\noindent
\(
+a^2 \bigl(q^3+q^4+2 q^5+2 q^6+2 q^7+q^8+q^9+q^4 t
+2 q^5 t+4 q^6 t+6 q^7 t+8 q^8 t+7 q^9 t+6 q^{10} t
+3 q^{11} t+q^{12} t+q^5 t^2+2 q^6 t^2+5 q^7 t^2
+8 q^8 t^2+13 q^9 t^2+15 q^{10} t^2+15 q^{11} t^2
+10 q^{12} t^2+5 q^{13} t^2+q^{14} t^2+q^6 t^3
+2 q^7 t^3+5 q^8 t^3+9 q^9 t^3+15 q^{10} t^3
+20 q^{11} t^3+24 q^{12} t^3+19 q^{13} t^3
+11 q^{14} t^3+4 q^{15} t^3+q^{16} t^3+q^7 t^4+2 q^8 t^4
+5 q^9 t^4+9 q^{10} t^4+16 q^{11} t^4+22 q^{12} t^4
+29 q^{13} t^4+27 q^{14} t^4+16 q^{15} t^4+6 q^{16} t^4
+q^{17} t^4+q^8 t^5+2 q^9 t^5+5 q^{10} t^5+9 q^{11} t^5
+16 q^{12} t^5+23 q^{13} t^5+31 q^{14} t^5+30 q^{15} t^5
+19 q^{16} t^5+6 q^{17} t^5+q^{18} t^5+q^9 t^6
+2 q^{10} t^6+5 q^{11} t^6+9 q^{12} t^6+16 q^{13} t^6
+23 q^{14} t^6+32 q^{15} t^6+30 q^{16} t^6+16 q^{17} t^6
+4 q^{18} t^6+q^{10} t^7+2 q^{11} t^7+5 q^{12} t^7
+9 q^{13} t^7+16 q^{14} t^7+23 q^{15} t^7+31 q^{16} t^7
+27 q^{17} t^7+11 q^{18} t^7+q^{19} t^7+q^{11} t^8
+2 q^{12} t^8+5 q^{13} t^8+9 q^{14} t^8+16 q^{15} t^8
+23 q^{16} t^8+29 q^{17} t^8+19 q^{18} t^8+5 q^{19} t^8
+q^{12} t^9+2 q^{13} t^9+5 q^{14} t^9+9 q^{15} t^9
+16 q^{16} t^9+22 q^{17} t^9+24 q^{18} t^9+10 q^{19} t^9
+q^{20} t^9+q^{13} t^{10}+2 q^{14} t^{10}+5 q^{15} t^{10}
+9 q^{16} t^{10}+16 q^{17} t^{10}+20 q^{18} t^{10}
+15 q^{19} t^{10}+3 q^{20} t^{10}+q^{14} t^{11}
+2 q^{15} t^{11}+5 q^{16} t^{11}+9 q^{17} t^{11}
+15 q^{18} t^{11}+15 q^{19} t^{11}+6 q^{20} t^{11}
+q^{15} t^{12}+2 q^{16} t^{12}+5 q^{17} t^{12}
+9 q^{18} t^{12}+13 q^{19} t^{12}+7 q^{20} t^{12}
+q^{21} t^{12}+q^{16} t^{13}+2 q^{17} t^{13}+5 q^{18} t^{13}
+8 q^{19} t^{13}+8 q^{20} t^{13}+q^{21} t^{13}+q^{17} t^{14}
+2 q^{18} t^{14}+5 q^{19} t^{14}+6 q^{20} t^{14}
+2 q^{21} t^{14}+q^{18} t^{15}+2 q^{19} t^{15}
+4 q^{20} t^{15}+2 q^{21} t^{15}+q^{19} t^{16}
+2 q^{20} t^{16}+2 q^{21} t^{16}+q^{20} t^{17}
+q^{21} t^{17}+q^{21} t^{18}\bigr)
\)

\vfil
%\smallskip
\noindent
\(
+a \bigl(q+q^2+q^3+q^4+q^5+q^2 t+2 q^3 t+3 q^4 t
+4 q^5 t+5 q^6 t+4 q^7 t+2 q^8 t+q^9 t+q^3 t^2
+2 q^4 t^2+4 q^5 t^2+6 q^6 t^2+9 q^7 t^2+10 q^8 t^2
+9 q^9 t^2+5 q^{10} t^2+2 q^{11} t^2+q^4 t^3
+2 q^5 t^3+4 q^6 t^3+7 q^7 t^3+11 q^8 t^3+14 q^9 t^3
+16 q^{10} t^3+13 q^{11} t^3+6 q^{12} t^3+2 q^{13} t^3
+q^5 t^4+2 q^6 t^4+4 q^7 t^4+7 q^8 t^4+12 q^9 t^4
+16 q^{10} t^4+20 q^{11} t^4+20 q^{12} t^4+13 q^{13} t^4
+4 q^{14} t^4+q^{15} t^4+q^6 t^5+2 q^7 t^5+4 q^8 t^5
+7 q^9 t^5+12 q^{10} t^5+17 q^{11} t^5+22 q^{12} t^5
+24 q^{13} t^5+17 q^{14} t^5+7 q^{15} t^5+q^{16} t^5
+q^7 t^6+2 q^8 t^6+4 q^9 t^6+7 q^{10} t^6+12 q^{11} t^6
+17 q^{12} t^6+23 q^{13} t^6+26 q^{14} t^6+18 q^{15} t^6
+7 q^{16} t^6+q^{17} t^6+q^8 t^7+2 q^9 t^7+4 q^{10} t^7
+7 q^{11} t^7+12 q^{12} t^7+17 q^{13} t^7+23 q^{14} t^7
+26 q^{15} t^7+17 q^{16} t^7+4 q^{17} t^7+q^9 t^8
+2 q^{10} t^8+4 q^{11} t^8+7 q^{12} t^8+12 q^{13} t^8
+17 q^{14} t^8+23 q^{15} t^8+24 q^{16} t^8+13 q^{17} t^8
+2 q^{18} t^8+q^{10} t^9+2 q^{11} t^9+4 q^{12} t^9
+7 q^{13} t^9+12 q^{14} t^9+17 q^{15} t^9+22 q^{16} t^9
+20 q^{17} t^9+6 q^{18} t^9+q^{11} t^{10}+2 q^{12} t^{10}
+4 q^{13} t^{10}+7 q^{14} t^{10}+12 q^{15} t^{10}
+17 q^{16} t^{10}+20 q^{17} t^{10}+13 q^{18} t^{10}
+2 q^{19} t^{10}+q^{12} t^{11}+2 q^{13} t^{11}
+4 q^{14} t^{11}+7 q^{15} t^{11}+12 q^{16} t^{11}
+16 q^{17} t^{11}+16 q^{18} t^{11}+5 q^{19} t^{11}
+q^{13} t^{12}+2 q^{14} t^{12}+4 q^{15} t^{12}
+7 q^{16} t^{12}+12 q^{17} t^{12}+14 q^{18} t^{12}
+9 q^{19} t^{12}+q^{20} t^{12}+q^{14} t^{13}
+2 q^{15} t^{13}+4 q^{16} t^{13}+7 q^{17} t^{13}
+11 q^{18} t^{13}+10 q^{19} t^{13}+2 q^{20} t^{13}
+q^{15} t^{14}+2 q^{16} t^{14}+4 q^{17} t^{14}
+7 q^{18} t^{14}+9 q^{19} t^{14}+4 q^{20} t^{14}
+q^{16} t^{15}+2 q^{17} t^{15}+4 q^{18} t^{15}
+6 q^{19} t^{15}+5 q^{20} t^{15}+q^{17} t^{16}
+2 q^{18} t^{16}+4 q^{19} t^{16}+4 q^{20} t^{16}
+q^{21} t^{16}+q^{18} t^{17}+2 q^{19} t^{17}
+3 q^{20} t^{17}+q^{21} t^{17}+q^{19} t^{18}
+2 q^{20} t^{18}+q^{21} t^{18}+q^{20} t^{19}
+q^{21} t^{19}+q^{21} t^{20}\bigr).
\)
}
\renewcommand{\baselinestretch}{1.2} 
\vfil
\smallskip

{\sf  Betti numbers} (conjecturally):
\begin{align} \label{bett3}
&\h(q\!=\!1,t,a\!=\!0)=
1 + 5 t + 12 t^2 + 20 t^3 + 28 t^4 + 34 t^5 
+ 37 t^6 \\
&+ 37 t^7
+ 36 t^8 + 33 t^9 + 29 t^{10}
+ 25 t^{11} + 21 t^{12} 
+ 17 t^{13} + 13 t^{14}\notag\\
&+ 10 t^{15} + 7 t^{16} + 5 t^{17} + 3 t^{18} + 2 t^{19} 
+ t^{20} + t^{21}.\notag
\end{align}
```

### 4.6 (2,25)-cable of $T(4,3)$, $\vec r=\{4,2\},\vec s=\{3,1\}$ — CD14 (4)

`1408.4348/src/cable.tex:2193-2359`
```latex
{\mathversion{bold}$\h_{\{4,2\},\{3,1\}}$.}
The next example will be of the same deg$_a=5$.
Recall that deg${}_a\le (\ss_1\rr_2-1)$ 
in the uncolored case with $2$ iterations; conjecturally
they coincide for algebraic knots. Generally, this
must be multiplied by \,ord$(\la_b)$.
$$
(\mathbf 4):\ \vec\rr=\{4,2\},\, \vec\ss=\{3,1\},\,
\t= C\!ab(25,2)T(4,3);\ 
\h_{\,\vec\rr,\,\vec\ss}\,(\square\,;\,\,q,t,a)=
$$

\vfil
\renewcommand{\baselinestretch}{0.5} 
{\small
\(
1+q t+q^2 t+q^3 t+q^4 t+q^5 t+q^2 t^2+q^3 t^2+2 q^4 t^2
+2 q^5 t^2+3 q^6 t^2+q^7 t^2+q^8 t^2+q^3 t^3+q^4 t^3
+2 q^5 t^3+3 q^6 t^3+4 q^7 t^3+3 q^8 t^3+3 q^9 t^3
+q^4 t^4+q^5 t^4+2 q^6 t^4+3 q^7 t^4+5 q^8 t^4+4 q^9 t^4
+5 q^{10} t^4+q^{11} t^4+q^5 t^5+q^6 t^5+2 q^7 t^5
+3 q^8 t^5+5 q^9 t^5+5 q^{10} t^5+6 q^{11} t^5
+q^{12} t^5+q^6 t^6+q^7 t^6+2 q^8 t^6+3 q^9 t^6
+5 q^{10} t^6+5 q^{11} t^6+7 q^{12} t^6+q^{13} t^6
+q^7 t^7+q^8 t^7+2 q^9 t^7+3 q^{10} t^7+5 q^{11} t^7
+5 q^{12} t^7+6 q^{13} t^7+q^{14} t^7+q^8 t^8+q^9 t^8
+2 q^{10} t^8+3 q^{11} t^8+5 q^{12} t^8+5 q^{13} t^8
+5 q^{14} t^8+q^9 t^9+q^{10} t^9+2 q^{11} t^9+3 q^{12} t^9
+5 q^{13} t^9+4 q^{14} t^9+3 q^{15} t^9+q^{10} t^{10}
+q^{11} t^{10}+2 q^{12} t^{10}+3 q^{13} t^{10}
+5 q^{14} t^{10}+3 q^{15} t^{10}+q^{16} t^{10}
+q^{11} t^{11}+q^{12} t^{11}+2 q^{13} t^{11}
+3 q^{14} t^{11}+4 q^{15} t^{11}+q^{16} t^{11}
+q^{12} t^{12}+q^{13} t^{12}+2 q^{14} t^{12}
+3 q^{15} t^{12}+3 q^{16} t^{12}+q^{13} t^{13}
+q^{14} t^{13}+2 q^{15} t^{13}+2 q^{16} t^{13}
+q^{17} t^{13}+q^{14} t^{14}+q^{15} t^{14}
+2 q^{16} t^{14}+q^{17} t^{14}+q^{15} t^{15}
+q^{16} t^{15}+q^{17} t^{15}+q^{16} t^{16}
+q^{17} t^{16}+q^{17} t^{17}+q^{18} t^{18}
\)

\vfil
\noindent
\(+a^5 \bigl(q^{15}+q^{16} t+q^{17} t^2+q^{18} t^3\bigr)
\)

\vfil
\noindent
\(
+a^4 \bigl(q^{10}+q^{11}+q^{12}+q^{13}+q^{14}
+q^{11} t+2 q^{12} t+2 q^{13} t+2 q^{14} t
+2 q^{15} t+q^{12} t^2+2 q^{13} t^2+3 q^{14} t^2
+3 q^{15} t^2+2 q^{16} t^2+q^{13} t^3+2 q^{14} t^3
+3 q^{15} t^3+3 q^{16} t^3+2 q^{17} t^3+q^{14} t^4
+2 q^{15} t^4+3 q^{16} t^4+2 q^{17} t^4+q^{18} t^4
+q^{15} t^5+2 q^{16} t^5+2 q^{17} t^5+q^{18} t^5
+q^{16} t^6+2 q^{17} t^6+q^{18} t^6+q^{17} t^7
+q^{18} t^7+q^{18} t^8\bigr)
\)

\vfil
\noindent
\(
+a^3 \bigl(q^6+q^7
+2 q^8+2 q^9+2 q^{10}+q^{11}+q^{12}+q^7 t+2 q^8 t
+4 q^9 t+5 q^{10} t+6 q^{11} t+4 q^{12} t+3 q^{13} t
+q^{14} t+q^8 t^2+2 q^9 t^2+5 q^{10} t^2
+7 q^{11} t^2+10 q^{12} t^2+7 q^{13} t^2
+5 q^{14} t^2+q^{15} t^2+q^9 t^3+2 q^{10} t^3
+5 q^{11} t^3+8 q^{12} t^3+12 q^{13} t^3
+9 q^{14} t^3+6 q^{15} t^3+q^{16} t^3+q^{10} t^4
+2 q^{11} t^4+5 q^{12} t^4+8 q^{13} t^4+13 q^{14} t^4
+9 q^{15} t^4+5 q^{16} t^4+q^{17} t^4+q^{11} t^5
+2 q^{12} t^5+5 q^{13} t^5+8 q^{14} t^5+12 q^{15} t^5
+7 q^{16} t^5+3 q^{17} t^5+q^{12} t^6+2 q^{13} t^6
+5 q^{14} t^6+8 q^{15} t^6+10 q^{16} t^6+4 q^{17} t^6
+q^{18} t^6+q^{13} t^7+2 q^{14} t^7+5 q^{15} t^7
+7 q^{16} t^7+6 q^{17} t^7+q^{18} t^7+q^{14} t^8
+2 q^{15} t^8+5 q^{16} t^8+5 q^{17} t^8+2 q^{18} t^8
+q^{15} t^9+2 q^{16} t^9+4 q^{17} t^9+2 q^{18} t^9
+q^{16} t^{10}+2 q^{17} t^{10}+2 q^{18} t^{10}
+q^{17} t^{11}+q^{18} t^{11}+q^{18} t^{12}\bigr)
\)

\vfil
\noindent
\(
+a^2 \bigl(q^3+q^4+2 q^5+2 q^6+2 q^7+q^8+q^9+q^4 t
+2 q^5 t+4 q^6 t+6 q^7 t+7 q^8 t+6 q^9 t+5 q^{10} t
+2 q^{11} t+q^{12} t+q^5 t^2+2 q^6 t^2+5 q^7 t^2
+8 q^8 t^2+12 q^9 t^2+12 q^{10} t^2+11 q^{11} t^2
+5 q^{12} t^2+2 q^{13} t^2+q^6 t^3+2 q^7 t^3+5 q^8 t^3
+9 q^9 t^3+14 q^{10} t^3+17 q^{11} t^3+17 q^{12} t^3
+8 q^{13} t^3+3 q^{14} t^3+q^7 t^4+2 q^8 t^4+5 q^9 t^4
+9 q^{10} t^4+15 q^{11} t^4+19 q^{12} t^4+20 q^{13} t^4
+9 q^{14} t^4+3 q^{15} t^4+q^8 t^5+2 q^9 t^5+5 q^{10} t^5
+9 q^{11} t^5+15 q^{12} t^5+20 q^{13} t^5+20 q^{14} t^5
+8 q^{15} t^5+2 q^{16} t^5+q^9 t^6+2 q^{10} t^6
+5 q^{11} t^6+9 q^{12} t^6+15 q^{13} t^6+19 q^{14} t^6
+17 q^{15} t^6+5 q^{16} t^6+q^{17} t^6+q^{10} t^7
+2 q^{11} t^7+5 q^{12} t^7+9 q^{13} t^7+15 q^{14} t^7
+17 q^{15} t^7+11 q^{16} t^7+2 q^{17} t^7+q^{11} t^8
+2 q^{12} t^8+5 q^{13} t^8+9 q^{14} t^8+14 q^{15} t^8
+12 q^{16} t^8+5 q^{17} t^8+q^{12} t^9+2 q^{13} t^9
+5 q^{14} t^9+9 q^{15} t^9+12 q^{16} t^9+6 q^{17} t^9
+q^{18} t^9+q^{13} t^{10}+2 q^{14} t^{10}+5 q^{15} t^{10}
+8 q^{16} t^{10}+7 q^{17} t^{10}+q^{18} t^{10}
+q^{14} t^{11}+2 q^{15} t^{11}+5 q^{16} t^{11}
+6 q^{17} t^{11}+2 q^{18} t^{11}+q^{15} t^{12}
+2 q^{16} t^{12}+4 q^{17} t^{12}+2 q^{18} t^{12}
+q^{16} t^{13}+2 q^{17} t^{13}+2 q^{18} t^{13}
+q^{17} t^{14}+q^{18} t^{14}+q^{18} t^{15}\bigr)
\)

\vfil
\noindent
\(
+a \bigl(q+q^2+q^3+q^4+q^5+q^2 t+2 q^3 t+3 q^4 t
+4 q^5 t+5 q^6 t+3 q^7 t+2 q^8 t+q^9 t+q^3 t^2
+2 q^4 t^2+4 q^5 t^2+6 q^6 t^2+9 q^7 t^2+8 q^8 t^2
+7 q^9 t^2+3 q^{10} t^2+q^{11} t^2+q^4 t^3+2 q^5 t^3
+4 q^6 t^3+7 q^7 t^3+11 q^8 t^3+12 q^9 t^3
+13 q^{10} t^3+7 q^{11} t^3+2 q^{12} t^3
+q^5 t^4+2 q^6 t^4+4 q^7 t^4+7 q^8 t^4
+12 q^9 t^4+14 q^{10} t^4+17 q^{11} t^4
+10 q^{12} t^4+3 q^{13} t^4+q^6 t^5+2 q^7 t^5
+4 q^8 t^5+7 q^9 t^5+12 q^{10} t^5+15 q^{11} t^5
+19 q^{12} t^5+11 q^{13} t^5+3 q^{14} t^5
+q^7 t^6+2 q^8 t^6+4 q^9 t^6+7 q^{10} t^6
+12 q^{11} t^6+15 q^{12} t^6+19 q^{13} t^6
+10 q^{14} t^6+2 q^{15} t^6+q^8 t^7+2 q^9 t^7
+4 q^{10} t^7+7 q^{11} t^7+12 q^{12} t^7
+15 q^{13} t^7+17 q^{14} t^7+7 q^{15} t^7
+q^{16} t^7+q^9 t^8+2 q^{10} t^8+4 q^{11} t^8
+7 q^{12} t^8+12 q^{13} t^8+14 q^{14} t^8
+13 q^{15} t^8+3 q^{16} t^8+q^{10} t^9
+2 q^{11} t^9+4 q^{12} t^9+7 q^{13} t^9
+12 q^{14} t^9+12 q^{15} t^9+7 q^{16} t^9
+q^{17} t^9+q^{11} t^{10}+2 q^{12} t^{10}
+4 q^{13} t^{10}+7 q^{14} t^{10}+11 q^{15} t^{10}
+8 q^{16} t^{10}+2 q^{17} t^{10}+q^{12} t^{11}
+2 q^{13} t^{11}+4 q^{14} t^{11}+7 q^{15} t^{11}
+9 q^{16} t^{11}+3 q^{17} t^{11}+q^{13} t^{12}
+2 q^{14} t^{12}+4 q^{15} t^{12}+6 q^{16} t^{12}
+5 q^{17} t^{12}+q^{14} t^{13}+2 q^{15} t^{13}
+4 q^{16} t^{13}+4 q^{17} t^{13}+q^{18} t^{13}
+q^{15} t^{14}+2 q^{16} t^{14}+3 q^{17} t^{14}
+q^{18} t^{14}+q^{16} t^{15}+2 q^{17} t^{15}
+q^{18} t^{15}+q^{17} t^{16}+q^{18} t^{16}
+q^{18} t^{17}\bigr).
\)
}
\renewcommand{\baselinestretch}{1.2} 
\smallskip


{\sf Betti numbers} (conjecturally):
\begin{align} \label{bett4}
&\h(q\!=\!1,t,a\!=\!0)\!=\!
1\!+5 t\!+\!11 t^2\!+\!17 t^3\!+\!22 t^4\!+\!24 t^5\!+
\!25 t^6\!+\!24 t^7\\
+22 t^8\!&+\!19 t^9\!+\!16 t^{10}
\!+\!12 t^{11}\!+\!10 t^{12}
\!+\!7 t^{13}\!+\!5 t^{14}\!+\!3 t^{15}\!+\!2 t^{16}
\!+\!t^{17}\!+\!t^{18}.\notag
\end{align}
```

### 4.7 DAHA-Betti polynomials only (other algebraic cables) — CD14

`1408.4348/src/cable.tex:2423-2551`
```latex
{\sf Further Betti numbers.}
First of all, our DAHA-Betti polynomials (as far as we
calculated them) fully match
the table before Conjecture 23 in \cite{Pi}. For instance,
the coefficients of $\h_{\{5,2\},\{2,1\}}(q\!=\!1,t,a\!=\!0)$ 
for $C\!ab(21,2)T(5,2)$ 
constitute the line there for the Puiseux exponents
$(4,10,11)$. Without 
posting the (known) formulas for the superpolynomials, let us
provide below only the DAHA-Betti polynomials,
which conjecturally give the corresponding Betti numbers.
\smallskip

{\mathversion{bold}$\h_{\{4,2\},\{3,1\}}$.} In this example,
the knot, the ring of the singularity, 
the valuation semigroup and the arithmetic genus
are as follows:
$$ 
\t=C\!ab(20,3)T(3,2),\  \r=C[[t^6,t^9+t^{11}]],\
\Gamma=\lan 6,9,20\ran,\ \de=22. 
$$
The counterpart of the DAHA-Betti polynomial (\ref{bett3}) reads:
\begin{align} \label{bett3x}
&\h(q\!=\!1,t,a\!=\!0)=
1 + 5 t + 12 t^2 + 20 t^3 + 28 t^4 + 35 t^5 + 39 t^6\\
&+ 40 t^7 + 
39 t^8 +  37 t^9 + 33 t^{10} + 29 t^{11} + 25 t^{12} + 21 t^{13} 
+ 17 t^{14}\notag\\
&+ 13 t^{15} + 
 10 t^{16} + 7 t^{17} + 5 t^{18} + 3 t^{19} + 2 t^{20} + 
t^{21} + t^{22}.\notag
\end{align}
The arithmetic genus of $\c=\c_{\{3,3\},\{2,2\}}$
equals $\de=22$ and $e(J(\c))=423$ in this case; the 
Puiseux exponents are $(6,9,11)$, deg$_a$=5.
\medskip

{\mathversion{bold}$\h_{\{5,2\},\{3,2\}}$.} The Puiseux exponents
are $(6,15,17)$, deg$_a=5$,
$$
\t=C\!ab(32,3)T(5,2),\  \r=\C[[t^6,t^{15}+t^{17}]],\
\Gamma=\lan 6,15,32\ran,\ \de=37, 
$$
%and the counterpart of (\ref{bett3}) is
\begin{align} \label{bett3xx}
&\h(q\!=\!1,t,a\!=\!0)=
1+5 t+15 t^2+32 t^3+55 t^4+81 t^5+108 t^6\\
&+134 t^7+157 t^8+175 t^9+186 t^{10}+192 t^{11}+192 t^{12}
+189 t^{13}\notag\\
&+181 t^{14}+172 t^{15}+159 t^{16}+147 t^{17}+132 t^{18}
+120 t^{19}+105 t^{20}\notag\\
&+93 t^{21}+79 t^{22}+68 t^{23}+56 t^{24}+47 t^{25}+37 t^{26}
+30 t^{27}+23 t^{28}\notag\\
&+18 t^{29}+13 t^{30}+10 t^{31}+7 t^{32}
+5 t^{33}+3 t^{34}+2 t^{35}+t^{36}+t^{37}.\notag
\end{align}
The Euler number of $J(\c_{\{5,2\},\{3,2\}})$ 
is $3031$.
\medskip


{\mathversion{bold}$\h_{\{3,4\},\{2,1\}}$.} In this
case, deg$_a=7$ and the Puiseux exponents are $(8,12,13)$.
 One has: 
$$
\t=C\!ab(25,4)T(3,2),\  \r=\C[[t^8,t^{12}+t^{13}]],\
\Gamma=\lan 8,12,25\ran,\ \de=40, 
$$
and the DAHA-Bette polynomial is
\begin{align} \label{bett3xxx}
&\h(q\!=\!1,t,a\!=\!0)=
1+7 t+24 t^2+56 t^3+104 t^4+165 t^5+232 t^6\\
&+297 t^7+355 t^8
+402 t^9+435 t^{10}+454 t^{11}+461 t^{12}+456 t^{13}\notag\\
&+442 t^{14}+420 t^{15}+394 t^{16}
+362 t^{17}+330 t^{18}+295 t^{19}\notag\\
&+262 t^{20}+229 t^{21}+199 t^{22}
+168 t^{23}+143 t^{24}+118 t^{25}\notag\\
&+97 t^{26}+78 t^{27}+63 t^{28}+48 t^{29}
+38 t^{30}+28 t^{31}+21 t^{32}\notag\\
&+15 t^{33}+11 t^{34}
+7 t^{35}+5 t^{36}+3 t^{37}+2 t^{38}+t^{39}+t^{40}.\notag
\end{align}
The Euler number of $J(\c_{\{3,4\},\{2,1\}})$ 
is $7229$. The latter follows from \cite{ObS,Ma};
the Betti numbers are conjectured.

\smallskip
{\mathversion{bold} {\em Double iteration\,}: $\vec\rr=\{3,2,2\},
\vec\ss=\{2,1,1\}$}. The corresponding knot is
$\t=C\!ab(53,2)C\!ab(13,2)T(3,2)$ and deg$_a=7$.
In this case, 
$$
%\t=C\!ab(53,2)C\!ab(13,2)T(3,2),\  
\r=\C[[t^8,t^{12}+t^{14}+t^{15}]],\
\Gamma=\lan 8,12,26,53\ran,\ \de=42, 
$$
and the DAHA-Bette polynomial is
\begin{align} \label{bett3xxxx}
&\h(q\!=\!1,t,a\!=\!0)=
1+7 t+24 t^2+56 t^3+104 t^4+166 t^5\\
&+236 t^6+306 t^7
+370 t^8+424 t^9+465 t^{10}+492 t^{11}+507 t^{12}\notag\\
&+510 t^{13}+504 t^{14}+488 t^{15}+466 t^{16}+437 t^{17}
+406 t^{18}+370 t^{19}\notag\\
&+335 t^{20}+298 t^{21}+264 t^{22}
+230 t^{23}+199 t^{24}+168 t^{25}+143 t^{26}\notag\\
&+118 t^{27}
+97 t^{28}+78 t^{29}+63 t^{30}+48 t^{31}+38 t^{32}+28 t^{33}\notag\\
&+21 t^{34}+15 t^{35}+11 t^{36}
+7 t^{37}+5 t^{38}+3 t^{39}+2 t^{40}+t^{41}+t^{42}.\notag
\end{align}
\medskip
The Euler number
of $J(\c_{\{3,2,2\},\{2,1,1\}})$ is $8512$. This case
is listed in \cite{Pi} (right before Section 4) as an example
beyond the technique there.
\smallskip 

We omit the complete DAHA superpolynomials and the DAHA-Khovanov
polynomials in this and the previous examples (they are all
known).
Quite a few interesting relations for the (conjectural) 
Bette numbers can be seen from these and other examples 
of the DAHA-Bette polynomials. For instance, we expect 
that {\em $b_2$ of the Jacobian factors 
coincide with deg$_a$ of the corresponding DAHA 
superpolynomials\,} (though do not see algebraic reasons 
for this).  
```

### 4.8 Non-algebraic (pseudo-algebraic) (2,9) and (2,11) cables of trefoil — CD14 (6),(7)

`1408.4348/src/cable.tex:2885-3023`
```latex
\subsection{\bf Pseudo-algebraic knots}
This section is a theoretical challenge.
It indicates that some non-algebraic iterated
knots are very similar to the
algebraic ones. The positivity
of their DAHA-superpolynomials is the key 
``algebraic" feature.
Recall that we call such knots {\em pseudo-algebraic\,},
which conjecturally includes all algebraic knots
and their mirror images for any rectangle Young diagrams. 

The two uncolored knots posses such a positivity and all 
other DAHA properties of algebraic knots, for instance
the existence of the tilde- normalization from
Theorem \ref{JONITER} for any $A_N$ 
and the strict equality for deg$_a$ in (\ref{tilde-sup}).
We give below some other examples of pseudo-algebraic 
knots that do not have these two features; see
(\ref{nac7},\ref{nac-3}).
 
In spite of the presence of negative $\ss_i$ there, we think that
the coefficients of the DAHA-Bette polynomials in 
(\ref{bett6},\ref{bett7}) can be expected to be
the Bette numbers of certain Jacobian factors for the corresponding
curves (\ref{yxcurve}), however this is unknown at the moment.


$$
(\mathbf 6):\ \vec\rr=\{3 ,2\},\,
\vec\ss=\{2,-3\},\,
\t= C\!ab(9,2)T(3,2);\ 
\h_{\,\vec\rr,\,\vec\ss}\,(\square\,;\,\,q,t,a)=
$$
\renewcommand{\baselinestretch}{0.5} 
{\small
\(
1+a^3 q^6+q t+q^2 t+q^3 t+q^2 t^2+q^3 t^2+2 q^4 t^2
+q^3 t^3+q^4 t^3+q^5 t^3+q^4 t^4+q^5 t^4+q^5 t^5
+q^6 t^6 +a^2 \bigl(q^3+q^4+q^5+q^4 t+2 q^5 t
+q^6 t+q^5 t^2+q^6 t^2+q^6 t^3\bigr)
\)

\noindent
\(+a \bigl(q+q^2+q^3+q^2 t+2 q^3 t+3 q^4 t
+q^5 t+q^3 t^2+2 q^4 t^2+3 q^5 t^2+q^4 t^3
+2 q^5 t^3+q^6 t^3+q^5 t^4+q^6 t^4+q^6 t^5\bigr).
\)
}
\renewcommand{\baselinestretch}{1.2} 
\smallskip

{\sf Betti numbers} (presumably):
\begin{align} \label{bett6}
\h(q\!=\!1,t,a\!=\!0)=
1+3 t+4 t^2+3 t^3+2 t^4+t^5+t^6.
\end{align}

This sequence of coefficients and the one in the next example
do not come from algebraic knots. However they fit the table  
of Betti numbers from \cite{Pi} for the Puiseux
exponents $(4,6,v)$. The first two entries there
(our examples $(1),(2)$) are for $(4,6,7),(4,6,9)$;  
the knots from $(6),(7)$  {\em formally\,}
correspond to non-algebraic $(4,6,3),(4,6,5)$. 
\smallskip

{\sf Khovanov polynomial:}
$$(6):\ D\!A\!H\!A'K\!h(C\!ab(9,2)T(3,2))=$$
\renewcommand{\baselinestretch}{0.5} 
{\small
\(
 1 + q^4t^2 + q^6t^3 + q^6t^4 + q^{10}t^5 + q^8t^6 + q^{12}t^7 + 
  2q^{12}t^8 + 2q^{14}t^9 + q^{18}t^{11} + q^{20}t^{12}.
\)
}
\renewcommand{\baselinestretch}{1.2}
 
This polynomial, as well as the next $D\!A\!H\!A'K\!h$\~polynomial,
coincide with the corresponding (reduced) Khovanov polynomials.
This confirms the Connection Conjecture $(ii)$ and also demonstrates
that the differentials necessary to extract the Khovanov polynomials
from the stable ones satisfy the assumptions from Section 3.6 
in \cite{CJ}. 
\medskip

{\mathversion{bold}$\h_{\{3,2\},\{2,-1\}}$.}
This knot is pseudo-algebraic but non-algebraic as well. In a
sense, it is the greatest non-algebraic, since 
$C\!ab(13,2)T(3,2)$ from $(1)$ above is the first algebraic one
in this series.

$$
(\mathbf 7):\ \vec\rr=\{3, 2\},\,
\vec\ss=\{2,-1\},\,
\t= C\!ab(11,2)T(3,2);\ 
\h_{\,\vec\rr,\,\vec\ss}\,(\square\,;\,\,q,t,a)=
$$
\renewcommand{\baselinestretch}{0.5} 
{\small
\(
1+q t+q^2 t+q^3 t+q^2 t^2+q^3 t^2+2 q^4 t^2
+q^3 t^3+q^4 t^3+2 q^5 t^3+q^4 t^4+q^5 t^4
+q^6 t^4+q^5 t^5+q^6 t^5+q^6 t^6+q^7 t^7
+a^3 \bigl(q^6+q^7 t\bigr)
\)

\noindent
\(
+a^2 \bigl(q^3+q^4+q^5+q^4 t+2 q^5 t
+2 q^6 t+q^5 t^2+2 q^6 t^2+q^7 t^2
+q^6 t^3+q^7 t^3+q^7 t^4\bigr)
\)

\noindent
\(
+a \bigl(q+q^2+q^3+q^2 t+2 q^3 t+3 q^4 t
+q^5 t+q^3 t^2+2 q^4 t^2+4 q^5 t^2
+q^6 t^2+q^4 t^3+2 q^5 t^3+3 q^6 t^3
+q^5 t^4+2 q^6 t^4+q^7 t^4+q^6 t^5
+q^7 t^5+q^7 t^6\bigr).
\)
}
\renewcommand{\baselinestretch}{1.2} 
\smallskip

{\sf Betti numbers} (presumably):
\begin{align} \label{bett7}
\h(q\!=\!1,t,a\!=\!0)=
1+3 t+4 t^2+4 t^3+3 t^4+2 t^5+t^6+t^7.
\end{align}

The general formula is as follows:
\begin{align}\label{bet-gen-pos}
&\h_{\{3,2\},\{2,2n-15\}}(q\!=\!1,t,a\!=\!0)=
1+3 t+4\bigl(t^2\!+\!\ldots\!+\!t^{n-4}\bigr)\\
+ 3 &t^{n-3}\!+2 t^{n-2}+t^{n-1}+t^n 
\for C\!ab(2n\!-\!3,2)T(3,2),\ n\!\ge\! 6.\notag
\end{align}

```

### 4.9 Non-algebraic family $Cab(2m+1,2)T(3,2)$, $m\le 3$: (2,7),(2,5),**(2,3)**,**(2,1)**,(2,-1),(2,-3),(2,-5) — CD14 (a)–(g)

The (2,3)-cable is (c), the (2,1)-cable is (d).

`1408.4348/src/cable.tex:3064-3381`
```latex
\section{\sc Further aspects}
\subsection{\bf Non-algebraic knots}
Continuing the previous section, let
us discuss non-algebraic knots that 
are not {\em pseudo-algebraic\,}. 

The first knot we will consider is actually transitional.
It is {\em pseudo-algebraic\,}, but really
different from knots $({\bf 6,7})$
(see (\ref{bett6},\ref{bett7})).
For instance, the corresponding 
\,deg$_a$ is smaller than $(\ss_1\rr_2-1)=3$; 
the coincidence holds for  $({\bf 6,7})$
(and is conjectured for any algebraic ones). 
%The (formal) DAHA-Bette polynomials are
%``chaotic"  for the next knots in this family until 
%$C\!ab(-3,2)T(3,2)$.


$$
({\bf a})\ \vec\rr=\{3,-2\},\,
\vec\ss\!=\!\{2,5\},\,
\t\!=\! C\!ab(7,2)T(3,2);\, 
\h_{\vec\rr,\vec\ss}(\square;q,t,a)=
$$
\renewcommand{\baselinestretch}{0.5} 
{\small
\(
1+q t+q^2 t+q^3 t+q^2 t^2+q^3 t^2+q^4 t^2+q^3 t^3+q^4 t^3+q^4 t^4
+q^5 t^5
+a^2 (q^3+q^4+q^4 t+q^5 t+q^5 t^2)
+a (q+q^2+q^3+q^2 t +2 q^3 t+2 q^4 t
+q^3 t^2+2 q^4 t^2+q^5 t^2 +q^4 t^3+q^5 t^3+q^5 t^4).
\)
}
\renewcommand{\baselinestretch}{1.2} 
\smallskip
\begin{align}\label{nac7}
&\h(q\!=\!1,t,a\!=\!0)\!=\!
1+3 t+3 t^2+2 t^3+ t^4+t^5,\\
&D\!A\!H\!A'K\!h(C\!ab(7,2)T(3,2))\!=\!
1+q^4 t^2+q^6 t^3+q^6 t^4+q^{10} t^5\notag\\
&+q^8 t^6+q^{12} t^7+q^{12} t^8
+q^{14} t^9+q^{18} t^{11}+q^{20} t^{12}.
\notag
\end{align}

The first formula matches (\ref{bet-gen-pos}).
The DAHA-Khovanov polynomial here coincides with
the actual Khovanov one in this case (reduced and
under the tilde normalization); recall
that we use the standard parameters in the
formulas for DAHA-Khovanov polynomials.
\smallskip


$$
({\bf b})\ \vec\rr=\{3,-2\},\,
\vec\ss\!=\!\{2,7\},\,
\t\!=\! C\!ab(5,2)T(3,2);\, 
\h_{\vec\rr,\vec\ss}(\square;q,t,a)=
$$
\renewcommand{\baselinestretch}{0.5} 
{\small
\(
1-\frac{a^3 q^5}{t}+q t+q^2 t+q^2 t^2+q^3 t^2+q^3 t^3
+q^4 t^4+a^2 \left(q^3-q^5-\frac{q^4}{t}+q^4 t\right)\\
+a \left(q+q^2-q^4+q^2 t+2 q^3 t+q^3 t^2+q^4 t^2+q^4 t^3\right).
\)
}
\renewcommand{\baselinestretch}{1.2} 
\smallskip
\begin{align}\label{nac5}
&\h(q\!=\!1,t,a\!=\!0)\!=\!
1+2 t+2 t^2+t^3+t^4,\\
&D\!A\!H\!A'K\!h(C\!ab(5,2)T(3,2))\!=\!
1+q^4 t^2+q^6 t^3+q^6 t^4+q^{10} t^5\notag\\
&+q^{10} t^6+q^{12} t^7
-q^{12} t^9-q^{14} t^{10}-q^{18} t^{12}-q^{20} t^{13}.
\notag
\end{align}

The corresponding DAHA-Jones polynomial for $A_1$ is
\begin{align}\label{jon5-2}
&\tilde{J\!D}_{\{3,-2\},\,\{2,7\}}^{A_1}(\om_1\,;\,q,t)=
\h(q,t,a\!=\!-t^2) =\\
1+q t+q^2 t-&q t^2+q^3 t^2+q^4 t^2
-q^2 t^3-q^3 t^3-q^4 t^3-q^5 t^4+q^5 t^5.\notag
\end{align}

Its value at $t=q$ is the (reduced) {\em Jones 
polynomial\,} $\tilde{J}=1+q^2-q^7-q^9+q^{10}$ 
(for our $q$ and under the tilde normalization); the 
(reduced) one 
produced by the KnotTheory package
is $q^4+q^6-q^{11}-q^{13}+q^{14}$. 

We note that
(\ref{jon5-2}) has little to do with 
the polynomial obtained in \cite{Sam}, Section 5.2
(for the same knot).
We cannot comment on this discrepancy and the general approach 
used in \cite{Sam}. The latter stimulated our one, but the
construction there remains unclear to us. 
\smallskip

The knot $C\!ab(5,2)T(3,2)$ is the first with non-positive 
$\h_{\vec\rr,\vec\ss}(\square;q,t,a)$ in the
considered family. In this and 
further non-positive examples below the actual Khovanov 
polynomials can be obtained from DAHA-Khovanov by the 
substitution $-1\mapsto 1/t$, though the tilde-normalization
creates here some problems and $D\!A\!H\!A'K\!h(C\!ab(1,2)T(3,2))$ 
contains  the term $q^{10}(t^4 - t^6)$, which has no counterpart
in the corresponding $\tilde{K\!h}$, presumably due to the 
(topological) differential. Examples show that such a simple 
substitution is not generally sufficient.
%\smallskip

$$
({\bf c})\ \vec\rr=\{3,-2\},\,
\vec\ss\!=\!\{2,9\},\,
\t\!=\! C\!ab(3,2)T(3,2);\, 
\h_{\vec\rr,\vec\ss}(\square;q,t,a)=
$$
\renewcommand{\baselinestretch}{0.5} 
{\small
\(
1-q^2+q t+q^2 t
-q^3 t+q^2 t^2+q^3 t^3+a^3 \left(-\frac{q^4}{t^2}
-\frac{q^5}{t}\right)+
a^2 \bigl(q^3-q^4-q^5-\frac{q^3}{t^2}-\frac{q^3}{t}
-\frac{2 q^4}{t}\bigr)
+a \left(q+q^2-2 q^3-q^4-\frac{q^2}{t}
-\frac{q^3}{t}+q^2 t+q^3 t-q^4 t+q^3 t^2\right).
\)
}
\renewcommand{\baselinestretch}{1.2} 
\smallskip
\begin{align}\label{nac3}
&\h(q\!=\!1,t,a\!=\!0)\!=\!
t+t^2+t^3,\\
&D\!A\!H\!A'K\!h(C\!ab(3,2)T(3,2))\!=\!
1+q^4 t^2+q^6 t^3-q^4 t^4+q^6 t^4\notag\\
&-q^6 t^5+q^{10} t^5-q^{10} t^7
-q^{12} t^8-q^{12} t^9-q^{14} t^{10}-q^{18} t^{12}-q^{20} t^{13}.
\notag
\end{align}
%\smallskip
$$
({\bf d})\ \vec\rr=\{3,-2\},\,
\vec\ss\!=\!\{2,11\},\,
\t\!=\! C\!ab(1,2)T(3,2);\, 
\h_{\vec\rr,\vec\ss}(\square;q,t,a)=
$$
\renewcommand{\baselinestretch}{0.5} 
{\small
\(
-q+t-2 q^2 t+q t^2-q^3 t^2+q^2 t^3+a^3 \bigl(-q^5-\frac{q^3}{t^2}
-\frac{q^4}{t}\bigr)+a^2 \bigl(-2 q^3-2 q^4
-\frac{q^2}{t^2}-\frac{q^2}{t}-\frac{2 q^3}{t}-q^4 t-q^5 t\bigr)
+a \bigl(-3 q^2-q^3-\frac{q}{t}
-\frac{q^2}{t}+q t-3 q^3 t-q^4 t+q^2 t^2-q^4 t^2\bigr).
\)
}
\renewcommand{\baselinestretch}{1.2} 

%\smallskip
\begin{align}\label{nac1}
&\h(q\!=\!1,t,a\!=\!0)\!=\!
-1-t+t^3,\\
&D\!A\!H\!A'K\!h(C\!ab(1,2)T(3,2))\!=\!
q^2-q^2 t^2+q^6 t^2-q^4 t^3\notag\\
&+q^8 t^3-q^6 t^4
+q^{10} t^4-2 q^8 t^5+q^{12} t^5-q^{10} t^6-q^{12} t^7\notag\\
&-q^{14} t^8
-q^{14} t^9-q^{16} t^{10}-q^{20} t^{12}-q^{22} t^{13}.
\notag
\end{align}
%\smallskip
$$
({\bf e})\ \vec\rr=\{3,-2\},\,
\vec\ss\!=\!\{2,13\},\,
\t\!=\! C\!ab(-1,2)T(3,2);\, 
\h_{\vec\rr,\vec\ss}(\square;q,t,a)=
$$
\renewcommand{\baselinestretch}{0.5} 
{\small
\(
1+2 q t-t^2+q t^2+2 q^2 t^2-q t^3+q^3 t^3
+a^3 \left(q^4+\frac{q^2}{t^2}+\frac{q^3}{t}+q^5 t\right)\\
+a^2 \left(2 q^2+2 q^3+\frac{q}{t^2}+\frac{q}{t}
+\frac{2 q^2}{t}+q^2 t+2 q^3 t+2 q^4 t+q^4 t^2+q^5 t^2\right)\\
+a \left(3 q+q^2+\frac{1}{t}+\frac{q}{t}+q t+4 q^2 t
+q^3 t-q t^2+q^2 t^2+3 q^3 t^2+q^4 t^2+q^4 t^3\right).
\)
}
\renewcommand{\baselinestretch}{1.2} 
\smallskip

Note that this is the first case of nontrivial
$\h(q=0,t,a=0)=1+a/t$; it was $1$ (as for algebraic
knots) for $C\!ab(2\mm+1\ge 3,2)T(3,2)$ and then $t$ before 
this example.
For all further negative $\mm$, 
$\h(q=0,t,a=0)=1+a/t$  in the family $C\!ab(2\mm+1,2)T(3,2)$.


\begin{align}\label{nac-1}
&\h(q\!=\!1,t,a\!=\!0)\!=\!
1+2 t+2 t^2,\\
&D\!A\!H\!A'K\!h(C\!ab(-1,2)T(3,2))\!=\!
1 - q^4 + q^2 t + q^4 t^2 \notag\\
&- q^8 t^2+ q^6 t^3 - q^{10} t^3 + q^8 t^4 + 
2 q^{10} t^5 + q^{14} t^6 + q^{14} t^7 \notag\\
&+ q^{16} t^8 + q^{16} t^9 
+ q^{18} t^{10} + 
q^{22} t^{12} + q^{24} t^{13}.
\notag
\end{align}
%\smallskip
$$
({\bf f})\ \vec\rr=\{3,-2\},\,
\vec\ss\!=\!\{2,15\},\,
\t\!=\! C\!ab(-3,2)T(3,2);\, 
\h_{\vec\rr,\vec\ss}(\square;q,t,a)=
$$
\renewcommand{\baselinestretch}{0.5} 
{\small
\(
1+2 q t+q t^2+2 q^2 t^2-q t^3+q^2 t^3+2 q^3 t^3+q^4 t^4
+a^3 \left(q^4+\frac{q^2}{t^2}+\frac{q^3}{t}+q^5 t+q^6 t^2\right)
\\
+a^2 \bigl(2 q^2+2 q^3+\frac{q}{t^2}+\frac{q}{t}
+\frac{2 q^2}{t}+q^2 t+2 q^3 t+2 q^4 t+q^3 t^2
+2 q^4 t^2+2 q^5 t^2\\
+q^5 t^3+q^6 t^3\bigr)
+a \bigl(3 q+q^2+\frac{1}{t}+\frac{q}{t}+q t
+4 q^2 t
+q^3 t+2 q^2 t^2+4 q^3 t^2
+q^4 t^2\\+q^3 t^3+3 q^4 t^3+q^5 t^3+q^5 t^4\bigr).
\)
}
\renewcommand{\baselinestretch}{1.2} 

%\smallskip
\begin{align}\label{nac-3}
&\h(q\!=\!1,t,a\!=\!0)\!=\!
1+2 t+3 t^2+2 t^3+t^4,\\
&D\!A\!H\!A'K\!h(C\!ab(-3,2)T(3,2))\!=\!
1+q^2 t+q^4 t^2-q^8 t^2+q^6 t^3\notag\\
&+q^8 t^4+q^{10} t^5
+q^{12} t^6+q^{14} t^6+q^{12} t^7+2 q^{14} t^7
+q^{18} t^8\notag\\
&+q^{18} t^9+q^{20} t^{10}
+q^{20} t^{11}+q^{22} t^{12}+q^{26} t^{14}+q^{28} t^{15}.
\notag
\end{align}
%\smallskip
$$
({\bf g})\ \vec\rr=\{3,-2\},\,
\vec\ss\!=\!\{2,17\},\,
\t\!=\! C\!ab(-5,2)T(3,2);\, 
\h_{\vec\rr,\vec\ss}(\square;q,t,a)=
$$
\renewcommand{\baselinestretch}{0.5} 
{\small
\(
1+2 q t+q t^2+2 q^2 t^2+q^2 t^3+2 q^3 t^3+q^3 t^4
+2 q^4 t^4+q^5 t^5+a^3 \bigl(q^4+\frac{q^2}{t^2}
+\frac{q^3}{t}+q^5 t+q^6 t^2+q^7 t^3\bigr)
+a^2 \bigl(2 q^2+2 q^3+\frac{q}{t^2}
+\frac{q}{t}+\frac{2 q^2}{t}+q^2 t
+2 q^3 t+2 q^4 t+q^3 t^2+2 q^4 t^2
+2 q^5 t^2+q^4 t^3+2 q^5 t^3+2 q^6 t^3
+q^6 t^4+q^7 t^4\bigr)+a \bigl(3 q+q^2
+\frac{1}{t}+\frac{q}{t}+q t+4 q^2 t+q^3 t
+2 q^2 t^2+4 q^3 t^2+q^4 t^2+q^2 t^3
+2 q^3 t^3+4 q^4 t^3+q^5 t^3
+q^4 t^4+3 q^5 t^4+q^6 t^4+q^6 t^5\bigr).
\)
}
\renewcommand{\baselinestretch}{1.2} 

%\smallskip
\begin{align}\label{nac-5}
&\h(q\!=\!1,t,a\!=\!0)\!=\!
1 + 2 t + 3 t^2 + 3 t^3 + 3 t^4 + t^5,\\
&D\!A\!H\!A'K\!h(C\!ab(-5,2)T(3,2))\!=\!
1+q^2 t+q^4 t^2+q^6 t^3+q^8 t^4\notag\\
&+q^{10} t^5+q^{12} t^6
+q^{12} t^7+q^{14} t^7+q^{16} t^8
+q^{18} t^8+q^{16} t^9+2 q^{18} t^9
\notag\\
&+q^{22} t^{10}+q^{22} t^{11}+q^{24} t^{12}
+q^{24} t^{13}+q^{26} t^{14}+q^{30} t^{16}+q^{32} t^{17}.
\notag
\end{align}

One has that $\tilde{K\!h}-D\!A\!H\!A'K\!h=
q^{12}(t^5 - t^7) + q^{16}(t^7 - t^9),$ as can be
expected, though the first difference is 
obviously not of ``top" $q$\~degree.

This is a pseudo-algebraic knot, as well as all
further knots in the considered family, however 
numerically quite different from those in 
(\ref{bett6},\ref{bett7}). The polynomial 
$\h(q\!=\!1,t,a\!=\!0)$  becomes
$1+2 t+3 t^2+3 t^3+4 t^4+3 t^5+t^6$ for $C\!ab(-7,2)T(3,2)$
and 
\begin{align}\label{bet-gen-neg}
&\h_{\{3,-2\},\{2,2n+7\}}(q\!=\!1,t,a\!=\!0)=
1+2 t+3 t^2+3 t^3\\
+4\bigl(t^4\!+&\!\ldots\!+\!t^{n-2}\bigr)+
3 t^{n-1}\!+\!t^n \for C\!ab(-2n\!+\!5,2)T(3,2),\, n\!\ge\! 5.\notag
\end{align}
\medskip
```

### 4.10 (2,17)-cable of $T(4,3)$, $\vec r=\{4,-2\},\vec s=\{3,7\}$ — CD14 (h)

`1408.4348/src/cable.tex:3383-3472`
```latex
Let us consider a more involved
{\em pseudo-algebraic\,} knot:

$$
({\bf h})\ \vec\rr=\{4,-2\},\,
\vec\ss\!=\!\{3,7\},\,
\t\!=\! C\!ab(17,2)T(4,3);\, 
\h_{\vec\rr,\vec\ss}(\square;q,t,a)=
$$
\renewcommand{\baselinestretch}{0.5} 
{\small
\(
1+q t+q^2 t+q^3 t+q^4 t+q^5 t+q^2 t^2+q^3 t^2+2 q^4 t^2
+2 q^5 t^2+3 q^6 t^2+q^7 t^2+q^8 t^2+q^3 t^3+q^4 t^3+2 q^5 t^3
+3 q^6 t^3+4 q^7 t^3+3 q^8 t^3+2 q^9 t^3+q^4 t^4+q^5 t^4+2 q^6 t^4
+3 q^7 t^4+5 q^8 t^4+4 q^9 t^4+3 q^{10} t^4+q^5 t^5+q^6 t^5
+2 q^7 t^5+3 q^8 t^5+5 q^9 t^5+4 q^{10} t^5+2 q^{11} t^5+q^6 t^6
+q^7 t^6+2 q^8 t^6+3 q^9 t^6+5 q^{10} t^6+3 q^{11} t^6+q^{12} t^6
+q^7 t^7+q^8 t^7+2 q^9 t^7+3 q^{10} t^7+4 q^{11} t^7+q^{12} t^7
+q^8 t^8+q^9 t^8+2 q^{10} t^8+3 q^{11} t^8+3 q^{12} t^8+q^9 t^9
+q^{10} t^9+2 q^{11} t^9+2 q^{12} t^9+q^{13} t^9+q^{10} t^{10}
+q^{11} t^{10}+2 q^{12} t^{10}+q^{13} t^{10}+q^{11} t^{11}
+q^{12} t^{11}+q^{13} t^{11}+q^{12} t^{12}+q^{13} t^{12}
+q^{13} t^{13}+q^{14} t^{14}+a^4 \bigl(q^{10}+q^{11}+q^{12}
+q^{11} t+2 q^{12} t+q^{13} t+q^{12} t^2+2 q^{13} t^2+q^{14} t^2
+q^{13} t^3+q^{14} t^3+q^{14} t^4\bigr)+a^3 \bigl(q^6+q^7+2 q^8
+2 q^9+2 q^{10}+q^{11}+q^7 t+2 q^8 t+4 q^9 t+5 q^{10} t+5 q^{11} t
+2 q^{12} t+q^8 t^2+2 q^9 t^2+5 q^{10} t^2+7 q^{11} t^2
+7 q^{12} t^2+2 q^{13} t^2+q^9 t^3+2 q^{10} t^3+5 q^{11} t^3
+7 q^{12} t^3+5 q^{13} t^3+q^{14} t^3+q^{10} t^4+2 q^{11} t^4
+5 q^{12} t^4+5 q^{13} t^4+2 q^{14} t^4+q^{11} t^5+2 q^{12} t^5
+4 q^{13} t^5+2 q^{14} t^5+q^{12} t^6+2 q^{13} t^6+2 q^{14} t^6
+q^{13} t^7+q^{14} t^7+q^{14} t^8\bigr)+a^2 \bigl(q^3+q^4+2 q^5
+2 q^6+2 q^7+q^8+q^9+q^4 t+2 q^5 t+4 q^6 t+6 q^7 t+7 q^8 t+6 q^9 t
+4 q^{10} t+q^{11} t+q^5 t^2+2 q^6 t^2+5 q^7 t^2+8 q^8 t^2
+12 q^9 t^2+11 q^{10} t^2+7 q^{11} t^2+q^{12} t^2+q^6 t^3+2 q^7 t^3
+5 q^8 t^3+9 q^9 t^3+14 q^{10} t^3+14 q^{11} t^3+7 q^{12} t^3
+q^{13} t^3+q^7 t^4+2 q^8 t^4+5 q^9 t^4+9 q^{10} t^4+14 q^{11} t^4
+11 q^{12} t^4+4 q^{13} t^4+q^8 t^5+2 q^9 t^5+5 q^{10} t^5
+9 q^{11} t^5+12 q^{12} t^5+6 q^{13} t^5+q^{14} t^5+q^9 t^6
+2 q^{10} t^6+5 q^{11} t^6+8 q^{12} t^6+7 q^{13} t^6+q^{14} t^6
+q^{10} t^7+2 q^{11} t^7+5 q^{12} t^7+6 q^{13} t^7+2 q^{14} t^7
+q^{11} t^8+2 q^{12} t^8+4 q^{13} t^8+2 q^{14} t^8+q^{12} t^9
+2 q^{13} t^9+2 q^{14} t^9+q^{13} t^{10}+q^{14} t^{10}
+q^{14} t^{11}\bigr)+a \bigl(q+q^2+q^3+q^4+q^5+q^2 t+2 q^3 t
+3 q^4 t+4 q^5 t+5 q^6 t+3 q^7 t+2 q^8 t+q^9 t+q^3 t^2+2 q^4 t^2
+4 q^5 t^2+6 q^6 t^2+9 q^7 t^2+8 q^8 t^2+6 q^9 t^2+2 q^{10} t^2
+q^4 t^3+2 q^5 t^3+4 q^6 t^3+7 q^7 t^3+11 q^8 t^3+12 q^9 t^3
+10 q^{10} t^3+3 q^{11} t^3+q^5 t^4+2 q^6 t^4+4 q^7 t^4+7 q^8 t^4
+12 q^9 t^4+13 q^{10} t^4+10 q^{11} t^4+2 q^{12} t^4+q^6 t^5
+2 q^7 t^5+4 q^8 t^5+7 q^9 t^5+12 q^{10} t^5+12 q^{11} t^5
+6 q^{12} t^5+q^{13} t^5+q^7 t^6+2 q^8 t^6+4 q^9 t^6+7 q^{10} t^6
+11 q^{11} t^6+8 q^{12} t^6+2 q^{13} t^6+q^8 t^7+2 q^9 t^7
+4 q^{10} t^7+7 q^{11} t^7+9 q^{12} t^7+3 q^{13} t^7+q^9 t^8
+2 q^{10} t^8+4 q^{11} t^8+6 q^{12} t^8+5 q^{13} t^8+q^{10} t^9
+2 q^{11} t^9+4 q^{12} t^9+4 q^{13} t^9+q^{14} t^9+q^{11} t^{10}
+2 q^{12} t^{10}+3 q^{13} t^{10}+q^{14} t^{10}+q^{12} t^{11}
+2 q^{13} t^{11}
+q^{14} t^{11}+q^{13} t^{12}+q^{14} t^{12}+q^{14} t^{13}\bigr).
\)
}
\renewcommand{\baselinestretch}{1.2} 
\smallskip

The corresponding $D\!A\!H\!A'K\!h(C\!ab(17,2)T(4,3))$
equals
\renewcommand{\baselinestretch}{0.5} 
{\small
\(
1+q^4 t^2+q^6 t^3+q^6 t^4+q^{10} t^5+q^8 t^6
+q^{12} t^7+q^{10} t^8+q^{12} t^8+2 q^{14} t^9+q^{12} t^{10}
+q^{14} t^{10}+2 q^{16} t^{11}+q^{18} t^{11}+3 q^{16} t^{12}
+q^{20} t^{12}+2 q^{18} t^{13}+2 q^{20} t^{13}+q^{18} t^{14}
+q^{20} t^{14}+q^{22} t^{14}+4 q^{22} t^{15}+q^{20} t^{16}
+q^{22} t^{16}+q^{24} t^{16}+q^{26} t^{16}+2 q^{24} t^{17}
+q^{26} t^{17}+q^{24} t^{18}+q^{28} t^{18}+q^{26} t^{19}
+q^{28} t^{19}+q^{28} t^{20}+q^{30} t^{20}
+2 q^{30} t^{21}+q^{32} t^{22}+q^{36} t^{24}+q^{34} t^{25}.
\)
}
\renewcommand{\baselinestretch}{1.2} 
%\smallskip

One has:
$\tilde{K\!h}-D\!A\!H\!A'K\!h=q^{34}(t^{23} - t^{25})$.
The corresponding (reduced and tilde-normalized) 
HOMFLYPT polynomial coincides with   
$\h(\square;q,q,-a)$ as well as in all examples
we considered.

```

(CD14:3474-3531 contains a further example "(i)" inside `\comment{}` labelled $\vec r=\{4,-2\},\vec s=\{3,7\}$, $Cab(7,2)T(4,3)$; the
label is inconsistent ($\{4,-2\},\{3,7\}$ is (h)); not transcribed.)

### 4.11 A_1 DAHA-Jones polynomial of $Cab(5,2)T(3,2)$ — CD14 [checked: $=\mathcal H_{(b)}(a=-t^2)$]

`1408.4348/src/cable.tex:3145-3160`
```latex

The corresponding DAHA-Jones polynomial for $A_1$ is
\begin{align}\label{jon5-2}
&\tilde{J\!D}_{\{3,-2\},\,\{2,7\}}^{A_1}(\om_1\,;\,q,t)=
\h(q,t,a\!=\!-t^2) =\\
1+q t+q^2 t-&q t^2+q^3 t^2+q^4 t^2
-q^2 t^3-q^3 t^3-q^4 t^3-q^5 t^4+q^5 t^5.\notag
\end{align}

Its value at $t=q$ is the (reduced) {\em Jones 
polynomial\,} $\tilde{J}=1+q^2-q^7-q^9+q^{10}$ 
(for our $q$ and under the tilde normalization); the 
(reduced) one 
produced by the KnotTheory package
is $q^4+q^6-q^{11}-q^{13}+q^{14}$. 

```

### 4.12 C24 cable examples (non-algebraic, related to Gorenstein non-planar rings)

$Cab(5,2)T(3,2)$ (identical to CD14 (b) [checked]) and $Cab(7,2)T(3,2)$ at $a=-t/q$:

`2410.14703/src/superpolynoms.tex:7830-7854`
```latex
Let $\r=\F[[z^4,z^5,z^6]]\subset \o=\F[[z]]$, which is a
Gorenstein ring with $\de=4, \cc=z^8\o$. Then $L(q,t,a=-1/q)$ coincides with
(uncolored) $\h^{daha}_K(q, t, a=-t/q)$ upon the  substitution $q\mapsto q t$
for the cable $K=C\!ab(5,2)T(3,2)$. Recall that the smallest algebraic
knot in this family is $C\!ab(13,2)T(3,2)$. The general formula is due
to \cite{ChD1} and Proposition 4.2 in \cite{ChW}:\, $\h_K=$
\renewcommand{\baselinestretch}{0.5} 
{\small
\(
-\frac{a^3 q^5}{t}+a^2 (-q^5+q^4 t-\frac{q^4}{t}+q^3)+
a (q^4 t^3+q^4 t^2-q^4+q^3 t^2+2 q^3 t+q^2 t+q^2+q)+
q^4 t^4+q^3 t^3+q^3 t^2+q^2 t^2+q^2 t+q t+1.
\)
}
\renewcommand{\baselinestretch}{1.2} 

One has: $\h^{daha}_K(q,t,a=-t/q)=
1-t+q^3 t-q^3 t^4+q^4 t^4$, which coincides with
$L(q/t,t,a=-t/q)$ and with
$\h^{mot}(q,t,a=-t/q)$. 

The next cable in this family is $K'=C\!ab(7,2)T(3,2)$, which
corresponds to $\r'=\F_q[[z^4,z^6,z^7]]$. One has:
$\h^{daha}_{K'}(q,t)=1-t+q^3 t-q^3 t^2+q^4 t^2-q^4 t^5+q^5 t^5=
L(q/t,t)=\h^{mot}(q,t)$ for $a=-t/q$,
```

$Cab(9,2)T(3,2)$ (identical to CD14 (6) [checked]) and $Cab(11,2)T(4,3)$ (full $\mathcal H$), $Cab(7,2)T(4,3)$ at $a=-t/q$:

`2410.14703/src/superpolynoms.tex:8051-8186`
```latex
{\sf Non-algebraic cables.}
The relation between Gorenstein $\r$ 
and non-algebraic $K$ seems  of sufficiently
general nature, but we mostly deal with examples by now. 
For instance, let
$\r=\F_q[[z^{\upsilon \rr}, z^{\upsilon \ss}, z^{\tt}]]$, where 
$\tt=\upsilon\rr\ss+p$ for $p\in \Z$, 
$\rr,\ss,\tt, \upsilon>0$, gcd$(\rr,\ss)=1$, gcd$(|p|,\nu)=1$,
and the ring is assumed Gorenstein. Then $\Ga=\lan \rr,\ss,\tt\ran$,
$G=\Z_+\setminus \Ga$; compare with (\ref{Ga-cab}).
The corresponding
cable will be $K=C\!ab(\tt,\upsilon)T(\rr,\ss)$. We will use
$q$ from $\h$ in this section.

Given $\rr,\ss$, we expect that there exists $p_0\ge 0$ such that
 for $p\ge p_0$, 
\begin{align}\label{conjKG}
\h^{daha}_K(q,t,a\!=\!-\frac{t}{q})\!=\!L(q/t,t,a\!=\!-\frac{t}{q})\!=\!
(1\!-\!q)(1\!-\!t)\varrho(q,t)+ q^\de
\end{align}
for $\varrho(q,t)\!=\!
\!\sum_{G\ni x>y\in \Ga}q^{g(x)}t^{v(y)-1}=(q t)^{\de-1}
\varrho(t^{-1},q^{-1})$
from Section \ref{sec:rhocab}, where 
$v(x)=|\{\nu\!\in\! \Ga \mid 0\!\le\! \nu\!\le\! x\}|$ and 
$g(x)=|\{g\!\in\! G \mid 0\!\le\! g\! <\! x\}|$. 
%Generally, this can be up to some factor $t^\bullet q^\bullet$.
Let us provide a couple of examples; we mostly follow \cite{ChW, ChS}.
\vskip 0.2cm

For $\rr=3,\ss=2,\upsilon=2,p=-4$, $K=C\!ab(9,2)T(3,2)$, and
$\h_K^{daha}=$
{\small $1+a^3 q^6+q t+q^2 t+q^3 t+q^2 t^2
+q^3 t^2+2 q^4 t^2+q^3 t^3+q^4 t^3+q^5 t^3
+q^4 t^4+q^5 t^4+q^5 t^5+q^6 t^6
+a^2 (q^3+q^4+q^5+q^4 t+2 q^5 t+q^6 t+q^5 t^2+q^6 t^2+q^6 t^3)
+a (q+q^2+q^3+q^2 t+2 q^3 t+3 q^4 t+q^5 t+q^3 t^2+2 q^4 t^2
+3 q^5 t^2+q^4 t^3+2 q^5 t^3+q^6 t^3+q^5 t^4+q^6 t^4+q^6 t^5).$
}

Then $\varrho=$ {\small $1+q+q^2+q^3+q^4+q^5+q^3 t+q^4 t+q^5 t+q^4 t^2
+q^5 t^2+q^5 t^3+q^5 t^4+q^5 t^5$} in this case, and
 $\h_K^{daha}(q,t,a=-t/q)=$
{\small $1-t+q^3 t-q^3 t^2+q^4 t^2-q^4 t^3+q^5 t^3-q^5 t^6+q^6 t^6$},
which 
coincides with $L(q/t,t,a=-t/q)$. 
\vskip 0.2cm

For $\rr=4,\ss=3,\upsilon=2,p=-13$, the ring $\r=\F_q[[z^6,z^8,z^{11}]]$
is still Gorenstein, $K=C\!ab(11,2)T(4,3)$, and the coincidence holds.
One has:
$\varrho=$ {\small 
$1+q+q^2+q^3+q^4+q^5+q^6+q^7+q^8
+q^9+q^10+q^5 t+q^6 t+q^7 t+q^8 t+q^9 t+q^10 t+q^6 t^2+q^7 t^2+q^8 t^2
+q^9 t^2+q^{10} t^2+q^8 t^3+q^9 t^3+q^{10} t^3+q^8 t^4+q^9 t^4
+q^{10} t^4
+q^9 t^5+q^{10} t^5+q^{10} t^6+q^{10} t^7+q^{10} t^8+q^{10} t^9
+q^{10} t^{10},$
}
and $L(q/t,t,a\!=\!-t/q)=$
{\small $1-t+q^5 t-q^5 t^2+q^6 t^2-q^6 t^3+q^8 t^3-q^8 t^5+q^9 t^5
-q^9 t^6+q^{10} t^6-q^10 t^11+q^{11} t^{11}$.}
\vfil

%\comment{
Let us provide the DAHA superpolynomial for the  $K$ above:\, 
$\h_K^{daha}=$
\renewcommand{\baselinestretch}{1.} 
{\footnotesize
\( 
t^{11} q^{11}-t^5 q^{11}+t^{10} q^{10}+t^9 q^{10}+t^8 q^{10}+
t^7 q^{10}+t^6 q^{10}-2 t^4 q^{10}-t^3 q^{10}+t^9 q^9+t^8 q^9
+2 t^7 q^9+2 t^6 q^9+3 t^5 q^9-3 t^3 q^9-t^2 q^9+t^8 q^8+t^7 q^8
+2 t^6 q^8+3 t^5 q^8+4 t^4 q^8+t^3 q^8-3 t^2 q^8-t q^8+t^7 q^7
+t^6 q^7+2 t^5 q^7+3 t^4 q^7+4 t^3 q^7-2 t q^7+t^6 q^6+t^5 q^6
+2 t^4 q^6+3 t^3 q^6+3 t^2 q^6-q^6+t^5 q^5+t^4 q^5+2 t^3 q^5
+2 t^2 q^5+t q^5+t^4 q^4+t^3 q^4+2 t^2 q^4+t q^4+t^3 q^3+t^2 q^3
+t q^3+t^2 q^2+t q^2+t q+a^5 \bigl(-\frac{q^{14}}{t}
-\frac{q^{13}}{t^2}-\frac{q^{12}}{t^3}\bigr)+a^4 \bigl(-t q^{14}
-q^{14}-t q^{13}-\frac{2 q^{13}}{t}-2 q^{13}-\frac{3 q^{12}}{t}
-\frac{2 q^{12}}{t^2}-2 q^{12}+t q^{11}-\frac{2 q^{11}}{t}
-\frac{2 q^{11}}{t^2}-\frac{q^{11}}{t^3}-\frac{q^{10}}{t^2}
-\frac{q^{10}}{t^3}+q^{10}\bigr)+a^3 \bigl(-t^2 q^{14}-t^3 q^{13}
-2 t^2 q^{13}-3 t q^{13}-q^{13}-3 t^2 q^{12}-5 t q^{12}
-\frac{q^{12}}{t}-5 q^{12}+t^5 q^{11}+t^4 q^{11}+2 t^3 q^{11}
+t^2 q^{11}-5 t q^{11}-\frac{5 q^{11}}{t}-\frac{q^{11}}{t^2}
-7 q^{11}+t^4 q^{10}+2 t^3 q^{10}+4 t^2 q^{10}+2 t q^{10}
-\frac{5 q^{10}}{t}-\frac{3 q^{10}}{t^2}-5 q^{10}+t^3 q^9+2 t^2 q^9
+4 t q^9-\frac{3 q^9}{t}-\frac{2 q^9}{t^2}-\frac{q^9}{t^3}+q^9
+t^2 q^8+2 t q^8-\frac{q^8}{t^2}+2 q^8+t q^7+q^7+q^6\bigr)
+a^2 \bigl(-t^4 q^{13}-t^3 q^{13}-t^2 q^{13}-t^4 q^{12}-4 t^3 q^{12}
-4 t^2 q^{12}-2 t q^{12}+t^8 q^{11}+t^7 q^{11}+2 t^6 q^{11}
+2 t^5 q^{11}+2 t^4 q^{11}-2 t^3 q^{11}-9 t^2 q^{11}-7 t q^{11}
-3 q^{11}+t^7 q^{10}+2 t^6 q^{10}+4 t^5 q^{10}+6 t^4 q^{10}
+6 t^3 q^{10}-2 t^2 q^{10}-11 t q^{10}-\frac{2 q^{10}}{t}-7 q^{10}
+t^6 q^9+2 t^5 q^9+5 t^4 q^9+8 t^3 q^9+9 t^2 q^9-2 t q^9
-\frac{4 q^9}{t}-\frac{q^9}{t^2}-9 q^9+t^5 q^8+2 t^4 q^8+5 t^3 q^8
+8 t^2 q^8+6 t q^8-\frac{4 q^8}{t}-\frac{q^8}{t^2}-2 q^8+t^4 q^7
+2 t^3 q^7+5 t^2 q^7+6 t q^7-\frac{q^7}{t}-\frac{q^7}{t^2}+2 q^7
+t^3 q^6+2 t^2 q^6+4 t q^6+2 q^6+t^2 q^5+2 t q^5+2 q^5+t q^4+q^4
+q^3\bigr)+a \bigl(-t^5 q^{12}-t^4 q^{12}-t^3 q^{12}+t^{10} q^{11}
+t^9 q^{11}+t^8 q^{11}+t^7 q^{11}+t^6 q^{11}-3 t^4 q^{11}-4 t^3 q^{11}
-2 t^2 q^{11}+t^9 q^{10}+2 t^8 q^{10}+3 t^7 q^{10}+4 t^6 q^{10}
+5 t^5 q^{10}+2 t^4 q^{10}-5 t^3 q^{10}-7 t^2 q^{10}-3 t q^{10}
+t^8 q^9+2 t^7 q^9+4 t^6 q^9+6 t^5 q^9+9 t^4 q^9+5 t^3 q^9-6 t^2 q^9
-7 t q^9-2 q^9+t^7 q^8+2 t^6 q^8+4 t^5 q^8+7 t^4 q^8+10 t^3 q^8
+5 t^2 q^8-5 t q^8-\frac{q^8}{t}-4 q^8+t^6 q^7+2 t^5 q^7+4 t^4 q^7
+7 t^3 q^7+9 t^2 q^7+2 t q^7-\frac{q^7}{t}-3 q^7+t^5 q^6+2 t^4 q^6
+4 t^3 q^6+6 t^2 q^6+5 t q^6-\frac{q^6}{t}+t^4 q^5+2 t^3 q^5
+4 t^2 q^5+4 t q^5+q^5+t^3 q^4+2 t^2 q^4+3 t q^4+q^4+t^2 q^3+2 t q^3
+q^3+t q^2+q^2+q\bigr)+1.
\)
}
\renewcommand{\baselinestretch}{1.2}

%}

\vskip 0.2cm
However, there is no relation between $\h_K^{daha}(q,t,a)$
for (non-algebraic) $K=C\!ab(7,2)T(4,3)$ 
and $L(q/t,t,a)$ evaluated at $a=-t/q$  for
the corresponding Gorenstein 
ring, which is $\r=\F_q[[z^6,z^7,z^{8}]]$
for $\rr=4,\ss=3,\upsilon=2, p=-17$. They are 
{\small $1-t+q^5 t+q^6 t-q^7 t-q^5 t^2+q^6 t^2+q^7 t^2-q^8 t^2-
2 q^6 t^3+q^7 t^3+q^8 t^3-q^7 t^4+q^8 t^4-q^8 t^9+q^9 t^9$},
and {\small $1 - t + q^5 t - q^5 t^4 + q^8 t^4 - q^8 t^9 + 
q^9 t^9$} correspondingly. Both satisfy the superduality
$(q t)^\de \h(1/t,1/q,a)=\h(q,t,a)$.
The consideration
of sufficiently 
general Gorenstein rings (if they are non-planar) is needed for 
such a coincidence. The same reservation
may be needed to ensure that  
$\h^{daha}(q,t,a)=L(qt,t,a)$ at $a=-t/q$, which is not clear by now.

```

(Typo in the source: `q^10` instead of `q^{10}` inside $\varrho$ (C24:8104, 8112); the full $\mathcal H$ at C24:8121-8161 parses fine.)

C24 iteration formula relating $Cab(2rs+1,2)T(r,s)$, $T(2r,2s)$ and $T(r,s)$ colored by $2\omega_1$ (motivic; DAHA version in [ChW]):

`2410.14703/src/superpolynoms.tex:5903-5927`
```latex
\begin{proposition}\label{thm:h-iter}
{\sf (a)} Let $\h^{mot}_{\rr,\ss;2}$ be the uncolored motivic 
superpolynomial for 
$\r'=\F_q[[x'\!=\!z^{2\rr}, y'\!=\!z^{2\ss}\! +
u\,\!z^{2\ss+1}]]$, where $u\neq 0$, 
gcd$(\rr,\ss)=1$ and the knot is $C\!ab(2\rr\ss+1,2)
T(\rr,\ss)$. {\sf (b)} Let $\h^{mot}_{2\rr,2\ss}$ be that for 
$\tilde{\r}=\F_q[[\tilde{x}=uz_1^\rr+z_2^\rr,\,
\tilde{y}=z_1^\ss+z_2^\ss]]
\subset \tilde{\o}=\F_q[[e_1,e_2,z_1,z_2]]$ as above, where 
$u^{\ss}\neq 0,1$. 
This is for $T(2\rr,2\ss)$.
{\sf (c)} Let $\h^{mot}_{\rr,\ss}(2)$ be
the motivic superpolynomial 
for $T(\rr,\ss)$ colored
by $2\om_1$, i.e. for 
 $\r_\circ=\F_q[[x=z^\rr,y=z^\ss]]$ acting diagonally in
$\Om=\F_q[[z]]\vep_1\oplus \F_q[[z]]\vep_2$. The following 
counterpart of the corresponding $\h^{daha}$-identity (for $m=0$) 
from \cite{ChW} holds:
\begin{align}\label{h-dtor}
(1+a q)\h^{mot}_{\rr,\ss}(2)+(q-1)\h^{mot}_{\rr,\ss;2}=
q\h^{mot}_{2\rr,2\ss}.
\end{align}
\end{proposition}
```

`2410.14703/src/superpolynoms.tex:6004-6022`
```latex
The simplest version of this formula is the following 
identity
$$(q-1)\Bigl(\h^{mot}_{1,1;2}-(1+a q)\h^{mot}_{1,1}(2)\Bigr)=
q\Bigl(\h^{mot}_{2,2}-(1+a q)\h^{mot}_{1,1}(2)\Bigr).
$$
Here the iterated knot $C\!ab(3,2)T(1,1)$ is isotopic to
$T(3,2)$. Also, $\h^{mot}_{1,1}(2)=1$ since $T(1,1)$ is unknot, 
The identity becomes
$$(q-1)\Bigl(\h^{mot}_{3,2}-(1+aq)\Bigr)=
q\Bigl(\h^{mot}_{2,2}-(1+aq)\Bigr).$$

It clarifies very well the similarity of
$\h^{mot}_{3,2}=qt+(1+aq)$
and $\h^{mot}_{2,2}=(q-1)t+(1+aq)$, not immediate in
the DAHA approach. Recall that  $\h^{daha}_{1,1;2}=\h^{daha}_{3,2}$ 
corresponding to $C\!ab(3,2)T(1,1)=T(3,2)$ follows from the DAHA
identity 
$\{\tau_-\tau_+\tau_-(P_\la)\}_{ev}=
\{\tau_-\bigl(\tau_+\tau_-(P_\la)(1)\bigr)\}_{ev}$.
```

### 4.13 CP16 cable examples

$Cab(15,2)T(3,2)$ (same as CD14 (2)):

`1605.00978/src/jacobian-fin-full.tex:3131-3175`
```latex
For $\r=\C[[ t^4, t^6+t^9]]$ corresponding to the
$\Ga=\lan 4,6,15\ran$ and
 cable $C\!ab(15,2)T(3,2)$ the situation is very similar. We
checked that:
\begin{align*}
&\sum_{m=0}^\infty \,
\sum_{\{D_0=D,\ldots,D_m\}} q^{|D|+m\,} 
t^{\de-\hbox{\tiny dim}
(J^m_{\r}[D])\,} a^m=
\h_{\,\{3,2\},\{2,3\}}(\square\,;\,\,q,t,a)=
\end{align*}

\renewcommand{\baselinestretch}{0.5} 
{\small
\(
1 + q t + q^9 t^9 + q^2 \bigl(t + t^2\bigr) + q^3 \bigl(t + t^2 
+ t^3\bigr) + 
 a^3 \bigl(q^6 + q^7 t + q^8 t^2 + q^9 t^3\bigr) 
+ q^4 \bigl(2 t^2 + t^3 + t^4\bigr) + 
 q^5 \bigl(2 t^3 + t^4 + t^5\bigr) + q^6 \bigl(2 t^4 + t^5 
+ t^6\bigr) + 
 q^7 \bigl(2 t^5 + t^6 + t^7\bigr) + q^8 \bigl(t^6 + t^7 
+ t^8\bigr) + 
 a^2 \bigl(q^3 + q^4 \bigl(1 + t\bigr) + q^5 \bigl(1 + 2 t 
+ t^2\bigr) + 
    q^6 \bigl(2 t + 2 t^2 + t^3\bigr) + q^7 \bigl(2 t^2 + 2 t^3 
+ t^4\bigr) + 
    q^8 \bigl(2 t^3 + 2 t^4 + t^5\bigr) + q^9 \bigl(t^4 + t^5 
+ t^6\bigr)\bigr) + 
 a \bigl(q + q^2 \bigl(1 + t\bigr) + q^3 \bigl(1 + 2 t + t^2\bigr)
 + q^4 \bigl(3 t + 2 t^2 + t^3\bigr) + 
    q^5 \bigl(t + 4 t^2 + 2 t^3 + t^4\bigr) + q^6 \bigl(t^2 
+ 4 t^3 + 2 t^4 + t^5\bigr) + 
    q^7 \bigl(t^3 + 4 t^4 + 2 t^5 + t^6\bigr) + 
    q^8 \bigl(t^4 + 3 t^5 + 2 t^6 + t^7\bigr) + q^9 \bigl(t^6 
+ t^7 + t^8\bigr)\bigr).
\)
}
\renewcommand{\baselinestretch}{1.2} 

\smallskip
\noindent
See the same section of \cite{ChD1}. Note that formulas
(3.1) and (3.2) there for the DAHA-Betti polynomials are
obtained from the superpolynomials as $a=0,q=1$. 
```

$Cab(31,2)T(7,2)$, $\vec r=\{7,2\},\vec s=\{2,3\}$:

`1605.00978/src/jacobian-fin-full.tex:3206-3320`
```latex
Compare them with the coefficients of $a^3$ of the corresponding
superpolynomial from \cite{ChD1}, which we provide:

\vskip -0.5cm
$$\r=\C[[ z^4,z^{14}+z^{17}]],\ \t=C\!ab(31,2)T(7,2),\ 
\h_{\{7,2\},\{2,3\}}(\square;q,t,a)=$$
%\smallskip

\renewcommand{\baselinestretch}{0.5} 
{\small
\(
1+q t+q^2 t+q^3 t+q^2 t^2+q^3 t^2+2 q^4 t^2+q^5 t^2+q^6 t^2+q^3 t^3
+q^4 t^3+2 q^5 t^3+2 q^6 t^3+2 q^7 t^3+q^8 t^3
+q^9 t^3+q^4 t^4+q^5 t^4+2 q^6 t^4+2 q^7 t^4+3 q^8 t^4+2 q^9 t^4
+2 q^{10} t^4+q^5 t^5+q^6 t^5+2 q^7 t^5+2 q^8 t^5+3 q^9 t^5
+3 q^{10} t^5+3 q^{11} t^5+q^6 t^6+q^7 t^6+2 q^8 t^6+2 q^9 t^6
+3 q^{10} t^6+3 q^{11} t^6+4 q^{12} t^6+q^7 t^7+q^8 t^7+2 q^9 t^7
+2 q^{10} t^7+3 q^{11} t^7+3 q^{12} t^7+4 q^{13} t^7+q^8 t^8
+q^9 t^8+2 q^{10} t^8+2 q^{11} t^8+3 q^{12} t^8+3 q^{13} t^8
+4 q^{14} t^8+q^9 t^9+q^{10} t^9+2 q^{11} t^9+2 q^{12} t^9
+3 q^{13} t^9+3 q^{14} t^9+4 q^{15} t^9+q^{10} t^{10}+q^{11} t^{10}
+2 q^{12} t^{10}+2 q^{13} t^{10}+3 q^{14} t^{10}+3 q^{15} t^{10}
+3 q^{16} t^{10}+q^{11} t^{11}+q^{12} t^{11}+2 q^{13} t^{11}
+2 q^{14} t^{11}+3 q^{15} t^{11}+3 q^{16} t^{11}+2 q^{17} t^{11}
+q^{12} t^{12}+q^{13} t^{12}+2 q^{14} t^{12}+2 q^{15} t^{12}
+3 q^{16} t^{12}+2 q^{17} t^{12}+q^{18} t^{12}+q^{13} t^{13}
+q^{14} t^{13}+2 q^{15} t^{13}+2 q^{16} t^{13}+3 q^{17} t^{13}
+q^{18} t^{13}+q^{14} t^{14}+q^{15} t^{14}+2 q^{16} t^{14}
+2 q^{17} t^{14}+2 q^{18} t^{14}+q^{15} t^{15}+q^{16} t^{15}
+2 q^{17} t^{15}+2 q^{18} t^{15}+q^{19} t^{15}+q^{16} t^{16}
+q^{17} t^{16}+2 q^{18} t^{16}+q^{19} t^{16}+q^{17} t^{17}
+q^{18} t^{17}+2 q^{19} t^{17}+q^{18} t^{18}+q^{19} t^{18}
+q^{20} t^{18}+q^{19} t^{19}+q^{20} t^{19}+q^{20} t^{20}
+q^{21} t^{21}
\)

%\smallskip
\noindent
\(
+a^3 \bigl(q^6+q^7 t+q^8 t+q^9 t+q^8 t^2+q^9 t^2
+2 q^{10} t^2+q^{11} t^2+q^{12} t^2+q^9 t^3+q^{10} t^3+2 q^{11} t^3
+2 q^{12} t^3+2 q^{13} t^3+q^{10} t^4+q^{11} t^4+2 q^{12} t^4
+2 q^{13} t^4+3 q^{14} t^4+q^{11} t^5+q^{12} t^5+2 q^{13} t^5
+2 q^{14} t^5+3 q^{15} t^5+q^{12} t^6+q^{13} t^6+2 q^{14} t^6
+2 q^{15} t^6+3 q^{16} t^6+q^{13} t^7+q^{14} t^7+2 q^{15} t^7
+2 q^{16} t^7+3 q^{17} t^7+q^{14} t^8+q^{15} t^8+2 q^{16} t^8
+2 q^{17} t^8+2 q^{18} t^8+q^{15} t^9+q^{16} t^9+2 q^{17} t^9
+2 q^{18} t^9+q^{19} t^9+q^{16} t^{10}+q^{17} t^{10}
+2 q^{18} t^{10}+q^{19} t^{10}+q^{17} t^{11}+q^{18} t^{11}
+2 q^{19} t^{11}+q^{18} t^{12}+q^{19} t^{12}+q^{20} t^{12}
+q^{19} t^{13}+q^{20} t^{13}+q^{20} t^{14}+q^{21} t^{15}\bigr)
\)

%\smallskip
\noindent
\(
+a^2 \bigl(q^3+q^4+q^5+q^4 t+2 q^5 t+3 q^6 t
+2 q^7 t+q^8 t+q^5 t^2+2 q^6 t^2+4 q^7 t^2+4 q^8 t^2+4 q^9 t^2
+2 q^{10} t^2+q^{11} t^2+q^6 t^3+2 q^7 t^3+4 q^8 t^3+5 q^9 t^3
+6 q^{10} t^3+5 q^{11} t^3+3 q^{12} t^3+q^7 t^4+2 q^8 t^4
+4 q^9 t^4+5 q^{10} t^4+7 q^{11} t^4+7 q^{12} t^4+5 q^{13} t^4
+q^8 t^5+2 q^9 t^5+4 q^{10} t^5+5 q^{11} t^5+7 q^{12} t^5
+8 q^{13} t^5+6 q^{14} t^5+q^9 t^6+2 q^{10} t^6+4 q^{11} t^6
+5 q^{12} t^6+7 q^{13} t^6+8 q^{14} t^6+6 q^{15} t^6+q^{10} t^7
+2 q^{11} t^7+4 q^{12} t^7+5 q^{13} t^7+7 q^{14} t^7+8 q^{15} t^7
+6 q^{16} t^7+q^{11} t^8+2 q^{12} t^8+4 q^{13} t^8+5 q^{14} t^8
+7 q^{15} t^8+8 q^{16} t^8+5 q^{17} t^8+q^{12} t^9+2 q^{13} t^9
+4 q^{14} t^9+5 q^{15} t^9+7 q^{16} t^9+7 q^{17} t^9+3 q^{18} t^9
+q^{13} t^{10}+2 q^{14} t^{10}+4 q^{15} t^{10}+5 q^{16} t^{10}
+7 q^{17} t^{10}+5 q^{18} t^{10}+q^{19} t^{10}+q^{14} t^{11}
+2 q^{15} t^{11}+4 q^{16} t^{11}+5 q^{17} t^{11}+6 q^{18} t^{11}
+2 q^{19} t^{11}+q^{15} t^{12}+2 q^{16} t^{12}+4 q^{17} t^{12}
+5 q^{18} t^{12}+4 q^{19} t^{12}+q^{16} t^{13}+2 q^{17} t^{13}
+4 q^{18} t^{13}+4 q^{19} t^{13}+q^{20} t^{13}+q^{17} t^{14}
+2 q^{18} t^{14}+4 q^{19} t^{14}+2 q^{20} t^{14}+q^{18} t^{15}
+2 q^{19} t^{15}+3 q^{20} t^{15}+q^{19} t^{16}+2 q^{20} t^{16}
+q^{21} t^{16}+q^{20} t^{17}+q^{21} t^{17}+q^{21} t^{18}\bigr)
\)

%\smallskip
\noindent
\(
+a \bigl(q+q^2+q^3+q^2 t+2 q^3 t+3 q^4 t+2 q^5 t+q^6 t+q^3 t^2
+2 q^4 t^2+4 q^5 t^2+4 q^6 t^2+4 q^7 t^2+2 q^8 t^2+q^9 t^2+q^4 t^3
+2 q^5 t^3+4 q^6 t^3+5 q^7 t^3+6 q^8 t^3+5 q^9 t^3+4 q^{10} t^3
+q^{11} t^3+q^5 t^4+2 q^6 t^4+4 q^7 t^4+5 q^8 t^4+7 q^9 t^4
+7 q^{10} t^4+7 q^{11} t^4+2 q^{12} t^4+q^6 t^5+2 q^7 t^5
+4 q^8 t^5
+5 q^9 t^5+7 q^{10} t^5+8 q^{11} t^5+9 q^{12} t^5+3 q^{13} t^5
+q^7 t^6+2 q^8 t^6+4 q^9 t^6+5 q^{10} t^6+7 q^{11} t^6
+8 q^{12} t^6
+10 q^{13} t^6+3 q^{14} t^6+q^8 t^7+2 q^9 t^7+4 q^{10} t^7
+5 q^{11} t^7+7 q^{12} t^7+8 q^{13} t^7+10 q^{14} t^7+3 q^{15} t^7
+q^9 t^8+2 q^{10} t^8+4 q^{11} t^8+5 q^{12} t^8+7 q^{13} t^8
+8 q^{14} t^8+10 q^{15} t^8+3 q^{16} t^8+q^{10} t^9+2 q^{11} t^9
+4 q^{12} t^9+5 q^{13} t^9+7 q^{14} t^9+8 q^{15} t^9+9 q^{16} t^9
+2 q^{17} t^9+q^{11} t^{10}+2 q^{12} t^{10}+4 q^{13} t^{10}
+5 q^{14} t^{10}+7 q^{15} t^{10}+8 q^{16} t^{10}+7 q^{17} t^{10}
+q^{18} t^{10}+q^{12} t^{11}+2 q^{13} t^{11}+4 q^{14} t^{11}
+5 q^{15} t^{11}+7 q^{16} t^{11}+7 q^{17} t^{11}
+4 q^{18} t^{11}+q^{13} t^{12}+2 q^{14} t^{12}
+4 q^{15} t^{12}+5 q^{16} t^{12}+7 q^{17} t^{12}
+5 q^{18} t^{12}+q^{19} t^{12}+q^{14} t^{13}
+2 q^{15} t^{13}+4 q^{16} t^{13}+5 q^{17} t^{13}
+6 q^{18} t^{13}+2 q^{19} t^{13}+q^{15} t^{14}
+2 q^{16} t^{14}+4 q^{17} t^{14}+5 q^{18} t^{14}
+4 q^{19} t^{14}+q^{16} t^{15}+2 q^{17} t^{15}
+4 q^{18} t^{15}+4 q^{19} t^{15}+q^{20} t^{15}
+q^{17} t^{16}+2 q^{18} t^{16}+4 q^{19} t^{16}
+2 q^{20} t^{16}+q^{18} t^{17}+2 q^{19} t^{17}
+3 q^{20} t^{17}+q^{19} t^{18}+2 q^{20} t^{18}
+q^{21} t^{18}+q^{20} t^{19}+q^{21} t^{19}+q^{21} t^{20}\bigr).
\)
}
\renewcommand{\baselinestretch}{1.2} 
```

$Cab(22,3)T(3,2)$, $\vec r=\{3,3\},\vec s=\{2,4\}$:

`1605.00978/src/jacobian-fin-full.tex:5204-5429`
```latex
\subsection{\bf Superpolynomial}\label{sec:sup6-9-22}
Let us provide the corresponding DAHA superpolynomial,
which is for $\vec\rr=\{3,3\},\,
\vec\ss=\{2,4\}$. Recall that
the DAHA construction requires the matrices from
$PSL_2(\Z)$ with the first columns $(\rr_i,\ss_i)^{tr}$,
namely:
{\small $
\begin{pmatrix}3,*\\2,*\end{pmatrix}\,=\,
\tau_+\tau_-^2\,,\ \ 
 \begin{pmatrix}3,*\\4,*\end{pmatrix}\,=\,
\tau_-\tau_+^2\tau_-\,.
$} Then:
$$
\vec\rr=\{3,3\},\,
\vec\ss=\{2,4\},\,
\t= C\!ab(22,3)T(3,2);\ 
\h_{\,\vec\rr,\,\vec\ss}\,(\square\,;\,\,q,t,a)=
$$
\renewcommand{\baselinestretch}{0.5} 
{\small
\(
1+q t+q^2 t+q^3 t+q^4 t+q^5 t+q^2 t^2+q^3 t^2+2 q^4 t^2+2 q^5 t^2
+3 q^6 t^2+2 q^7 t^2+q^8 t^2+q^3 t^3+q^4 t^3+2 q^5 t^3+3 q^6 t^3
+4 q^7 t^3+4 q^8 t^3+4 q^9 t^3+q^{10} t^3+q^4 t^4+q^5 t^4
+2 q^6 t^4+3 q^7 t^4+5 q^8 t^4+5 q^9 t^4+6 q^{10} t^4
+4 q^{11} t^4+q^{12} t^4+q^5 t^5+q^6 t^5+2 q^7 t^5+3 q^8 t^5
+5 q^9 t^5+6 q^{10} t^5+7 q^{11} t^5+6 q^{12} t^5+4 q^{13} t^5
+q^{14} t^5+q^6 t^6+q^7 t^6+2 q^8 t^6+3 q^9 t^6+5 q^{10} t^6
+6 q^{11} t^6+8 q^{12} t^6+7 q^{13} t^6+6 q^{14} t^6+3 q^{15} t^6
+q^7 t^7+q^8 t^7+2 q^9 t^7+3 q^{10} t^7+5 q^{11} t^7+6 q^{12} t^7
+8 q^{13} t^7+8 q^{14} t^7+7 q^{15} t^7+3 q^{16} t^7+q^{17} t^7
+q^8 t^8+q^9 t^8+2 q^{10} t^8+3 q^{11} t^8+5 q^{12} t^8
+6 q^{13} t^8+8 q^{14} t^8+8 q^{15} t^8+8 q^{16} t^8+3 q^{17} t^8
+q^9 t^9+q^{10} t^9+2 q^{11} t^9+3 q^{12} t^9+5 q^{13} t^9
+6 q^{14} t^9+8 q^{15} t^9+8 q^{16} t^9+7 q^{17} t^9+3 q^{18} t^9
+q^{10} t^{10}+q^{11} t^{10}+2 q^{12} t^{10}+3 q^{13} t^{10}
+5 q^{14} t^{10}+6 q^{15} t^{10}+8 q^{16} t^{10}+8 q^{17} t^{10}
+6 q^{18} t^{10}+q^{19} t^{10}+q^{11} t^{11}+q^{12} t^{11}
+2 q^{13} t^{11}+3 q^{14} t^{11}+5 q^{15} t^{11}+6 q^{16} t^{11}
+8 q^{17} t^{11}+7 q^{18} t^{11}+4 q^{19} t^{11}+q^{12} t^{12}
+q^{13} t^{12}+2 q^{14} t^{12}+3 q^{15} t^{12}+5 q^{16} t^{12}
+6 q^{17} t^{12}+8 q^{18} t^{12}+6 q^{19} t^{12}+q^{20} t^{12}
+q^{13} t^{13}+q^{14} t^{13}+2 q^{15} t^{13}+3 q^{16} t^{13}
+5 q^{17} t^{13}+6 q^{18} t^{13}+7 q^{19} t^{13}+4 q^{20} t^{13}
+q^{14} t^{14}+q^{15} t^{14}+2 q^{16} t^{14}+3 q^{17} t^{14}
+5 q^{18} t^{14}+6 q^{19} t^{14}+6 q^{20} t^{14}+q^{21} t^{14}
+q^{15} t^{15}+q^{16} t^{15}+2 q^{17} t^{15}+3 q^{18} t^{15}
+5 q^{19} t^{15}+5 q^{20} t^{15}+4 q^{21} t^{15}+q^{16} t^{16}
+q^{17} t^{16}+2 q^{18} t^{16}+3 q^{19} t^{16}+5 q^{20} t^{16}
+4 q^{21} t^{16}+q^{22} t^{16}+q^{17} t^{17}+q^{18} t^{17}
+2 q^{19} t^{17}+3 q^{20} t^{17}+4 q^{21} t^{17}+2 q^{22} t^{17}
+q^{18} t^{18}+q^{19} t^{18}+2 q^{20} t^{18}+3 q^{21} t^{18}
+3 q^{22} t^{18}+q^{19} t^{19}+q^{20} t^{19}+2 q^{21} t^{19}
+2 q^{22} t^{19}+q^{23} t^{19}+q^{20} t^{20}+q^{21} t^{20}
+2 q^{22} t^{20}+q^{23} t^{20}+q^{21} t^{21}+q^{22} t^{21}
+q^{23} t^{21}+q^{22} t^{22}+q^{23} t^{22}+q^{23} t^{23}
+q^{24} t^{24}
\)

\vfil\smallskip
\noindent
\(
+a^5 \bigl(q^{15}+q^{16} t+q^{17} t+q^{17} t^2
+q^{18} t^2+q^{19} t^2+q^{18} t^3+q^{19} t^3+q^{20} t^3
+q^{21} t^3+q^{19} t^4+q^{20} t^4+q^{21} t^4+q^{20} t^5
+q^{21} t^5+q^{22} t^5+q^{21} t^6+q^{22} t^6+q^{22} t^7
+q^{23} t^7+q^{23} t^8+q^{24} t^9\bigr)
\)

%\vfil%\smallskip
\noindent
\(
+a^4 \bigl(q^{10}
+q^{11}+q^{12}+q^{13}+q^{14}+q^{11} t+2 q^{12} t+3 q^{13} t
+3 q^{14} t+3 q^{15} t+q^{16} t+q^{12} t^2+2 q^{13} t^2
+4 q^{14} t^2+5 q^{15} t^2+5 q^{16} t^2+3 q^{17} t^2+q^{18} t^2
+q^{13} t^3+2 q^{14} t^3+4 q^{15} t^3+6 q^{16} t^3+7 q^{17} t^3
+5 q^{18} t^3+3 q^{19} t^3+q^{20} t^3+q^{14} t^4+2 q^{15} t^4
+4 q^{16} t^4+6 q^{17} t^4+8 q^{18} t^4+6 q^{19} t^4+3 q^{20} t^4
+q^{21} t^4+q^{15} t^5+2 q^{16} t^5+4 q^{17} t^5+6 q^{18} t^5
+8 q^{19} t^5+6 q^{20} t^5+3 q^{21} t^5+q^{16} t^6+2 q^{17} t^6
+4 q^{18} t^6+6 q^{19} t^6+8 q^{20} t^6+5 q^{21} t^6+q^{22} t^6
+q^{17} t^7+2 q^{18} t^7+4 q^{19} t^7+6 q^{20} t^7+7 q^{21} t^7
+3 q^{22} t^7+q^{18} t^8+2 q^{19} t^8+4 q^{20} t^8+6 q^{21} t^8
+5 q^{22} t^8+q^{23} t^8+q^{19} t^9+2 q^{20} t^9+4 q^{21} t^9
+5 q^{22} t^9+3 q^{23} t^9+q^{20} t^{10}+2 q^{21} t^{10}
+4 q^{22} t^{10}+3 q^{23} t^{10}+q^{24} t^{10}+q^{21} t^{11}
+2 q^{22} t^{11}+3 q^{23} t^{11}+q^{24} t^{11}+q^{22} t^{12}
+2 q^{23} t^{12}+q^{24} t^{12}+q^{23} t^{13}+q^{24} t^{13}
+q^{24} t^{14}\bigr)
\)

\vfil%\smallskip
\noindent
\(
+a^3 \bigl(q^6+q^7+2 q^8+2 q^9+2 q^{10}
+q^{11}+q^{12}+q^7 t+2 q^8 t+4 q^9 t+6 q^{10} t+7 q^{11} t
+6 q^{12} t+4 q^{13} t+2 q^{14} t+q^8 t^2+2 q^9 t^2+5 q^{10} t^2
+8 q^{11} t^2+12 q^{12} t^2+12 q^{13} t^2+10 q^{14} t^2
+5 q^{15} t^2+2 q^{16} t^2+q^9 t^3+2 q^{10} t^3+5 q^{11} t^3
+9 q^{12} t^3+14 q^{13} t^3+17 q^{14} t^3+16 q^{15} t^3
+11 q^{16} t^3+5 q^{17} t^3+2 q^{18} t^3+q^{10} t^4+2 q^{11} t^4
+5 q^{12} t^4+9 q^{13} t^4+15 q^{14} t^4+19 q^{15} t^4
+21 q^{16} t^4+16 q^{17} t^4+9 q^{18} t^4+3 q^{19} t^4+q^{20} t^4
+q^{11} t^5+2 q^{12} t^5+5 q^{13} t^5+9 q^{14} t^5+15 q^{15} t^5
+20 q^{16} t^5+23 q^{17} t^5+19 q^{18} t^5+10 q^{19} t^5
+3 q^{20} t^5+q^{12} t^6+2 q^{13} t^6+5 q^{14} t^6+9 q^{15} t^6
+15 q^{16} t^6+20 q^{17} t^6+24 q^{18} t^6+19 q^{19} t^6
+9 q^{20} t^6+2 q^{21} t^6+q^{13} t^7+2 q^{14} t^7+5 q^{15} t^7
+9 q^{16} t^7+15 q^{17} t^7+20 q^{18} t^7+23 q^{19} t^7
+16 q^{20} t^7+5 q^{21} t^7+q^{14} t^8+2 q^{15} t^8+5 q^{16} t^8
+9 q^{17} t^8+15 q^{18} t^8+20 q^{19} t^8+21 q^{20} t^8
+11 q^{21} t^8+2 q^{22} t^8+q^{15} t^9+2 q^{16} t^9+5 q^{17} t^9
+9 q^{18} t^9+15 q^{19} t^9+19 q^{20} t^9+16 q^{21} t^9
+5 q^{22} t^9+q^{16} t^{10}+2 q^{17} t^{10}+5 q^{18} t^{10}
+9 q^{19} t^{10}+15 q^{20} t^{10}+17 q^{21} t^{10}
+10 q^{22} t^{10}+2 q^{23} t^{10}+q^{17} t^{11}+2 q^{18} t^{11}
+5 q^{19} t^{11}+9 q^{20} t^{11}+14 q^{21} t^{11}
+12 q^{22} t^{11}+4 q^{23} t^{11}+q^{18} t^{12}+2 q^{19} t^{12}
+5 q^{20} t^{12}+9 q^{21} t^{12}+12 q^{22} t^{12}+6 q^{23} t^{12}
+q^{24} t^{12}+q^{19} t^{13}+2 q^{20} t^{13}+5 q^{21} t^{13}
+8 q^{22} t^{13}+7 q^{23} t^{13}+q^{24} t^{13}+q^{20} t^{14}
+2 q^{21} t^{14}+5 q^{22} t^{14}+6 q^{23} t^{14}+2 q^{24} t^{14}
+q^{21} t^{15}+2 q^{22} t^{15}+4 q^{23} t^{15}+2 q^{24} t^{15}
+q^{22} t^{16}+2 q^{23} t^{16}+2 q^{24} t^{16}+q^{23} t^{17}
+q^{24} t^{17}+q^{24} t^{18}\bigr)
\)

\vfil%\smallskip
\noindent
\(
+a^2 \bigl(q^3+q^4+2 q^5
+2 q^6+2 q^7+q^8+q^9+q^4 t+2 q^5 t+4 q^6 t+6 q^7 t+8 q^8 t
+7 q^9 t+6 q^{10} t+3 q^{11} t+q^{12} t+q^5 t^2+2 q^6 t^2
+5 q^7 t^2+8 q^8 t^2+13 q^9 t^2+15 q^{10} t^2+15 q^{11} t^2
+10 q^{12} t^2+5 q^{13} t^2+q^{14} t^2+q^6 t^3+2 q^7 t^3
+5 q^8 t^3+9 q^9 t^3+15 q^{10} t^3+20 q^{11} t^3+24 q^{12} t^3
+19 q^{13} t^3+12 q^{14} t^3+5 q^{15} t^3+q^{16} t^3+q^7 t^4
+2 q^8 t^4+5 q^9 t^4+9 q^{10} t^4+16 q^{11} t^4+22 q^{12} t^4
+29 q^{13} t^4+28 q^{14} t^4+21 q^{15} t^4+11 q^{16} t^4
+4 q^{17} t^4+q^{18} t^4+q^8 t^5+2 q^9 t^5+5 q^{10} t^5
+9 q^{11} t^5+16 q^{12} t^5+23 q^{13} t^5+31 q^{14} t^5
+33 q^{15} t^5+29 q^{16} t^5+16 q^{17} t^5+6 q^{18} t^5
+q^{19} t^5+q^9 t^6+2 q^{10} t^6+5 q^{11} t^6+9 q^{12} t^6
+16 q^{13} t^6+23 q^{14} t^6+32 q^{15} t^6+35 q^{16} t^6
+32 q^{17} t^6+19 q^{18} t^6+6 q^{19} t^6+q^{20} t^6+q^{10} t^7
+2 q^{11} t^7+5 q^{12} t^7+9 q^{13} t^7+16 q^{14} t^7
+23 q^{15} t^7+32 q^{16} t^7+36 q^{17} t^7+32 q^{18} t^7
+16 q^{19} t^7+4 q^{20} t^7+q^{11} t^8+2 q^{12} t^8+5 q^{13} t^8
+9 q^{14} t^8+16 q^{15} t^8+23 q^{16} t^8+32 q^{17} t^8
+35 q^{18} t^8+29 q^{19} t^8+11 q^{20} t^8+q^{21} t^8+q^{12} t^9
+2 q^{13} t^9+5 q^{14} t^9+9 q^{15} t^9+16 q^{16} t^9
+23 q^{17} t^9+32 q^{18} t^9+33 q^{19} t^9+21 q^{20} t^9
+5 q^{21} t^9+q^{13} t^{10}+2 q^{14} t^{10}+5 q^{15} t^{10}
+9 q^{16} t^{10}+16 q^{17} t^{10}+23 q^{18} t^{10}
+31 q^{19} t^{10}+28 q^{20} t^{10}+12 q^{21} t^{10}+q^{22} t^{10}
+q^{14} t^{11}+2 q^{15} t^{11}+5 q^{16} t^{11}+9 q^{17} t^{11}
+16 q^{18} t^{11}+23 q^{19} t^{11}+29 q^{20} t^{11}
+19 q^{21} t^{11}+5 q^{22} t^{11}+q^{15} t^{12}+2 q^{16} t^{12}
+5 q^{17} t^{12}+9 q^{18} t^{12}+16 q^{19} t^{12}
+22 q^{20} t^{12}+24 q^{21} t^{12}+10 q^{22} t^{12}
+q^{23} t^{12}+q^{16} t^{13}+2 q^{17} t^{13}+5 q^{18} t^{13}
+9 q^{19} t^{13}+16 q^{20} t^{13}+20 q^{21} t^{13}
+15 q^{22} t^{13}+3 q^{23} t^{13}+q^{17} t^{14}+2 q^{18} t^{14}
+5 q^{19} t^{14}+9 q^{20} t^{14}+15 q^{21} t^{14}
+15 q^{22} t^{14}+6 q^{23} t^{14}+q^{18} t^{15}+2 q^{19} t^{15}
+5 q^{20} t^{15}+9 q^{21} t^{15}+13 q^{22} t^{15}+7 q^{23} t^{15}
+q^{24} t^{15}+q^{19} t^{16}+2 q^{20} t^{16}+5 q^{21} t^{16}
+8 q^{22} t^{16}+8 q^{23} t^{16}+q^{24} t^{16}+q^{20} t^{17}
+2 q^{21} t^{17}+5 q^{22} t^{17}+6 q^{23} t^{17}+2 q^{24} t^{17}
+q^{21} t^{18}+2 q^{22} t^{18}+4 q^{23} t^{18}+2 q^{24} t^{18}
+q^{22} t^{19}+2 q^{23} t^{19}+2 q^{24} t^{19}+q^{23} t^{20}
+q^{24} t^{20}+q^{24} t^{21}\bigr)
\)

\vskip 0.1cm
%\vfil
\noindent
\(+a \bigl(q+q^2+q^3+q^4+q^5
+q^2 t+2 q^3 t+3 q^4 t+4 q^5 t+5 q^6 t+4 q^7 t+2 q^8 t+q^9 t
+q^3 t^2+2 q^4 t^2+4 q^5 t^2+6 q^6 t^2+9 q^7 t^2+10 q^8 t^2
+9 q^9 t^2+5 q^{10} t^2+2 q^{11} t^2+q^4 t^3+2 q^5 t^3+4 q^6 t^3
+7 q^7 t^3+11 q^8 t^3+14 q^9 t^3+16 q^{10} t^3+13 q^{11} t^3
+6 q^{12} t^3+2 q^{13} t^3+q^5 t^4+2 q^6 t^4+4 q^7 t^4+7 q^8 t^4
+12 q^9 t^4+16 q^{10} t^4+20 q^{11} t^4+20 q^{12} t^4
+14 q^{13} t^4+6 q^{14} t^4+2 q^{15} t^4+q^6 t^5+2 q^7 t^5
+4 q^8 t^5+7 q^9 t^5+12 q^{10} t^5+17 q^{11} t^5+22 q^{12} t^5
+24 q^{13} t^5+21 q^{14} t^5+13 q^{15} t^5+4 q^{16} t^5
+q^{17} t^5+q^7 t^6+2 q^8 t^6+4 q^9 t^6+7 q^{10} t^6
+12 q^{11} t^6+17 q^{12} t^6+23 q^{13} t^6+26 q^{14} t^6
+25 q^{15} t^6+17 q^{16} t^6+7 q^{17} t^6+q^{18} t^6+q^8 t^7
+2 q^9 t^7+4 q^{10} t^7+7 q^{11} t^7+12 q^{12} t^7+17 q^{13} t^7
+23 q^{14} t^7+27 q^{15} t^7+27 q^{16} t^7+18 q^{17} t^7
+7 q^{18} t^7+q^{19} t^7+q^9 t^8+2 q^{10} t^8+4 q^{11} t^8
+7 q^{12} t^8+12 q^{13} t^8+17 q^{14} t^8+23 q^{15} t^8
+27 q^{16} t^8+27 q^{17} t^8+17 q^{18} t^8+4 q^{19} t^8
+q^{10} t^9+2 q^{11} t^9+4 q^{12} t^9+7 q^{13} t^9+12 q^{14} t^9
+17 q^{15} t^9+23 q^{16} t^9+27 q^{17} t^9+25 q^{18} t^9
+13 q^{19} t^9+2 q^{20} t^9+q^{11} t^{10}+2 q^{12} t^{10}
+4 q^{13} t^{10}+7 q^{14} t^{10}+12 q^{15} t^{10}
+17 q^{16} t^{10}+23 q^{17} t^{10}+26 q^{18} t^{10}
+21 q^{19} t^{10}+6 q^{20} t^{10}+q^{12} t^{11}+2 q^{13} t^{11}
+4 q^{14} t^{11}+7 q^{15} t^{11}+12 q^{16} t^{11}
+17 q^{17} t^{11}+23 q^{18} t^{11}+24 q^{19} t^{11}
+14 q^{20} t^{11}+2 q^{21} t^{11}+q^{13} t^{12}+2 q^{14} t^{12}
+4 q^{15} t^{12}+7 q^{16} t^{12}+12 q^{17} t^{12}
+17 q^{18} t^{12}+22 q^{19} t^{12}+20 q^{20} t^{12}
+6 q^{21} t^{12}+q^{14} t^{13}+2 q^{15} t^{13}+4 q^{16} t^{13}
+7 q^{17} t^{13}+12 q^{18} t^{13}+17 q^{19} t^{13}
+20 q^{20} t^{13}+13 q^{21} t^{13}+2 q^{22} t^{13}+q^{15} t^{14}
+2 q^{16} t^{14}+4 q^{17} t^{14}+7 q^{18} t^{14}+12 q^{19} t^{14}
+16 q^{20} t^{14}+16 q^{21} t^{14}+5 q^{22} t^{14}+q^{16} t^{15}
+2 q^{17} t^{15}+4 q^{18} t^{15}+7 q^{19} t^{15}+12 q^{20} t^{15}
+14 q^{21} t^{15}+9 q^{22} t^{15}+q^{23} t^{15}+q^{17} t^{16}
+2 q^{18} t^{16}+4 q^{19} t^{16}+7 q^{20} t^{16}+11 q^{21} t^{16}
+10 q^{22} t^{16}+2 q^{23} t^{16}+q^{18} t^{17}+2 q^{19} t^{17}
+4 q^{20} t^{17}+7 q^{21} t^{17}+9 q^{22} t^{17}+4 q^{23} t^{17}
+q^{19} t^{18}+2 q^{20} t^{18}+4 q^{21} t^{18}+6 q^{22} t^{18}
+5 q^{23} t^{18}+q^{20} t^{19}+2 q^{21} t^{19}+4 q^{22} t^{19}
+4 q^{23} t^{19}+q^{24} t^{19}+q^{21} t^{20}+2 q^{22} t^{20}
+3 q^{23} t^{20}+q^{24} t^{20}+q^{22} t^{21}+2 q^{23} t^{21}
+q^{24} t^{21}+q^{23} t^{22}+q^{24} t^{22}+q^{24} t^{23}\bigr).
\)
}
\renewcommand{\baselinestretch}{1.2} 
```

### 4.14 Links containing these cables (not knots; only pointers)

CD15 §11 "Double Cab(13,2)T(3,2)" (`CD15:13328-14436`) gives $\hat{\mathcal H}^{min}$ of the 2-component link $Cab(13,2)T(6,4)$, label
$\{3,2\},\{2,1\}$ with two arrowheads; its $q=1$ value is $(1+a)(1+a+t)\times\mathcal H_{Cab(13,2)T(3,2)}(q=1)$ (CD15:14420-14435).
CD15 also has $T(3,2)$ colored by $(2,1)$ only inside 2-links (`CD15:9466-9613`); there is no knot colored by a non-rectangular diagram among the cable examples.

---------------------------------------------------------------------------------------------------

## 5. Positivity; algebraic vs non-algebraic

CD14 (uncolored and rectangle-colored algebraic knots; pseudo-algebraic knots; KhR):

`1408.4348/src/cable.tex:1055-1091`
```latex
{\sf What DAHA provides.}
Conjecturally, the coefficient of the DAHA superpolynomial
$\h_{\{\vec\rr,\vec\ss\,\}}$ (to be defined below)
at $a=0$ and for $ q=1$ in the uncolored case 
is $\sum_{k=0}^{\de}h_{2k}t^k$ for {\em any\,} unibranch 
$\c_{\vec\rr,\vec\ss}$,
which matches the tables/conjecture from \cite{Pi}. 
This includes the {\em van Straten-Warmt conjecture\,}, 
which claims that all odd Betti numbers 
of $J(\c_{\vec\rr,\vec\ss})$ vanish, since
the DAHA superpolynomials contain only integral powers of $q,t,a$.    
In the opposite direction, this would
imply the positivity of $\h_{\{\vec\rr,\vec\ss\,\}}$ at
$a=0,q=1$. The whole (uncolored) $\h_{\{\vec\rr,\vec\ss\,\}}$
describes the geometry of Hilbert schemes of 
$\,\c_{\vec\rr,\vec\ss}\,$ due to 
Conjecture 2 from \cite{ORS}.

The general positivity conjecture claims that the
coefficient of the whole $\h_{\{\vec\rr,\vec\ss\,\}}$ 
are positive for any rectangle Young diagrams and algebraic 
knots. Assuming such a positivity
(the corresponding knots will be called
{\em pseudo-algebraic\,}), $\h_{\{\vec\rr,\vec\ss\,\}}$  
is expected to coincide with the {\em stable reduced 
Khovanov-Rozansky polynomial\,} of
$\t(\vec\rr,\vec\ss)$, producing (by definition) the 
Khovanov-Rozansky polynomials for $sl_n$ upon $a=-t^{n+1}$
for sufficiently large $n$. The substitution $a=-1$ here
is expected to result in the {\em Heegard-Floer homology\,}.
This is stated in the DAHA parameters $q,t,a$;  
see Parts $(ii,iii)$  of Conjecture \ref{CONCONJ}. 
We restrict ourselves there to the uncolored case, 
because this is needed in \cite{ORS} and in the definition  
of the Khovanov-Rozansky triply-graded homology. 
Though see \cite{Web} and other
papers on the {\em categorification\,} (via Quantum Group).
```

`1408.4348/src/cable.tex:1445-1452`
```latex
After papers \cite{GoN,CJJ},
the last unproven ``DAHA-intrinsic" 
conjecture from \cite{CJ} is the positivity of 
all $q,t,a$\~coefficients of  (reduced) 
$\,\h_{\rr,\ss}(b\,;\,q,t,a)\,$ for any rectangle $\la_b$. 
We expect this 
to hold for $\h_{\,\vec\rr,\,\vec\ss}\,(b\,;\,q,t,a)$ for 
any {\em algebraic\,} knots.
```

`1408.4348/src/cable.tex:2885-2903`
```latex
\subsection{\bf Pseudo-algebraic knots}
This section is a theoretical challenge.
It indicates that some non-algebraic iterated
knots are very similar to the
algebraic ones. The positivity
of their DAHA-superpolynomials is the key 
``algebraic" feature.
Recall that we call such knots {\em pseudo-algebraic\,},
which conjecturally includes all algebraic knots
and their mirror images for any rectangle Young diagrams. 

The two uncolored knots posses such a positivity and all 
other DAHA properties of algebraic knots, for instance
the existence of the tilde- normalization from
Theorem \ref{JONITER} for any $A_N$ 
and the strict equality for deg$_a$ in (\ref{tilde-sup}).
We give below some other examples of pseudo-algebraic 
knots that do not have these two features; see
(\ref{nac7},\ref{nac-3}).
```

`1408.4348/src/cable.tex:3168-3180`
```latex
\smallskip

The knot $C\!ab(5,2)T(3,2)$ is the first with non-positive 
$\h_{\vec\rr,\vec\ss}(\square;q,t,a)$ in the
considered family. In this and 
further non-positive examples below the actual Khovanov 
polynomials can be obtained from DAHA-Khovanov by the 
substitution $-1\mapsto 1/t$, though the tilde-normalization
creates here some problems and $D\!A\!H\!A'K\!h(C\!ab(1,2)T(3,2))$ 
contains  the term $q^{10}(t^4 - t^6)$, which has no counterpart
in the corresponding $\tilde{K\!h}$, presumably due to the 
(topological) differential. Examples show that such a simple 
substitution is not generally sufficient.
```

CD15 (links; non-rectangle diagrams fail; positivity after dividing by $(1-t)^\bullet$):

`1509.08351/src/links.tex:6245-6265`
```latex
(ii) Now let $\{\l,\,'\!\l^\vee\}$ be a pair of trees 
such that $\rr_i^j,\,\ss_i^j,\,'\rr_i^j,\,'\ss_i^j>0$, and 
$\,'\ss_1\ss_1>\,'\rr_1\rr_1$, where $1\le j\le \kappa$, $ 1\le 
i\le \ell^j$ and $1\le j\le \,'\!\kappa$, $ 1\le i\le 
\,'\!\ell^j$. Here and in $(iii)$ only 
$\,'\!\l^\vee$ with $\vee$ (can be empty) is considered and the 
integral form is $\hat{P}_\la=J_\la$. Then formula 
(\ref{deg-a-jj}) for deg$_a$ is conjectured to be exact. In the 
uncolored case, we also expect the positivity of the following
series: 
%diagrams $(\la^j,\,'\!\la^j)$, we also expect the positivity:
\begin{align}\label{posit-claim}
\hat{\h}^{min}_{\,\l,\,'\!\l^\vee}\,(q,t,a)/(1-t)^
{\kappa+'\!\kappa-1}\in
\Z_+[[q,t,a]]
\end{align}
upon the natural $t$\~expansion of this ratio. For any diagrams,
it is conjectured to hold for sufficiently
large powers of $(1-t)(1-q)$\, (provided $\rr_i^j,\,\ss_i^j,\,
'\rr_i^j,\,'\ss_i^j,\,
'\ss_1\ss_1\!-\!'\rr_1\rr_1>0$ and with the usage of $\vee$).
```

`1509.08351/src/links.tex:6367-6394`
```latex
One of the reasons for such an uncertainty
is that the positivity of the polynomials 
$\hat{\h}^{min}_{\,\l,\,'\!\l^\vee}\,(q,t,a)$
does not hold for links (including uncolored ones) and for knots
if the corresponding Young diagrams are non-rectangle. We
address it in the positivity conjecture from Part $(ii)$,
which however corresponds to the {\em non-reduced theory}. 
Recall that (\ref{posit-claim}) is for any  
{\em positive pairs\,}  $\{\l,\,'\!\l^\vee\}$ and can be 
extended to any weights. 
 The experiments 
show that the positivity of the series there almost always  
fails when the inequalities from Part $(ii)$  are not satisfied 
(even in the uncolored case). This is in contrast to
the Khovanov-Rozansky theory, where knots are arbitrary.  

For quite a few 
colored algebraic links, the division by $(1-t)^\bullet$ is 
sufficient for the positivity in Part $(ii)$. However 
we have examples when the positivity holds with 
$(1-q)^\bullet$ and fails for the division by
any powers of $(1-t)$.
There are also examples when only proper powers of
$(1-q)(1-t)$ ensure the positivity. We did 
not reach any conjectures concerning the occurrence of
$t$ or $q$ here
and the minimal powers of these corresponding binomials.
\smallskip
```

C24:

`2410.14703/src/superpolynoms.tex:673-686`
```latex
Adding $a$ and the colors is beyond the Hasse-Weil-Deligne theory 
\cite{Del1,Del2}, as well as the following feature.
All coefficients of DAHA superpolynomials are presumably 
positive integers for rectangle Young diagrams and
algebraic knots. This is the last 
unproven {\em intrinsic} conjecture from the list
in \cite{CJ}, extended to any algebraic knots 
in  \cite{ChD1}.  For links and non-rectangle 
Young diagrams, the positivity conjecturally holds only 
upon the division by
some powers of $(1-t), (1-q)$. For instance, the division by
$(1-t)^{\kap-1}$ is presumably sufficient
for uncolored algebraic links with $\kap$ components.
\vskip 0.2cm
```

CP16:

`1605.00978/src/jacobian-fin-full.tex:236-244`
```latex
The case
of rectangle Young diagrams is exceptional due to 
the conjectured positivity of the corresponding
{\em reduced\,} DAHA 
superpolynomials for algebraic knots \cite{CJ,ChD1}.
The switch from  the rank-one torsion free modules
in the definition of compactified Jacobians 
to arbitrary ranks is expected here (among other modifications),
which is in progress.
```

`1605.00978/src/jacobian-fin-full.tex:400-410`
```latex
Importantly, the positivity was conjectured in
\cite{CJ,ChD1} for 
the {\em reduced\,} superpolynomials of algebraic 
knots colored by any {\em rectangle\,} Young diagrams,
i.e. not only in the uncolored case. 
We are going to incorporate such diagrams in our approach. 
Generally, adding colors requires using links, but 
we see some direct ways for such diagrams. 
%symmetric/wedge powers of the fundamental representation 
%there are more direct ways.

```

`1605.00978/src/jacobian-fin-full.tex:4016-4026`
```latex
\subsection{\bf Adding colors}
In contrast to \cite{ChD2} and previous paper, we
restrict this one to algebraic uncolored knots and
only type $A$ is addressed in our conjecture. Let
us briefly comment on this. Adding colors is
expected via the curves in \cite{Ma}, so algebraic 
{\em links\,} are/seem necessary for this. Apart
from  rectangle Young diagrams, the coefficients
of DAHA-superpolynomials are non-positive, which is
an obvious challenge for the geometric interpretation.
This is the same for links (included uncolored ones). 
```

Summary of what the examples show [checked by coefficient signs]:

* algebraic cables (all $r_i,s_i>0$): all coefficients positive (examples (1)-(5), CP16, $2\omega_1$ included).
* trefoil family $Cab(2m+1,2)T(3,2)$: positive for $2m+1\ge 7$ (with (a) $Cab(7,2)$ having deg$_a=2<3$), **non-positive** for
  $2m+1=5,3,1,-1,-3$ (examples (b)-(f)), positive again for $2m+1=-5$ (g) "and all further knots in the considered family".
  In particular the **(2,3)- and (2,1)-cables of the trefoil have non-positive $\mathcal H$** (and $t^{-1}$ powers).
* pseudo-algebraic = positive uncolored $\mathcal H$; conjecturally contains algebraic knots and their mirrors (CD14:2892-2894).
* deg$_a$: conjecturally $=(s_1r_2\cdots r_\ell-1)|\lambda|$ for algebraic knots ($s_1<r_1$); this is $\le$ in general (CD14 Thm 2.3(iv)).

---------------------------------------------------------------------------------------------------

## 6. My numerical checks (scratch, not from the sources)

Parser: LaTeX -> Laurent polynomial (scratchpad `lpoly.py`, `check2.py`, `check3.py`). Results:

* $\mathcal H(q,t\mapsto q,a\mapsto -1)$ equals $\pm q^\bullet\,\Delta_{T(r_1,s_1)}(q^{r_2})\,\Delta_{T(r_2,|a_2|)}(q)$ for **all** knot examples transcribed
  above: CD14 (1),(2),(3),(4),(6),(7),(a)-(h); CP16 $Cab(15,2)T(3,2)$, $Cab(31,2)T(7,2)$, $Cab(22,3)T(3,2)$; C24 $Cab(5,2)T(3,2)$,
  $Cab(9,2)T(3,2)$, $Cab(11,2)T(4,3)$. This confirms both the transcription and the dictionary (winding number $=r_2$, slope $a_2$).
  (Alexander cannot distinguish $a_2=\pm1$, i.e. (d) vs (e); rely on the source labels there.)
* $\mathcal H(q=1,a=0)$ reproduces every DAHA-Betti polynomial printed in CD14 (bett1-bett7, nac7, nac5, nac3, nac1, nac-1, nac-3, nac-5).
* CD14 (5) ($2\omega_1$): positive, deg$_a=6$, and $\mathcal H(2\omega_1;q=1)=\mathcal H(\omega_1;q=1)^2$.
* CP16 $Cab(22,3)T(3,2)$: superduality $\mathcal H(q,t,a)=(qt)^{24}\mathcal H(t^{-1},q^{-1},a)$ holds exactly.
* CD14 (b) at $a=-t^2$ equals (jon5-2).
* Standard variables (from `conv.py`): 
  * (2,13)-cable: $\mathcal H_{st}=1+q^{4}t^{2}+q^{6}t^{4}+q^{8}t^{4}+\dots+q^{32}t^{16}+a^2(q^{2}t^{3}+\dots)+a^4(\dots)+a^6(q^{12}t^{15}+q^{16}t^{17}+q^{20}t^{19})$;
  * (2,3)-cable (c): $\mathcal H_{st}=1+q^{4}t^{2}-q^{4}t^{4}+q^{6}t^{4}+q^{8}t^{4}-q^{8}t^{6}+q^{12}t^{6}+a^2(q^{2}t^{3}-q^{2}t^{5}+q^{4}t^{5}-q^{4}t^{7}+q^{6}t^{5}-2q^{6}t^{7}+q^{8}t^{7}-q^{8}t^{9}+q^{10}t^{7}-q^{10}t^{9})+a^4(-q^{2}t^{8}-q^{4}t^{8}+q^{6}t^{8}-2q^{6}t^{10}-q^{8}t^{10}-q^{10}t^{12})+a^6(-q^{4}t^{11}-q^{8}t^{13})$;
  * (2,1)-cable (d): $\mathcal H_{st}=q^{2}-q^{2}t^{2}+q^{6}t^{2}-2q^{6}t^{4}+q^{10}t^{4}-q^{10}t^{6}+a^2(-t^{3}-q^{2}t^{5}+q^{4}t^{3}-3q^{4}t^{5}-q^{6}t^{7}+q^{8}t^{5}-3q^{8}t^{7}-q^{10}t^{9}-q^{12}t^{9})+a^4(-t^{6}-q^{2}t^{6}-2q^{4}t^{8}-2q^{6}t^{8}-2q^{8}t^{10}-q^{10}t^{10}-q^{12}t^{12})+a^6(-q^{2}t^{9}-q^{6}t^{11}-q^{10}t^{13})$.
