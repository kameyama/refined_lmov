# 論文の記号表（2026-10-09、執筆前の準備）

論文で使う記号と規約を、conventions.md から論文用に抜き出したもの。根拠の節番号は conventions.md のもの。
執筆中に記号を変えるときは、この表と conventions.md を同時に直す。

## 1. 変数
| 記号 | 意味 | 根拠 |
|---|---|---|
| (a, q, t) | KN の変数。q = e^{iε₁}, t = e^{iε₂}、a = t^N。Ω 背景 (z₁, z₂) ↦ (q z₁, t^{−1} z₂) | §1.1、KN :407-409 |
| (A, Q, T) = (a, q, t)^{1/2} | 実装の半整数冪の変数。論文では使わず、必要な箇所だけ脚注で触れる | §1.1 |
| (𝐚, 𝐪, 𝐭) | 太字の変数（DGR の変数）。a = −𝐚²𝐭、q^{1/2} = −𝐪𝐭、t^{1/2} = 𝐪（KN (cov)）。t = q は 𝐭 = −1 | §1.1 |
| 標準の変数（CD14） | t = q_st²、q = (q_st t_st)²、a = a_st² t_st。(q_st, t_st, a_st) は太字の変数と同じ（符号の分枝を除く） | §1.2、claims §7 の 5 |
| a_DAHA | a_DAHA = −a_KN | §1.2 |
| 符号の分枝 | √q = −𝐪𝐭（KN (cov) に従う）。全体の符号は最低次の項を正にする向き | §8 |

## 2. 対称関数と不変量
| 記号 | 意味 | 根拠 |
|---|---|---|
| P_λ(x; q, t) | monic な Macdonald 多項式 | §2 |
| g_λ(q,t) | (q/t)^{\|λ\|}/⟨P_λ, P_λ⟩_{q,t} | §2、KN |
| rCS_λ(K; a, q, t) | reduced の refined CS 不変量（framing 0、正規化 (N1)(N2')）。U(N) と SU(N) の違いは無視 | §4, §5 |
| T_λ | q^{n(λ')} t^{−n(λ)}。refined の framing の固有値（(N2') の T_λ） | §5 |
| n(λ), κ_λ | n(λ) = Σ(i−1)λ_i、n(λ') − n(λ) = Σ content。MV の κ_R = 2Σ content | §17、MV :647 |
| f^q_μ, f^{\bar t}_μ | KN の (refined-GT1)(GT2) の reformulated invariants。f^q は付録 B の尺度 | §8 |
| \hat f_ρ | f^q = M(t)\hat f、f^{\bar t} = M(q^{−1})\hat f で決まる量 | §8 |
| \hat{𝐍}_{ρ;i,j,k} | G 展開の前の係数：a 次数をずらした \hat f_ρ を太字にした Laurent 多項式の係数 | §8 |
| \hat N_{ρ,g,β,J_r} | G 展開の後の refined LMOV（refined BPS）不変量。G = (𝐪𝐭 − 𝐪⁻¹𝐭⁻¹)(𝐪 − 𝐪⁻¹) の冪で展開 | §8 |
| a 次数のずらし | nontorus 規約：\|ρ\| = 1 で a^{−1/2}、\|ρ\| = 3 で a^{+1/2} | §8 |
| 正値性 | 全体の符号を一つ選ぶと、\hat{𝐍} と \hat N の全ての項が非負整数 | §8 |
| σ | DGR の対称性 (e_𝐚, e_𝐪, e_𝐭) ↦ (e_𝐚, −e_𝐪, e_𝐭 − e_𝐪)。G 展開が閉じる条件 | REPORT_knot_table §2 |

## 3. framing
| 記号 | 意味 | 根拠 |
|---|---|---|
| θ_λ(q,t)（準備稿の記号、2026-10-10 に framing[λ] から変更） | (−(q/t)^{1/2})^{\|λ\|} q^{n(λ')} t^{−n(λ)} = (−(q/t)^{1/2})^{\|λ\|} T_λ。Taki の f_λ(t,q) の逆数。太字では 𝐭^{\|λ\|+2n(λ')} 𝐪^{2(n(λ')−n(λ))} | §17 |
| x | framing の整数。x = 0 が KN の場合（Seifert framing を基準にした (N2')） | §17 |
| rCS^{(x)}_λ | framing[λ]^x rCS_λ。Adams の写しでは framing[λ](q^k, t^k)^x（符号は k によらない） | §17 |
| \hat f^{(x)}_ρ, \hat{𝐍}^{(x)}, \hat N^{(x)} | rCS^{(x)} から同じ手順で作る量 | §17 |
| F | (N2') の framing の整数。ケーブルでは p₂q₂（q₂ は companion の Seifert framing で測る） | §5, §13 |

## 4. 結び目
| 記号 | 結び目 | 根拠 |
|---|---|---|
| T(m, n) | トーラス結び目。DGR の向き（三葉結び目 T(2,3) は負の交差） | §1.2 |
| C_{p,q}K | K の (p, q)-ケーブル。q は companion の Seifert framing で測る。DAHA の Newton 対は r = (p₁, p₂), s = (q₁, q₂ − p₁q₁p₂) | §5 |
| K1, K3, K6a, K9 | T(2,3) の (2,13)、T(2,3) の (2,15)、T(2,5) の (2,21)、T(2,7) の (2,31)-ケーブル（代数的） | RESULTS.md §1 |
| K8, K7, K10 | T(3,4) の (2,25)、T(2,3) の (3,19)、K1 の (2,53)-ケーブル（代数的） | RESULTS.md §1 |
| K2, K4, K5a, K5b, K5c, K6b | 非代数的なケーブル（付録） | REPORT.md |
| K_p = K_{p,1}, K_{m,n} | ツイスト結び目、ダブルツイスト結び目（KNTZ 1902.02275 の規約。GNSSS の (2m−1, 1, 2n−1)-pretzel K_{s,t} と同じ） | conventions.md §15 |
| δ | 種数型の量。T(m,n) は (m−1)(n−1)/2、ケーブルは δ ↦ pδ + (p−1)(q−1)/2 | §13 |
| b | ブレイドの本数 Π p_i | §13 |

## 5. 正規化
| 記号 | 意味 | 根拠 |
|---|---|---|
| (N1) | t = q で framing 0 の reduced 彩色 HOMFLY-PT（DGR の規約）に一致 | §5 |
| (N2) | [1] の最低 𝐭 次数が 0（正値性に影響しない (q/t)^c の自由度を決めるだけ） | §5 |
| (N2') | M_λ = M_□^{\|λ\|} T_λ^F | §5 |
| M_□ | (−q/t)^{b−1}(a(q/t)^{1/2})^δ | §13 |
| tilde 正規化 | CD14 の正規化。本研究の [1] の太字の P = 𝐚^{2δ}𝐪^{−2δ} × (tilde) | §13、claims §7 の 5 |
| H_sat | refined 衛星公式（仮説） | §13.1 |

## 6. 論文での表記の決めごと（執筆時に守る）
- **ヤング図は KN と同じく youngtab の `\yng(2,1)` で書く**（箱 4pt、2026-10-09、ユーザーの指示）：角括弧の [1], [2,1] は引用の番号 [n] と紛らわしいので本文では使わない。言葉を添えるときは「基本表現 ρ = \yng(1)」のように書く。`\yng` は見出しと caption の中では壊れるので、見出しは言葉で書き、caption には短い版を付ける。作業用の文書（claims.md など）では従来どおり [2,1] と書いてよい。
- **framing の因子は $\theta_\lambda(q,t)$**（2026-10-10、ユーザーの判断）：彩色 HOMFLY-PT の framing の固有値 θ_λ の refined 版として。Taki の $f_\mu(t,q)$ とは $f_\mu(t,q)=\theta_\mu(q,t)^{-1}$。Taki の記号は KN の $f^q_\mu$, $\widehat f_\rho$ とぶつかるので使わない。作業用の文書・コードでは従来どおり framing[λ] と書いてよい。
- **colored の訳語**（2026-10-10、ユーザーの指示で確認）：「彩色」ではなく「色付き」と書く（色付き HOMFLY-PT 多項式・ホモロジー、色付き Jones 多項式）。日本語の文献の用法：東京大学のトポロジー火曜セミナー（2013 年）の講演要旨「色付きHOMFLYホモロジー」「色付きHOMFLY多項式」、村上斉の講演資料と川添の博士論文（明治大学）の「色付きJones多項式」。「彩色」は結び目理論では Fox の n 彩色など別の概念に使われる。表現そのものを指すときは「表現」（例：表現ごとの単項式）、uncolored は「色付きでない」。
- **refined LMOV 不変量の二つの表示**（2026-10-10、ユーザーの指示）：genus 展開の前の $\widehat{\mathbf N}_{\rho;i,j,k}$ を **spin 表示**、後の KN の (bold-form) の係数 $\widehat{\mathbf N}_{\rho,g,\beta,F}$ を **genus 表示** と呼ぶ（2026-10-10 に記号を変更：$\widehat{\mathbf N}_{\rho,g,\beta,F}=\sum_{2(J_r-S_R)=F}(-1)^{2S_R}\widehat N_{\rho,g,\beta,J_r,S_R}$ で、$\widehat N_{\rho,g,\beta,J_r}$ と一致するのは $S_R$ が自明に作用するという予想の下だけ。$\widehat N_{\rho,g,\beta,J_r}$ は KN の予想の文脈でだけ使う）。根拠（KN の原典で確認）：genus 表示の $g$ は M2-brane の曲線 $\Sigma_{g,h}$ の bulk の genus で、$G^g$ は Jacobian $H^*((T^2)^g)$ の寄与（KN :735、:793-796）。spin 表示の単項式は refined index の $q^{J_1-S_R}t^{S_R-J_2}=\pm\mathbf q^{2(J_1-J_2)}\mathbf t^{2(J_1-S_R)}$ で、$J_1-J_2$ は OV の $N_{R,Q,s}$ の spin $s$（非 refined の index $q^{J_1-J_2}$、KN :665）。spin 表示は固有値（重み）でラベルされた符号付きの数。knots–quivers の文献の「LMOV invariants」（genus 展開しない型）とは呼び方が逆なので注意。
- **展開の用語**（2026-10-10、ユーザーの指示）：KN（:900）にならい、(hat-f-2) の展開は「（M2-brane の）genus 展開」と呼ぶ。準備稿では「G 展開」と書かない。$G=(\mathbf{qt}-\mathbf{q}^{-1}\mathbf{t}^{-1})(\mathbf q-\mathbf q^{-1})$ は略記として残す。展開の前の係数は KN の「素朴な変数変換 (cov)」の係数。作業用の文書・コードでは従来どおり G 展開と書いてよい。
- **rCS のバー**（2026-10-10、ユーザーの指示）：KN と同じく、不変量そのもの・生成関数・framing の定義・計算の対象は unreduced の $\overline{\rm rCS}_\lambda$。reduced の多項式そのものを指す箇所（[1] の関係式 $\widehat f_{[1]}=(a^{1/2}-a^{-1/2}){\rm rCS}_{[1]}$、転置対称性、2 箱の公式の入力、ホモロジーの Poincaré 多項式、係数の負、ケーブルの公式、正規化）は KN の脚注どおりバーなしの ${\rm rCS}_\lambda$。
- **KN に合わせる記号**（2026-10-09）：トーラス結び目は $T_{m,n}$。unreduced を $\overline{\rm rCS}_\lambda$、reduced を ${\rm rCS}_\lambda$（$\overline{\rm rCS}_\lambda(K)=\overline{\rm rCS}_\lambda(\bigcirc)\,{\rm rCS}_\lambda(K)$）。生成関数は unreduced で書く。G 展開の前の係数は KN の $\widehat{\mathbf N}_{\rho;i,j,k}$（$\mathbf a^{2i}\mathbf q^{2j}\mathbf t^k$）、後は $\widehat N_{\rho,g,\beta,J_r}$（KN (hat-f-2)）。予想は KN と同じく枠（tcolorbox）で囲む。
- superpolynomial は、ホモロジーの Poincaré 多項式（DGR などの意味）にだけ使う。rCS や DAHA の値は superpolynomial と呼ばない（§8 の用語）。
- 「refined LMOV 不変量」は \hat{𝐍} と \hat N だけを指す。f^q, f^{\bar t} の係数は付録の補足。
- Zhu の "refined LMOV conjecture"（非 refined の概念）と区別することを、序論で一度書く。
- 主張には、定理・命題・計算結果・予想・仮説の区別をつける（claims.md）。
