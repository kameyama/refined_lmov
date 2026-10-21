# 進捗記録（第 3 の研究：refined Kauffman の large N duality）

## 現在の状況（再開用、2026-10-21 更新）

新しいセッションでは、この節を読めば作業を再開できるようにしてある。日々の詳細はこの下の日付ごとの記録（続き 1〜52）にある。

### 1. 環境
- **worktree**：`/Users/masaya/projects/ai_research/refined_cs/.claude/worktrees/dgr-quadruply-graded-lomv-2a7fb6`
  - ブランチは `claude/refined-kauffman-large-n-00d58e`。
  - main への直接の commit とローカルでのマージは禁止。
- **PR**：[kameyama/refined_lmov#11](https://github.com/kameyama/refined_lmov/pull/11)（OPEN）
  - PR の head は 3de70c5。それ以降の **54 commit（最新 a7c5ac9、2026-10-21）は push していない**。
  - push は外部送信なので、毎回ユーザーの許可を得る。
  - SSH の `git push` はサンドボックスで通らない。前回は、許可を得たうえで gh api（Git Data API）で同じ SHA の commit を作った。
    - 手順：blobs → trees（base_tree は親）→ commits（作者・日時・メッセージを同じに）→ refs。
    - 各段で SHA がローカルと一致することを確かめる。
- **Python**：main checkout の `.venv` を使う（`/Users/masaya/projects/ai_research/refined_cs/.venv/bin/python`）。
  - worktree 直下から相対で書くと `../../../.venv/bin/python`。
  - スクリプトは worktree の `src` と `scripts` を sys.path の先頭に入れている。
- **テスト**：`../../../.venv/bin/pytest -q tests/test_rk_*.py`
  - 既定の 58 件（約 72 秒）。
  - `-m results` の 6 件（約 437 秒）。
  - 2026-10-20 に全件通過。ログは logs/refined_kauffman/paper_tests_{default,results}.log。
- **TeX**：`TEXMFVAR=$TMPDIR/texmfvar lualatex main.tex` を 2 回実行する。
  - ページを画像にするには gs を使う（pdftoppm は無い）。
  - .aux/.log/.out/.toc は git で追跡しない。
- **文献**：`references/<arXiv id>/src/`。取得は `just fetch-refs`（justfile の papers に一覧がある）、または `curl`（allowed_domains: arxiv.org）。書誌の確認は arXiv API（export.arxiv.org）で行う。
- **Dropbox の 2018 年のノートブック**は読むだけで、書き換えない。

### 2. 現在の段階：論文の草稿（日本語の準備稿）
- 場所：`docs/refined_kauffman/paper/`
  - `outline.md`：論文の軸、決定事項、章立て
  - `claims.md`：主張の格付け表
  - `ja/main.tex`、`ja/sections/*.tex`、`ja/main.pdf`（30 ページ）
- **決定事項（2026-10-20、ユーザーが確認済み）**：
  - 本研究の framed knot の論文（docs/paper/）とは別の 1 本にする。
  - 日本語で準備し、最終版は英語にする。
  - 著者は亀山昌也の単著。2018 年の試みは「著者が以前に行った計算」と書く。
  - 正値性は「未解決」とし、障害の定式化を成果の一つにする。
  - 題：「refined Kauffman 不変量の large N 双対：orientifold の refined 指数とクロスキャップの構造」。
- **草稿の状態**：別のエージェントによる記録との照合（誤り 7 件、要確認 14 件）を反映済み（続き 42）。
- 付録の式は `docs/refined_kauffman/tex/gen/*.tex`（`scripts/rk_report_tex.py` の生成物）をそのまま `\input` している。
- 参考文献は `ja/sections/bibliography.tex`。書誌は arXiv API で確認済みで、掲載誌の無いものは「掲載誌未確認」と書いてある。
- **次にやること（ユーザーの指示待ち）**：
  1. 日本語の本文の推敲（読み手向けの説明、記号の統一）
  2. 英語版の作成
  3. PR #11 の更新（push の許可が要る）

### 3. 主な結果（論文の主張。詳細は claims.md と theory_orientifold_index.md）
| 主張 | 種別 | 場所 |
|---|---|---|
| 双対 Cauchy 核（t-brane）、一様な規格化 (t/q)^{δ\|μ\|/2}、reformulation は Mariño の M(q^{−1}) | 定式化＋計算 | report.md §4・§6、conventions.md |
| 全 d の構造 g = f − ½f^{cov} = ½E/(1+t) + C（対・不変・凍結の分類） | 導出 | theory §2、草稿 §3 |
| E と C の分離は指数からは決まらない（g から決まるのは E の 2(1+t) を法とした類だけ） | 導出 | 草稿 prop:EC（theory §3 は 10-20 に訂正済み） |
| 閉弦：D_N の積表示と AgSc §6 は N の 1 次式を除いて一致（N = 2..7）。AM の Klein 瓶は凍結した状態と不変な状態の並べ替え | 計算＋導出 | theory §8、scripts/rk_b_closed_orientifold.py |
| unknot・Hopf（a, b ≤ 4）の整数性、ν~ = ν(q/t)^{1/2}、E_22・E_24 の一意性 | 計算 | report.md §4・§7 |
| T(2,2m+1) の [1]（m ≤ 5）：(N)(P)、c=2 の閉じた形 Π_m、c=1 の G 展開、判定式 | 計算＋予想 | report.md §8 |
| DAHA-Jones（D_N, GL_N、一般の t）で、trefoil の 2・3 箱と T(2,5) の 2 箱を一般の N で計算。大きな N のクロスキャップは整数、(1+t) は p₂ 方向だけ、t = q で Mariño を再現 | 計算 | report.md §9、tex/report.tex §11・§12 |
| 普遍因子：ĝ − ĝ\|_{t=q}(ν → νw^k) は k = +1 でだけ (q − t)(ν² − t) で割り切れる | 計算 | theory §8.3 |
| 正値性：[1] は判定式が成立。多箱は商の符号が混在。障害は凍結した状態の多重項の基底と E・C の分離 | 計算＋未解決 | theory §4〜§10、草稿 §7 |

### 4. 試して失敗した方法（同じことを繰り返さないため）
| 方法 | 結果 | 記録 |
|---|---|---|
| Mariño 型（bosonic 核、½、Σ_{d 奇数}）の素朴な refined 化 | unknot の 2 箱で U の極が残る | report.md §3 |
| 多箱の c=1 部分の G 展開 | σ 対称でなく展開できない。T_± を使っても同じ | report.tex §11、theory §8.2 |
| 自由エネルギーの右辺の基底の取り替え、二次量子化の指数、多重項の変数の探索 | 符号はそろわない | 続き 28〜31 |
| E = P − 2Q（状態の部分和）の整数計画 | trefoil で解なし（T(2,5) は未実行） | 続き 34 |
| Adams 型の E、t = −1 の差→和の規則、G の捻り | 条件を満たさない。当てはめになるので打ち切った | 続き 35・36 |
| 凍結した状態の展開の仮定（Z[G,w] と X×Z[G,w]） | [1] の c=1 だけ成立 | 続き 37 |
| 多重項ごとの捻ったトレースの整数計画 | 結び目の境界の因子は z_q^+。trefoil は一意に解けるが、T(2,5) は解なし | 続き 39、theory §9 |
| 被覆の随伴 R_adj を正の Poincaré 多項式として読む | どの符号規則でも正にならない（Cherednik–Elliot の記述と一致） | 続き 40、theory §10 |

### 5. 未解決（論文の §8 にも記載）
- 凍結した状態（向き付けられない面）の refined の取っ手の因子と、多重項の基底。
- E と C の分離。固定セクターの状態の空間への Ĩ の作用を、指数ではなく幾何から求める必要がある。Hopf の (4,4) の E も未決定。
- 普遍因子のずれの向き（k = +1）と、閉弦の Kähler の三重化との関係。開弦の比 w³ と閉弦の比 w^{−1} が合わず、保留している。
- 非トーラス結び目の大きな N の値（入力が無い）。T(2,5) の 3 箱は、m = 2 の composite と D_N の 3 箱が要り、数時間〜半日の見込み。
- 環状の HOMFLY ホモロジー（逆向きの 2 本のケーブル）と refined CS の composite の関係。未着手で、見込みは低い。

### 6. 主なファイル
- **文書**：
  - report.md（K0〜K4' と多箱の詳細）
  - tex/report.tex・report.pdf（計算結果の報告書。§11〜§12 に 10-20 までの結果）
  - theory_orientifold_index.md（指数の理論 §1〜§10）
  - literature.md（§4 に非 U(N) の rCS の文献）
  - conventions.md
- **コード**：
  - src/rcsknots/：kauffman_*.py、koornwinder.py、daha_dn.py、daha_fast.py、flint_rf*.py
  - scripts/：rk_*.py。rk_c_* は一般の N、rk_b_* は構造と正値性の解析、rk_k* は固定 N と Hopf
- **データ**：
  - data/refined_kauffman/daha/：R_U_m*_*.pkl（SO の一般の N）、Radj_A_m*.pkl、Rcomp_A_m1_*.pkl、largeN_2box_m{1,2}.pkl、largeN_3box_m1.pkl
  - data/refined_kauffman/：hopf_*.pkl

## 2026-10-14
- ディレクトリ docs/refined_kauffman/ を作成（README.md、progress.md）。
- 過去のノートブック 42 本（2018-07-10〜10-02）を scripts/nb2txt.py でテキスト化して読解を開始。同梱の PDF は Cherednik 1111.6195（DAHA-Jones、B/C 型の hyperpolynomial）と Cherednik–Elliot 1503.01441（refined composite invariants of torus knots via DAHA）（arXiv API で題名を確認）。
- 文献の取得：justfile の papers に 0904.1088（Mariño）、hep-th/0405083・hep-th/0411227（Bouchard–Florea–Mariño）、hep-th/0505008（Borhade–Ramadevi）、1003.5282（Paul–Borhade–Ramadevi）、1003.2861（Stevan）、hep-th/0512298（Gukov–Walcher）、1202.4456（Aganagic–Schaeffer）、1111.6195、1503.01441、1310.2240（Nawata–Ramadevi–Zodinmawia）、1007.1656（Chen–Chen）、hep-th/0108064（Mariño–Vafa）を追加して `just fetch-refs` を実行（題名・内容は原典で確認してから literature.md に書く）。
- **unrefined の構造の整理**：literature.md §1（Mariño の g_R, h_R, M_{RS}, 予想、証拠、弦の解釈）、conventions.md（Mariño の規約表）。
- **refined の既知の構成要素**：literature.md §2（AS の SO(2N) refined CS と large N の自由エネルギー、GW の Kauffman 超多項式、Cherednik–Elliot の composite DAHA 超多項式、NRZ の対称表現の refined Kauffman と 4_1 の unrefined の表）。
- **2018 年の refined SO unknot の検算**【計算】（scripts/rk_so_unknot_check.py）：ノートブック check kauf の refined SO(N) unknot（nakami[λ,μ] = (−1)^{(|μ|−rank μ)} P_{λ/μ} を p_n → (a^{n/2}−a^{−n/2})/(t^{n/2}−t^{−n/2}) で評価し、μ を自己共役な分割で足したもの）を、AS の SO(2N) refined CS の unknot（D_N 型 Macdonald 多項式の主特殊化 M_λ(t^ρ)、a = t^{2N−1}）と比べた。
  - AS (5.5) の t = q^β の積と、一般 t の評価式 t^{−⟨λ,ρ⟩} Π_{α>0} Π_{i=0}^{⟨λ,α⟩−1} (1 − q^i t^{⟨ρ,α⟩+1})/(1 − q^i t^{⟨ρ,α⟩}) は N = 3, 4、β = 2, 3、λ = [1], [2], [1,1] で一致。
  - [1]：ノートブックと AS は一致（N = 5, 6 で記号的に 0）。
  - [2]：差 = (q − t)(q − t^{N−1})/((qt − 1)(q t^{N−1} − 1))。[1,1]：差 = −(q − t)(1 + t + … + t^{2N−3} の因子)/(q t^{2N−3} − 1)。いずれも t = q で消える。
  - したがって、2018 年の試行の「refined SO(N) unknot」は、AS の SO(2N) refined CS の unknot と 2 箱で異なる。試行の入力（AS の trefoil (6.15) など）が AS の理論なら、unknot の正規化が整合していなかったことになる（ノートブック全体の読解で確認する）。
- **unrefined の再現**【計算】（scripts/rk_unrefined_check.py、sympy で厳密）：
  - trefoil：colored HOMFLY は src/rcsknots/homfly.py の Rosso–Jones（KN 変数 a = ν², Q = t、sign = −1）、colored Kauffman [2] は NRZ 1310.2240 (6.11) の t = −1、[1,1] は対称性 F_{[1,1]}(ν,t) = F_{[2]}(−ν, t^{−1})、composite ([1],[1]) は Mariño (4.4)。\hat g_[1,1] は Mariño (4.8) と全係数一致。\hat g_[2] は z^9 の項だけ不一致：Mariño の (4.8) の z^9 の項「−11ν^11 + 78ν^9 − 79ν^7 + 12ν^5」は z^8 の項の複写（奇数 z の項は偶数 ν 冪のはず）で、計算値は z^9(11ν^12 − 78ν^10 + 78ν^8 − 11ν^6)。**(4.8) の z^9 は原論文の誤植と判断**（他の全項は一致）。
  - Mariño の本文の P(3_1) = 2ν² − ν⁴ − z²ν² は量子群の標準値 2ν² − ν⁴ + z²ν²（= ν²(t² + t^{−2}) − ν⁴）と z² の符号が違う。(4.8) の再現には標準値を使った（Mariño の F(3_1) は NRZ の r = 1 と一致）。本文の P の符号は誤植の可能性（未確認）。
  - 4_1：\hat g_[2], \hat g_[1,1] は NRZ (C.10) と全係数一致。\hat h_[2], \hat h_[1,1] の z 展開から読んだ N^{c=0} の表は NRZ の表と完全一致（NRZ (C.9) の表示式は計算値と合わず、表示の誤植と判断）。
  - 2018 年 7 月のノートブック（h_hat_fig8, g_hat fig8）はこの NRZ の計算をなぞったもので、数値も一致する。
- **refined unknot の極**【計算】（scripts/rk_refined_unknot_poles.py）：AS の SO(2N) refined unknot（D_N Macdonald の主特殊化）は [2] で分母 (1 − q t^{N−1})、[1,1] で (1 − q t^{2N−3})、[3] で (1 − q² t^{2N−3}) などの **N に依存する分母**をもち、a = t^{2N−1} の Laurent 多項式にならない（t = q では約分される）。U(N) の composite ([1],[1]) = (2,1^{N−2}) の A 型 Macdonald 次元も分母 (1 − q t^{N−1})（a = t^N で (t − qa)）をもつ。2018 年 9 月 30 日の composite for adjoint ノートブックが見つけた「真の refined composite unknot の分母 (t − qa)」はこれと一致【補助エージェントの読解＋本計算で確認】。
- **AS §7 の閉弦の構造から読める示唆**【導出】：AS の自由エネルギーの oriented 部分 Σ_d (1/d) Q^d q^{−d/2}/((q^{d/2} − q^{−d/2})(t^d − t^{−d})) は、分母が (t^d − t^{−d}) = ((t²)^{d/2} − (t²)^{−d/2}) なので「(q, t²) の refined conifold」の形をしている（½ は現れない）。½F_cov ＋ unoriented の形に書き直すと、偶数 d の引き算に (1 + q^d)/(1 + t^d) の因子が付く（単純な d → 2d の置換ではない）。2018 年の試行は Σ_{d 奇数} と ½ を unrefined のまま置いたので、この構造を含んでいない。

## 2026-10-14（K2、途中）
- **閉じた式**【計算】（src/rcsknots/kauffman_refined.py）：AS の SO(2N) refined unknot dim_λ と計量 g_λ（D_N、AS (B.4)）を U = t^N の有理式として実装（dim は N = 5, 6 の積と一致を確認）。U(M) の composite [λ,μ] の Macdonald 次元（対称正規化）を、長方形の部分の積を列の対で伸縮させて A = t^M の有理式にし、M = 8, 9 の直接計算と 13 例で一致。
  - dim_[2] の分母 (Uq − t)、dim_[1,1] の分母 (U²q − t³)；g_[2] は (Uq − t)²、g_[1,1] は (U²q − t³)² を分母にもつ。composite ([1],[1]) の分母 (Aq − t)、([2],[1]) は (Aq² − t)、([1],[1,1]) は (Aq − t²)。
  - 訂正（2026-10-14 夕）：[3] の分母は (Uq² − t)（q²t^{N−1} = 1）であって (1 − q²t^{2N−3}) ではない。[2,1] は (U²q² − t³)、[1,1,1] は (U²q − t⁵)、[4] は (Uq² − t)(Uq³ − t)、[2,2] は (Uq − t²)(U²q² − t³)(U²q³ − t³)。
- **極の打ち消しの不可能性**【導出】：Z_G = Σ_λ dim_λ w_λ B_λ(v) で w_λ が U に依らなければ、2 次の [2] 成分の極 (Uq − t) と [1,1] 成分の極 (U²q − t³) は、composite ([1],[1]) の極 (Aq − t) 一つ（A は一つの値）では同時に打ち消せない。w_λ = 1/g_λ でも、(Uq² − t²), (Uq − t³) が composite の 2 箱の項には現れない。したがって「D_N の次元 × U に依らない重み」と「U(M) の composite 次元」の組み合わせでは、unknot の 2 箱で Laurent 多項式にならない。
- **積の形の検査**【計算】（scripts/rk_k2_littlewood.py、有理点 2 点）：w ∈ {1, b^A_λ, (q/t)^{|λ|/2} b^A_λ, 1/g_λ, b^A_λ/g_λ}、基底 ∈ {Schur, P_λ(v;q,t), P_{λ^T}(−v;t,q)} の 15 通りすべてで、log Z_G の 3 次の p_1³, p_1p_2 の係数が 0 でない（積の形でない）。t = q では w = 1・Schur で 0（Littlewood の恒等式）になることを確認（検査の妥当性）。
- **文献**：Rains–Warnaar 1506.02755（Memoirs AMS 270 (2021) 1317）と Lee–Rains–Warnaar 2007.03174（SIGMA 16 (2020) 142）を取得。
  - (q,t)-Littlewood 恒等式（Macdonald の 4 つと Kawanaka の 1 つ）、有界版（長方形の D_n, B_n, C_n Macdonald 多項式を A 型 Macdonald で展開）。
  - 辞書：P^{(D_n,D_n)}_λ(x;q,t) = P^{(B_n,B_n)}_λ(x;q,t,1)（ℓ(λ) < n）= Koornwinder K_λ(x;q,t; −1, −q^{1/2}, 1, q^{1/2})。AS の D_N Macdonald はこれ。
  - LRW の GL(2n) ⊃ O(2n) 型の (q,t)-分岐則（長方形の λ）は K_μ(x;q,t; 1, −1, t^{1/2}, −t^{1/2}) を使う：**AS の D_N（±1, ±q^{1/2}）と境界パラメータの q ↔ t が違う**（t = q で一致）。
  - 普遍指標の (q,t)-版は lifted Koornwinder 多項式 \tilde K_λ(x;q,t,T;t̄)（Rains 2005）。T = t^n で K_λ に戻る。T に極をもつ（Eq (2.?)：t_0t_1t_2t_3 = q^{1−â(s)} t^{l̂(s)+1} T^{−2}）。Cauchy 恒等式は双対型（Σ (−1)^{|λ|} \hat K_λ(x;t,q,T) \tilde K_{λ'}(y;q,t,T) = Π(1 − x_iy_j)）。
- 次：Koornwinder 多項式の主特殊化（評価式）を両パラメータ（±1,±q^{1/2} と ±1,±t^{1/2}）で求め、AS の D_N の式を再現することを確認した上で、極の構造を composite 側と比べる。
- **Koornwinder の主特殊化**【計算】（kauffman_refined.koornwinder_eval、van Diejen–Sahi の式）：パラメータ (1, −1, q^{1/2}, −q^{1/2})・評価点 (1, t, …, t^{n−1}) で AS の D_N の unknot と 6 つの λ で一致（辞書 P^{(D_n,D_n)} = K(±1, ±q^{1/2}) と実装の両方を確認）。LRW 型 (1, −1, t^{1/2}, −t^{1/2}) では [1^k] は極を持たず、[2] は qQ = 1、[3] は q³Q = 1、[2,1] は qQ = t に極（Q = U²/t）。
- **2 箱の極の打ち消し**【計算、scripts/rk_k2_poles.py】：composite ([1],[1]) の寄与は v の p_1² 方向だけにあるので、打ち消しには Kauffman 側の [2] と [1,1] が同じ極を残数の比 f^{(1,1)}_{(1),(1)} でもつ必要がある。D_N では極の位置が違い（U = t/q と U² = t³/q）、LRW 型では [1,1] に極が無い。重み 1/g_λ（D_N の計量）と A = U/t では (Uq − t²) が共有されるが、(Uq² − t²), (Uq − t³) が残り、これらは composite の 3 箱の極の位置に一致する（次数が合わない）。**結論：AS の D_N 理論でも LRW 型でも、U(M) の composite の Macdonald 次元を v の積の基底で組み合わせる限り、unknot の 2 箱で Laurent 多項式にならない**【導出＋計算】。
- K2 の停止点に到達（研究計画 §2）。

## 2026-10-14（夕）deep research の報告の読解
- ユーザーが deep research で得た文献調査の報告を docs/refined_kauffman/literature_deep_research_2026-10-14.md に保存（外部の要約、未検証）。要点：refined かつ開弦の orientifold の整数性を定式化した文献は無い（AS 自身が第一量子化は未定義と明言）。U(N) では KN（1703.05408）と Cheng–Sułkowski（2104.00713）が分母を Macdonald 基底・計量・refined holonomy に吸収する方式を確立。Avetisyan–Mkrtchyan（2205.12832、閉弦）はクロスキャップの数ごとに Kähler 変数の量子シフトが違う「三重化」を提案し、Klein 瓶の項の分母 (t^n + t^{−n}) を与える。Cherednik–Danilenko 1509.08351 は B/C/D の DAHA-hyperpolynomial の存在と対称性を未解決としている。Koornwinder 型（C^∨C_n）の refined CS は文献に無い。
- **D_N の unknot の N 依存の分母は Koornwinder の評価式の C⁺ 因子（lifted Koornwinder の T の極）に由来する**【計算】：C⁺_λ(t₀t₁t₂t₃ t^{2n−2}/q) = Π_{(i,j)∈λ}(1 − q^{λ_i+j−1} t^{2n−i−λ'_j})（t₀t₁t₂t₃ = q）の因子のうち分子 (t^n, t₀t_r t^{n−1}; q,t)_λ で約分されないものが、[1]〜[4], [2,1], [1,1,1], [2,2] の分母と一致する（上の訂正を含む）。報告の「Q3 の推測」は [3] の分母の誤記（本セッションの以前の記録に由来）を除いて正しい。
- 報告の推奨のうち、本研究で採用する順：(1) 整数性の判定基準を AS (3.13) 型の第二量子化（積表示）に移す、(2) Kauffman 側の母関数の基底を lifted Koornwinder（T を Kähler 変数）にして、Rains の lifted Cauchy 恒等式 Σ_λ (−1)^{|λ|} \hat K_λ(x;t,q,T) \tilde K_{λ'}(y;q,t,T) = Π(1 − x_i y_j) の上で組み直す、(3) クロスキャップごとの Kähler 変数のずれ（Avetisyan–Mkrtchyan）を開弦に持ち込む。いずれも【仮説】として研究計画 K2' に追加する。

## 2026-10-14（夜）K2'
- Koornwinder 多項式を q-差分作用素の固有関数として実装（src/rcsknots/koornwinder.py、n = 2, 3）。van Diejen–Sahi の評価式と一致、Mimachi の双対 Cauchy 恒等式を (n,m) ≤ (3,2) で確認【計算】。tests/test_rk_koornwinder.py。
- 双対 Cauchy 型の refined orientifold OV 核の定式化：AS 型（ホロノミー側）と LRW 型（probe 側）のパラメータが両側に現れること、unknot が第二量子化の積表示を自動的に満たすこと、Hopf 絡み目の明示式（scripts/rk_k2p_hopf.py、n = 2, 3 で直接の和と一致）。y^n z^n Z は n = 2, 3 で既約（t = q でも既約）。report.md §4–5 に整理。
- 停止点：K2' の判定（unknot の 2 箱で整合する定式化）は「双対 Cauchy 型なら積表示として整合、Laurent 多項式性の意味では未定義」。次の方向は report.md §5。

## 2026-10-14（深夜）K2'(a)：U を形式変数とした Hopf 絡み目
- 訂正：Askey–Wilson 漸化式（src/rcsknots/askey_wilson.py）の monic 化が誤っていた（KS (3.1.4) の漸化式は p̃_n = p_n/(a^{−n}(ab,ac,ad;Q)_n) に対するもの）。修正後、固有関数として解いた Koornwinder と n ≤ 4 で一致。昨夜の report.md §4 の「係数は t^n の有理式にならない」はこの誤りによるもので、撤回（§4 に訂正を記した）。
- probe 側の上位係数 c_{2j}(n) は U_n = t^n の有理式（有理補間、次数 (2j,2j)）。閉じた形（report.md §4.1）を 4 つの有理点で j ≤ 5、n ≤ 30 で確認【計算】。記号的な証明は未（【予想】扱い）。
- Hopf 絡み目（probe 各 1 個）を virtual 正規化 Z~ = y^N z^N Z で Q(q,t,U) の上で計算（src/rcsknots/kauffman_hopf.py、scripts/rk_k2p_hopf_U.py、a,b ≤ 4、約 1 分）。U = t^N（N = 3..8）の直接計算と一致【計算】。
- 結果【計算】：f_{10} = −1 − (ν−ν^{−1})/(t^{1/2}−t^{−1/2})、f_{a0} = 0 (a ≥ 2)；f_{11} = −f_{10} z_q(ν~−ν~^{−1})、ν~ = U q^{1/2}/t（t = q で Mariño の \hat g_{[1],[1]} = z(ν−ν^{−1}) と unknot 因子・符号を除いて一致）；f_{1b} = f_{11} χ_{b−1}(U/t)；全 f_{ab} は (1−q)/(1−t)（(偶,偶) は /(1−t²)）× 整数係数 Laurent 多項式。t = q では分母なし。
- テスト：tests/test_rk_askey_wilson.py（3 件）、tests/test_rk_hopf.py（2 件、A = B = 2、N = 3, 4 との照合と f_{10}, f_{11}, f_{12} の式）。
- 停止点：K2'(a) の結果（Hopf 絡み目で整数性が成り立つ）。次の候補：(i) probe 2 個（2 変数 Koornwinder、Mimachi の (n,m) = (N,2)）で 2 箱の表現を分離し Mariño の [2],[1,1] と比較、(ii) trefoil（AS の T 行列）、(iii) 閉じた形 c_{2j}(n) の証明（Rains の二項公式）。
- ユーザーの指示（2026-10-14）：(i) → (ii) の順で進める。

## 2026-10-14（深夜 2）K2'(b)：成分 1 に probe 無限個
- 訂正：K2'(a) の報告で「t = q の 1 変数 probe は Sp(2) 指標で混ざる」と書いたのは誤り。SO(2N) では Λ^k V（k < N）が既約なので、双対 Cauchy から K_{(m)}(y;q,q;t̄) = y^m + y^{−m}（D_1 指標）になり、virtual 正規化では y^k。Askey–Wilson の直接計算（K_{(2)} = y² + y^{−2}）と閉じた形（t = q で c_{2j} = 0）とも整合。report.md §4.1 を訂正。
- S 行列の対称性で成分 1 の和を縮約すると、probe の数によらず 1 変数の \hat K_{(k)}(z) と冪和だけで Z~(P;z) が書ける（多変数 Koornwinder 不要）。|R| ≤ 4、b ≤ 2 で f_{R,b} を Schur 基底と Macdonald 基底で計算（約 20 秒）。
- 結果【計算】：Schur 基底では全 f_{R,b} が (1−q)/(1−t)（b=2 かつ |R| 偶数は /(1−t²)）× 整数係数 Laurent 多項式。Macdonald 基底 P_R(y;t,q) では (1−qt), (1−q²t), (1−qt²) などの分母が出て整数性が破れる。f_{(a),b} = 1 変数の f_{ab}、f_{R,0} = 0（|R| ≥ 2）。report.md §4.2。
- テスト：tests/test_rk_hopf.py::test_multi_probe_consistency（A = B = 2、約 3 秒）。
- 次：(ii) trefoil。AS 1202.4456 §6 の (6.15)（T(2,3)、[1]、unknot で正規化）と T(2,2m+1) の一般式を使う。t = q で AS の値は Mariño の F_{3_1}(z,ν) の鏡像（ν → 1/ν、z → −z）× ν^{−6}（framing 6 単位）に一致【計算】。

- 未解決：f_{ab} の (偶,偶) の分母 (1−t²) の意味（refined LMOV の分母の構造との対応）。正値性の判定基準（refined BPS の抽出式）は未定義。

## 2026-10-14（深夜 3）K3'：trefoil の [1]
- AS 1202.4456 §6 の T(2,2m+1) の一般式を実装（src/rcsknots/kauffman_trefoil.py）。m = 1 で (6.15) に一致。t = q で Mariño の F_{3_1} の鏡像 × ν^{−6}【計算】。
- 古典的 Rosso–Jones（D_N 指標）は t = q で AS の値に係数 1 で一致するが、素朴な refined Rosso–Jones（Koornwinder 展開、T_μ = q^{(μ,μ)/2}t^{(μ,ρ)}）は AS と単項式倍でも一致しない【計算、n = 2, 3】。p_2 の Koornwinder 展開：K_{(1,1)} の係数は A 型と同じ (1+q)(1−t)/(1−qt)、K_{(1)} の係数 0、定数項は N 依存（n = 2 で −a11²）。
- refined の h/g 分解（仮説）：g^{ref}_[1] = W_SO − ½h^{ref}（½h^{ref} は DGR の超多項式の鏡像に 𝐚_U² = ν²(t/q)^{1/2} を入れたもの × (ν−ν^{−1})/(𝐪−𝐪^{−1})）。t = 1 の極の打ち消しから (6.15) に掛ける単項式は r^{−5} s^{d} U^{−6}（d 奇数）に限られ、d = 11（T_[1]^{−6}(q/t)^{1/2}）で最短の 11 項。整数係数。t = q で (ν²−ν⁴)(𝐪²+1+𝐪^{−2}) + ν⁶（z の偶数冪のみ、N^{c=2}_{[1]} = 0 と整合）。report.md §6。
- テスト：tests/test_rk_trefoil.py（3 件、約 1 秒）。
- 停止点：trefoil [1] の整数性は仮説の単項式のもとで成立。未解決：余分な (q/t)^{1/2} の由来、refined Rosso–Jones が D 型で AS と合わない原因、1^k（k ≥ 2）の refined SO 不変量の計算法（SL(2,Z)、有限レベル）。

## 2026-10-14（深夜 4）K4'：refined のクロスキャップ展開
- ユーザーの指示：refined のクロスキャップ展開（c = 1, 2）を unknot・Hopf で定義し、正値性と DGR 変数での正整数性を見る項目を計画に追加（research_plan.md K4'）。効果の高い順に (1) unknot・Hopf [1],[1]、(2) trefoil [1]、(3) Hopf (a),(b)、(4) 結び目の 1^k。
- 被覆（composite 表現の U(N) refined Hopf）を実装（kauffman_cover.py）。Pieri 係数の形を数値で特定、記号版（A = ν²）は有限 N と一致。unknot・[1],[1] の解析式と一致【計算】。
- Hopf の g_{ab} = f_{ab} − ½f^{cov}_{ab} は a,b ≤ 4 で ν の Laurent 多項式（有理補間）【計算】。
- 記号計算：(1,b) は g_{1b} = z_q(ν~−ν~^{−1})χ_{b−1}(ν~/q^{1/2})、t は ν~ 経由のみ（refined クロスキャップ = unrefined with ν → ν~）。(2,3),(3,3) などは Laurent・整数だが t 依存。(2,2),(2,4),(4,4) は 2(1+t), 4(1+t)(1+t²) の分母が残る（AS §7 の (t^{d/2}+t^{−d/2}) の開弦版）。閉弦の (1+q^d)/(1+t^d) 規則では解消しない。report.md §7。
- クロスキャップ展開（G = z_q z_t の冪 + 奇数因子 X）を unknot・Hopf [1],[1]・trefoil [1] に適用（scripts/rk_k4_crosscap.py）：unknot は c=1 の −1、Hopf [1],[1] は X = z_q の c=2（±1）、trefoil の奇数部分は G 環に入らない（仮説依存）。
- テスト：tests/test_rk_cover.py（5 件、約 6 秒）。
- 未解決：refined の偶数 d の引き算の正しい形（(2,2) の留数を記録）、trefoil の単項式、結び目の 1^k。

## 2026-10-17　K4'(a)：refined の偶数 d の構造
- AS 1202.4456 §3（M 理論の軌道体 z₂ → −z₂、oriented 状態は偶数モード Π_{l}(1 − q^{m₁+l₁}t^{m₂+2l₂}Q)、unoriented 状態は z₂ = 0 に凍結 Π_{l₁}(1 − q^{m₁+l₁}t^{m₂}Q)）と §7（SO(2N) の large N：F = Σ_d (1/d)[Q^d q^{−d/2}/((q^{d/2}−q^{−d/2})(t^d−t^{−d})) + Q^{d/2}t^{d/2}q^{−d/2}/(q^{d/2}−q^{−d/2})]）を読み直した。
- 【導出】開弦（probe は t-brane、開弦の分母 (1−t)）では z₂ → −z₂ が probe に沿うので、全 d の自由エネルギーの係数は f_{ab} = ½f^{cov}_{ab} + ½E_{ab}/(1+t) + C_{ab}（E：Ĩ 固定の oriented 状態、(偶,偶) のみ；C：凍結した unoriented 状態）。閉弦では E = (1−t)f^{cov}（全状態が固定）で AS §7 の (t^{d/2}+t^{−d/2}) の項を再現する。
- 【計算】（scripts/rk_k4_hopf_alld.py、a,b ≤ 4）全 d の g_{ab} = f_{ab} − ½f^{cov}_{ab} = g^{odd}_{ab} − ½g^{odd}_{a/2,b/2}(x²) の分母は (偶,偶) で 2(1+t) のみ（(4,4) に (1+t²) は出ない）、それ以外は整数 Laurent。2(1+t)g_{ab} は全セクターで整数 Laurent。unknot の (2,0) = −½ は全 d では 0。E_{22}(−1) = z_q(z_q^+)²(ν~−ν~^{−1})(q^{−1/2}ν~ − q^{1/2}ν~^{−1})、E_{24}(−1) = E_{22}(−1)(q^{−1}ν~² + qν~^{−2})。
- 【未解決】E と C の分離。候補 (A) C = ĝ^{unref}(ν→ν~)、(B) E = 積の形、はどちらも整数 Laurent だが (4,4) の留数を予測しない。決めるには被覆の composite セクター別の Ĩ の作用（Koike の composite 指標のプレシズム）か、結び目のデータが要る。
- 先に試した閉弦の規則（引き算 × (1+q^d)/(1+t^d)）が失敗した理由：それは E_{22} = −(1+q)g_{11}(x²) を仮定するが、実際の E_{22}(−1) とは因子 (q − ν~²)/(q(1+ν~²)) だけ違う（開弦では固定される状態が倍化したクロスキャップと一致しない）。
- テスト：tests/test_rk_cover.py::test_hopf_crosscap_all_d_even_sector（記号計算、約 1 秒）。docs：report.md §7、conventions.md、research_plan.md を更新。
- 次：(b) trefoil の 1^k（SL(2,Z) の有限レベル実装）。(a) の結論により、結び目の判定基準は「全 d の f_a が (1−t^{e_a})f_a 整数 Laurent（e_a = 2 for a 偶数、1 for a 奇数）、g_a = f_a − ½f^{cov}_a が奇数 a で整数 Laurent、偶数 a で 2(1+t)g_a 整数 Laurent」。

## 2026-10-17　(b) trefoil の 1^k：方法の検討
- AS 1105.5117 §7.1：SU(N) の trefoil でも超多項式との比に √(q/t) が余る（AS 自身が「open flat coordinate のずれ」と注記）。§6 の仮説の単項式 T^{−6}(q/t)^{1/2} と同じ構造【文献】。
- conventions.md §7 の記録：A 型でも素朴な refined Rosso–Jones（q,t 不変の Adams × T^{n/m}）は不成立。K3' の D 型の不一致は普遍的な現象で、A 型の正しい構成は GN の楕円 Hall 代数【確認】。
- 【計算】SO(4) = SU(2)²、SO(6) = SU(4) で GN により AS (6.15) を厳密に再現（単項式 q^{−1}、q^{−1}t^{−4}）。(6.15)(N=2) = [SU(N) 式(N=2) の鏡像 × unknot]² も記号的に確認。tests/test_rk_trefoil.py に追加（4 件）。
- GN の正規化（conventions.md (N2')）では det との積が q^{2mn|λ|}（U(1) の framing、N に依らない）だけずれ、双対表現も単項式倍で一致【計算、N = 4, 5】。composite（λ − 1^N）の trefoil を GN で作るときはこの規則で補正する。K_{111}|D₃ = P_[2] + P_[2,2,2]|A₃（10 ⊕ 10̄）【計算】。
- GN を sympy の field（q,t 記号）で動かすと [2,1,1] で HeuristicGCDFailed。記号的には点ごとの厳密計算＋補間が要る。
- **(b) の実行可能性の結論**：1^k（k ≥ 2）の ν 依存性（整数性・正値性の判定に必要）には D_N の SL(2,Z)（有限レベル、可積分表現すべて・spinor 含む）か DAHA（一般の t では nonsymmetric Koornwinder の機械が必要）を N ≥ 6 程度まで動かす必要があり、本セッションの範囲を超える。N = 2, 3 の厳密値（GN）で固定 N の構造（全 d の分母、framing の走査）を見ることは数時間で可能だが、ν の構造は見えない。停止点として報告。

## 2026-10-18　(a) の残り：E と C の分離（ユーザーの指示「4へ進め」）
- 被覆をセクター別に分けた（kauffman_cover.cover_f_sym4：probe を x 側 y1 と 1/x 側 y2 に、成分 2 も同様。A = B = 4 は約 6 分、scripts/rk_k4_sector_cover.py、data/refined_kauffman/hopf_fcov_sector_A{2,4}.pkl）。セクター和 = 1 変数の f^{cov}、Ĩ 対称性を確認【計算】。
- 最初の仮説「固定セクターの全状態が固有値 +1」（閉弦と同じ E = (1−t)f^{cov}_sec）は t = −1 の条件を満たさず否定【計算】。
- E = P − 2Q（Q は固定セクターの状態の部分和）＋ t = −1 の条件：(2,2) で 16 解、(2,4) で 590 万解、(4,4) で ~10^{110}（いずれも矛盾は無い：t = −1 の留数は状態の台の上にある）。E_24 = E_22χ_1 の倍化の仮定は 16 解すべてと両立（選別しない）。
- t = q で C = ĝ（Mariño）を加えて整数線形計画（scripts/rk_k4_fixed_sector.py、scipy の milp、約 10 分）：**(2,2) と (2,4) は一意**で積の形（E_22 = 候補 B、E_24 = E_22χ_1(ν~²/q)）。**(4,4) は不定**（164/180 変数が自由、Σ|E| 最小でも 74 が自由、倍化の仮定 4 種は不成立）。data/refined_kauffman/hopf_E_fixed.pkl（(4,4) は MILP の一つの解で、一意ではないので使わないこと）。
- テスト：tests/test_rk_cover.py に 2 件（セクター別の被覆と P_22 の閉じた形、E_22 の一意性：約 15 秒）。
- 停止点：(4,4) の E を決めるには Ĩ の作用を表現論的に（composite の (λ,λ) 型セクターでの swap の固有値、Koike の composite 指標のプレシズム）決める必要がある。次の候補：(i) その表現論的な構成、(ii) (b) の選択肢 1（N = 2, 3 の trefoil 1^k）、(iii) 区切りとして push/PR（許可が必要）。

## 2026-10-18（続き）選択肢 1：Ĩ の作用の表現論的な検査（ユーザーの指示「1へ進め」）
- 成分 1 に probe 無限個の被覆（kauffman_cover.cover_f_multi）を実装。1 変数への特殊化を確認（tests/test_rk_cover.py::test_multi_probe_cover_and_adams_b2）。
- 【計算】Adams の予言（E_{R,b} = Σ c_{2;ρ}^R ε_{ρ,b}）が |R| ≤ 4、b ≤ 4 で成立。ε_{(1,1),(1)} は一意で積の形。ρ ⊢ 2、b = 2 の規則 P = P_22 X_ρ²、ε = ±E_22 ψ²X_ρ。
- 【計算】b = 4 では ε_{(2),(2)}、ε_{(1,1),(2)} とも一意でない。Adams は E_44 を決める情報を加えない【導出】。
- 実行の記録：cover_f_multi(4,4) 542 秒、SO multi(4,4) 757 秒（data/refined_kauffman/hopf_fcov_multi_A4B4.pkl、hopf_f_multi_A4B4.pkl に保存、logs/rk_k4_adams_b4.log）。
  **見積もりの誤り**：計算部分だけを測り、後処理（全 59 鍵の g を sympy の式の cancel で約分）を測らなかったため、後処理が 10 時間止まった（CPU 95%、出力なし）。プロセスを停止し、必要な 11 鍵だけを有理関数体 Q(r,s,ν) で処理するように書き換えて 1 分で完了。今後、後処理も小さい例で時間を測ってから実行する。
- 停止点：E_44 の決定には、Ĩ の単一粒子 BPS 状態への作用を直接構成する必要がある（Hopf のデータ・Adams・t = q の条件では足りない）。次の候補：(i) N = 2, 3 の trefoil 1^k（GN）、(ii) 区切りとして push/PR（許可が必要）。

## 2026-10-18（続き 2）(b) trefoil の 1^k を N = 2, 3 で（ユーザーの指示「2へ進め」）
- 見積もり：GN の再構成は [2,1,1] 1 秒、[2,2,2] 8 秒（data/refined_kauffman/torus、scripts/rk_k5_trefoil_colors.py）。固定 N の自由エネルギーは N = 3、a ≤ 6 で 34 秒。
- 【計算】data は DGR の向き（AS の鏡像）。A = a^{1/2} = s^N（最初 s^{2N} と誤って代入し、極の検査が全部失敗した。修正済み）。
- 【導出】固定 N では Z~ の Laurent 性と f の整数 Laurent 性が同値で、内容は N 依存の極の打ち消し ⇔ reduced 不変量が極の上で単項式倍。
- 【計算】SO(4)（[2] が q t = 1 で単項式）、SO(6)（[2,1,1] が q t³ = 1 で単項式、[1,1]・[2]・[2,2,2] が q t = 1 で互いに単項式倍）。条件を満たす framing で f_a は整数 Laurent（N = 2 は a ≤ 4、N = 3 は a ≤ 6、data/refined_kauffman/trefoil_fixedN_f.pkl）。m_1 = K3' の仮説の値。
- 走査の初版は sympy の式の展開で 10 分以上かかったので止め、有理関数体に替えて 1 秒（見積もりを先に測る規則を守れていなかった。今回は止めてから測った）。
- テスト：tests/test_rk_trefoil.py::test_trefoil_fixed_N_pole_cancellation（約 1 秒）。
- 結論と停止点：K2' の枠組みは trefoil の 1^k（k ≤ N）で N = 2, 3 の必要条件（極の打ち消し）を満たす。固定 N では framing・ν の構造・正値性は決まらない。次の候補：(i) 一般の N（D_N の SL(2,Z) か DAHA の実装、規模大）、(ii) (4,4) の E（Ĩ の単一粒子状態への作用）、(iii) 区切りとして push/PR（許可が必要）。

## 2026-10-18（続き 3）用語と PR
- 用語：「交差帽」（crosscap の直訳）を、数学で普通に使う「クロスキャップ」に置き換えた（docs/refined_kauffman の .md 全部、deep research の報告を含む）。コードとテストのコメントには「交差帽」が残る（scripts/rk_k4_*.py、tests/test_rk_cover.py）。
- PR #11 を作成（gh api で 27 commit を同じ SHA で作り、gh pr create）。Auto-fix はアプリ側で有効。
- 文献の確認：1202.4456 の著者は Aganagic–Schaeffer（TeX の \author と謝辞の K.S. で確認）。チャットで一度「Aganagic–Shakirov」と書いた誤りを訂正（docs には誤りなし）。

## 2026-10-18（続き 4）タスク A：T(2,2m+1) の [1]
- PR #11 に 4a1a8a6 を gh api で push（ユーザーの指示「ghでpushしてAへ進め」）。
- 【計算】m = 1..5 で g^{ref}_[1] を Laurent にする単項式を走査（scripts/rk_k3_torus_family.py、1 つの m に数十秒）：符号 +、U^{−2(2m+1)}、r^{−(3m+2)}、s の冪は偶奇だけ。
- 【計算】s の冪は t = q の Kauffman 型の条件（r → −1/r 不変）で d = 7m+4 に決まる。単項式は T_[1]^{−2(2m+1)}(q/t)^{δ/2}（δ = m）で、AS 1105.5117 §7 の SU(N) の因子と同じ形【文献で確認】。K3' の「最短」で選んだ d = 11 はこの条件と一致。
- 【計算】g^{ref} は整数 Laurent、t = q で z の奇数冪なし。クロスキャップ展開：c = 1 は G 環、c = 2 は (w − 1/w) で割れるが商は G 環に入らない。符号は混在。
- 最初、t = q の対称性の判定を expand で書いて全 m で失敗（unknot の因子の分母のため）。cancel に直した。
- テスト：tests/test_rk_trefoil.py::test_torus_family_framing_monomial（約 1 秒）。
- 【計算・予想】refined でだけ現れる部分（𝐚 の奇数冪）は (−1)^{m+1}(𝐚²−𝐪²)(𝐭²−1)·単項式·Π_m、Π_m は係数 +1 の m(m+1)/2 項の閉じた形（m ≤ 3 から推測、m = 4, 5 で確認）。tests/test_rk_trefoil.py::test_torus_family_refined_only_part。
- 次：B（refined のクロスキャップ BPS 数の定義）。A のデータ（T(2,2m+1) の g^{ref}、c = 1 の G 展開、c = 2 の (𝐚²−𝐪²)(𝐭²−1)Π_m の形）を試験台にする。

## 2026-10-18（続き 5）タスク B：refined のクロスキャップ BPS 数の定義
- PR #11 に 9e87c2e, dea8f25 を gh api で push（ユーザーの指示「ghでpushしてBへ進め」）。
- 【計算】T(2,2m+1)（m ≤ 5）の c = 1 部分を ν^β ごとに分類：正の単項式（β = −(4m+2)）、(𝐭²−1)×正（m − 1 個、refined だけ）、負、正。混在は無い。G 展開の前は m ≥ 2 で混在（G 自身の符号の混在のため）。
- 定義と予想を report.md §8、conventions.md に記録（【仮説】・【予想】）。tests/test_rk_trefoil.py::test_refined_crosscap_c1_blocks（約 3 秒）。
- 停止点：判定式の候補はできたが、根拠はトーラス結び目の [1] だけ。次の候補：(i) 偶数セクター（2(1+t)g、E ≠ 0）の Hopf で同じ分類が成り立つか（(2,2),(2,4) の C は決まっている）、(ii) (𝐭²−1) の意味、(iii) 一般の N（C、規模大）。

## 2026-10-18（続き 6）B の続き：(1) 偶数セクターの Hopf、(2) (𝐭² − 1) の意味（ユーザーの指示「1、２へ進め」、push の指示は無し）
- 汎用の分類器 scripts/rk_b_classify.py（式を c = 1/c = 2 に分け、ν^β ごとに coherent／tfactor／mixed）。
- (1)【計算】Hopf：(0,1), (1,1) は分類どおり。2 箱以上は 1 変数 probe でも多 probe の Schur 基底でも G で書けない。t = q でも z の多項式になるのは自己共役な R だけ（q^{κ_R/2} の framing のずれ）。判定式の検査には多 probe の Hopf に refined の framing を入れる必要がある（未実施、数時間規模）。偶数セクターの C, E も同じ理由で検査できない。
- (2)【計算】(𝐚² − 𝐪²)：SO(2) で g^{ref} = (−𝐭)^m = (q/t)^{δ/2}（U(1) の自明さ、§6 の framing 因子そのもの）。(𝐭² − 1)：g^{ref} の 𝐭 の偶奇が m と同じ（g(−𝐭) = (−1)^m g(𝐭)）なので、𝐭 = −1 で消える部分は 𝐭 = ±1 で消える。m = 1..5 で確認。
- テスト：tests/test_rk_trefoil.py::test_refined_crosscap_special_points（約 1 秒）。
- 停止点：普遍的な因子は説明がついた。残る内容は商の符号がそろうこと（正値性）で、根拠はトーラス結び目の [1] だけ。次の候補：(i) 多 probe の Hopf に refined の framing を入れて 2 箱以上のセクターで判定式を検査（数時間）、(ii) 一般の N（C）。

## 2026-10-18（続き 7）Hopf の hatted 不変量（ユーザーの指示「ghでpushして1へ進め」）
- PR #11 に 10afd1c, e287486 を gh api で push。
- 【訂正】前回の「2 箱以上は framing のずれで検査できない」は誤り。原因は Mariño の reformulation ĝ = M^{−1}g を飛ばしていたこと（0904.1088 で確認）。M は KN の M(q^{−1})（= Mariño の M_{RS}）。
- 【計算】ĝ_{R,[1]}（|R| ≤ 4）：hook で (−1)^d z_q(ν~−ν~^{−1})ν~^{2d−n+1}、非 hook で 0。t は ν~ だけを通る。判定式どおり（c = 2、符号がそろう）。t = q で Mariño の Hopf の値（R ↔ R^T、符号 −1）。scripts/rk_b_hopf_hatted.py、tests/test_rk_cover.py::test_hopf_hatted_crosscap_invariants。
- 残り：偶数×偶数（E ≠ 0）の hatted 不変量には成分 2 も多 probe の g_{R,S} が要る（被覆 cover_f_multi と SO 側の両方を 2 成分とも多 probe にする。見積もりは未）。

## 2026-10-18（続き 8）見積もり：偶数×偶数の hatted 不変量（両成分とも多 probe）
- PR #11 に be2183d を gh api で push（ユーザーの指示「ghでpushして見積もりを出せ」）。
- 目的：E ≠ 0 のセクターで ĝ_{R,S}（|R|,|S| 偶数）を作り、判定式（report.md §8）を検査する。hatted の基底での符号のそろい方が E_44 を絞る可能性もある。
- 必要なもの：成分 2 も多 probe にした g_{R,S}（S = [2], [1,1] の両方。1 変数 probe では g_{R,[1,1]} が見えない）。
- 測定：既存の Koornwinder の解法（記号の x で q 差分作用素を展開）は n = 4 の K_(1,1) で 532 秒、K_(2) は止めた。n について補間する用途には使えない。
- 段階と見積もり（|R|, |S| ≤ 2 の最小版）：
  | 段階 | 内容 | 開発 | 計算 | 危険 |
  |---|---|---|---|---|
  | S1 | SO 側：成分 2 の核の 2 次（½p₁(z)²m₁(x)² − ½p₂(z)m₂(x)）を D_N の K_(1), K_(1,1), K_(2) で展開する係数を、U = t^N の閉じた形で求める（数値の x で作用素の行列 d_{μν} を当てはめ、n = 2..12 から U について補間）。成分 1 は x = t^ρ q^λ での exp の式 | 0.5〜1 日 | 数分〜1 時間 | 中（U の次数が高いと n を増やす必要） |
  | S2 | 被覆側：成分 2 の 2 変数（brane と image）で次数 ≤ 2 の composite の積（e_k e_l(x)、e_k(x)e_l(1/x)、k + l ≤ 2）を A の記号のまま Pieri で展開 | 0.5 日 | 数秒 | 低 |
  | S3 | 組み立て：両成分の plethystic 反転、g = f − ½f^cov、ĝ = M^{−1}(q^{−1})⊗M^{−1}(q^{−1})、固定セクター（(ρ,ρ)×(σ,σ)）の E と C の分離、分類 | 0.5 日 | 数分 | 中（E の分離が一意でない可能性） |
  合計：開発 1.5〜2.5 日、計算は 1 時間以内。得られるもの：ĝ_{[2],[2]}, ĝ_{[2],[1,1]}, ĝ_{[1,1],[1,1]} と hatted の E（E ≠ 0 のセクターでの判定式の最初の検査）。
- 拡張（|R|, |S| ≤ 4、E_44 を hatted の基底で）：K_λ（|λ| ≤ 4、12 の形）の閉じた形と、被覆の次数 4 の composite（Pieri の反復）が要る。開発 3〜5 日追加、危険大。計算は数時間の見込み（実行前に再開・分割の準備を書く）。
- 位置づけ：判定式の検査の範囲を E ≠ 0 に広げるもので有効性は中。主要な制約（refined SO の結び目の入力がトーラス結び目の [1] しかない）は解消しない（それは C：一般の N の 1^k）。

## 2026-10-18（続き 9）最小版：偶数×偶数の hatted 不変量（ユーザーの指示「ghでpushして最小版へ進め」）
- PR #11 に 0904903 を gh api で push。
- 実装：src/rcsknots/kauffman_hopf2.py（SO 側 so_f2：K_(1), K_(1,1), K_(2) の展開を評価の対称性で U の形式変数のまま；被覆側 cover_f2：成分 2 の brane・image 次数 ≤ 2）。計算は両方とも数秒（見積もりの「数分〜1 時間」より速い：遅い Koornwinder の解法を使わずに済んだ）。開発は見積もり（1.5〜2.5 日）より短い。
- 【計算】成分 2 を 1 変数にすると既存の値に一致、g_{R,S} = g_{S,R}。Adams の予言（E_{RS} = ±ε）で偶数×偶数の C が Laurent。hatted の Ê は ∓z_q(ν~−ν~^{−1})(tq^{−1/2}ν~ + q^{1/2}ν~^{−1})、非対角 Ĉ は ν~ だけの関数、対角 Ĉ の refined だけの部分は (𝐚²−𝐪²)(𝐭²−1) × 冪ごとに符号のそろう多項式。
- テスト：tests/test_rk_cover.py::test_two_multi_probe_hopf_even_sector（約 7 秒）。データ：data/refined_kauffman/hopf_f_two_multi.pkl, hopf_fcov_two_multi.pkl, hopf_g_two_multi.pkl。
- 判定式の述べ直し（report.md §8）：「t = q の部分 ＋ 普遍的な因子 (𝐚²−𝐪²)(𝐭²−1) または (𝐭²−1) × 冪ごとに符号のそろう多項式」。
- 次の候補：(i) 拡張版（|R|,|S| ≤ 4、E_44 の hatted 版）、(ii) 一般の N（C）。

## 2026-10-18（続き 10）拡張版（|R|, |S| ≤ 4）：S1 の準備（実行前に記録）
- PR #11 に 3de70c5 を gh api で push（ユーザーの指示「ghでpushして1へ進め」）。
- S1（D_N の K_λ、|λ| ≤ 4、U 形式）：評価の対称性で未定係数を決める（kauffman_hopf2.koornwinder_general）。x → −x の偶奇で未知数を減らした（K_λ には |μ| ≡ |λ| の S_μ だけ）。
- 測定（1 コア）：|λ| ≤ 3 は数秒。sympy の体の Gauss の消去では (2,2)（未知数 5）が 465 秒、flint の Bareiss（分母を払って Cramer）でも 316 秒。(3,1)（未知数 6）は 1 コアで 30 分を超えても終わらない。最初の素朴な版（偶奇を使わない）は nmax = 4 で 30 分以内に終わらなかった（停止）。
- 見積もり：行列式は Cramer の n + 1 個が独立なので 8 並列（コア 8、メモリ 8 GB、使用メモリは数十 MB）。(2,2) 約 1 分、(3,1) 10〜30 分、(4)（7×7）1〜3 時間。合計 1〜4 時間（24 時間の上限より十分小さい）。
- 再開：λ ごとに data/refined_kauffman/koornwinder_DN/K_<λ>.pkl に保存し、再開時は読み込んで飛ばす。
- 実行：scripts/rk_b_koornwinder_DN.py を nohup で独立に実行、ログ logs/rk_b_koornwinder_DN.log、完了の目印は最後の行「保存 完了」。検算：最小版の係数と一致、全 66 組の評価の対称性。
- 並行して S2（被覆側の次数 4）を実装する。

## 2026-10-18（続き 11）拡張版の途中経過
- S1：K_λ を |λ| ≤ 3 まで数秒、(2,2) 119 秒、(3,1) 1063 秒（8 並列）。(4)（7×7 の行列式 8 個）は 1.5 時間を超えて計算中（data/refined_kauffman/koornwinder_DN/ に λ ごとに保存、再開可）。見積もり（1〜4 時間）の範囲内。
- S2：被覆 Σ N N f^{cov}（|R|,|S| ≤ 4）は 3745 秒で完了（data/refined_kauffman/hopf_fcovLR_A4.pkl）。
- S3 の一部（|R| ≤ 4, |S| ≤ 2、K_(4) を使わない）：SO 側 103 秒。Adams の予言は 10 セクターすべてで成立、hatted Ĉ の refined だけの部分は 10/10 で (ν~²−q)(q−t) で割り切れ、符号のそろい方は 8/10（転置の対 2 成分で破れ）。report.md §8。
- テスト：tests/test_rk_cover.py::test_extension_building_blocks_evaluation_symmetry（約 15 秒）、results の印の test_saved_42_sectors_adams_and_hatted_factor（約 9 分）。
- 残り：K_(4) の完了後に |R| = |S| = 4（ε_{ρσ}、ρ,σ ⊢ 2 が未知、E_44 を含む）。

## 2026-10-18（続き 12）拡張版の完了（|R|, |S| ≤ 4）
- S1：K_(4)（7×7 の行列式 8 個、8 並列）は 18435 秒（約 5.1 時間）。見積もり（1〜4 時間）を超えた（行列式の大きさが 1 段で約 9 倍という見込みより伸びた：各プロセスの CPU 時間 158〜172 分）。全 66 組で評価の対称性が成立、最小版と一致。data/refined_kauffman/koornwinder_DN/。
- S3：SO 側 so_f_general(4,4) は 3268 秒。**実行前に見積もりを書かなかった**（(3,3) の 109 秒から軽く見た。鍵の数 49 → 144 で、事後の見積もりは 15〜75 分）。規則の違反として記録する。
- 【計算】g_{R,S} = g_{S,R}、奇数のセクターは Laurent、1 変数の値が hopf_g_alld（別の経路）と一致。4 箱×4 箱の 25 組の留数が Adams（4 つの ε_{ρσ}(−1)）で整合、ε_{(2),(2)}(−1) = E_44(−1)（1 変数）。data/refined_kauffman/hopf_g_two_multi_44.pkl、hopf_eps44_m1.pkl。
- E_44 は依然として t = −1・t = q 以外では決まらない。テスト：results の印の test_saved_44_block_symmetry_and_adams_at_minus1。
- 次の候補：(i) 4 方向の被覆（固定セクター）で ε_{ρσ} を部分和として絞る（被覆の (4,4) を 4 方向でやると鍵が 38×38、数時間〜の見込み、実行前に見積もる）、(ii) 一般の N（C）。

## 2026-10-18（続き 13）固定 N の 2 箱（ユーザーの指示「Nを固定して2箱を小さく検証せよ。」）
- 見積もり：1 回の自由エネルギーは約 1 秒、m_(2) の走査（±、|α| ≤ 12、|β| ≤ 50、s を 2 点に固定）は N ごとに約 1 分。全体で 127 秒（logs/refined_kauffman/trefoil_fixedN_2box.log）。メモリはサンドボックスで `time -l` の一部が失敗して未計測。
- 方法：probe 2 個の Z~(Y) を D_N の K_(1), K_(1,1), K_(2) で展開した。W_λ = m_λ R_λ K_λ(t^ρ)、R は data の鏡像（AS の向き）。
- 【導出】s_(2)(Y) ↔ Λ²V、s_(1,1)(Y) ↔ Sym²V なので、f_(2) は W_(2) を含まない。新しい検査は f_(1,1)（D_N の (2)）。
- 【計算】
  - f_(1), f_(2) は 1 probe の計算と全係数一致した。
  - f_(1,1) を Laurent にする m_(2) は一意だった（N=2：s^8、N=3：r^4 s^20）。SO(4) では m_(1,1)² に等しく、SU(2)² の積と整合する。
  - f はすべて整数係数の Laurent 多項式で、符号は混在する。
  - M^{-1} f（M(q^{±1}), M(t^{±1})）は Laurent でない。p_2 成分が 1+q でも 1+t でも割り切れず、t = q では割り切れる。
  - 保存先：data/refined_kauffman/trefoil_fixedN_2box.pkl。
- 走査の初版（−A0/A1 を極の上で sp.solve）は解が空か nan になり、判定に使えなかった。s を固定した Laurent 性の篩に替えた。
- テスト：tests/test_rk_trefoil.py::test_trefoil_fixed_N_two_probe（N = 2、約 4 秒）。
- 結論：固定 N の 2 箱でも、極の打ち消しと整数性は成立する。ただし固定 N の f だけでは、hatted の検査（M^{-1} の (1+q)）も正値性も判定できない。
- 次の候補：(i) 固定 N の被覆（U(3), U(5) の composite の trefoil、[2,1]・[2,1,1,1] など）で g = f − ½f^{cov} と ĝ を作り、ĝ が Laurent かを見る（GN の再構成が要る。実行前に見積もる）、(ii) 一般の N（C）。

## 2026-10-19（続き 14）固定 N の被覆と hatted（ユーザーの指示「1」：前回の候補 (i)）
- 見積もり：GN の [2,1,1,1]（SU(5) の随伴）は 2.4 秒。被覆の自由エネルギーは N ごとに約 1 秒。走査は N ごとに 1 分程度。いずれも 30 分の目安を大きく下回る。
- 方法：U(N_A)（N_A = 3, 5）の核の |Y| ≤ 2 を Macdonald で展開した（Pieri の c, κ を確認）。荷電共役で W_μ̄ = W_μ と置いた（data の R_μ̄/R_μ が単項式であることを 4 組で確認）。
- 【計算】
  - unknot は g = (−1, 0, 0)。g_(1) は K3' の −g^{ref}_[1] に一致。
  - 被覆の極は随伴だけ。単項式の族のうち t = q で 1 のもの (t/q)^{N_A/2} を採った。
  - **framing の規則を追加**：極に加えて「t = q で 1」。SO 側は m_(1,1) = (t/q)^{N−1} になり、このとき m_(2) は一意に (t/q)^2（t = q で自動的に 1）。前回の m_(1,1) = s^4・s^16 はこの規則を満たさない。conventions.md「固定 N の framing」に記録した。
  - 2g は整数 Laurent で、g_(2) + g_(1,1) は整数。半整数は p_2 方向だけ（Adams の予言と整合）。
  - **hatted**：K3' の m_(1) = (t/q)^{1/2} では、g の p_2 成分が 1+q でも 1+t でも割り切れない。被覆の p_2 成分は割り切れるので、原因は SO 側である。
  - m_(1) を走査すると、t = q で 1 かつ割り切れるのは t/q だけだった。このとき ĝ は z を除いて Laurent になり、SO 側の単項式はすべて (t/q)^{|μ|/2} の一様な形になる。
  - ところが m_(1) = t/q では、大きな N の g^{ref}_[1] に 1/(1 − t) の極が残る。2 つの要請は両立しない（未解決、report.md §6）。
- 検査の欠陥の修正：整数性の判定（scripts/rk_k5_trefoil_2box.py の stats、scripts/rk_k5_trefoil_fixedN.py）が分母の定数（2 など）を見ていなかった。修正後に前回の結論（f は整数）を再確認した。誤りは g の判定で見つかった（g_(2) は半整数）。
- データ：data/refined_kauffman/trefoil_fixedN_2box.pkl（新しい framing で作り直し）、trefoil_fixedN_2box_g.pkl（2 案の g, ĝ）、trefoil_fixedN_f.pkl（1 probe、新しい framing）。ログ：logs/refined_kauffman/trefoil_fixedN_{2box,2box_cover,1probe}.log、gn_T2_3_2111.log。
- テスト：tests/test_rk_trefoil.py::test_trefoil_fixed_N_two_probe（更新）、::test_trefoil_fixed_N_cover_and_hatted（新規、約 3 秒）。
- 次の候補：(i) 3 箱（|R| = 3）で同じ検査をして、m_(1) の 2 案のどちらが残るかを見る。SO 側は D_N の K_λ（|λ| ≤ 3）が既にある。被覆は composite（|α|+|β| ≤ 3）の trefoil が要り、[3,1^{N_A−2}]・[2,2,1^{N_A−3}] などを GN で作る。実行前に見積もる。(ii) T(2,5) の固定 N の 2 箱で同じ緊張が出るか。(iii) 一般の N（C）。

## 2026-10-19（続き 15）3 箱（ユーザーの指示「1」：前回の候補 (i)）
- 見積もり：GN の 5 彩色は計 28 秒（logs/refined_kauffman/gn_T2_3_3box.log）。次数 6 の Macdonald P は 0.2 秒。3 箱の本体は N = 2 で 13 秒、N = 3 で 42 秒（m_(1) の 9 通りの走査を含む）。30 分の目安を大きく下回る。
- 方法：
  - SO 側は型 A の経路（D_3 = A_3、D_2 = A_1 × A_1）にした。大きな N の安定な K_λ に U = t^N を入れると、ℓ(λ) > N の項が下の K に崩れて扱いが微妙になるためである。
  - 被覆は composite の縮約 d_α を主特殊化から求め、分解を数値で確認した。
  - sympy の疎多項式の gcd が次数 6 で失敗する（HeuristicGCDFailed、代わりの方法に切り替えない）ので、flint の fmpz_mpoly の上の小さい有理関数クラス src/rcsknots/flint_rf.py を作った。
- 【計算】
  - 型 A の経路は 2 箱で D_N の経路と全係数一致した。unknot の f_R（|R| = 2, 3）は 0。
  - framing は「t = q で 1」の族で一意だった。被覆の (2;1), (11;1) は (t/q)^{N_A/2}。SO の新しい表現の単項式は m_(1) に連動する。
  - m_(1) = t/q のとき、SO 側の単項式は全部 (t/q)^{|μ|/2} になる（予言どおり）。
  - SO 側の f̂ = M^{-1}(q^{-1}) f（|R| ≤ 3）は、m_(1) = t/q のときだけ整数 Laurent だった（β1 ∈ [−2, 6] を走査、N = 2, 3）。2 箱は D_N の経路でも確認した。M(t) では閉じない。
  - g（|R| = 3）はどの m_(1) でも整数で、奇数のセクターの E = 0 と整合する。
  - ĝ は p_3 方向では m_(1) = t/q のときだけ Laurent。混ざった方向（p_2p_1、p_1³、2 箱の p_1²）は被覆側が原因で閉じず、t = q では閉じる。
- 修正した誤り：
  - 被覆の随伴 (1;1) の分割を (∅;1) と取り違えていた（unknot の検査で発見）。
  - U(3) の [1,1,1] を SU(3) の自明表現（W = 1）に落としていた。被覆の framing に解が無くなったことで発見し、data の U(3) の値に直した。
- データ：data/refined_kauffman/trefoil_fixedN_3box_g.pkl、torus/T2_3_{33,321,31,3111,2211}.json。ログ：logs/refined_kauffman/trefoil_fixedN_3box.log。
- テスト：tests/test_rk_trefoil.py::test_flint_rf_basic、::test_trefoil_fixed_N_three_box（N = 2、計約 2 秒）。
- 結論：固定 N の 3 箱でも m_(1) = t/q が一意に選ばれる。大きな N の [1]（K3'）との食い違いは未解決のまま残る。
- 次の候補：
  - (i) 大きな N の [1] の側を見直す。被覆の [1] の規格化 𝐚_U² = A(t/q)^{1/2} と、AS (6.15) と data の対応を、m_(1) = t/q の前提で調べる。計算は軽い。
  - (ii) 被覆の composite の hatted（M ⊗ M）を定義し、混ざった方向が閉じるかを見る。
  - (iii) T(2,5) の固定 N で同じ選択（t/q）が起きるかを確かめる。
  - (iv) 一般の N（C）。

## 2026-10-19（続き 16）大きな N の [1] の見直しと被覆の composite の hatted（ユーザーの指示「1と２へすすめ」）
- 見積もり：GN の 3 彩色（U(5) の [1,1,1,1]・[1,1,1,1,1]・[2,2,2,2]）は計 32 秒。解析は N = 2 で 3 秒、N = 3 で 7 秒。
- 着眼：被覆の普通の表現の framing は固定 N の極で決まらず、[1] は K3' に合わせて 1 と置いていた。K3' が大きな N で決めていたのは SO 側と被覆側の比だけである。
- 【計算】
  - 一様な規格化 W_μ = (t/q)^{δ|μ|/2} R_μ P_μ(t^ρ) で、共役・行列式のずらしの比が正確に 1 になった。trefoil で 16 組（U(2)〜U(5)）、T(2,5)・T(2,7) で 3 組ずつ。data の行列式表現は (q/t)^{δ n/2}。
  - 一様な規格化では、固定 N の framing の単項式は SO 側・被覆側とも全部 1 になる。走査で決まる単項式は予言と一致した（SO の新しい表現、被覆の (2;1)・(11;1)）。
  - 大きな N の [1]：一様な規格化での g は (t/q)^{1/2} g^{ref}（K3'）で、Laurent のまま。SO 側の単項式はちょうど T_[1]^{−6} で、固定 N = 2, 3 の data と全係数一致した。T(2,5)・T(2,7) でも T_[1]^{−2(2m+1)} で同様（タスク F の解決）。
  - 被覆の 2 変数版（brane に Y、image に Y′）を作った（cover_Z2、plog2）。Y′ = Y に制限すると 1 変数版と一致した。
  - 一様な規格化では、M ⊗ M の f̂_{ρ,σ}（|ρ| + |σ| ≤ 3）、対角の f̂^{cov}、SO の f̂、ĝ が全方向で Laurent になった。前回の混ざった方向の破れは、被覆の framing を 1 と置いたことによる。
  - 整数性：SO の f̂_R・対角の f̂^{cov}_R・ĝ_R（|R| = 3）は整数、ĝ_R（|R| = 2）は 2ĝ が整数（E の ½）。
- 文書：conventions.md に「一様な規格化」を既定として記録した。report.md §6 の K3'・タスク A・食い違いの項に解決を書いた。
- テスト：tests/test_rk_trefoil.py::test_trefoil_uniform_normalization、::test_torus_family_uniform_normalization（計約 1 秒）。
- 次の候補：
  - (i) 一様な規格化で、トーラス結び目の [1] と 2 箱・3 箱の ĝ を G 展開し、正値性の判定式（タスク B）を作り直す。
  - (ii) T(2,5) の固定 N の 2 箱・3 箱で同じ構造（hatted の整数性）を確かめる。data（T2_5 の彩色）は一部ある。
  - (iii) 2 変数の f̂_{R,S} の Schur 基底での整数性。
  - (iv) 一般の N（C）。

## 2026-10-19（続き 17）判定式を一様な規格化で作り直す（ユーザーの指示「１」）
- 見積もり：[1] の族（m ≤ 5）は計 6 秒、固定 N の偶奇の比較は N = 2 で 5 秒、N = 3 で 14 秒（logs/refined_kauffman/b_uniform.log）。
- 【計算】T(2,2m+1) の [1]（m = 1..5）、一様な規格化 g^{true} = (t/q)^{δ/2} g^{ref}：
  - SO(2) の点で g^{true} = 1。𝐭 について偶。
  - c = 2 部分 = −(𝐚² − 𝐪²)(𝐭² − 1) × 単項式 × Π_m で、符号は m によらない。
  - c = 1 の G 展開の分類は以前と同じ。
  - G 展開の係数は w² ↔ w^{−2} で対称でなく、SU(2)_R の指標には分解しない。
- 【計算】固定 N の 2 箱・3 箱は、ν が q, t と混ざるので G 展開できない。代わりに偶奇を調べた。
  - g_R（|R| ≤ 3）は、一様な規格化のときだけ r → −r で不変だった（ほかの 3 通りは 2/6）。
  - ĝ_R の偶奇は (−1)^{|R|−1}（M の正規化から）。
- 文書：report.md §8 に「一様な規格化での判定式」、conventions.md に判定式の述べ直しを書いた。
- テスト：tests/test_rk_trefoil.py::test_crosscap_criterion_uniform_normalization（約 5 秒）。
- 結論：一様な規格化では、判定式のうち規格化に関わる部分（SO(2) の点、𝐭 の偶奇）が単純な形 (1, 偶) になる。正値性の内容（c = 1 の G 展開、c = 2 の普遍的な因子 × 正）はそのまま残る。多箱の結び目の正値性は、大きな N の入力が無いので判定できない。
- 次の候補：
  - (i) 一般の N（C）：多箱の ν 依存性と正値性を見る唯一の道。規模は大きい。
  - (ii) T(2,5) の固定 N の 2 箱・3 箱で、一様な規格化の hatted の整数性と偶奇を確かめる。data の彩色の追加が要る。
  - (iii) 2 変数の f̂_{R,S} の Schur 基底での整数性。

## 2026-10-19（続き 18）T(2,5) の固定 N の 2 箱・3 箱（ユーザーの指示「２」）
- 見積もり：GN の T(2,5) は 1 彩色 2 秒（[2,2], [2,1,1] で計測）から、残り 7 彩色で数分と見積もった。実績は計 2.2 分（最大 [3,1,1,1] 56 秒、logs/refined_kauffman/gn_T2_5_3box.log）。解析は N = 2 で 16 秒、N = 3 で 54 秒。
- 準備：
  - scripts/rk_k5_trefoil_2box.py の data の結び目を KNOT で、scripts/rk_k5_trefoil_3box.py の一様な規格化の δ を DELTA で選べるようにした。
  - so_uniform は核に現れる表現から W^{δ|μ|} を作る形にした（trefoil の値と一致を確認）。
  - 被覆の普通の表現の単項式も上書きできるようにした（mo）。
  - 解析の本体は scripts/rk_k5_torus_fixedN_hatted.py（引数：結び目と δ）。trefoil で回して以前の結果を再現した。
- 【計算】T(2,5)、N = 2, 3、|R| ≤ 3：
  - 一様な規格化の比が 13 組すべてで 1。
  - framing の走査で、全表現の単項式が一意に (t/q)^{δ|μ|/2}。
  - hatted は全方向で Laurent。整数性と偶奇は trefoil と同じ形。
- データ：data/refined_kauffman/torus/T2_5_{22,211,31,222,2111,2211,3111,33,321}.json、data/refined_kauffman/{T2_3,T2_5}_fixedN_3box_g.pkl。ログ：logs/refined_kauffman/fixedN_hatted_{T2_3,T2_5}.log。
- テスト：tests/test_rk_trefoil.py::test_T25_fixed_N_uniform_hatted（約 10 秒、結び目に依存するキャッシュを前後でクリア）。
- 結論：一様な規格化の構造（framing が全部 1、全方向の hatted、偶奇）は δ = 2 でも成り立つ。
- 次の候補：
  - (i) 一般の N（C）。
  - (ii) T(2,7) でも同じ確認。data の彩色の追加で数分〜十数分の見込み。
  - (iii) 2 変数の f̂_{R,S} の Schur 基底での整数性。
  - (iv) 非トーラス結び目（例えば 4_1）の固定 N の入力を A 型の refined 不変量から作れるかの検討。

## 2026-10-19（続き 19）T(2,7) と非トーラス結び目（ユーザーの指示「2,3,1の順で実行せよ」の 2 と 3）
- **2：T(2,7)**
  - 見積もり：GN の 1 彩色 3〜4 秒（[2,2], [3,1] で計測）。残り 7 彩色は数分の見込みで、実績は計 1.8 分（最大 [3,2,1] 31 秒）。解析は N = 2 で 56 秒、N = 3 で 217 秒。
  - 【計算】T(2,5) と同じ項目がすべて成り立った（δ = 3）。logs/refined_kauffman/fixedN_hatted_T2_7.log。
- **3：非トーラス結び目の検討**
  - 入力：長方形の彩色は KNTZ の refined 公式（主プロジェクトの scripts/x1_check.triply）から作れる。trefoil で GN の data と全係数一致した。
  - 名前と K_{m,n} の対応は KNTZ の TeX（arXiv から取得、Cyclotomic-Rectangular.tex :521-531）の図 1 で確認した。
  - 可否の結論：
    - N = 2 の SO 側だけの検査（3 箱まで）は、任意の二重ツイスト結び目でできる。
    - 被覆は 4_1 の 2 箱だけ（[2,1] は KN の候補）。
    - N = 3 と大きな N はできない。Gukov–Walcher の Kauffman superpolynomial は、AS によると refined SO の値と別物である。
  - 【計算】4_1, 5_2, 6_1, 7_2, 7_4, 8_1, 8_3：
    - data のずれの指数 c は [1] の superpolynomial の δ 次数に等しい（0 か 1、種数はすべて 1）。
    - N = 2 の SO 側で、framing の走査が一意に β = c|μ|。f は整数で、SO の f̂ は全方向で整数 Laurent。
    - 4_1 の被覆の 2 箱で、KN の [2,1] の候補は極の打ち消しと hatted の Laurent 性を満たした。
    - 偶奇は (r, s) → (−r, −s) だけで成り立つ（r → −r 単独はトーラス結び目に固有）。
  - データ：data/refined_kauffman/torus/{4_1,5_2,6_1,7_2,7_4,8_1,8_3}_*.json。4_1_21.json は KN の候補で、meta に確立していないと明記した。7_4 の [2,2,2], [1,1,1,1], [3,3] と 8_x の [2,2,2] 以降は作っていない（SMALL=1、N = 2 の検査に不要）。
  - ログ：logs/refined_kauffman/{f41_fixedN,dt_fixedN_a,dt_fixedN_b,dt_fixedN_c}.log。
  - テスト：tests/test_rk_trefoil.py::test_double_twist_offset_and_41_fixed_N（約 5 秒）。

## 2026-10-19（続き 20）1：一般の N（タスク C）の実行前の見積もりと計画【停止点：確認待ち】
- 目的：refined SO(2N) の torus 結び目 T(2,2m+1) の多箱（|R| ≤ 2〜3）の不変量を、ν を形式変数として得る。そうすれば ĝ_R を G 展開し、多箱の正値性（判定式の (c = 1)(c = 2)）を検査できる。固定 N（N = 2, 3）は A 型との同一視で済んでいるが、N ≥ 4（D_4 以上）は A 型にならない。
- 経路の候補：
  - (α) **Cherednik の DAHA-Jones 多項式を D_N で実装する**（1111.6195 の定義は任意の根系で成り立つ）。
    - 同じ論文の Prop. SH-JONES（jones.tex :1056）の Gauss 変換による計算は k_ν = 0, 1（t = 1, q）に限られる。一般の t には、多項式表現（Demazure–Lusztig の T_i、π、X、Y）の上で τ± の作用を組み、γ̂(P_λ) を評価する必要がある。
    - 段階：
      1. A_1, A_2 で一般の t の DAHA-Jones を実装し、trefoil の [1], [2] を GN の data と照合する（1〜2 日）。
      2. D_2, D_3 で本研究の固定 N の data（一様な規格化）と照合する（2〜3 日）。
      3. D_4〜D_6 で λ ∈ {(1), (1,1), (2)} を計算する（計算は数時間〜数日の見込み、6 変数の Laurent 多項式と q, t の有理式）。
      4. U = t^N で補間し、[1] を AS (6.15) と照合する（1〜2 日）。
    - 危険：τ+ の作用（X の Gauss 函数による共役）の多項式表現での扱い、D_N の DAHA と Koornwinder (±1, ±q^{1/2}) の対応の規約、N についての安定性（D 型では未証明、AS の [1] の式が N に関して安定なのは支え）。
  - (β) **有限レベルの refined SL(2,Z)（AS の方法）**：t = q^β（β は正の整数）、q は 1 の冪根でしか定義されない。多くの (レベル, β) で計算して q, t の有理式を再構成する。spinor を含む可積分表現の全体が要る。実装と補間で 1〜2 週間、危険は (α) より大きい。
- 見積もり：(α) で合計 1〜3 週間（実装・検証が主で、計算そのものは各段階で数時間〜数日）。**24 時間を大きく超えるので、CLAUDE.md の規則に従い、実行前に確認を求める。**
- 準備として決めておくこと（承認後）：
  - 段階ごとの停止点：1 と 2 の照合が通らなければ止める。
  - 中間結果の保存：data/refined_kauffman/daha/ に (根系, N, λ, 結び目) ごとに保存し、保存済みは飛ばす。
  - 並列：N と λ ごとに独立に 8 並列で回す。
  - 長時間の実行：nohup で回し、ログと完了の目印を logs/refined_kauffman/daha_*.log に置く。
- 代わりの小さな候補（承認が無い場合）：
  - 4_1 の被覆の 3 箱の極の条件から、未知の [3,1] に一次の条件を課す（主プロジェクトの f41_box4 と同じ問題を、orientifold の側から絞る）。
  - 2 変数の f̂_{R,S} の Schur 基底での整数性。

## 2026-10-19（続き 21）一般の N：DAHA の実装と 2 箱（ユーザーの指示「一般のNへすすめ。ただし計画を立てて実行に時間がかかる(1時間以上)ところで確認せよ。」）
- 方法の決定：Cherednik 1111.6195 の Prop. SH-JONES は t = 1, q に限られるので、一般の t は多項式表現で直接計算する。
  - D_N のベクトル表現 ω_1 は minuscule なので、Y_{ε_1} = π_1 T_{u_1}（l(u_1) = 2N − 2、有限の Demazure–Lusztig だけ）と書ける。
  - τ_−(X_{ε_1}) = q^{1/2} Y X、τ_+(Y_{ε_1}) = q^{−1/2} X Y、他の X_{ε_i}, Y_{ε_i} は T_i の共役。T(2m+1, 2) は γ = τ_+^m τ_−^2。
  - 評価は {H} = (H·1)(q^{−ρ_k})。src/rcsknots/daha_dn.py（D_N と GL_N）。
- 【計算】検証：
  - D_N（N = 2, 3, 4）と GL_N（N = 2, 3, 4）で、2 次関係、組紐関係、交差関係、Y の可換性が成り立つ。
  - 保存済みの D_N の Macdonald 多項式（U = t^N）は、DAHA の Macdonald 作用素 Σ(Y_{ε_i} + Y_{ε_i}^{−1}) の固有関数（N = 2, 3, 4、|λ| ≤ 2）。
  - D_N の DAHA-Jones = T_λ^{−6} × 鏡像（固定 N の一様な規格化の W_λ）（N = 2, 3、λ = (1), (1,1), (2)）。N = 4, 5 の [1] は鏡像の AS (6.15) にちょうど一致する。
  - GL_N の DAHA-Jones = 単項式 × 鏡像（KN の data の一様な規格化）。単項式は q^{−3(b,b)} t^{−6(b,ρ)+2|b|(N−1)}（[1], [2], [1,1] で読み、[3], [2,1], [1,1,1] と随伴 (1;1)（U(3), U(5) の固定 N の値）で予言どおり）。
  - GL_N の π_1 は q^{(ω_{1*}, b)} の ω_{1*} = ω_{N−1} ≡ −ε_N（最初 ε_1 として Y が可換にならなかった）。
- 【計算】一般の N：
  - D_N の DAHA-Jones を N = 2〜9 で計算した（scripts/rk_c_daha_dn.py、N = 9 の 2 箱で各約 3 分、最大 400 MB、data/refined_kauffman/daha/）。
  - reduced な R_λ = T_λ^{−6} × 鏡像(Z)/K_λ(t^ρ) を U = t^N で補間した（scripts/rk_c_interp.py、Newton の差分商）。
    - [1] は U^{−5}〜U^{−2} の 4 項で、N = 2〜5 から N = 6, 7 を予言し、AS (6.15) × T^{−6} に一致する。
    - (1,1), (2) は U^{−10}〜U^{−4} の 7 項で、N = 2〜8 から N = 9 を、N = 3〜9 から N = 2 を予言した。
  - 被覆の随伴を GL_N の DAHA で N_A = 2〜11 で計算し（scripts/rk_c_daha_cover.py、各 11 秒以内）、A^{1/2} の冪 −10〜−4 の 7 項で補間した。N_A = 2〜8 から 9〜11 を、5〜11 から 2〜4 を予言した。
  - 大きな N の 2 箱のクロスキャップ部分 g = f − ½f^{cov}（scripts/rk_c_2box_largeN.py、変数 (r, s, U)、src/rcsknots/flint_rf3.py）：
    - g_(1) = −(t/q)^{1/2} g^{ref}_[1]（K3'）。
    - g_(2), g_(1,1) の分母は 2 × 単項式 × (1 + t) だけ（N 依存の極はすべて消える）。
    - U = t², t³ で固定 N の g と全係数一致。
    - g_(2) + g_(1,1) と (1 + t)(g_(2) − g_(1,1)) は整数 Laurent（Adams の予言：(1 + t) は p_2 方向だけ）。
- 実行時間：どれも数分以内で、1 時間を超えるものは無かった（確認の要る計算は無し）。
- 次：hatted の ĝ（大きな N）、E と C の分離、G 展開と正値性の判定。t = q での Mariño の ĝ_[2], ĝ_[1,1]（0904.1088 :1174-1196）との照合（規約の対応が要る）。

## 2026-10-19（続き 22）大きな N の hatted と、長い計算の見積もり
- 【計算】大きな N の hatted ĝ（2 箱）：p_1² 成分は 2 × 単項式の分母だけ、p_2 成分は 2 × 単項式 × (1 + t)。data/refined_kauffman/daha/largeN_2box_m1_ghat.pkl。
- 【計算】G 展開の対称性（V ↔ 1/V、q ↔ t^{−1}）は 2 箱の c = 1 部分で成り立たない（ĝ_(2) + ĝ_(1,1)、R ↔ R^T の入れ替え、U 固定・ν 固定のいずれでも）。多箱の BPS 数の定義は未解決（report.md §9）。
- 【計算】T(2,5)（m = 2）：D_N の DAHA-Jones は N = 2, 3 で W = T_λ^{−10} × 鏡像(Z) の関係を満たし、一様な規格化の data と一致した。
- 計算時間（多項式表現、有理関数の係数）：
  | 計算 | N = 3 | N = 4 | N = 5 | 外挿 N = 8 |
  |---|---|---|---|---|
  | trefoil の 3 箱（各 λ） | 2.5 秒 | 19 秒 | 120 秒 | 約 10 時間 |
  | T(2,5) の 2 箱（各 λ） | 1.4 秒 | 10 秒 | 70 秒 | 約 7 時間 |
  補間の係数 c_j の s の幅は小さい（trefoil の 2 箱で 4〜10）ので、大きな N を一つ取れば c_j を直接読める（2 箱は N ≥ 8）。それでも 3 箱・T(2,5) には N = 8〜12 程度が要り、このままでは 1 時間を大きく超える。
- 次：高速化（DAHA の作用には分母が出ないので、γ̂(S_μ)·1 を整数係数の Laurent 多項式で計算し、最後に K_λ の係数を掛ける）。速度を測ってから、1 時間を超える計算はユーザーに確認する。
- テスト：tests/test_rk_daha.py（5 件、約 4 秒）。

## 2026-10-19（続き 23）高速化の試みと停止点（ユーザーの指示「実行に時間がかかる(1時間以上)ところで確認せよ」）
- 試み：DAHA の作用に分母が出ないことを使い、多項式を {(X の指数, r の指数, s の指数): 整数} の辞書で持つ版を作った。結果は従来版と一致したが、約 10 倍遅かった（(r, s) の冪を別々の項にすると項の数が増える）。採用せず削除した。
- 測定（cProfile、N = 4 の 3 箱 λ = (2,1)、24 秒）：時間の大半は Python の層での有理関数オブジェクトの生成と和（約 200 万回）で、flint の演算そのものではない。
- 1 時間を超える計算（確認待ち）：
  | 計算 | 必要な N | 現在の方式での見積もり |
  |---|---|---|
  | trefoil の 3 箱を一般の N で（(3), (2,1), (1,1,1)） | 8（窓の分離で直接読む）と検証用の 9 | 各 λ で N = 8 が約 10 時間、N = 9 が約 3 日 |
  | T(2,5) の 2 箱を一般の N で | 10〜12 程度（U の冪が増える） | 各 λ で数日以上 |
- 選択肢：
  - (a) 高速な実装を作る（X の多項式全体を一つの flint の多変数多項式として持ち、差分商を C 側の割り算にする、または q, t を数値にして素数を法とする評価から再構成する）。開発は数時間〜1 日の見込みで、速度は作ってから測る。
  - (b) 重い計算はせず、2 箱の結果の解析を進める。t = q での Mariño の ĝ_[2], ĝ_[1,1] との規約の対応を決め、多箱の refined クロスキャップ BPS 数の定義（G 展開の対称性が 2 箱で成り立たない問題）を検討する。
  - (c) 現在の方式のまま、trefoil の 3 箱の N = 8 を 1 つの λ で回す（約 10 時間、nohup・保存・再開あり）。


## 2026-10-19（続き 24）(b) 2 箱の結果の解析（ユーザーの指示「(b) を先に進め、並行して (a) を作れ」）
- (a) は別セッション「DAHA-Jones の高速な実装（タスク C の (a)）」（daha-fast、新しいファイルだけを作る）に任せた。報告を受けてから結果と差分を確認する。
- 【計算】KN の二つの枠の整合（scripts/rk_b_frames.py、3 秒、logs/refined_kauffman/b_frames.log）：
  - q 枠（Macdonald の Cauchy 核、因子 (Q − 1/Q)(T/Q)^{|μ|}、M(t)）と t 枠（双対 Cauchy 核、因子 (1/T − T)、M(q^{−1})）の hatted 不変量が、大きな N の 2 箱の SO・被覆・クロスキャップ部分で一致した。
- 【計算】t = q での Mariño 0904.1088 との辞書（scripts/rk_b_marino_tq.py、7 秒）：
  - SO の W_λ と被覆の W_μ（一様な規格化）は、Mariño の G_λ, H_μ（標準 framing）の鏡像（ν → 1/ν, t → 1/t）に比 1 で一致する。
  - g_(1) = −鏡像(g^Mar_[1])。2 箱では、全 d の g を奇数 d に直した g^odd が g^odd_(2) = 鏡像(g^Mar_[1,1]) − ½、g^odd_(1,1) = 鏡像(g^Mar_[2]) + ½。±½ は unknot の g_(1) = −1 の d = 2 の寄与（核の規約の違い、結び目によらない）。
  - unknot の値を引いて M(q^{−1})^{−1} で hatted にすると、Mariño の ĝ_R の鏡像の −1 倍（R はそのまま）。(4.8) の印刷された ĝ_[1,1] とも直接一致した（ĝ_[2] は z⁹ の項の誤植を訂正した値と一致）。
  - 途中で M(q^{−1}) の行列の向きを取り違えて M(q) を使っていた（bps.M_matrix(2, x) = [[1/x, −x], [−x, 1/x]]）。訂正後の値で上の結論を得た。和の方向 (1,1) の固有値は符号しか違わないので、和についての判定は訂正の前後で同じだった。
- 【計算】大きな N の 2 箱の hatted ĝ の構造（scripts/rk_b_2box_structure.py、6 秒、logs/refined_kauffman/b_2box_structure.log）：
  - ĝ_(2) + ĝ_(1,1) は整数 Laurent。(1 + t) の極は p_2 方向（差）だけにある。
  - Mariño の奇数 d の規約を素朴な Adams（q, t, U → q², t², U²）で refined に持ち上げると (1 + q) の極が加わる。refined では全 d の規約が自然である（Hopf の結果と同じ）。
  - (1 + t) の留数（t = −1 での (t^{1/2} + t^{−1/2})(ĝ_(2) − ĝ_(1,1))/2）は −i(U² − 1)²(qU² − 1)(q³U⁶ + q³U⁴ − qU² + 1)/(2q^{9/2}U^{12})。g_(1) の Adams 型の候補とは一致しない。
  - G 展開の対称性（σ：q → t^{−1}, t → q^{−1}、ν 固定）：c = 1 部分（V の偶奇で分ける。ν の偶奇と完全に連動する）は [1] で対称、2 箱では ĝ_(2)、ĝ_(1,1)、和、奇数 d 版のどれでも非対称。
    - 和の非対称の部分は (q − t)(ν²q − t)(qt − 1) × 多項式（(qt − 1) は σ 反対称な関数に必ず付く因子）。
    - σ で固定する Kähler 変数を U(rs)^c と変えても、[1] と和の両方で対称になる c は無い（U 固定では和は全体が c = 2 になり、c = 1 部分は 0）。
  - SO の簡約した不変量 R_λ には、σ・q ↔ t・鏡像のいずれに ν → ±ν^{±1}(rs)^{c1}(r/s)^{c2}（|c1|, |c2| ≤ 2）と R ↔ R^T を組み合わせた写像でも、比が単項式になる対称性は無い（[1] だけでも無い）。type A の super-duality に当たるものが SO 側には無く、SO の f̂ が σ 対称でないことと整合する。被覆の f̂ は R ごとに σ 対称。
  - **述べ直しの判定式（report §8、Hopf から得たもの）は 2 箱で成り立つ**：refined だけの部分 X = ĝ − ĝ|_{t=q}(ν → ν~)（ν~ = ν(q/t)^{1/2}）は、ĝ_(1)、ĝ_(2)、ĝ_(1,1)、和のすべてで (q − t)(ν² − t)（太字で (𝐭² − 1)(𝐚² − 𝐪²) と単項式の積）で割り切れる。(ν² − t) で割れることは「SO(2) の点 U = t で ĝ が t によらない」ことと同じ。ν~ を使わず ν のままでは割り切れない。
  - 商の符号（太字、𝐚 の冪ごと）：[1] は 3 項で冪ごとにそろう。2 箱の和では c = 2（𝐚 の偶数冪）の 𝐚^{−10}, 𝐚^{−8}, 𝐚^{−6}, 𝐚^{−4} でそろい（+, −, +, −）、𝐚^{−12} と c = 1（奇数冪）は混在する。商自体も σ 対称でない。正値性は未判定。
- テスト：tests/test_rk_daha.py に 3 件（Mariño の辞書、KN の枠、構造の判定は results の印）。既定 7 件、約 11 秒。
- 未解決：多箱のクロスキャップの c = 1 部分の展開の変数（G 展開は 2 箱で使えない）、c = 1 部分の符号の規則、(1 + t) の留数の意味。
- 次：daha-fast の報告を確認する。3 箱・T(2,5) の一般の N が得られれば、述べ直しの判定式（(𝐭² − 1)(𝐚² − 𝐪²) で割り切れること）を 3 箱と別の結び目で確かめる。

## 2026-10-19（続き 25）TeX の報告書（ユーザーの指示「これまで計算した結果を確認したい。日本語のtexでまとめて報告せよ。実際に計算した結果も載せよ。」）
- docs/refined_kauffman/tex/report.tex（lualatex、17 ページ、report.pdf）に K0〜タスク C (b) の結果をまとめた。
- 付録の式は scripts/rk_report_tex.py が保存済みのデータとコードから生成する（docs/refined_kauffman/tex/gen/、15 秒）。生成時に Mariño との辞書、ĝ_(1) = −g^{true}_[1]、[1] の (N)(P)(c=2)（m ≤ 3）、(1+t) の留数の式を assert で確かめる。
- 載せた計算結果：Mariño の trefoil の ĝ_[2], ĝ_[1,1]（z⁹ の誤植を訂正した値）、T(2,2m+1)（m ≤ 3）の g^{true}_[1]・c = 1 の G 展開・Π_m、DAHA の R_(1), R_(2), R_(1,1)（U の係数）と被覆の随伴、大きな N の 2 箱の ĝ_(2) + ĝ_(1,1)・(1+t)(ĝ_(2) − ĝ_(1,1))・商と符号の表。
- テスト：tests/test_rk_*.py の既定 55 件が通る（約 75 秒）。daha-fast の結果は未確認なので載せていない。
- 再生成：`.venv/bin/python scripts/rk_report_tex.py` のあと docs/refined_kauffman/tex で lualatex を 2 回（TEXMFVAR を書き込める場所にする）。

## 2026-10-19（続き 26）(a) 高速版の確認（daha-fast セッションの commit 3ed3257）
- daha-fast は新しいファイルの追加だけ（src/rcsknots/daha_fast.py、scripts/rk_c_fast_daha_dn.py、tests/test_rk_daha_fast.py、docs/refined_kauffman/daha_fast.md、logs、data の新規 67 件）で、既存のファイルは変更していない。内容は docs/refined_kauffman/daha_fast.md。
- 報告：trefoil の 3 箱は N = 8 で 24 秒、N = 9 で 43 秒（現行版は外挿で約 1.4 日、約 10 日）。trefoil の 3 箱を N = 2〜12、T(2,5) の 2 箱を N ≤ 14 で保存した。
- このセッションでの独立な確認【計算】：
  - tests/test_rk_daha_fast.py の 3 件が通る（5 秒）。
  - 高速版で保存した 3 箱を現行版で計算し直して全係数一致：N = 4 の (1,1,1)（22 秒）、N = 5 の (2,1)（118 秒）。Z と K の両方。
  - N ≤ 9 の現行版のデータから作った一般の N の式（R_(1), R_(2), R_(1,1)、data/refined_kauffman/daha/R_U_m1_2box.pkl）が、高速版の N = 10, 11, 12 の値を 9 件すべて全係数で予言した。
- 次（未着手、ユーザーの確認待ち）：trefoil の 3 箱と T(2,5) の 2 箱の一般の N の補間、大きな N のクロスキャップ部分、述べ直しの判定式（(𝐭² − 1)(𝐚² − 𝐪²) で割り切れること）の 3 箱・T(2,5) での確認。どれも 1 件 30 分未満の見込み。
- 追記（ユーザーの質問「refinedのunknotとhopfの結果はどこ？」）：TeX の付録 A.4 に unknot と Hopf の計算結果を加えた（25 ページ）。1 probe の全 d の 2(1+t)g_22、g_23、g_33、C_22（ν~ の冪ごと）、多 probe の hatted Ĉ_{R,S}（|R| ≤ 4, |S| ≤ 2）の t = q の部分と refined だけの部分の商、§6.3 に符号の表（14 成分中、refined だけの部分で混在は Ĉ_{(3,1),(1,1)}, Ĉ_{(2,1,1),(2)}）。生成時に unknot の g、g_{1b}（b ≤ 4）、E_22、E_24 の閉じた形を保存済みのデータと照合した（scripts/rk_report_tex.py の gen_hopf、約 2 分）。

## 2026-10-19（続き 27）正値性の現状の整理（ユーザーの質問：U(N) 以外の rCS の正整数性、refined での符号の混在、スピンごとに分けると正値性が見えるか）
- 【計算】SO(2N) の rCS そのもの（一般の N の trefoil、太字 r = −𝐭𝐪, s = 𝐪, U = 𝐚𝐪）：R_(1) は 9 項（正 5・負 4）、R_(2) は 55 項（29・26）、R_(1,1) は 44 項（24・20）。(𝐚, 𝐪, 𝐭) の符号の取り替えと全体の符号のどれでも同符号にならない（[1] でも 𝐭⁴𝐪𝐚³ と −𝐭⁴𝐪^{−1}𝐚³ のように 𝐪² だけ違う項が逆符号）。U(N) の [1] のトーラス結び目（superpolynomial、正）と違い、SO では [1] でも rCS は正でない（AgSc も「一般に負の符号を含む」と書く）。
- 【計算】クロスキャップ部分の展開の可能性（大きな N の trefoil）：
  - z_q の多項式（q^{1/2} → −q^{−1/2} で不変、t と ν_k = ν(q/t)^{k/2}（|k| ≤ 2）を固定）には、[1]（m = 1, 2）も 2 箱（ĝ_(2), ĝ_(1,1), 和）も書けない。
  - U(1)×U(1) の電荷（単項式）で分け、符号を ±(−1)^{α e_𝐪 + γ e_𝐭} で説明できるか：2 箱の和の商で符号がそろわない 𝐚 ブロック（𝐚^{−12}, 𝐚^{−11}, 𝐚^{−9}, 𝐚^{−7}, 𝐚^{−5}）は、どの (α, γ) でもそろわない。
- 【文献、AgSc 1202.4456 の TeX】orientifold は SU(2)_L × SU(2)_R を U(1)_1 × U(1)_2 に壊すので、unoriented の BPS は完全なスピン多重項をなさなくてよい（§7）。二次量子化の整数 M, M~ は整数性を明らかにするが、スピン多重項への整理を見えなくする（§4）。
- 【文献、KN 1703.05408 の TeX】unrefined の LMOV 不変量は \hat N_{μ,g,β} = Σ (−1)^{2J_r} \hat N_{μ,g,β,J_r,S_R}（符号付きの和）。トーラス結び目では ρ, g, β ごとに 2J_r の偶奇がそろい、\hat N_{ρ,g,β} = ±Σ_{J_r} \hat N_{ρ,g,β,J_r}。4_1 では unrefined で打ち消しが起きる。
- 次の候補（未着手）：結び目の固定セクターの被覆（Hopf で行ったもの）で E を P と Q（どちらも非負であるべき）に分ける。C を二次量子化の積の指数（AgSc (3.13) の M~ の開弦版）に直し、電荷ごとの符号を見る。高速版で得た 3 箱・T(2,5) で同じ検査をする。

## 2026-10-19（続き 28）自由エネルギーの右辺の基底の取り替え（ユーザーの質問「logZの右辺のschurを別の基底に置き換えるパターンは試した？」）
- これまで試したもの：K2 の Mariño 型（Schur、P_λ(v; q,t)、P_{λ^T}(−v; t,q) × 重み 5 通り、unknot の 2 箱で極が打ち消さない）、K2'(b) の Hopf の多 probe（Schur と P_R(y; t,q)、後者は (1 − qt) などの分母で整数性が破れる）、hatted の M(x)（Hopf で M(q^{−1}), M(t), M(t^{−1}), M(q)）。結び目では Schur + M(q^{−1}) だけで、正値性を目的に基底を替えたことは無かった。
- 【計算】大きな N の trefoil の 2 箱の g（全 d）を基底を替えて展開した（scripts/rk_b_2box_basis_den.py）：
  - Schur・HL(t)・HL(q)・P(q,t)：両方の係数に (1 + t) の極（P(q,t) は P_(1,1) の係数に (1 − qt) も）。
  - **P_λ(Y; t, q)（双対 Cauchy の probe 側の Macdonald 基底）と P(t^{−1}, q^{−1})：(1 + t) の極は P_(2) の係数だけに集まり、P_(1,1) の係数 x_(1,1) は (1 − qt) の分母だけ。** x_(1,1) = g_(1,1) + (q − t)/(1 − qt) g_(2) で、E_(2) = −E_(1,1) から (1 + t) が打ち消す【導出】。
  - Q(q,t)・Q(t,q)・H̃：(1 − t)、(1 − q)、(q − t) などの分母が加わる。
- 【計算】E を含まない方向で、冪和の定数倍型の基底（÷ φ）を走査した（scripts/rk_b_2box_basis_sign.py、logs/refined_kauffman/b_2box_basis_sign.log、24 秒）：
  - X = g_(2) + g_(1,1) を割り切る φ は z_q（= 1 − q と単項式）だけ（z_t, 1 − t, 1 − qt などでは割り切れない）。どの φ でも c = 1 部分は σ で単項式倍を除いても対称にならない。
  - 普遍因子 (q − t)(ν² − t) で割った商の 𝐚 ブロックごとの符号：Schur の ĝ_sum = X/z_q は 5/10、**P(Y; t, q) の (1 − qt) x_(1,1)/z_q は 8/11** でそろう。後者では c = 1（𝐚 の奇数冪）の 𝐚^{−11}, 𝐚^{−9}, 𝐚^{−7}, 𝐚^{−5}, 𝐚^{−3} がすべてそろい（+, −, +, −, −）、混在は 𝐚^{−12}, 𝐚^{−10}, 𝐚^{−8}（c = 2 の 3 ブロック）だけ。(1 − qt) は Laurent にするための因子で、⟨P_(1,1), P_(1,1)⟩_{t,q} = (1 − qt)(1 − t)/((1 − q²)(1 − q)) に含まれる（Q 基底の規格化に当たる）。
  - 位置づけ：符号のそろい方は基底に依存し、probe 側の自然な基底で改善した【計算、trefoil の 2 箱だけ】。正値性が見えたとは言えない（c = 2 の 3 ブロックが混在）。
- 次の候補：P(Y; t, q) 基底で Hopf の Ĉ_{R,S}（符号が混在した Ĉ_{(3,1),(1,1)}, Ĉ_{(2,1,1),(2)}）と固定 N の 3 箱（E = 0 のセクター）を見直す。

## 2026-10-20（続き 29）P_λ(Y; t, q) 基底での見直し（ユーザーの指示「次にできることを実行せよ」）
スクリプト：scripts/rk_b_pbasis.py。ログ：logs/refined_kauffman/b_pbasis.log、b_pbasis_hopf_K.log、b_pbasis_fixedN3.log。
- 【計算】trefoil の 2 箱で順序を替える：hatted（M(q^{−1})^{−1}）にしてから P 基底に替えると、P_(1,1) の係数の商は 0/11 でそろう（悪化）。続き 28 の 8/11 は「P 基底の係数を取ってから z_q で割る」順序に限る。
- 【計算】Hopf の 2 箱 × 2 箱（E を引いた C）：
  - hatted にしてから両成分を P 基底に替える：t = q の部分は 4/4 でそろうが、refined だけの部分は (2),(2) 以外の 3 成分で普遍因子が (q − t) だけに減り、符号も混在する。
  - P 基底の係数を取ってから z_q² で割る：(2),(2) に (1 − q) の分母が残り、t = q の部分の符号も 4 成分すべてで混在する。
  - Schur 基底の hatted Ĉ は 4 成分とも符号がそろっていた（続き 26 の表）ので、Hopf では P 基底のほうが悪い。|R| = 4 の成分は sympy の 4 次の基底変換が 20 分を超えたので打ち切った（2 × 2 の結果から改善は見込めない）。
- 【計算】固定 N = 2, 3 の 3 箱（E = 0 のセクター）：Schur 基底の g_R は整数 Laurent。P(Y; t, q) 基底では P_(3) の係数は整数 Laurent だが、P_(2,1) に (1 − q t²)、P_(1,1,1) に (1 − qt)(1 − q² t) の分母が出る（ĝ でも同じ）。これは Macdonald の内積 ⟨P_λ, P_λ⟩_{t,q} の型の分母で、K2'(b) の Hopf と同じ現象。
- 結論：P(Y; t, q) 基底は、trefoil の 2 箱の一つの方向・一つの順序でだけ符号のそろい方を改善し、Hopf では悪化し、3 箱では整数性を壊す。**基底の取り替えだけで正値性が見えるという証拠にはならない。** 続き 28 の 8/11 は、特定の組み合わせでの偶然の可能性を否定できない。
- 次の候補（未着手）：二次量子化の形（E を P と Q に分ける、C を凍結した状態の積の指数に直す）で、電荷ごとの符号を見る（続き 27）。

## 2026-10-20（続き 30）二次量子化の形での電荷ごとの符号（ユーザーの指示「二次量子化の形で電荷ごとの符号を調べよ」）
スクリプト：scripts/rk_b_second_quant.py。ログ：logs/refined_kauffman/b_second_quant.log、b_second_quant_2.log、b_second_quant_kn.log（どれも数秒）。
- 定義：Z~ = Π (1 − (電荷の単項式)·y^α)^{−c}（AgSc (3.13) 型）の指数 c。c_μ = Σ_R K_{Rμ} f_R（Y の単項対称式の係数）で、f の分母 (1 − t^p) を t の冪級数（モードの塔）に展開した係数が、電荷 ν^β t^{k/2} q^{a/2} ごとの「状態数 × (−1)^F」。E と C の分け方の任意性によらず決まる。
- 【計算】SO の f の分母は単項式 × (1 − t)（[1]）・(1 − t²)（2 箱）だけで、展開した指数はすべて整数。(1 − t²) は AgSc の「Ĩ で固定される状態は t の偶数か奇数のモードだけ」と合う。被覆は (1 − t)。
- 【計算】[1]（M の変換が無いので c がそのまま状態数）：
  - T(2,3)：SO の f の (β, k) ブロック 20 個がすべて同符号で、ν^β ごとに符号が一つに決まる（G 展開なしで ε_β の型の正値性）。
  - T(2,5)：33 中 3 ブロック（β = −8, −7, −6）、T(2,7)：48 中 10 ブロックで混在。混在するのは §8 の判定式で (𝐭² − 1) が掛かるブロック。ν^β ごとに (q − t) の冪をくくり出すと 25 中 2、35 中 6 に減るが残る。
  - 被覆の [1] は全ブロックでそろう。
- 【計算】2 箱：
  - 単項式の基底の指数（hatted にしない）は SO も被覆もほぼ全ブロックで混在（SO の m_(1,1) は 113/113）。被覆でも崩れるので、多箱の開弦ではこの指数は状態数ではない（境界の巻き付きの構造 M を外す必要がある）。
  - hatted（M(q^{−1})^{−1}）にしてから展開：SO の f̂_(2) は 123 中 16、f̂_(1,1) は 110 中 15 ブロックで混在。混在の多くは β が奇数のブロック（被覆の f̂^cov は β が偶数だけなので、クロスキャップだけから来る部分）。被覆の f̂^cov は 42 中 2、9 ブロックで混在。
  - 検算：KN の t 枠 (t^{−1/2} − t^{1/2}) f̂^cov を G 展開すると、(2) は全 β でそろい、(1,1) は β = −10, −8 で混在が残る（被覆の brane と image をまたぐ寄与には対角の M でなく M ⊗ M が要る、という続き 13〜16 の結果と整合）。SO の f̂ は σ 対称でなく G 展開できない。
- 結論：電荷ごとの数え直しだけでは、正値性が見えるのは trefoil の [1] まで。残る混在は (a) [1] の (𝐭² − 1) を含むブロック（t = q で打ち消す refined だけのボソン・フェルミオンの対）、(b) 2 箱の β が奇数のクロスキャップのブロック、(c) 被覆の composite のセクター（M ⊗ M が要る）にある。
- 次の候補：(c) 被覆と SO を brane・image の 2 変数に分けた生成関数（固定 N の cover_Z2 の大きな N 版）で M ⊗ M にする。(a)(b) クロスキャップの多重項の変数を w = (q/t)^{1/2} の (w − w^{−1})^k 型の因子として探す。

## 2026-10-20（続き 31）次の候補の実行：被覆の brane・image の分割と M ⊗ M、多重項の変数（ユーザーの指示「次の候補を実行せよ」）
スクリプト：scripts/rk_b_cover_split.py、scripts/rk_b_multiplet_var.py。ログ：logs/refined_kauffman/b_cover_split*.log、b_multiplet_var.log。データ：data/refined_kauffman/daha/largeN_2box_m1_cover_split.pkl。
- 【計算】(c) 大きな N の被覆を brane（Y）と image（Y'）に分けた生成関数：Z = 1 − W_1(p_1 + p_1') + [W_(2) e_2(Y) + W_(1,1) P_(2)(Y; t,q) + (Y → Y')] + (W_adj + κ) p_1(Y)p_1(Y')。Y = Y' で cover_f に一致（p_1, p_1², p_2）。
  - **brane–brane のセクター（f̂_(2), f̂_(1,1)、M ⊗ M）は、KN の t 枠 (t^{−1/2} − t^{1/2}) f̂ で G 展開でき、全 β で符号がそろう**（KN の正値性が成立）。続き 30 の対角の M での (1,1) の混在は、brane–image をまたぐ寄与から来ていた。
  - brane–image のセクター x = W_adj + κ − W_1²：(t^{−1/2} − t^{1/2}) x は Laurent だが σ 対称でない。**(t^{−1/2} − t^{1/2}) x / z_q は整数 Laurent で σ 対称**（(t^{−1/2} − t^{1/2})² x も同じ）。z_t、z_q ± z_t では割り切れない。M ⊗ M の φ_(1)φ_(1) = 1 ではなく z_q が要る（つまり対角の M の割り方と同じ）。
  - その G 展開：ν^{−12}, ν^{−6}, ν^{−4}, ν^0 はそろう、**ν^{−10} は −(W⁴ − 1) × 正（W⁴ − 1 = 𝐭² − 1、クロスキャップと同じ普遍因子）**、ν^{−8}, ν^{−2} は混在。向き付けられた被覆の brane–image のセクターが、クロスキャップと同じ (𝐭² − 1) の構造をもつ。
  - brane–image の部分をクロスキャップ側に戻した g' = f̂_SO − ½(brane–brane + image–image)：(1 − t) の分母が残り、(t^{−1/2} − t^{1/2}) を掛けても σ 対称にならない。
- 【計算】(a)(b) 多重項の変数：ν^β ごとのブロックを (r, w)（w = (q/t)^{1/2}）で書き、w → 1/w、w → −1/w、r → −1/r、r → 1/r の対称性（単項式倍まで）を調べた。2 箱の ĝ_(2) + ĝ_(1,1) とその商、[1] の T(2,3), T(2,5) のどれでも、項の少ない端のブロック以外は対称性が無い。[1] の G 展開（V と W の両方を含む）に代わる単一の変数（w ± 1/w、z_q、q^{1/2} + q^{−1/2}）の多重項は見つからない。
- 結論：
  - 被覆（向き付けられた部分）の符号の破れは、brane–image をまたぐセクターに限られる。brane–brane は KN の正値性を満たす。
  - brane–image のセクターとクロスキャップは、同じ普遍因子 (𝐭² − 1) をもつ「Ĩ に関わる」部分で、どちらにも混在するブロックが残る。
  - 正値性の判定の対象は、(brane–brane の被覆) と (brane–image ＋ クロスキャップ) に分けて考えるのが自然だが、後者の多重項の構造（G に代わる展開）はまだ見つからない（未解決）。
- 追記（2026-10-20、ユーザーの指示「texに今回の結果を追記して報告せよ」）：TeX の §11「正値性の検討」に続き 27〜31（符号の解釈、SO の rCS の符号、基底の取り替え、二次量子化の指数、被覆の brane・image の分割、多重項の変数）を加え、§12 の表と §13 の未解決を更新した（27 ページ）。

## 2026-10-20（続き 32）① 一般の N の補間（3 箱・T(2,5)）と ② T(2,5) の 2 箱のクロスキャップ（ユーザーの指示「1,2を実行せよ」）
スクリプト：scripts/rk_c_read_largeN.py（窓を s の次数の傾きから決めて Newton 補間）、scripts/rk_c_fast_adj.py、scripts/rk_c_2box_largeN_m.py、scripts/rk_b_2box_structure_m.py。ログ：logs/refined_kauffman/c_*.log、b_2box_structure_m2.log、b_cover_split_m2.log。
- 【計算】① 窓の決め方：R_λ(N) の s の次数の最小・最大は大きな N で 2N·j + 定数なので、隣り合う 3 つの N から [j_min, j_max] を決め、大きな N をノードにして補間し、残りの N で照合する。trefoil の [1], (2), (1,1) は既知の窓（U^{−5}〜U^{−2}、U^{−10}〜U^{−4}）を自動で再現し、補間に使わない N でも一致。
  - T(2,5) の [1]：U^{−9}〜U^{−4} の 6 項（N = 9〜14 で補間、N = 2〜8 で一致）。(2), (1,1)：U^{−18}〜U^{−8} の 11 項（N = 4〜14 で補間、N = 2, 3 で一致）。data/refined_kauffman/daha/R_U_m2_{1,2,11}.pkl。
  - 実装の修正：rk_c_interp.s_range が固定のずらし（400）で s^{−400} 未満を扱えなかったので、分母の単項式を引く形に直した。
  - **trefoil の 3 箱**：(3), (2,1), (1,1,1) とも U^{−15}〜U^{−6} の 10 項（N = 4〜13 で補間、(3), (2,1) は N = 2, 3、(1,1,1) は N = 3 で一致。D_N で ℓ(λ) ≤ N のため (1,1,1) は N ≥ 3）。照合点を増やすため高速版で N = 13 を追加計算した（4.9 分）。補間は sympy の約分が重く (3) で 24 分。data/refined_kauffman/daha/R_U_m1_{3,21,111}.pkl。
  - 【計算】t = q での照合（scripts/rk_c_check_tq.py）：R_λ|_{t=q, U = ν q^{1/2}} は NRZ の [r] 彩色の Kauffman 多項式の鏡像に比 1（[1], [2], [3]）、[1,1], [1,1,1] は F_{[r]}(−ν, t^{−1}) の鏡像に比 1。**refined SO(2N) の trefoil の 3 箱の不変量が一般の N で得られた**（(2,1) は t = q の独立な照合の相手が無い）。
  - テスト：tests/test_rk_daha.py に results の印の 1 件（3 箱の t = q、T(2,5) の普遍因子、約 30 秒）。
- 3 箱のクロスキャップ部分（大きな N）には被覆の composite (2;1), (11;1) の一般の N が要る（未着手、判断待ち）。
- 【計算】② GL の DAHA と一様な規格化の単項式の規則（T(2,5)）：M_b = q^{−5(b,b)} t^{−10(b,ρ) + 3|b|(N−1)}（普通の表現 9 例、N = 2〜4）。随伴でも U(3), U(5) の data と比 1。m = 1 と合わせて M_b = q^{−n(b,b)} t^{−2n(b,ρ) + (m+1)|b|(N−1)}（n = 2m + 1）。
  - 随伴を N_A = 2〜24 で計算（各 8 秒以内）、A^{1/2} の冪 −18〜−8 で補間（N_A = 14〜24 で補間、2〜13 で一致）。data/refined_kauffman/daha/Radj_A_m2.pkl。
  - 大きな N の T(2,5) の 2 箱のクロスキャップ部分（data/refined_kauffman/daha/largeN_2box_m2.pkl）：g_(2), g_(1,1) の分母は 2 × 単項式 × (1 + t)、g_(2) + g_(1,1) と (1 + t)(g_(2) − g_(1,1)) は整数 Laurent、U = t², t³ で固定 N の g と全係数一致。trefoil と同じ構造。
  - 検査（trefoil と比較）：
    - **述べ直しの判定式の普遍因子 (q − t)(ν² − t)（太字で (𝐭² − 1)(𝐚² − 𝐪²)）は、ĝ_(1), ĝ_(2), ĝ_(1,1), 和のすべてで割り切れる**（trefoil に固有でない）。
    - c = 1 部分は 2 箱で σ 対称でない（trefoil と同じ）。
    - 商の 𝐚 ブロックの符号：ĝ_(1) は 5/5 でそろう、和は 6/18（trefoil は 5/10）。
    - 二次量子化の指数（hatted、t の冪級数）：SO の f̂_(2) は 296 中 102、f̂_(1,1) は 277 中 92 ブロックで混在。被覆の f̂^cov は 99 中 29。[1] は DAHA から作っても AS の式と同じく 33 中 3。
    - 被覆の brane・image の分割：brane–brane は G 展開の後に全 β でそろう（KN の正値性）。brane–image は z_q で割ると σ 対称で、G 展開の ν^β のブロックは「そろう」か「(𝐭² − 1) × そろう」のどちらか（ν^{−20}, −18, −8, −6, −4, −2, 0）で、中央の ν^{−16}〜ν^{−10} の 4 ブロックだけ混在。

## 2026-10-20（続き 33）3 箱のクロスキャップ部分（大きな N、trefoil）（ユーザーの指示「3 箱のクロスキャップへ進め」）
スクリプト：scripts/rk_c_fast_composite.py、scripts/rk_c_3box_largeN.py、scripts/rk_b_3box_structure.py。ログ：logs/refined_kauffman/c_fast_composite_m1.log、c_3box_largeN.log、b_3box_structure.log。
- 【計算】被覆の composite (α;1)（α = (2), (1,1)）の一般の N：
  - 高速版の拡張：X^{−1} 側（β = −ε_1）の対称化で γ̂'(p_1(Y))·1 を作り、X 側の連鎖の出発点にして {γ̂(P_α(X) p_1(X^{−1}))} を計算（gl_adj_fast の一般化）。
  - 既約への分解 P_α(X) p_1(X^{−1}) = P_b + d_α P_(1)（b = α − ε_N）。ev は乗法的で、ev(P_λ) = P_λ の主特殊化（a = t^{−N}, q^{−1}, t^{−1}）（N = 3, 4 の 7 例で多項式と一致）から d_α, K_b を閉じた形で求め、N ≤ 5 では多項式の分解と照合した。
  - U(3) の [3,1], [2,2]、U(5) の [3,1,1,1], [2,2,1,1]（KN の data、一様な規格化）と比 1（単項式の規則 M_b は m = 1 のまま）。
  - N_A = 2〜20 で計算（N = 20 で 1 件 3 分、合計 19 分）、A^{1/2} の冪 −14〜−6 の 9 項で補間（N_A = 12〜20 で補間、2〜11 で一致）。data/refined_kauffman/daha/Rcomp_*.pkl。
- 【計算】大きな N の 3 箱の組み立て（data/refined_kauffman/daha/largeN_3box_m1.pkl）：
  - SO：Y の 3 次までの核を D_N の K_λ（|λ| ≤ 3、U は形式変数）の逆行列で K 基底に直し、⟨K_λ⟩ = R_λ(U) K_λ(t^ρ)。被覆：⟨P_λ(x) P_μ(x^{−1})⟩ を普通の表現・随伴・composite から。
  - U = t², t³ で SO の f、被覆の f^cov、g_R（|R| = 3）が固定 N の値と全係数一致。
  - **g_(3), g_(2,1), g_(1,1,1) は整数 Laurent で (1 + t) の極が無い**（奇数のセクターに E は無いという全 d の構造の予言どおり）。
- 【計算】3 箱の検査（E を含まないセクター）：
  - **g_R と ĝ_R（|R| = 3）の refined だけの部分は、3 つとも普遍因子 (q − t)(ν² − t)（太字で (𝐭² − 1)(𝐚² − 𝐪²)）で割り切れる。**述べ直しの判定式の普遍因子は trefoil の 2 箱・3 箱、T(2,5) の 2 箱のすべてで成立。
  - c = 1 部分は σ 対称でない（2 箱と同じ）。
  - 普遍因子で割った商の 𝐚 ブロックの符号：ĝ_(3) 2/15、ĝ_(2,1) 3/14、ĝ_(1,1,1) 2/13 でそろう（大部分が混在）。
  - 二次量子化の指数（hatted の SO の f̂、t の冪級数）：326 中 124、313 中 130、276 中 121 のブロックで混在。
- テスト：tests/test_rk_daha.py に results の印の 1 件（3 箱の固定 N との一致、整数性、普遍因子、約 90 秒）。
- 結論：3 箱でも、整数性・全 d の構造・普遍因子は成り立つが、符号のそろい方（正値性）は商の段階でも二次量子化の段階でも見えない。正値性の判定には、クロスキャップ（と brane–image）の多重項の構造が要る（未解決）。

## 2026-10-20（続き 34）多重項の構造（タスク 1）：E を固定セクターの状態に分ける試み（ユーザーの指示「1」、方針 (b)）
スクリプト：scripts/rk_b_knot_fixed_sector.py、rk_b_knot_fixed_sector_diag.py、rk_b_knot_fixed_sector_resid.py。ログ：logs/refined_kauffman/b_knot_fixed_sector*.log。
- 方針：Hopf（rk_k4_fixed_sector.py）と同じく、Ĩ で固定される状態の数を E = P − 2Q（P は固定セクターの全状態、Q は Ĩ の固有値 −1 の状態、0 ≤ Q_m/P_m ≤ 1）と書く。結び目の 2 箱の固定セクターは被覆の brane–image のセクター x で、P = (1 − t)x。条件 (i) C = g − E/(2(1+t)) が整数 Laurent、(ii) t = q で E = −(1 + q)(g_(1)(q², ν²) + 1)（Mariño との辞書から導いた値、Hopf の E_22 = −(1 + q)ĝ_11(x²) と同じ形）。
- 【計算】trefoil：整数計画に解が無い。E_(2) = −E（フェルミオン的な符号）にしても無い。切り分け：(i) だけでも 0 ≤ Q/P ≤ 1 の範囲では解が無い（範囲を外せば整数解はある）。(ii) だけなら解がある。(i) と (ii) を両方課すと、範囲を外しても解が無い。
- 【計算】t = −1 での比較：
  P(−1) = (ν² − 1)(ν² + 1)(q − 1)(qν² − 1)(q³ν⁶ − q³ν⁴ + qν² + 1)/(ν^{12}q^5)、
  R(−1) = (ν² + 1)²(q + 1)(qν² + 1)(q³ν⁶ − q³ν⁴ − qν² − 1)/(ν^{12}q^5)。
  R は P の因子 (x − 1) を (x + 1) に替えた形（x = ν², q, qν²、最後の因子でも (qν² + 1) の符号が反転）。Hopf の P_22 → E_22（z_q² → −(z_q^+)²、第 2 因子の符号の反転）と同じ型。
  - Hopf では E_22 = i z_q [P_22/z_q](q^{1/2} → i q^{1/2}, ν → ±iν) と書けるが、この規則（と ±1, ±i の変数の代入の総当たり）は結び目では (i) も (ii) も満たさない。
- 解釈【予想】：P は係数の符号が混在する（25 正・25 負）指数なので、「|E_m| ≤ |P_m|」（P_m が状態の数であることを前提）は使えない。t = −1 の対応は、フェルミオンの振動子の因子 (1 − x) が (−1)^N 型の挿入で (1 + x) になる形であり、Ĩ は固定セクターの Fock 空間の因子ごとにフェルミオン数の符号として働く。E は「P の部分和」ではなく「P の積の構造に Ĩ を挿入したトレース」として求める必要がある。
- 次の候補：(1) 固定セクター x を一粒子の Fock の因子の積（あるいは brane と image の円板の対）に分解する式を、DAHA の側（brane–image の composite (1;1) の構造）から作り、Ĩ の挿入をその上で定義する。(2) Hopf の ε_{ρ,(1)} = ±E_22 ψ²[X_ρ] のような Adams 型の規則を結び目の x に当てる。

## 2026-10-20（続き 35）(2) Adams 型の規則、(1) Ĩ の挿入の構成の検討（ユーザーの指示「(2) (1)と進め」）
スクリプト：scripts/rk_b_knot_adams_E.py。ログ：logs/refined_kauffman/b_knot_adams_E.log、b_knot_adams_E_resid.log。
- 【計算】(2)：E = F · (ψ²[g_(1)] + 1)（ψ²：q, t, ν → q², t², ν²(q/t)^{j/2}、|j| ≤ 2、F は (1 + q)、(1 + t)、(q + t)、(1 + qt) などとその単項式倍、符号 ±）は、trefoil・T(2,5) とも条件 (i) を満たさない。t = −1 で R(−1) と ψ²[g_(1)] + 1 の比も単純な因子にならない。t = q での AgSc の強い整数性の形 E = −(1 + q)(g_(1)(x²) + 1) は、refined ではそのままでは持ち上がらない。
- 【導出】(1) の検討：Ĩ を挿入したトレースは Tr(Ĩ) = 2 f_SO − f^{cov} = 2g = E/(1 + t) + 2C（全 d）。E（Ĩ で固定される向き付けられた状態、モードの塔 1/(1 + t)）と C（凍結した状態）は、t の冪級数の大きな次数の振る舞いからは E(t = −1) しか決まらない（E → E + 2(1 + t)X の任意性、続き 13 の Hopf と同じ）。固定セクターの指数 P = (1 − t)x は、この任意性を決めない。
- 【計算】続き 34 の切り分けとあわせると、(i)(ii) を満たす E は P の単項式の台に収まらない（範囲を外した整数解も無い）。つまり E は P に現れない電荷の状態を含む。P は (−1)^F の付いた指数なので、ボソンとフェルミオンの打ち消しで状態が見えなくなっている。
- 結論【予想】：E を状態の和として決めるには、指数 P ではなく、打ち消しの前の状態の数（KN 型の正の \hat N_{J_r} のような、スピンで分けた数）が要る。被覆の brane–image のセクターの G 展開は ν^{−8}, ν^{−2} で符号が混在した（続き 31）ので、そこでの正の分解はまだ無い。E の決定と、brane–image のセクターの正の多重項の分解は同じ問題になる。
- 判断が要る点：E と C の分離の原理を、(a) 物理の要請（固定された状態の塔の始まり、AgSc の積表示の l_2 ≥ 1 の条件）から決めるか、(b) 「E と C がともに符号のそろう形になる」ことを分離の基準にするか。(b) は判定したい正値性そのものを仮定に使うので、検証にならない。

## 2026-10-20（続き 36）(a) E と C の分離を物理の要請から決める試み（ユーザーの指示「(a) を進めよ」）
スクリプト：scripts/rk_b_knot_fixed_twist.py、rk_b_knot_fixed_twistG.py。ログ：logs/refined_kauffman/b_knot_fixed_twist.log、b_knot_fixed_twistG.log。
- 要請【仮説】：固定セクターの状態は多重項をなし、指数 P の差の因子 (x^{1/2} − x^{−1/2}) はフェルミオンのゼロモードの寄与。Ĩ がゼロモードに (−1)^F として働けば、Ĩ を挿入したトレース E では差の因子が和の因子 (x^{1/2} + x^{−1/2}) に替わる。Hopf の P_22 → E_22（z_q² → −(z_q^+)²、境界の z_q を一つ残す）はこの形。
- 【計算】trefoil・T(2,5) とも、P の分子 = c · t (t − ν²)(ν² − 1)(q − 1)(qν² − 1) · L（L は trefoil 8 項、T(2,5) 40 項）。
  - (ν² − 1), (q − 1), (qν² − 1) を和に替え、普遍因子 (ν² − t) を残すと、t = −1 の R の因子の形 (ν² + 1)²(q + 1)(qν² + 1) が正しく出る（両結び目）。
  - 残りの L は、単項式ごとの偶奇の符号の付け替え（±(−1)^{αa + βb + γc} の 16 通り）では条件 (i)(ii) を満たさない。T(2,5) では t = −1 で必要な値の係数に 2 が現れ、符号の付け替えでは作れない。
  - L を多重項とみて、Pk = (t^{−1/2} − t^{1/2})x/z_q（σ 対称）の G 展開で G → G^+ = z_q^+ z_t^+、z_q → z_q, z_q^+、ν^β の符号 1, (−1)^{β/2}, (−1)^β, i^β を試したが、t = −1 で R に一致しない。G 展開の g = 0 の項に ν の差の因子 (ν² − 1)²、(qν²/t − 1)（W²ν の形）が残り、ν についての捻りも要る。
- 判断：これ以上の規則の探索は、既知の値（t = −1, t = q）に合わせる当てはめになる（研究の規則で禁止されている場当たり的な調整）。ここで止める。
- 結論【予想】：E と C の分離には、brane–image のセクター（被覆の環状の弦。orientifold で割ると境界を一つもつ Möbius の帯）の BPS 状態の多重項の構造を、幾何（または M 理論の指数）から直接導く必要がある。データからの推定では決まらない。AgSc の強い整数性 \hat N_{(□,□),g,β} ≡ \hat N^{c=1}_{□,g/2,β/2} mod 2 は、brane–image の状態と c = 1 の状態の対応を示唆しており、この導出の手がかりになる。

## 2026-10-20（続き 37）理論の導出：orientifold の開弦の refined 指数（ユーザーの指示「理論の導出に進め」）
文書：docs/refined_kauffman/theory_orientifold_index.md。スクリプト：scripts/rk_b_frozen_ansatz.py。ログ：logs/refined_kauffman/b_frozen_ansatz.log、b_frozen_ansatz_parity.log。
- 【導出】商の一粒子状態を、Ĩ で入れ替わる対（全モード）、Ĩ で不変な状態（(−1)^l ε = 1 のモードだけ、Σ = ½[1/(1 − t) + ε/(1 + t)]）、凍結した状態（モードなし）に分けると、g = f − ½f^{cov} = ½E/(1 + t) + C、E = Σ_{不変} ε m、C = Σ_{凍結} m。report §7 で Hopf から帰納した全 d の構造の導出。
- 【導出】E と C の分離は指数からは決まらない（凍結した状態 m ≡ 不変な状態 m と mt の組）。固定セクターの指数 P は (−1)^F 付きの和で打ち消しを含み、E はその台に収まらない（続き 34・35）。
- 【導出】E は p₂ 型の方向にしか入らないので、p₁ⁿ の方向（2 箱の ĝ_(2) + ĝ_(1,1)）と奇数のセクター（3 箱）は E によらない。そこでの符号の混在（5/10、6/18、2/15 など）は、E の不定性ではなく、凍結した状態の多重項の基底が正しくないことによる。
- 【計算】凍結した状態の展開の仮定（c = 1 は Z[G, w^±]、c = 2 は X × Z[G, w^±]、X ∈ {z_q, z_t, (z_q ± z_t)/2}）：[1] の c = 1 だけが成立。[1] の c = 2 はどの X でも割り切れず（余分の因子は普遍因子 (𝐭² − 1)(𝐚² − 𝐪²)）、2 箱・3 箱は c = 1 部分から G 展開できない。
- 手がかり：AgSc §4 の強い整数性の導出では、向き付けられない面の取っ手の因子は「差 × 和」(q^{1/2} − q^{−1/2})(q^{1/2} + q^{−1/2})（反正則な対合の捻ったトレース）。refined 版の因子を幾何から決めることが残る課題【未解決】。
- 追記（2026-10-20、ユーザーの指示「texに今回の結果を追記して報告せよ」）：TeX の §12「一般の N の 3 箱・T(2,5) と、orientifold の指数の理論」に続き 32〜37 を加え、要旨、§13 の表、§14 の未解決を更新した。

## 2026-10-20（続き 38）文献の閉弦の refined クロスキャップ因子との照合（ユーザーの指示「1を進めよ」）
文献（arXiv から取得、references/）：Avetisyan–Mkrtchyan 2205.12832、2107.08679、2302.14319、Kim–Wei 2201.12264、Piazzalunga–Uranga 1405.6019。読んだのは 2205.12832 §2〜§3、1405.6019 §4.2、2201.12264 の序論、AgSc 1202.4456 §3〜§4・§6（いずれも TeX）。
スクリプト：scripts/rk_b_closed_orientifold.py、scripts/rk_b_shift_direction.py。ログ：logs/refined_kauffman/b_shift_direction.log。テスト：tests/test_rk_closed_orientifold.py（3 件、1 秒）。文書：theory_orientifold_index.md §5 の訂正と §8。
- 【計算】AgSc §6 の閉弦の式（½ 向き付けられた項 ＋ RP² の項 ＋ (t^{1/2} + t^{−1/2}) の項）を、D_N の積表示から作った −log Z の d 番目と比べた。N = 2..7 で、差は N の 1 次式（古典的な項と定数写像）だけで、t^N を含む項はすべて一致する。
- 【計算】AM の D 型の式は、全 n の凍結した状態の項と不変な状態の項の和に、n ごとの有理式として恒等的に等しい。
  - AM の式：c = 1（奇数の n）＋「Klein 瓶」（c = 2、全 n、分子 q/t − t/q）。
  - 不変な状態の項：−½ e^{nτ₀}/((Qⁿ − Q⁻ⁿ)(Tⁿ + T⁻ⁿ))。½F_A と合わせると z₂ の偶数のモードだけになる。
  - 結論【導出＋計算】：閉弦の refined SO(2N) は g = ½E/(1 + t) + C の分類で厳密に書ける（P¹ は ε = +1 の不変な状態、RP² は凍結した状態）。AM の「Klein 瓶」は独立な多重項ではない。奇数の d で整理したときの −½ψ²C と ½E/(1+t) の差である。
- 訂正【文献の読み直し】：続き 37 と theory §5 の「向き付けられない面の取っ手の因子は差 × 和」は誤り。
  - AgSc §4 の (q − q^{−1}) は、類 β/2 の c = 1 を d を倍にして並べたときの ψ²z である。
  - PU §4.2 によれば、被覆の取っ手一つの捻ったトレースは同じ変数で z 一つ（−2 sinh）。
- 【導出】refined の取っ手の因子は G = z_q z_t = T₊T₋ と分解できる。
  - T₊ = (q^{1/4} − q^{−1/4})(t^{1/4} + t^{−1/4})、T₋ はその q ↔ t。t = q では z になり、σ で T₊ ↔ −T₋。
  - しかし T₊, T₋ の多項式環は G の環と (V − V^{−1}) から作る環と同じなので、σ 対称でない多箱の c = 1 部分を展開できるようにはならない。
- 【計算】述べ直しの判定式のずれ ν → νw^k（w = (q/t)^{1/2}）は、[1]、2 箱、3 箱のすべてで k = +1 に一意に決まる（k = −3..3 を走査）。
  - 閉弦の凍結した状態と向き付けられた状態の比 x²/Q0 = w^{−1}（AgSc の規約）は、開弦の ν~²/𝐚_U² = w³ と一致しない。
  - 開弦と閉弦の変数の対応が決まっていないので、この比較は保留する。
- PU の解釈【文献】：凍結した状態の数は orientifold の偶奇で重みを付けた差 ñ = n⁺ − n⁻ である。
  - 【予想】[1] の c = 2 の (𝐭² − 1) = −𝐭(w − w^{−1}) は、SU(2)_R の二重項の捻ったトレースの形をしている。正値性は ñ ではなく n⁺, n⁻ のそれぞれについて問うことになる。これは指数からは分けられない（E と C の分離と同じ形の不定性）。
- 結論：閉弦の照合で、E と C の分類（全 d の構造）は閉弦の厳密な結果でも確かめられた。しかし開弦の E と C の分離や、凍結した状態の多重項の基底を決める新しい規則は得られなかった。閉弦では状態が二つだけで、打ち消しがないため分離が一意になる。
- Kim–Wei 2201.12264 は ON 面の 5-brane web の refined 頂点で、向き付けられない面は数えないので、照合の対象外とした。O5 面の refined 頂点は無い（同論文の序論）。
- 次の候補：
  - (i) 被覆の Jacobian と SU(2)_R への Ĩ の作用を、具体的な曲線（Hopf の brane–image の環状の弦）で幾何的に決め、n⁺, n⁻ を分ける。
  - (ii) T(2,5) の 3 箱で一般性を確かめる。
  - (iii) 確立した結果のまとめ。

## 2026-10-20（続き 39）固定セクターの多重項ごとの捻ったトレース（ユーザーの指示「1を進めよ」）
スクリプト：scripts/rk_b_twisted_multiplet.py。データ：data/refined_kauffman/daha/twisted_multiplet_m{1,2}.pkl。ログ：logs/refined_kauffman/b_twisted_multiplet_{hopf,m1,m2,hidden}.log。テスト：tests/test_rk_twisted_multiplet.py（results、15 秒）。文書：theory_orientifold_index.md §9。
- 方針：PU の機構（不変な状態は内部の空間への Ĩ の捻ったトレースで寄与する）を、固定セクター P = s z_q Pk の多重項に当てた。
  - E = s X Σ e_m b̃_m とおき、|e_m| ≤ |c_m|、e_m ≡ c_m (mod 2) のもとで条件 (i)(ii) の整数計画を解いた。
  - 捻りの規則【仮説】：q 平面の π 回転、t 平面の π 回転、単項式ごとの符号。境界の因子 X：z_q または z_q^+。
- 【計算】Hopf：(2,2), (2,4) の既知の E は X = z_q と単項式ごとの符号で一意に再現される。(4,4) は (i)+(ii) を満たす解が多数ある。
- 【計算】結び目：X = z_q ではすべての規則で解なし（目標の値が t = q でも t = −1 でも (1 + q) を含むため）。
  - X = z_q^+ では、trefoil は単項式ごとの符号で一意に解ける（規則性は見えない）。
  - T(2,5) はどの規則でも解なし。上限を緩めると、隠れた対の最小数は単項式の符号で 2、q 平面の回転で 5。
- 結論【計算＋予想】：
  - 不変な環状の弦では、Ĩ が二つの境界を入れ替えるので、境界の因子は和 z_q^+ になる（Hopf では対が二つで (z_q^+)²）。
  - 指数の多重項を上限にした捻ったトレースでは、E は一般には決まらない（T(2,5)）。打ち消しの前の固定セクターの状態の空間を幾何から求める必要がある【未解決】。
  - 隠れた対の数を最小にする解を採るのは当てはめになるので、採らない。
- 次の候補：
  - (i) 固定セクターの状態の空間を、指数でなく数として得る方法を探す。例：被覆の composite (1;1) の refined の composite DAHA の超多項式を、DGR 型のホモロジーの Poincaré 多項式（正の係数）として読む。
  - (ii) ここで止め、確立した結果をまとめる。

## 2026-10-20（続き 40）固定セクターの状態の数を正の多項式から得る試み（ユーザーの指示「(i)を進めよ」）
スクリプト：scripts/rk_b_adj_positivity.py。ログ：logs/refined_kauffman/b_adj_positivity.log。文書：theory_orientifold_index.md §10、literature.md §4。
- 【計算】被覆の随伴 (1;1) の reduced な値 R_adj（大きな N、GL の DAHA）は、KN の太字変換を鏡像で当てても係数の符号が混在する。
  - trefoil は 19 項中 7 項、T(2,5) は 56 項中 25 項が負。
  - 単項式の偶奇の符号の規則（64 通り）では、どれも係数がそろわない。
- 【文献、TeX】Cherednik–Elliot 1503.01441 の序論に、composite の DAHA 超多項式は正値性をもたず、annulus のホモロジーとの関係は予想しないとある。
- 検索の範囲では、トーラス結び目の composite の環状 HOMFLY ホモロジー（Queffelec–Rose–Sartori 1802.04131 の構成）の計算例は見つからなかった。composite の正の表示は Hopf 絡み目（2003.07836）だけ。どちらも抄録のみで、未確認。
- 結論：refined CS の composite の側からは、打ち消しの前の状態の数は得られない。ホモロジーの側の候補は、refined CS の値との関係が確立していない。E と C の分離（固定セクターの状態の空間）は、現時点の文献と道具では決められない【未解決】。
- 次の候補：
  - (a) 確立した結果のまとめに移る。
  - (b) 逆向きの 2 本のケーブル（K ∪ K^r）の三重次数の KR ホモロジーを計算し、t = q で被覆の □⊗□̄ の値と合うかを見る。refined での対応は仮説になる。重く、見込みは低い。

## 2026-10-20（続き 41）論文の草稿（ユーザーの指示「(a)を進めよ」）
場所：docs/refined_kauffman/paper/（outline.md、claims.md、ja/main.tex と sections/、ja/main.pdf 26 ページ）。
- outline.md：論文の軸（定式化・構造・計算・普遍因子・正値性の障害）、決定事項（案）、章立て、図表。決定事項はユーザーの確認待ち。
  - 案の中身：1 本の論文、日本語で準備して英語で最終版、単著、正値性は未解決として障害を定式化する。
- claims.md：主張を 6 群（枠組み S、指数 O、Hopf H、[1] T、一般の N D、正値性 P）に分け、種別・前提・根拠（テスト名）・新規性を付けた。
- TeX 草稿の構成：§1 序論（結果の表）、§2 準備、§3 orientifold の指数（命題 3 つ、閉弦の確認）、§4 unknot と Hopf、§5 トーラス結び目の [1]（判定式）、§6 一般の N（DAHA、補間の表、普遍因子とずれの向き）、§7 正値性の現状と障害、§8 議論、付録（生成した式、コードとテスト）。
- 参考文献は arXiv API で書誌を確認した（journal_ref・DOI）。無いものは「掲載誌未確認」と書いた。KNTZ 1902.02275 は Kameyama–Nawata–Tao–Zhang。
- 付録の式は tex/gen/（scripts/rk_report_tex.py の生成物）をそのまま読み込む。
- 次：ユーザーの確認（決定事項、題）。その後、本文の推敲、英語版、必要なら図の追加。

## 2026-10-20（続き 42）論文の草稿の検証と訂正（ユーザーの指示「その方針で進めよ」）
- outline.md の決定事項を「確認済み」にした（1 本、日本語で準備して英語で最終版、単著、正値性は未解決として障害を定式化、題）。
- テスト：tests/test_rk_*.py の既定の 58 件（72 秒）と results の 6 件（437 秒）がすべて通った（logs/refined_kauffman/paper_tests_{default,results}.log）。
- 閉弦の照合（N = 2..7）のログを残した（logs/refined_kauffman/b_closed_orientifold.log。テストは N = 2..5）。
- 別のエージェントに、草稿と記録（report.md、report.tex、theory、progress、ログ、テスト）を照合させた。確実な誤り 7 件と要確認 14 件があり、記録で確かめてから草稿・claims.md・theory 文書を直した。
  - 確実な誤り：
    - T_± の環の記述（G の環より大きい。言えるのは「σ 対称でない c = 1 部分は展開できない」ことだけ）
    - AM の Klein 瓶の分子の変数
    - brane–image で混在するブロック（trefoil は ν^{−8}, ν^{−2}、T(2,5) は中央の 4 ブロック）
    - 状態の部分和の計算は trefoil だけ（T(2,5) は未実行）
    - 3 箱の (1,1,1) は N = 3 だけで照合
    - [1] の商の符号は T(2,3), T(2,5) だけで確認
    - Mimachi の 6 通りの列挙
  - 主な要確認の反映：
    - 命題 EC を「E の 2(1+t) を法とした類だけが決まる」に直した。状態の同値は f と f^{cov} の両方を保つ形にした。
    - theory §3 の「凍結 m ≡ 不変 m と mt の組」は係数が合わない誤りだったので訂正した。
    - prop:alld の後半（p₂ 方向、E_(2) = −E_(1,1)）は注意に移し、計算で確かめた性質と明記した。
    - その他：Hopf の 25 組の照合の数え方、E_22 の一意性の前提（P の係数を状態の数とみなす仮説）、過大な表現（「示す」「決まらない」）、ν~ を AM の量子シフトの開弦版とする記述（予想に弱めた）、GL の単項式の規則の範囲、二次量子化の割合、Cherednik の引用、二重ツイスト結び目の範囲。
- 草稿は 26 ページで、コンパイルは警告なし。

## 2026-10-21（続き 43）deep research の報告の配置と、文献の追加（ユーザーの指示）
- 報告（ユーザーが Downloads に置いたもの）を docs/refined_kauffman/literature_deep_research_2026-10-21_AB.md に配置した（指示は deep_research_prompt_AB.md）。要点と検証できる予言を literature.md §5 にまとめた。報告は外部の出力で、原典を確認するまでは未検証の要約として扱う。
- 文献の一覧に Witten 1101.3216（Fivebranes and Knots）を追加した（literature.md §4、justfile の papers、references/README.md）。references/1101.3216/ に PDF と TeX を取得した。本文は未読。
- 次の候補（判断待ち）：
  - 報告の予言 (A-i)（[2,1] だけ σ 対称）、(A-ii)（U = t で t = q の値に一致）、(A-iii)（c = 2 の商の回文性）、(B-i)（unknot の多箱で E が一意）は、保存済みのデータで短い計算で検証できる。
  - Young 1603.05401（orientifold DT、CoHM）の原典の確認。
- 追記（ユーザーの説明）：Witten 1101.3216 は、KN で BPS の構造を理論的に考える際に参考にした重要論文。KN の参考文献には入っていない（bbl で確認）。本研究の理論の側（Ĩ の作用、凍結した状態）で参照する。literature.md §4 に記録。

## 2026-10-21（続き 44）Witten → KN の流れの SO(2N) 版と、結果の再考察（ユーザーの指示「2を進めよ。Witten → KN の流れを SO(N) 版にするつもりで考察し、これまでの結果を再考察せよ」）
原典で読んだもの：Witten 1101.3216 §3.4・§4.1・§5.1.4–5.1.5・§5.4、KN 1703.05408 §2–§3、AgSc 1202.4456 付録 A、Young 1603.05401 序論・§3.4。文書：theory_orientifold_index.md §11。スクリプト：scripts/rk_b_boundary_mode_rule.py。ログ：logs/refined_kauffman/b_boundary_mode_rule.log。テスト：tests/test_rk_boundary_mode.py（2 件）。
- 【導出】SO(2N) 版の trace：D 型の (2,0) 理論は M5 を X/ℤ₂ に巻いたもの（Witten §5.1.4、AgSc 付録 A）。被覆の射影 P₊ = (1+Ĩ)/2 から g = ½ Tr_{H₁^{cov}}(Ĩ·) + C。E は Witten の (zelbow) の「ζ を挿入したトレース」の型、C は商の理論にしかない状態。E と C の分離の不定性は、一つの数の二つの表示への分け方の不定性。
- 【導出】境界ゼロモードへの Ĩ の作用：Ĩ は L と I(L) の境界を入れ替えるので、対ごとに和が偶、差が奇。Ĩ 不変な曲線（境界の対 n 組）では、E の境界因子は z_q^{n−1}(z_q^+)^n、挿入なしの P は z_q^{2n−1}。系：E は brane と image の境界を同数もつセクター（p₂ 型）にだけ現れる。
- 【計算】規則の検証（保存済みのデータ、すべて一致）：Hopf の E_22, E_24 は z_q¹(z_q^+)²、P_22, P_24 は z_q³(z_q^+)⁰（n = 2）。結び目 T(2,3), T(2,5) の P は z_q¹(z_q^+)⁰（n = 1）。続き 39 の「結び目の境界の因子は和 z_q^+」と、続き 13 の E_22 の (z_q^+)² は、この規則の帰結。
- 【導出＋予想】σ は q_ℓ = (qt)^{1/2} ↦ q_ℓ^{−1}、q_r 固定、すなわち Lefschetz SU(2)_ℓ の Weyl 反射。向き付けられない面には Lefschetz が無いので、凍結した状態 C に σ 対称性を期待する理由は無い（多箱の c = 1 部分が G 展開できないことの説明）。[1] の c = 1 が σ 対称であることは未説明【未解決】。(𝐭² − 1) = q_r² − 1 は SU(2)_r 側の捻ったトレースの形【予想】。
- 【文献＋解釈】Young の CoHM：M_Q が H_Q 上で自由でないことが E と C の不定性の代数的な形。loop quiver では決まる（unknot の多箱で E が一意になるはずという予言、未実施）。
- 再考察の表は theory §11.6。
- 次の候補：(a) 境界モードの規則を組み込んだ E の整数計画（モジュライ部分だけを未知にする）、(b) unknot の多箱で E が一意か、(c) [1] の c = 1 の σ 対称性の理由。

## 2026-10-21（続き 45）§11 の (a)(b)(c) を順に実行（ユーザーの指示「全て順番に進めて報告せよ」）
文書：theory_orientifold_index.md §11.7。スクリプト：rk_b_twisted_multiplet.py（取っ手の捻りを追加）、rk_b_framed_unknot.py、rk_b_strong_integrality_refined.py。ログ：b_twisted_multiplet_handle.log、b_framed_unknot.log、b_strong_integrality_refined.log。
- (a)【計算】境界の因子を z_q^+ に固定し、取っ手の捻り G → X（5 候補）で E の整数計画を解いたが、T(2,3)・T(2,5) とも (i)+(ii) の解なし。モジュライ部分への Ĩ の作用は指数の多重項の符号の付け替えでは表せない。未知は H*(M) への Ĩ の作用に絞られた。
- (b)【計算】unknot：framing 0 では E = C = 0 で予言は自明。framing x = 2 では (i)+(ii) の解が一意（E = (ν²+1)(q+1)(ν²q−t)/t、境界因子 z_q^+ は §11.3 どおり）。x = 1 は解なし（条件 (ii) の形と奇数の framing の符号が未確定、保留）。Young の「loop quiver では決まる」と整合。
- (c)【導出】refined の強い整数性は E ≡ P_inv (mod 2)（整数計画の制約そのもの）。【計算】AgSc の unrefined の形を素朴に当てた mod 2 の検査は t = q でも合わず、辞書が未確定なので判定できない。[1] の σ 対称性は未解明のまま（予想のみ）。
- 次の候補：(1) H*(M_{g,h,β}) への Ĩ の作用を、被覆の composite (1;1) のモジュライ（brane–image の環状の弦）の幾何から決める。(2) framing 1 の unknot の条件 (ii) を、Mariño–Vafa の framed unknot の unrefined の値から固定する。(3) 論文の草稿に §11 の内容（trace 表示、境界モードの規則）を加える。

## 2026-10-21（続き 46）§11.7 の次の候補 1〜3 を順に実行（ユーザーの指示「順に実行せよ」）
- 1・2【計算】framed unknot（x = 1..4、条件 (ii) の 3 つの形）：P(x) = (ν²−1)(q−1)Σ_k (ν²q)^{x−1−k}t^k の積の構造。偶数の x で E が一意（E(x) = (ν²+1)(q+1)Σ_k(−1)^k(ν²q)^{x−1−k}t^k）、奇数の x は解なし（(ii) の形は未解決）。Ĩ は塔の k 番目に (−1)^k で働く（e^{iπJ₂}）。
  - 符号の規則の検証（rk_b_sign_rule.py）：unknot・Hopf (2,2) で成り立つ「t と ν の指数の偶奇」の規則は、Hopf (2,4) と結び目では成り立たない。
  - 結論：E の決定には固定セクターの状態の t 電荷を J₂ と S_R に分ける情報が要り、指数からは得られない（theory §11.8）。
- 3：論文の草稿に §11 の内容を反映（orientifold_index.tex に trace 表示と境界ゼロモードの規則、positivity.tex に強い整数性 E ≡ P (mod 2)・framed unknot・符号の規則の試み、claims.md に O9〜O11・P7・P8）。

## 2026-10-21（続き 47）§11.8 の残り 3 件を実行（ユーザーの指示「残りも全て実行せよ」）
文書：theory_orientifold_index.md §11.9・§11.10。スクリプト：rk_b_framed_unknot.py（enumerate_E）、rk_b_combined_rule.py、rk_b_one_box_vs_brane_image.py。ログ：b_framed_unknot.log、b_combined_rule.log、b_one_box_vs_brane_image.log。テスト：tests/test_rk_framed_unknot.py（x = 1 の全数を既定のテストに追加）。
- (A)【計算】奇数の framing の unknot：(ii) の形を仮定せず「t = q で C が整数」だけを課すと、x = 1, 2, 3 で E は一意。E(x) = (−1)^x(ν²+1)(q+1)Σ_k(−1)^k(ν²q)^{x−1−k}t^k。奇数で解が無かったのは全体の符号 (−1)^x の違い。
- (B)【計算】符号 ε(β,j) と取っ手の捻り X を組み合わせた 36 通りは、T(2,3)・T(2,5) で (i) も (ii) も満たさない。KN の S_R 自明の仮定のもとでは Ĩ の符号は偶奇で決まるはずなので、結び目の brane–image のセクターでは S_R が非自明に働くと考えられる【予想】。
- (C)【計算】[1] の c = 1 の ν^β のブロックの G の次数は、brane–image の ν^{2β} のブロックの G の次数のちょうど半分（全 β）。[1] の凍結した状態は brane–image の環状の弦の商で、商の取っ手が G を担う。[1] の σ 対称性はこれで説明できる【予想】。多箱の非対称性は未説明。
- 論文：§3 の σ 対称性の小節に (C) を、§7 に (A)(B) を反映。claims.md に P7 の更新と P9。

## 2026-10-21（続き 48）(J₂, S_R) の分離のための案（ユーザーの指示）
- 文書：docs/refined_kauffman/proposals_SR_action.md。
- 【導出】指数の精密化では J₂ と S_R を分けられない（超電荷と可換な電荷は J₁−S_R と S_R−J₂ で尽きる）。必要なのは BPS 状態の空間とその上の対称性の作用。
- 【導出＋計算】framed unknot の固定セクターは x = κ(1 − T_adj^x)（T₁² = T_adj）で、E は κ の三つの因子（境界 (1−q)、β=1 の (1−A)、随伴の極 1/(1−qA/t) の塔）を「励起した状態に −1」で捻ったもの。
- 案：A（κ の規則を組み込み未知を R_adj に限る、今すぐ）、B（Hilbert scheme と実構造）、C（四重次数）、D（Sp(2N) で C_{c=2} を分離）、E（B 型との比較）、F（固定点への局所化）。推奨は A → D → B。

## 2026-10-21（続き 49）案 A（κ の規則を組み込んだ E の整数計画）（ユーザーの指示「推奨の順に進めよ」）
- 【計算】x = R_adj(qn² − κ) + κ − W₁² の積の形に unknot の捻りの規則を当て、未知を R_adj の係数に限った整数計画は、Laurent 性の条件 (0)（ν² = −t/q で R̃_adj = 1）と (0′)（t = 1 の留数）が両立しない（T(2,3)・T(2,5)）。Ĩ の捻りは composite の分解の積に分配しない。案 A は不成立（theory §11.11）。
- 次：案 D（Sp(2N) の refined CS）。AM 2205.12832 (refCS) の定義は、任意のルート系で k_α = y を一様にとり q_α = q^{(α,α)/2} とする Macdonald (R,R∨) 型（付録 B：C_n は長い根の長さの 2 乗を 4 に正規化）。C_N の DAHA-Jones の実装がいる。

## 2026-10-21（続き 50）案 D：C_N の DAHA と refined Sp(2N) の [1]（ユーザーの指示「推奨の順に進めよ」の続き）
- 【計算】C_N の DAHA の多項式表現を実装（src/rcsknots/daha_cn.py）。関係式、u = 1 で D_N と一致、Cherednik の ℋ^C_{3,2}(ω₁) と一致。
- 【計算】ℋ^C_{2m+1,2}(ω₁; q,t,u,a) を N = 2..8 の値から復元（T(5,2) は新規）。
- 【計算】AM の refined Sp(2N) = C_N(u = t²)。[1] のクロスキャップ部分は T(2,3), T(2,5) で g^{Sp}(𝐚) = −g^{SO}(−𝐚)。c = 1（𝐚 偶数冪）の符号反転が refined でも成立。theory §11.12。
- テスト：tests/test_rk_daha_cn.py（3 件、1 秒）。
- 次：Sp の 2 箱（C_N の Macdonald 多項式）、案 B（Hilbert scheme と実構造）。

## 2026-10-21（続き 51）案 B の計画（ユーザーの指示「推奨の順に進めよ」の続き）
- 文献を取得して原典で確認：OS 1003.1568、ORS 1201.2115、DSV 1111.6533（references/、書誌は arXiv API で確認）。
- 辞書（ORS の重みの次数 ↔ −J₂、偏屈フィルトレーション ↔ J₁、a 次数 ↔ β）と、orientifold の実構造の模型、framed unknot での較正の対象を theory §11.13 に書いた。実構造つきの Hilbert scheme の計算は文献に無く、本セッションでは着手しない。
- 案 A〜D の結果のまとめ：A は不成立（積に分配しない）、D は成功（refined Sp(2N) の [1]、g^{Sp}(𝐚) = −g^{SO}(−𝐚)）、B は計画のみ。

## 2026-10-21（続き 52）案 D の 2 箱：refined Sp(2N) の trefoil の 2 箱（ユーザーの指示「1を進めよ」「続けよ」）
- 【計算】C_N の Macdonald 多項式 P_(2), P_(1,1) を実装（daha_cn.py の two_box_P、評価公式の前の因子を q^{−(ρ_k,b)} に修正）。検証：固有関数（C_2, C_3）、u = 1 で D_N の 2 箱と一致、Cherednik の ℋ^C_{3,2}(ω₂) と単項式を除いて一致。
- 【計算】N = 2..9（m = 1、最長 704 s）を x = t^N の有理関数に補間し、Sp（u = t²、ν = t^{N+1/2}）の large N の g_(2), g_(1,1) を得た（largeN_2box_sp_m1.pkl）。整数性は SO と同じ全 d の構造。**E^{Sp}(t=−1) = E^{SO}(t=−1)**。hatted の p₁² 成分で Sp(𝐚) = +SO(−𝐚)、p₂ 成分と g_R には ± も転置も成り立たない。theory §11.14、論文 §sec:sp、claims D13、tests/test_rk_sp_two_box.py（2 件、7 秒）。
- 【計算】写像の走査（rk_d_sp_so_map.py）：{id, σ, 鏡像, q↔t} × ν の符号・逆数 × 転置 × 単項式倍で結べるのは、[1] の g_R・f_R の (id, −ν, −1) と hatted の p₁² の (id, −ν, +1) だけ。2 箱の g_R・f_R と p₂ 成分は結べない。sympy の cancel が σ の写像で止まる（p₂、500 秒超）ため、判定を項の形の比較に書き換えた。
- 未解決：p₂ 方向の SO–Sp の関係（Adams 項 ψ² が 𝐚 → −𝐚 と可換でない）。m = 2 の Sp の 2 箱は N = 2..4 のみ（補間に足りない。N = 9 までで 1 N あたり最大 ~20 分の見込み）。
- 次：案 B（Hilbert scheme と実構造）、または論文の仕上げ。

## 2026-10-21（続き 53）m=2（T(2,5)）の Sp の 2 箱の見積もり（ユーザーの指示「再開せよ」、実行前の見積もり）
- 必要な N：SO の m=2 の 2 箱は R の x の窓が 11 項（随伴 x^{−18}〜x^{−8}）で Z は N=14 まで使った。Sp も同じ窓なら N=2..13（11 ノード + 照合 1 点）が要る。前回の「N=9 まで」は m=1 の窓（7 項）を流用した誤り。
- 実測（rk_d_cn_two_box_data.py、m=2、u=t²、1 コア）：N=2 0.6 s、N=3 5.3 s、N=4 43 s、**N=5 208 s**（1 段の比 8.8, 8.2, 4.8。m=1 の比 3.0〜2.0 の約 1.5 倍）。
- 外挿（比を m=1 の 1.5 倍とする）：N=6 ≈ 16 分、N=7 ≈ 1 時間、N=8 ≈ 3 時間、N=9 ≈ 10 時間、N=10 ≈ 1 日、N=11 ≈ 3 日、N=12 ≈ 8 日、N=13 ≈ 3 週間。N=13 までの合計は約 1 か月（並列にしても最長の 1 件で決まる）。**24 時間を超えるので実行しない**（研究上のルール）。
- 代替案：(a) m=2 を見送り、Sp の 2 箱は m=1 の結果（E の一致、p₁² の +）までとする。(b) 高速化：DAHA の多項式表現を (q,t,u) の数値（有限体または有理数）に特殊化して N ごとに計算し、係数を再構成する（flint_rf3 の有理関数の約分が重い部分。効果は N=5 で測ってから。実装に 1〜2 日）。(c) Cherednik の評価公式と P_λ の係数は m によらないので、m 依存は γ̂^m の作用だけ。γ̂ の作用を N の有理関数として直接当てはめる定式化（規模不明）。
- データ：cn_twobox_N5_m2_sp.pkl を保存（照合点として使える）。
- **決定（ユーザー、2026-10-21）：(a) m=2 の Sp の 2 箱は見送る。** Sp の 2 箱は m=1（trefoil）の結果までとし、論文もその範囲で書く。(b) の高速化は、多箱の σ 非対称性などで必要になった時点で着手する。

## 2026-10-21（続き 54）SO と Sp の 2 箱の差の同定（決定 (a) の後の次の作業）
- 【計算】f^{Sp}_R(𝐚) − f^{SO}_R(−𝐚) = ±ψ²[δ]（(2) は +、(1,1) は −）、δ = even_𝐚(g^{SO}_1) + 1 = (𝐚²−𝐪²)(𝐚²𝐭²−1)(𝐚²𝐪²𝐭²−1)/(𝐚⁶𝐪²𝐭⁴)。g でも同じ（被覆は共通）。p₂ 方向に規則が無かった理由は、[1] の規則が c=1 の符号だけを変えるのに対し Adams 項 ψ²[½(g_1+1)] が 𝐚→−𝐚 を見ないこと。定数 1 は unknot のクロスキャップと解釈【予想】。theory §11.15、論文 §sec:sp、claims D14、tests/test_rk_sp_so_adams.py。
- 次：(i) 定数 1 の由来（unknot の 2 箱で確かめる）、(ii) 多箱の c=1 の σ 非対称性、(iii) 論文の仕上げ、(iv) push（許可）。
