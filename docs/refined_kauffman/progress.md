# 進捗記録（第 3 の研究：refined Kauffman の large N duality）

## 現在の状況（再開用、2026-10-31 更新）

新しいセッションでは、この節を読めば作業を再開できるようにしてある。日々の詳細はこの下の日付ごとの記録（続き 1〜99）にある。

### 1. 環境
- **worktree**：`/Users/masaya/projects/ai_research/refined_cs/.claude/worktrees/rename-literature-survey-doc-9cc5fc`
  - ブランチは `claude/kauffman-research-progress-c23e34`（続き 63 以降）。続き 62 までは worktree `dgr-quadruply-graded-lomv-2a7fb6`（ブランチ `claude/refined-kauffman-large-n-00d58e`）。
  - main への直接の commit とローカルでのマージは禁止。
  - この worktree が無くなっていたら、push 済みのブランチから作り直す（main checkout で `git -C /Users/masaya/projects/ai_research/refined_cs worktree add <場所> claude/kauffman-research-progress-c23e34`、サンドボックス外）。
  - data/cache/（追跡しない）の n2_quotient のキャッシュ（twisted_n2_quotient*.pkl）などは、無ければ初回の実行で作り直される（数分）。
- **PR**：[kameyama/refined_lmov#11](https://github.com/kameyama/refined_lmov/pull/11)（続き 62 まで、MERGED、2026-10-26）。
  - 続き 63〜86 の commit（28 個）は 2026-10-30 にブランチ claude/kauffman-research-progress-c23e34 として origin に push した（先頭 4ba396b）。その後の push の記録の commit（90dc666）は未 push。PR はユーザーの指示でまだ作っていない（作成は毎回許可を得る）。
  - 前回の push は `git push` をサンドボックス外で実行した（gh api は不要だった）。
  - push は外部送信なので、毎回ユーザーの許可を得る。
  - SSH の `git push` はサンドボックスで通らない。前回は、許可を得たうえで gh api（Git Data API）で同じ SHA の commit を作った。
    - 手順：blobs → trees（base_tree は親）→ commits（作者・日時・メッセージを同じに）→ refs。
    - 各段で SHA がローカルと一致することを確かめる。
- **Python**：main checkout の `.venv` を使う（`/Users/masaya/projects/ai_research/refined_cs/.venv/bin/python`）。
  - worktree 直下から相対で書くと `../../../.venv/bin/python`。
  - スクリプトは worktree の `src` と `scripts/refined_kauffman` を sys.path の先頭に入れている。
- **テスト**：`../../../.venv/bin/pytest -q tests/refined_kauffman`（2026-10-26 に main の整理（PR #10）に合わせて scripts/refined_kauffman/、tests/refined_kauffman/、src/rcsknots/refined_kauffman/ に移した。スクリプトは `PYTHONPATH=src:scripts/refined_kauffman`）
  - 既定の 110 件（約 315 秒。2026-10-30 に全件通過）。
  - `-m results` の 8 件（2026-10-20 に 6 件通過、logs/refined_kauffman/paper_tests_results.log。2026-10-30 に追加した test_cover_parts_sum は 21 分で通過、logs/refined_kauffman/test_cover_parts_sum.log）。
- **実行上の注意**（2026-10-29〜30 に分かったこと）：
  - git の書き込み（commit・push）はサンドボックス外で `git -C <worktree> ...` として実行する。push は `git push` をサンドボックス外で実行すれば通る。
  - `/usr/bin/time -l` はサンドボックス内では最大メモリの行を出さず、終了コードが 1 になる。完了の目印は "real" の行にする（"maximum resident" を待つと終わらない）。
  - `ps`・`pkill` はサンドボックス外で実行する。
  - 30 分を超える計算は、実行前に見積もりと再開の準備を progress に書く（CLAUDE.md）。
- **TeX**：`TEXMFVAR=$TMPDIR/texmfvar lualatex main.tex` を 2 回実行する。
  - ページを画像にするには gs を使う（pdftoppm は無い）。
  - .aux/.log/.out/.toc は git で追跡しない。
- **文献**：`references/<arXiv id>/src/`。取得は `just fetch-refs`（justfile の papers に一覧がある）、または `curl`（allowed_domains: arxiv.org）。書誌の確認は arXiv API（export.arxiv.org）で行う。
- **Dropbox の 2018 年のノートブック**は読むだけで、書き換えない。

### 2. 現在の段階：論文の草稿（日本語の準備稿）
- 場所：`docs/refined_kauffman/paper/`
  - `outline.md`：論文の軸、決定事項、章立て
  - `claims.md`：主張の格付け表
  - `ja/main.tex`、`ja/sections/*.tex`、`ja/main.pdf`（31 ページ）
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
  3. 続き 63 以降の新しい PR の作成（ブランチは push 済み。PR の作成はユーザーの許可が要る。PR #11 はマージ済み）

### 2′. 進行中の研究の線：E の決定（続き 65〜99、theory §12.6〜§12.31）
- 目標：多箱の自由エネルギーの分解 g = f − ½f^{cov} = ½E/(1+t) + C のうち、E を状態の数と Ĩ の作用から決める。
- 入力（計算済み）：
  - trefoil の composite (2;2) 型（被覆、大きな N、続き 65〜68）。
  - SO の 4 箱の R_λ(U)（続き 70〜71）。
  - 条件 (i)（t = −1、SO の 4 箱）と (ii)（t = q）。
- 経過：
  - 指数の台の上の整数計画では解が無かった（続き 72・73）。
  - 固定セクターの因数分解へのテンソル積の捻りで解が出た（続き 74）。
  - N の下の系列 Sym^k V（V = 1 + t + qt）に対合 diag(1,1,−1)、上の系列に双対の対合を課すと、n = 1 の E が決まった。trefoil は (i) だけで、T(2,5) は (i)+(ii) で一意（続き 75）。
  - 同じ規則は n = 2 の N_λ にはそのままでは合わない。t = −1 だけで決まるのは 84/950 単項式で、入れ替わる対が 2845 個以上要る。下の端は Sym W（W = 1 + q² + q²t²）（続き 76）。
  - N_λ には ν 方向の系列が無く、ν のブロックの中に重み q, t, qt のセルの塔（1/((1 − q)(1 − t)(1 − qt)) の打ち切り）がある。全ブロックが鏡映 (q, t) → (t^{−1}, q^{−1}) で対称。λ = (2) と q^{−2}·λ = (1,1) は q の指数の小さい側で一致し、差は打ち切りの側に集まる。W = ψ²V′（続き 76 の「V(q², t²)」を訂正）（続き 77）。
  - 塔のセルの対角な対合（0, q, t, q+t）を部分ごとに選ぶ規則でも、入れ替わる対は最小 2052 で、規則なしの 20 に遠い。鏡映で対称な規則だけ・一様な規則・qt の縞ごとの 0 か q+t では解なし（続き 78）。
  - 規則なしで入れ替わる対が最小（20）の解では、950 中 308 単項式の e が一意に決まる。下の端の ν^{−24} 全体（(1 + qt) の qt の項だけ反転）と、中間ブロックの低い q の領域（λ によらない）。中間ブロックの決まった符号は部分ごとの対角な対合の形にならない（続き 79）。
  - 条件 (i) は 154 個の行（λ、ν、q の冪）の交代和 S = E_λ(q, −1, ν) だけを決める（階数で確認）。そのため対の最小化は行ごとの部分和の問題に分かれ、動的計画法で 1 秒未満に解ける（続き 79 の約 1 時間の整数計画と完全に一致）。決まった e に固有の構造は無く、行の値 S と行の中の |c| の組合せで尽きる。低い q への偏りは行が小さいこと、λ によらないことは行の問題が同じことによる（続き 80）。
  - (ii)（t = q）を加えると、割った形の条件はブロックごとに閉じ、交点 q = t = −1 では (i) と一致する。しかし N の台と |e| ≤ |c| では ν^{−18}〜ν^{−24} の 4 ブロックで解が無い。ν^{−24} では (ii) は (i) から一意に決まる値のちょうど符号違いを求める（因子 (1 − qt) のため交点では見えない）。解のある上のブロックでも、入れ替わる対の最小が大きく増える（ν^{−12} で 4 → 61）。n = 2 の (ii) の形が破れているという予想を補強した（続き 81）。
  - (ii) の導出を見直した。AgSc から導けるのは mod 2 の関係だけで、等式は「凍結した状態の二重被覆がすべて同じ符号 ε で入る」という仮定による。n = 1 の F の c = 2 の部分は t = q で 0 なので、n = 1 は ε の c 依存を見ない。c = 2 の符号を反転した (ii′) E|_{t=q} = −(1 + q)ψ²[F^{c=1} − F^{c=2}]（ε = (−1)^c、n によらない）では、n = 1 の解は変わらず、n = 2 の ν^{−18}〜ν^{−24} の矛盾が ν^{−22} の見えない対 1 個を除いて解消した（ν^{−24} は (i) の解と一致）。続き 72 の「一様な規則ではない」を訂正した。上のブロックの対の数の増加（ν^{−12} で 4 → 61）は (ii′) でも変わらない（続き 82）。
  - 上のブロックの対の増加を切り分けた。(ii′) だけ・(i) だけの最小はどちらも小さく（ν^{−12} で 26・4）、増加は両者を同時に満たすことから来る（61）。台を箱に広げて見えない対を許すと一部減る（61 → 39）が主因ではない。ν^{−12} では t = q で N_λ = −X_λ²、(ii′) の目標 H′ = Σ_λ ψ²[X_λ]（X ⊗ X の入れ替えのトレース）で、入れ替えの描像は対 148 個、λ ごとの条件と (i) の両立で最小 94。対の最小の原理のほうが自然でない、と予想した。ν^{−22} の見えない対の候補 4 個はすべてブロックの下端の角にある（続き 83）。
  - 平方の構造は ν^{−12}・t = q に限られる。ν^{−14} 以下の N_λ(q, q) は平方でも X_λ の倍数でもなく、ν の級数全体の平方根は Laurent でない。一般の t では N_λ は既約である。X_λ は g_λ の t = q での ν^{−3} の係数と単項式の違いしかない（X_2 = (q + q^{−1})X_11）（続き 84）。
  - N を積の部分 −Y_λ(ν)Y_λ(−ν)（image の側に符号 (−1)^c）と残りに分けた。(ii′) の目標は [3]_q で割り切れ、(q + q^{−1}) : 1 で λ ごとに分けられる。積の部分は ν^{−12} を完全に、ν^{−14} を大部分説明するが、ν^{−16} で残りに奇数の係数が現れ、λ ごとの目標は偶奇で解なしになる（続き 85）。
  - 合成表現の相関関数の項（R_b K_b、c₁W_adj、c₀）と低い次数の積は、どれも単独では整数の状態の数にならず、4 項の和で初めて Laurent になる。大きな N の展開の最高次は c₀・低次（ν^{−6}）、c₁（ν^{−8}）、合成表現（ν^{−10}）の順。Koike の分解は指標の分解で、状態の分解としては読めない（続き 86）。
  - T(2,7)（m = 3）の 2 箱のデータを作り（約 3 時間、大部分は 15 項の補間の約分）、n = 1 の符号の規則を照合した。**続き 75 の規則（k ≤ 2m で (−1)^q、k > 2m で (−1)^t）は T(2,7) では (i) だけでも (i)+(ii) でも解を持たない**。ブロックと符号ごとに 0/q/t/qt を選ばせても (i)+(ii) では解なし。破れるのは上の系列の最初のブロック k = 2m+1 だけで、下の系列と移り目は q、k ≥ 2m+2 は t の規則を保つ。下の規則＋上は一般の e でも、(i)+(ii) で E は一意に決まらない（k = 7 の 7 変数、実現できる組 3 通り、対の最小 15）（続き 88）。
  - **n = 1 の E を状態の空間の同定から作った（続き 89、theory §12.23）**：下の系列 Sym^k V（ε = (−1)^{a+b}）、移り目 Sym^{2m}V − S ⊗ S（S は簡約 superpolynomial、Ĩ は入れ替え、トレース ψ²S）、上の系列 −t^{K−2j}[Sym^j V + (q − t)R_j]（双対の対合 ε = −(−1)^{b+c}、補正 R_m = t^{m−1}、R_{m+1} = t^{m−2}(S′ ⊗ S′ + t²) に入れ替え）。条件 (i)・(ii) を使わずに、T(2,3)・T(2,5) の一意な E を全単項式で再現し、T(2,7) では (i)・(ii) をともに満たす（予言の的中、対 16）。
  - n = 2 の trefoil の N_λ に同じ同定を当てた（続き 90、theory §12.24）：上の端は B_0 = (t − q)/q、B_1 = B_0·VV′/q²（brane ⊗ image）だが、2 段目から ν 方向の Fock 型でない（Fock 型の項を引いた残りが B_2 と同程度）。各ブロックは qt = −1 でほぼ消え（残差 1% 程度）、(1 + qt) の因子とブロックの中の塔が主である。
  - N を式から読んだ（続き 91、theory §12.25）：n = 1 では κ の項が Sym^k V の系列、W₁W₁′（brane × image の積）が移り目の −S ⊗ S、(R_adj − 1)K_adj と枠 w が上の系列と (q − t) の補正を担う。枠なしの unknot は n = 1 でも n = 2 でも N = 0。n = 2 の塔は枠付き unknot の参照 N^{univ}_λ（分母 (q∓1)³(t∓1)³(q²+1)²(t²+1)²、ν のセル {t, q⁻¹, t²q⁻¹, tq⁻²}、(√q − √t) の因子）に由来し、R − 1 の項が打ち切る。
- **再開の時点の状態（2026-10-31、続き 99 の後）**：ユーザーの指示「これから 24 時間は自己判断で研究せよ」（2026-10-31 15 時ごろ）のもとで進めている。n = 1 の E の原理が T(2,3)〜T(2,13) の全ブロックで検証された（続き 95〜99、theory §12.29〜§12.31、claims P17〜P22・D15）：補正 R_{m+j′} = t^{m−1−j′}F_{j′}（母関数、j′ ≤ 4 で検証）、Ĩ の対角の部分 (−1)^{m+j′+1} ψ²(Sym^{j′}V + t·Sym^{j′−1}V)、残りは対。T(2,7)〜T(2,13) の E は予言として保存（twisted_E_n1_m{3..6}.pkl）。動いている計算は無い。次の候補：(A) n = 1 の原理の物理（F・X の母関数の由来）と論文の §7 への反映、(B) n = 2 の ν^{−14} 以降（t = q で正、[3]_q で割り切れ、平方でない）の同定、(C) n = 2 の λ = (2) の補正の Ĩ（平方の分け方の変更）。push と PR は許可待ち。
- 次の候補（続き 86 の末尾）：
  - (2″) 低いブロックで解を選ぶ原理。積の構造は ν^{−14} まで（続き 85）、相関関数の項は状態の分解にならない（続き 86）ので、この線はいったん区切る。ν^{−22} の見えない対の原因（(ii′) の形か台の端か）。
  - (2‴) ε = (−1)^c の物理的な導出（Ĩ の H¹ のゼロモードへの持ち上げの符号）を詰める。T(2,5) の n = 2 での (ii′) の照合は、SO の 4 箱の m = 2 が現行の D_N の DAHA では届かない（続き 87）ので、計算法を変える（続き 70 の案 1：q を有限体の数値、t の一変数）場合だけ実行できる。
  - (3) 中間ブロックの分子（×(1 − q)(1 − t)(1 − qt) の後）と塔の打ち切りの規則を同定する。
  - (4)・(5) は済み（続き 88・89）：符号の規則は T(2,7) で破れるが、状態の空間の同定（theory §12.23）で書けば (i) も (ii) も要らず、T(2,7) の E も決まる。
  - (6) n = 2 の trefoil の N_λ の同定：ν 方向の Fock 系列は 2 段目で破れる（続き 90）。(6b) 式からの導出（続き 91）で、N_λ = 枠付き unknot の参照 N^{univ}_λ ＋ R − 1 の項、と分かった。(6c) 項ごとの分解（続き 92、theory §12.26）：Wilson loop なしの背景（κ と c₀）は閉じた形をもつ 6 セルの Fock 空間（n = 1 の Sym^k V の n = 2 版、塔なし）、結び目の部分は λ = (2) で ν^{−10} から −(q − t)(1 + qt)⊗(1 + qt)。Ĩ の仮説「背景のセルに対角（n = 1 と同じ符号）、積に入れ替え」は上端 2 ブロックで (i)+(ii) を満たし E を一意に決める。(6d) ν^{−12} への延長（続き 93、theory §12.27）：結び目の部分は λ = (1,1) で −Y⊗Y − (q − t)t(1 + qt)⊗(1 + qt)（Y = (1 + q)(1 + t)(1 + qt)）、λ = (2) で −Y₂⊗Y₂ − (q − t)R（Y₂ = (1 + qt)Y、R は 17 状態）。Ĩ は平方に全因子の入れ替えで、λ = (1,1) は (i) だけで 187 通り中 1 通りに決まり、λ = (2) も Y によらず入れ替えだけが可能。(q − t) の補正の符号は (ii) に見えず（t = q で対が消える）、(i) だけでは決まらない。下の端 ν^{−22} は W の Fock ではない。(6e) 補正の符号の原理（続き 94、theory §12.28）：条件 (i)・(ii) からは決まらない（深い段の行を加えても 17 中 3 個、(ii) は対を見ない）。n = 1 の「状態 t·m を q·m に移す」読みは 18 通りの Y すべてで多重度が足りず、入れ替えの下では状態の ε が定義できない。λ = (1,1) の答えは「(1 + q), (1 + t), (1 + qt) の積の平方の符号付きで最小の和（3 項）に入れ替え」から一意に出るが、λ = (2) では最小の分解（基底 4 通り、背景の符号 4 通り）がどれも (i) と両立しない（最小 22 に対し (i) 整合は 26）。(1 + qt) の因子に対角な Ĩ も不可能。**この線はいったん区切り、同じ型の原理の検証は n = 1 の T(2,9)（候補 (7)）に移す。**
  - (7) T(2,9) の 2 箱のデータを作り、R_{m+1} の形と j = m+2 の補正の構造を検証する。**データは作った（続き 94）**：D_N の 2 箱 N ≤ 26（約 1 時間）、随伴 N_A ≤ 34（約 1 時間）、R_λ(U) は rk_c_bounded_save.py（r の次数ごとに U の冪の範囲を絞る補間、flint の約分）で (2) は N = 3〜8、(1,1) は N = 4〜10 から一意に決まり、使わない N ≤ 26 の全部で照合（sympy の Newton 補間なら日の単位の見込みだったものが 40 分）。largeN_2box_m4.pkl まで作成済み（g の整数性の検査は通過、固定 N の g は無いので U = t², t³ の照合は飛ばした）。**検証した（続き 95、theory §12.29）**：R_{m+1} = t^{m−2}(V ⊗ V + t²) は V = 1 + t + qt 固定で、続き 89 の S′ = T(2,2m−3) の読みは T(2,9) で破れる（訂正済み）。j ≤ m+1 の原理は (i)+(ii) と整合する。j = m+2 の補正 R₆ = t·P（35 状態、P は V³ を含まない）の Ĩ は (i)・(ii) に決まらず（2/15）、一様・入れ替えの仮説は破れる。補正の状態の数は 1, 10, 35 = C(2j′+3, 3)。(7b) T(2,11)（続き 97、theory §12.30）：F₁, F₂ は m によらず、R_{m+2} の Ĩ は予言どおり (i)+(ii) と整合、F₃ は 84 状態。母関数の予想（correction_F）を立てた。(7c)（続き 98〜99、theory §12.31）：j′ ≥ 3 の Ĩ も決まった——対角の部分は (−1)^{m+j′+1} ψ²(Sym^{j′}V + t·Sym^{j′−1}V)、残りは対。**T(2,13) で F₄（165 状態）と E 全体（対 140）が予言どおり (i)+(ii) を満たした。** n = 1 の原理は T(2,3)〜T(2,13) で完成（claims P17〜P22）。残るのは物理的な由来と論文への反映。
  - (8) 補正 (q − t)R_j の由来（双対性の破れ）を理論的に詰める。
- コード：rk_f_twisted_tensor.py（n1・n2・rule・full・lower・blocks・scan・n1choice・n1upper・n1freeset・principle・n2fock・n1terms・n2rule・n2series・n2tower・n2choice・n2pattern・n2rows・n2ii・n2iic・n2split・n2h22・n2square・n2skew・n2parts。rule・full・lower・blocks・scan・n1choice・n1upper・n1freeset は m を引数に取る）、rk_f_twisted22_lam.py（parts）、rk_f_n2_univ.py（univ・w1：unknot の参照、続き 91）、rk_f_n2_rterms.py（terms・analyze・background・fock・hyptop・hyptop2：R − 1 の項ごとの分解と背景の Fock 空間、上端の Ĩ の仮説の照合、続き 92；hyp12・hyp12ii・hyp12b・hyp12bii・hyp12c・hyp12ci：ν^{−12} への延長、続き 93；hyp12d・hyp12e・hyp12f：補正の符号の原理の候補、続き 94）、rk_c_bounded_save.py（絞った補間で 2 箱の R_U を保存、続き 94）、rk_f_twisted_tensor.py の principle に extra_R（j ≥ m+2 の補正）と未同定ブロックの整数計画、principle_R（principleR：補正の対の符号だけ自由、続き 95）。
- データ：data/refined_kauffman/daha/twisted_E_n1_m2.pkl（T(2,5) の一意な E）、twisted_n2_pattern.pkl（n = 2 の対が最小の解と e の範囲）、twisted_E_n1_m3.pkl（T(2,7) の E、原理からの予言、続き 89）、T(2,7) の 2 箱（Z_N{2..20}_m3_*、Wadj_N{2..26}_m3、R_U_m3_*、Radj_A_m3、largeN_2box_m3.pkl、続き 88）、n2_N_rterms_b{0,1}_m1.pkl・n2_N_rterms_series_b{0,1}_m1.pkl・n2_N_bg_knot_m1.pkl（項ごとの分解と背景・結び目の部分のブロック、続き 92）など。

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
- 凍結した状態（向き付けられない面）の refined の取っ手の因子と、多重項の基底。対になる取っ手は G（[1]、framed unknot、Hopf の両端で確認、続き 63）。自分自身に移る取っ手の多重項は未決：σ 非対称部分は Lefschetz の重みごとに (𝐭²−1)×整数多項式で t の有限の塔を含み、trefoil では同符号、T(2,5) では混在（theory §12）。t の塔は有限個の z₂ モードをもつ Ĩ 不変な状態では説明できず（指数では E と C に吸収、p₁² 方向に寄与しない）、凍結した状態の内部の多重項（続き 64、theory §12.5）。
- E と C の分離。固定セクターの状態の空間への Ĩ の作用を、指数ではなく幾何から求める必要がある。
  - n = 1 は符号の規則で決まった（§2′）。
  - n = 2 と Hopf の (4,4) の E は未決定。
- 普遍因子のずれの向き（k = +1）と、閉弦の Kähler の三重化との関係。開弦の比 w³ と閉弦の比 w^{−1} が合わず、保留している。
- 非トーラス結び目の大きな N の値（入力が無い）。T(2,5) の 3 箱は、m = 2 の composite と D_N の 3 箱が要り、数時間〜半日の見込み。
- 環状の HOMFLY ホモロジー（逆向きの 2 本のケーブル）と refined CS の composite の関係。未着手で、見込みは低い。

### 6. 主なファイル
- **文書**：
  - report.md（K0〜K4' と多箱の詳細）
  - tex/report.tex・report.pdf（計算結果の報告書。§11〜§12 に 10-20 までの結果）
  - theory_orientifold_index.md（指数の理論 §1〜§12。§12 は凍結した状態の多重項の基底と t の塔）
  - literature.md（§4 に非 U(N) の rCS の文献）
  - conventions.md
- **コード**：
  - src/rcsknots/：kauffman_*.py、koornwinder.py、daha_dn.py、daha_fast.py、flint_rf*.py
  - scripts/refined_kauffman/：rk_*.py。rk_c_* は一般の N、rk_b_* は構造と正値性の解析、rk_k* は固定 N と Hopf、rk_d_* は Sp、rk_e_* は Krefl–Walcher、rk_f_* は凍結した状態の基底と E の決定（rk_f_twisted*）
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
- 【計算】冪和の成分：Z_1^{Sp}(𝐚) = −Z_1^{SO}(−𝐚)、Z_{11}^{Sp}(𝐚) = Z_{11}^{SO}(−𝐚)、Z_2^{Sp}(𝐚) = Z_2^{SO}(−𝐚) + 1（符号は (−1)^{|ρ|}、p₂ 方向にだけ定数 1）。unknot の [1] で even(g_1^{SO}) = +1、even(g_1^{Sp}) = −1。被覆の規約（Sp でも SO と同じ ν の関数、A_cov = ν²s²）を theory §11.15 に明記。テストに Z₂・Z₁₁ の規則を追加。
- 【計算】定数 1 の由来（rk_d_sp_so_W_rule.py、301 s）：W_λ^{Sp}(𝐚) = (−1)^{|λ|} W_λ^{SO}(−𝐚)（全 λ）、冪和の係数は co[S2][()] だけ −2 違う → Z₂ の +1。SO と Sp の 2 箱の違いは (i) W の (−1)^{|λ|}・𝐚→−𝐚、(ii) p₂ の定数項のずれ、(iii) Adams 項の非可換性、で尽くされる（theory §11.15、claims D14、論文）。

## 2026-10-21（続き 55）多箱の c = 1 の σ 非対称性（ユーザーの指示「多箱の σ 非対称性を進めよ」）
スクリプト：rk_b_multibox_sigma.py、rk_b_multibox_sigma_split.py。ログ：b_multibox_sigma.log、b_multibox_sigma_split.log。文書：theory §11.16。
- 【計算】t = q で Mariño の z 偶奇（c の定義）と ν 偶奇を比べた：[1]・p₁²・(3)・(1,1,1) では ν 偶 = c 1、p₂ 方向と (2,1) では逆。hatted の 2 箱では ν 奇が c = 1（[1] と逆）。p₂ 方向の ν 偶には unrefined で既に c = 2 が混ざり（閉じた形 −2zδ(q^{4n}−1)/(q−1)…）、これが「多箱の c = 1 の σ 非対称性」の一部。
- 【計算】正しく c = 1（hatted の ν 奇）を取っても refined では σ 対称でない。前の因子・積・転置・Adams の補正では直らない。σ 非対称部分は (q−t)(qt−1) = 𝐪²(𝐭²−1)(V⁴−1) で割り切れる → c = 1 = σ 偶部分 S ＋ (𝐭²−1)(V²−V⁻²)·K′（K′ は σ 対称）。
- 【計算】S = σ 偶部分、O = σ 奇部分 = (q−t)(qt−1)K：S, K は半整数。G 展開のブロックの符号は trefoil では O 側が全ブロックでそろい（β で交互）、S 側は中央 3 ブロックだけ。T(2,5) では両方とも中央のブロックが混在。σ 偶・奇の分解だけでは多箱の正値性は出ない（theory §11.16 (3)(4)）。
- 結論：多箱の σ 非対称性 = (i) ν の偶奇による c の誤同定（unrefined）＋ (ii) refined の σ 奇部分（(𝐭²−1)(V²−V⁻²) の因子、t=q で消える）。(ii) は S_R が非自明に働く状態の G 展開側の表れと解釈【予想】。
- 次：論文への反映、push（許可）。

## 2026-10-21（続き 56）3 箱の hatted の規格化（ユーザーの指示）
- 【訂正】続き 55 で「3 箱の hatted は t=q で有理関数」としたのは、largeN_3box_m1.pkl の g が冪和基底であるのを Schur と取り違えたため（rk_b_3box_structure.py と tests/test_rk_daha.py は正しく変換している。影響は続き 55 の 3 箱の行だけ）。theory §11.16 の表を訂正。
- 【計算】Schur に直すと ĝ_λ（|λ| = 3）は refined でも U 一般で整数 Laurent。t = q では 3 つとも「ν 偶 = c 1、ν 奇 = c 2」。hatted の c の規則は β ≡ |R| − 1 (mod 2)（[1] 偶、2 箱 奇、3 箱 偶）。冪和では p₂ を含む方向だけ逆。Mariño の M_{RS} と KN の M_matrix(n, q^{−1/2}) の一致を n = 2, 3 で確認。theory §11.17、tests/test_rk_multibox_sigma.py に 3 箱の偶奇を追加。

## 2026-10-21（続き 57）案 B の第一歩（ユーザーの指示「案 B（Hilbert scheme と実構造）」）
- 【計算】ORS のセル分解（Γ = ⟨2,3⟩）を実装し、T(2,3) の DGR の unreduced superpolynomial を l ≤ 8 で再現（辞書を数値的に固定）。rk_b_planB_cells.py、tests/test_rk_planB_cells.py。実装の注意：σ の巡回規則 a_{m+1} = a_1 + n、生成元の条件は g − h ∉ 𝔧（h ≤ g）。
- 【計算】trefoil の固定セクターの指数 P = (1−t)x は (ν²−1)(q−1)(ν²−t)(ν²q−1)·N₈ の形（framed unknot と同じ境界因子）。P は被覆の [1] の ψ² と単純な関係に無い。theory §11.18。
- 未実施：固定セクター P のセル模型、実構造（y → −y × 複素共役）の各セルへの符号。

## 2026-10-21（続き 58）論文の仕上げ（1）：日本語稿の整合（ユーザーの指示「論文の仕上げ」）
- 要旨（5 項目に。Sp と多箱の σ 分解を追加）、序論（結果の概要に §sec:sp と σ 分解、表に 3 箱の hatted と Sp の行）、§sec:index の σ 対称性の小節、§sec:generalN（3 箱の hatted の整数性、G 展開の破れの注記、Sp の被覆の規約）、§sec:positivity のまとめ、議論（多重項・Sp・B 型・Hilbert scheme の模型）、付録（daha_cn.py、rk_b_planB_cells.py）を更新。ORS を文献に追加（書誌確認済み）。印 [D12][D14][P11] を追加。はみ出し 0、未定義参照 0、31 ページ。
- 残り：英語版（決定事項「最終版は英語にする」）。日本語稿は正値性が未解決のままの準備稿なので、英語化の時期は判断待ち。

## 2026-10-26（続き 59）論文の変数の統一（ユーザーの指示「変数が複雑過ぎるので統一せよ。基本は KN、Kauffman は NRZ」）
- 論文の変数を KN（q, t, 太字 𝐪, 𝐭, q_ℓ, q_r, a = t^{N_A}）と NRZ（Kauffman 変数 λ、表現は R）に統一した。廃止：U, A, V, W, w, r, s, ν, ν̃（→ t^N, a, q_ℓ^{1/2}, q_r^{1/2}, q_r, q^{1/2}, t^{1/2}, λ, λ̃）。表現の λ, μ は R, S に、簡約した不変量 R_λ は 𝓡_R に、framed unknot の framing x は 𝖿 に、Cherednik の a は a_C に改名。SO 側の 𝐚 は AgSc の 𝐚 = λ のままとし、KN の被覆の 𝐚_KN = 𝐚 (t/q)^{1/4} で結ぶ（判断の理由は conventions.md）。
- 生成付録（rk_report_tex.py）も同じ記号で再生成。31 ページ、未定義参照 0。コード・理論文書の内部変数は変更せず、対応表を conventions.md と論文の変数表に置いた。

## 2026-10-26（続き 60）閉弦の refined orientifold の文献確認（ユーザーの指示「未確認の文献を確認せよ」）
- arXiv API で書誌を確認し、ソースを取得：Sinha–Vafa hep-th/0012136、Walcher 0712.2775（書誌のみ）、Krefl–Walcher 0902.0616・1007.0263（LMP 95 (2011) 67–88）・1010.2635、Kim–Yagi 1709.01928（PRD 97 (2018) 026011）、Hayashi–Zhu 2012.13303（JHEP 04 (2021) 292）。literature.md §6、references/README.md、justfile の papers に追加。
- 要点：Krefl–Walcher は「β = 2 の Ω 背景の理論 = β = 1 の理論の orientifold」（N_f = 0, 2, 4）を示し、位相的 orientifold を β = 2 の精密化に埋め込む予想を述べる。Hayashi–Zhu は O5 面の refined O-vertex を摂動部分で構成。開弦の refined orientifold の BPS 不変量の文献は無い。
- 検討の余地【予想】：本研究の refined SO(2N) のデータを t = q²（β = 2）に特殊化したとき、unrefined の U(N) 側の量と結び付くか（Krefl–Walcher の関係の結び目版）。未着手。

## 2026-10-26（続き 61）Krefl–Walcher の関係の結び目版（ユーザーの指示）
- 【計算】unrefined の Kauffman [1]（trefoil、T(2,5)）と refined U(N′) の [1]（DGR）の β = 1/2, 1, 2 への特殊化は、Wilson loop の段階でも自由エネルギーの形（g_[1] = P̄(β=2) − ½P̄(β=1)）でも、a = λ^k q^j の写像と単項式倍の範囲で一致しない（否定的）。theory §11.19、rk_e_krefl_walcher.py。
- 読み：KW の関係は分配関数（閉弦）の性質で、Wilson loop には持ち上がらない（SO の指標と Macdonald 多項式の違い）【予想】。

## 2026-10-26（続き 62）push、main の取り込み、研究別ディレクトリへの再配置
- ユーザーの許可を得て c8e1b1f（＋437fbce）を origin に push（`git push` をサンドボックス外で実行）。PR #11 が main と競合（.gitignore）していたので、app の同期ツールで origin/main（PR #10「リポジトリを framed_knot/ と dgr_lmov/ に整理」）を本ブランチに merge し、.gitignore は main 版＋refined_kauffman の行で解消（b3129a5）。
- main の構成規則に合わせて再配置（dd5c428）：scripts/rk_*.py → scripts/refined_kauffman/、tests/test_rk_*.py → tests/refined_kauffman/、src/rcsknots/{askey_wilson, daha_cn, daha_dn, daha_fast, flint_rf, flint_rf3, flint_solve, kauffman_*, koornwinder} → src/rcsknots/refined_kauffman/、logs/ 直下の rk_*.log → logs/refined_kauffman/。import（rcsknots.X → rcsknots.refined_kauffman.X、相対 import は ..）、sys.path（scripts/refined_kauffman、framed_knot のスクリプトを使う 3 本は scripts/framed_knot も）、データのパス（data/torus 等 → data/framed_knot/…）、テストの ROOT（parents[2]）を更新。CLAUDE.md の項目 5、README の構成、付録を更新。
- 検証：既定のテスト 70 件通過（logs/refined_kauffman/paper_tests_default_2026-10-26.log）、論文 31 ページ。merge（b3129a5、Auto-fix の権限）と再配置（dd5c428、468c168、ユーザーの許可）を push 済み。PR #11 は MERGEABLE。

## 2026-10-26（続き 63）凍結した状態の refined の取っ手の因子と多重項の基底（ユーザーの指示「凍結した状態の refined の取っ手の因子と、多重項の基底を進めよ」）
- 作業場所：worktree `rename-literature-survey-doc-9cc5fc`（ブランチ `claude/kauffman-research-progress-c23e34`）に origin/main（PR #11 マージ後）を取り込んで作業した。
- 【導出】凍結した状態（電荷 β）の被覆は Ĩ 不変な曲線（電荷 2β）で、N のゼロモードは被覆の Ĩ 不変なモード。入れ替わる取っ手の対 → G、自分自身に移る取っ手 → σ 非対称の源、境界モードは h − 1 個で hatted が除く。hatted の c の規則 β ≡ |R| − 1 (mod 2) は β + |R| + c ≡ 0 (mod 2) の言い換え（theory §12.1）。
- 【計算】framed unknot（𝖿 = 1..4）の 2 箱の p₁² 方向の凍結した状態（E を含まない）は、結び目の [1] と同じ判定式を全部満たす（c=1 は σ 対称で G 展開の符号がそろうか (𝐭²−1)×同符号、c=2 は −(𝐚²−𝐪²)(𝐭²−1)×正）。多箱で判定式を満たす最初の例。Hopf では C_11 ↔ P_22、C_22 ↔ P_44 の両端の β で P_{2β} = 単項式 × C_β²（平方根の構造）。theory §12.2、rk_f_frozen_basis.py、rk_f_hopf_sqrt.py。
- 【計算】結び目の多箱のブロックを Lefschetz の重み m で分け、σ 非対称部分 d_m = c_m − c_{−m} を調べた。d_m は全対象で (1 − w²) で割り切れ（t = q で消えることの帰結）、商の符号は trefoil の 2 箱で全ブロック同符号（β で交互）、σ 対称部分 A⁻ も同符号で、全ブロックが基底 {G^g w^j} ∪ {(𝐭²−1) q_ℓ^m w^j} で同符号に展開できる。trefoil の 3 箱は非対称部分だけ同符号（A⁻ は中間で混在）、T(2,5) の 2 箱は中間の 4 ブロックで非対称部分も混在。非対称部分は q 電荷を固定した t の有限の塔を含む（例：trefoil 2 箱 ν^{−11}：−(1−w²)w^{−9}(2t + … + 2t⁵ + t⁶)）。theory §12.3。
- 試して不採用：基底 {G^g w^j, z_q² G^g w^j, z_t² G^g w^j}（同種のモードの対）での同符号展開の整数計画（rk_f_frozen_basis.py --ilp）。結び目の σ 非対称なブロックではほぼ解なし（trefoil 2 箱 ν^{−5}、3 箱の両端だけ解あり）。
- 記録：theory §12、claims P13、tests/refined_kauffman/test_rk_frozen_basis.py（3 件、42 秒）、logs/refined_kauffman/f_frozen_basis.log、f_hopf_sqrt.log。references/1202.4456/ に AgSc のソースを取得（未読。Ĩ の時空への作用の確認用）。
- 次の候補：(1) σ 非対称部分の t の塔を「有限個の z₂ モードをもつ Ĩ 不変な状態」として全 d の構造に組み込めるか（E と C の中間の型）。(2) framed unknot の 4 箱の composite 固定セクター（(2;2) 型）を作り、2 箱の凍結した状態との β ↔ 2β の対応と平方根の構造を確かめる。(3) T(2,5) の中間ブロックの混在を、同じ β に寄与する曲線の分離（(J₂, S_R) の分解）で説明できるか。

## 2026-10-26（続き 64）σ 非対称部分の t の塔と、有限個の z₂ モードをもつ Ĩ 不変な状態（ユーザーの指示「σ 非対称部分の t の有限の塔を…全 d の構造に組み込めるか（E と C の中間の型）を進めて」「解析を続けよ」）
- 【導出】長さ L の z₂ の塔をもつ状態の g への寄与：対は 0、Ĩ 不変で J₂ の塔なら ½εm[L]_{−t}（= E の 2 項）、S_R の塔なら ½εm[L]_t（対で凍結した状態と区別できない）、凍結は m[L]_t。中間の型は指数では E と C に吸収され、Ĩ 不変な状態は p₂ 方向にしか現れないので p₁² 方向には寄与しない。p₁² 方向の t の塔は凍結した状態の内部の多重項でなければならない（theory §12.5 (a)）。
- 【計算】p₁² 方向の σ 非対称部分 D = (1 − q/t) × e で、e の q の冪ごとの列は trefoil で全部、T(2,5) で 31 列中 26 列が t の指標 t^k[L]_t の同符号の和。最低のブロックは D = −(1 − q/t) w^{−(8m+1)} t([2n]_t + [2n−1]_t)（m = 1, 2）。T(2,5) の混在 5 列はどれも t の次数の順で符号の変化が 1 回（二つの多重項が t の別の区間に逆符号で入る）。
- 【予想】凍結した状態の t の多重項は、実の部分の Z/2 係数の Poincaré 多項式の形（RP^{L−1} 型）。検証には被覆の composite (2;2)（4 箱、大きな N）が要る。
- 記録：theory §12.5、claims P14、scripts/refined_kauffman/rk_f_frozen_towers.py、logs/refined_kauffman/f_frozen_towers.log、tests/refined_kauffman/test_rk_frozen_basis.py に 2 件追加（計 5 件、15 秒）。
- 次の候補：(1) trefoil の composite (2;2) の大きな N（GL の DAHA、4 箱）の見積もり（実構造の予想の検証）。(2) T(2,5) の混在列の二つの多重項を、ν の隣のブロックや A⁻ の G 展開の項と組にして説明できるか。

## 2026-10-26（続き 65）trefoil の composite (2;2) 型（被覆、4 箱、大きな N）の見積もり（ユーザーの指示「計算時間とメモリを小さい N で測って見積もって」）
目的：§12.5 の予想（凍結した状態の t の多重項は実の部分の Z/2 Poincaré 多項式の形）の検証。被覆の Ĩ 不変な曲線（電荷 2β）の側で同じ多重項が t² の刻みで現れるかを見るには、被覆の固定セクター（brane 2・image 2）を作る composite (α;β)（|α| = |β| = 2：(2;2), (2;11), (11;2), (11;11)。(11;2) は荷電共役で (2;11) と同じ）の大きな N の値が要る。
- 実装：scripts/refined_kauffman/rk_c_fast_composite22.py。rk_c_fast_composite.comp_raw（(α;1)）を (α;β) に広げた（X^{−1} 側で P_β(X^{−1}) を p_d(X^{−1}) の積で展開して γ̂' を当て、X 側の連鎖の出発点にする）。|α| = |β| なので評価の換算因子は 1。
- 検算【計算】：N = 3, 4 で (1;1) が gl_adj_fast + κ に、(2;1)・(11;1) が comp_raw に全係数一致。
- 実測（1 プロセス、(2;2)・(2;11)・(11;11) の生の値の計、途中結果を共有。logs/refined_kauffman/c_composite22_timing.log）：
  | N | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
  |---|---|---|---|---|---|---|---|
  | 時間 | 2.2 s | 3.6 s | 7.8 s | 21 s | 72 s | 141 s | 302 s |
  | 最大 RSS | 0.08 GB | 0.09 GB | 0.07 GB | 0.10 GB | 0.18 GB | 0.28 GB | 0.40 GB |
  時間の約 80% は (2;2)。N ごとの比は時間 約 2.0〜2.1（N = 8→9→10）、メモリ 約 1.4。
- 重い部分（cProfile、N = 7）：最後の評価の B（Y_{ε_1} を含む語）の作用での有理関数の和。flint_rf.RF.__init__ の gcd による約分が全体の約 70%。
- 必要な N【見込み】：補間の窓は、trefoil の随伴 (1;1)（2 箱）が A^{1/2} の冪 −10〜−4 の 7 項、(α;1)（3 箱）が −14〜−6 の 9 項なので、4 箱は −18〜−8 の 11 項と見込む（外挿）。(α;1) では大きな N の補間が N = 2 からの値と一致したので、多項式の形は N ≥ ℓ(α) + ℓ(β) = 4 から成り立つと見て、N = 4〜14 で補間、N = 15, 16 で照合とする。窓は N = 4〜10 の簡約した値の s の次数の傾きで、本計算の前に確かめる。
- 外挿（実行していない）：
  | N | 比 1.8 のとき | 比 2.1 のとき | メモリ |
  |---|---|---|---|
  | 12 | 16 分 | 22 分 | 0.8 GB |
  | 14 | 53 分 | 1.6 時間 | 1.5 GB |
  | 16 | 2.9 時間 | 7.2 時間 | 3.0 GB |
  | 18 | 9.2 時間 | 32 時間 | 5.9 GB |
  N = 4〜16 の合計は 6〜14 時間（逐次）。機械は 8 コア・8 GB なので、N ごとに独立なプロセスで並列にすると、N = 16 の 1 本が律速で壁時計 3〜8 時間。同時に走らせるとメモリの和が 8 GB を超えるので、N = 11〜14 と N = 15, 16 の 2 回に分ける。N = 18 まで要る場合はメモリが 6 GB に近く、1 本ずつしか走らない。
- 高速化の候補（未実装）：(1) 最後の一段を評価の汎関数にして多項式を作らない、(2) 係数を gcd 付きの RF から共通分母の Laurent 多項式（fmpz_mpoly）に替える（daha_fast.md の見込みで 2 倍程度）、(3) 和の途中で約分しない。(2)(3) は RF の約分が 70% を占めることから効果が大きいと見込む。
- 本計算に要るもの（未実装）：既約の composite への分解 P_α(X)P_β(X^{−1}) = Σ [α′;β′]（(1;1)、(0;0) を含む）、主特殊化からの K と分解係数の閉じた形、一様な規格化の単項式。どれも (α;1) と同じ方法で、1 件あたりの時間は生の値より小さい見込み。
- 再開：N と (α;β) ごとに data/refined_kauffman/daha/Rcomp22_N{N}_m1_{α}_{β}.pkl に保存し、保存済みは飛ばす。

## 2026-10-26（続き 66）composite (2;2) 型の高速化と本計算の計画（ユーザーの指示「RF の約分を減らす高速化…を先に実装し、N = 8〜10 で速さと一致を確かめてから本計算の計画を立てて」）
- 実装：src/rcsknots/refined_kauffman/flint_lp.py（約分なしの Laurent 多項式）、daha_fast_lp.py（GL_N の高速版の LP 版、基底の値 V[μ; ν] と最後の一段の評価の汎関数）、scripts/refined_kauffman/rk_c_composite22_lp.py、tests/refined_kauffman/test_rk_daha_fast_lp.py（2 件、1 秒）。詳細は daha_fast.md の最後の節。
- 【計算】旧版（RF）と全係数一致：N = 3〜5（poly・functional の両方）、N = 8, 9, 10（functional）。
- 速度：N = 10 で 288 s → 16 s（約 18 倍）、N = 12 で 66 s・1.0 GB。N ごとに時間 約 2.0 倍、メモリ 約 1.45〜1.7 倍。
- 試したこと：評価の汎関数の点の数は最後に N(N+1)/2 個で小さい。最初の版は点ごとの値を LP の逐次の和で求めて遅かった（N = 8 の 73%）ので、ずらしを numpy の行列積で求め、fmpz_mpoly を二分木で足す形にした（N = 10 で 35 s → 16 s）。

### 本計算の計画（未実行）
目的：被覆の固定セクター（brane 2・image 2）の大きな N の値を作り、theory §12.5 の予想（凍結した状態の t の多重項が被覆の側で t² の刻みの多重項になる）を検証する。
1. **既約の composite への分解と規格化**（実装、計算は数分）：
   - P_α(X) P_β(X^{−1})（|α| = |β| = 2）= P_[α;β] + c₁ P_[1;1] + c₀（[1;1] は随伴、c は q, t, t^N の有理式）。(2;11) は [1;1] を含まない場合もあるので、係数は一般の形で解く。
   - 係数は二つの点での評価で決める：点 t^ρ（ev）と、Macdonald の対称性 P_b(t^ρ q^{ε_1})/P_b(t^ρ) = P_(1)(t^ρ q^b)/P_(1)(t^ρ) を使う点 t^ρ q^{ε_1}。どちらも閉じた形なので、一般の N で c₁, c₀ が求まる。N = 4〜6 で多項式の分解（solve_symmetric）と照合する（(α;1) の d_α と同じ検証）。
   - 一様な規格化の単項式 M_b（b = α − β の並び、(α;1) と同じ規則）と、主特殊化の K_b を掛けて簡約した R_[α;β](N) を作る。U(4), U(5) の KN の data（移した分割 [4,2,2,0]・[3,3,1,1] 型）があれば、(α;1) と同じく比 1 を確かめる（無ければ固定 N の DAHA の値との一致で代える）。
2. **窓の確認**：N = 4〜12 の簡約した値（計 2 分程度）の s の次数の傾きから、t^N（A^{1/2}）の冪の窓を決める（rk_c_read_largeN.py の方法）。見込みは −18〜−8 の 11 項（外挿）。
3. **本計算**：窓が 11 項なら N = 4〜14 で補間、N = 15, 16 で照合（窓が広ければ N の上限を窓 + 5 にする）。
   - 見込み（N = 8〜12 の実測から外挿）：N = 13 で 2.2 分・1.5 GB、14 で 4.5 分・2.2 GB、15 で 9 分・3.3 GB、16 で 18 分・5 GB 前後。合計約 35 分。
   - メモリが 8 GB の機械なので、N ≥ 14 は 1 本ずつ（逐次）実行する。N ≤ 13 は並列でよい。
   - N = 17 以上が要る場合はメモリが 8 GB を超える見込み（N = 17 で 7〜11 GB）。そのときは、基底の値を ν ごとに別のプロセスで計算し、途中の多項式（B のべき）を捨てるように直してから測り直す。
   - 再開：N ごとに V[μ; ν] を data/refined_kauffman/daha/V22_N{N}_m1.pkl に保存し、保存済みは飛ばす。実行は nohup で、ログは logs/refined_kauffman/c_composite22_run.log、完了の目印は最後の行「done」。
4. **組み立てと検証**：被覆の 4 probe（brane 2・image 2）の f^cov を大きな N で作り（rk_c_3box_largeN の組み立てを 4 箱の被覆の部分だけに広げる）、固定セクターの指数 P = (1 − t)·f^cov_{固定}を求める。§12.5 の予想（trefoil の 2 箱の凍結した状態 C の ν^β のブロック ↔ P の ν^{2β} のブロックで、t の多重項 [L]_t ↔ t² の刻み）と、Hopf で見た平方根の構造（両端の β で P_{2β} = 単項式 × C_β²）を調べる。

## 2026-10-26（続き 67）composite [2;2], [2;11], [11;11] の一般の N（計画の 1〜3）（ユーザーの指示「計画の 1 と 2 を実装し、窓が見込みどおりなら 3 の本計算まで進めて」）
- 実装：scripts/refined_kauffman/rk_c_composite22_run.py（check・run・window）。tests/refined_kauffman/test_rk_daha_fast_lp.py に 1 件追加（計 3 件、2 秒）。
- **1. 既約への分解と規格化**【導出＋計算】：
  - P_α(X)P_β(X^{−1}) = P_b + c₁ P_adj + c₀（b = (α, 0, …, 0, −β の逆順)）。係数は二つの評価から解いた：ev（点 x = (t^{−(N−1)}, …, 1)）と、Macdonald の作用素 D の固有値（e_b = Σ q^{b_i} t^{N−i}）。点 x では D の局所化の和が i = N の項だけになり、ev(D f) = [N]_t f(x′)、x′ = (t^{−(N−1)}, …, t^{−1}, q)。
  - N = 3, 4 で多項式の分解（solve_symmetric）と、c₁, c₀, K_b, ev(D P_b) = e_b K_b がすべて全係数一致。[2;11] は c₀ = 0（[2;2], [11;11] は c₀ ≠ 0）。
  - 一様な規格化（M_b は (α;1) と同じ rk_c_fast_composite.mono）で、KN の data と比 1：U(4) の [3,1]（= [2;11]）・[2,2]（= [11;11]）、U(6) の [3,1,1,1]・[2,2,1,1]。[2;2] に対応する data（U(N) の [4,2^{N−2}]）は無いので、[2;2] の規格化は同じ規則によるもので、data との照合はしていない。
- **2. 窓**【計算】：N = 4〜12 の s の次数の範囲の傾きは x = s^N で [2;2]・[11;11] が x^{−20}〜x^{−8}（13 項）、[2;11] が x^{−18}〜x^{−8}（11 項）で、見込み（11 項）より広かった。しかし x の**偶数冪だけが現れ**、U = t^N で書くと [2;2]・[11;11] が U^{−10}〜U^{−4}（7 項）、[2;11] が U^{−9}〜U^{−4}（6 項）になる。
- **3. 本計算**：U = t^N の窓では N = 14 までで足りるので、N = 13, 14 を追加で計算した（170 s・1.1 GB、312 s・1.5 GB）。N = 15, 16 は不要になった。
  - 補間：ノード N = 6〜12（[2;11] は 7〜12）、補間に使わない N = 4, 5（[2;11] は 4〜6）と、補間の外側の N = 13, 14 で全係数一致。
  - 保存：data/refined_kauffman/daha/V22_N{N}_m1.pkl（基底の値）、Rcomp22_N{N}_m1_{α}_{β}.pkl（簡約した値、N = 4〜14）、Rcomp22_A_m1_{α}_{β}.pkl（U = t^N の係数）。ログ：logs/refined_kauffman/c_composite22_run.log、c_composite22_window.log。
- 次：計画の 4（被覆の brane 2・image 2 の固定セクターの組み立てと、§12.5 の予想・平方根の構造の検証）。

## 2026-10-26（続き 68）被覆の brane 2・image 2 の固定セクター（計画の 4）（ユーザーの指示「計画の 4 に進んで」）
- 実装：scripts/refined_kauffman/rk_f_cover_fixed22.py（分解の係数と次元を A = t^{N_A} の閉じた形に、⟨P_λ P_μ^∨⟩（|λ|, |μ| ≤ 2）、二変数の冪和の対数、固定セクター P、analyze）。tests/refined_kauffman/test_rk_daha_fast_lp.py に 1 件追加（計 4 件、2 秒）。データ：data/refined_kauffman/daha/largeN_cover_fixed22_m1.pkl。ログ：logs/refined_kauffman/f_cover_fixed22.log。
- 検算【計算】：A の閉じた形の c₁, c₀, K_b が N ごとの式と N = 5, 7 で一致。(1),(1) の成分が既存の brane–image の x に一致。
- 結果【計算】（theory §12.6、claims P15）：
  - 4P/z_q³ の ν のブロックは 12 個すべて σ 対称で G 展開できる。符号は両端でそろい中間で混在。
  - 凍結した状態の σ 非対称な t の塔に当たるものは被覆の側に無い。§12.5 の予想（被覆に t² の刻みの多重項）は素朴な形では支持されない（否定的）。
  - 平方根の構造（P_{2β} = 単項式 × C_β²）は [1] の最低のブロックだけで成り立ち、2 箱ではどのブロックでも成り立たない。
  - G の次数の 2 倍の関係（凍結した状態の ℓ 重みの最大 × 2 = 被覆の G の次数）は 6 ブロック中 5 で成立（β = −9 で 1 ずれる）。
  - 最低のブロック = −s^{20} r^{−41} Σ_{i,j=0}^{11} q^{max(i,j)} t^{min(i,j)}（t = q で完全平方 [12]_q²）。
- 注意（別件）：scripts/refined_kauffman/rk_c_3box_largeN.py の cover_parts は `from rcsknots import kauffman_refined` のままで、再配置（続き 62）の後は import に失敗するはず（モジュールは rcsknots.refined_kauffman.kauffman_refined）。今回は使っていないので直していない。
- 次の候補：(1) 凍結した状態の t の塔を、被覆の状態への Ĩ の作用の分け方（固定セクターの状態の Ĩ 固有値）として、Hopf (4,4) と同じ整数計画で調べる（4P の両端のブロックは状態の数が少ない）。(2) 最低のブロックの Σ q^{max}t^{min} の幾何的な意味（Hilb² / Sym² 型の模型）。

## 2026-10-26（続き 69）凍結した状態の t の塔を、被覆の固定セクターへの Ĩ の作用の分け方として調べる（ユーザーの指示）
- 実装：scripts/refined_kauffman/rk_f_twisted22.py（§9 の整数計画を brane n・image n の固定セクターに）、rk_f_twisted22_maps.py（refined の ψ² 型の候補の写像とずれの走査）。tests/refined_kauffman/test_rk_frozen_basis.py に 1 件追加（22 秒）。ログ：logs/refined_kauffman/f_twisted22.log。
- 【計算】n = 1（(1),(1)）：(i)+(ii) で E が一意（§9 の 34 変数の結果を再現）。n = 2（(2,2)）：t = −1 の条件 (i) には SO の 4 箱が要り使えないので (ii) だけで解いた。解はあるがほとんど決まらない（ν^{−24} で 78 中 74 が自由）。
- 【計算】refined の ψ² 型の分け方（凍結した状態の t の塔が t² の刻みの塔に移る形）は、t = q で ψ² に戻る 9 通りの写像と (s/r)^j のずれのどれでも、ほとんどのブロックで偶奇・|e| ≤ |c|・台の条件を破る（n = 1 で 7 中 2〜3、n = 2 で 12 中 0〜1）。n = 1 の一意な E も ψ² 型ではない。結論：凍結した状態の t の塔は、この形の Ĩ の分け方からは得られない（否定的）。theory §12.7、claims P16。
- 次の候補：n = 2 で E を決めるための SO の 4 箱（D_N の DAHA、大きな N、p_(2,2) 方向）の見積もり。daha_fast.md の 4 箱の見積もり（RF 版で N = 8 に 15 分・2 GB）は、D_N の LP 版（約分なし）を作れば大きく縮む見込み。補間に要る N の数は窓の確認から決める。

## 2026-10-26（続き 70）SO の 4 箱（D_N の DAHA、trefoil、大きな N）の見積もり（ユーザーの指示「n = 2 で E を決めるために、SO の 4 箱…を、D_N の約分なしの版も含めて小さい N で測って見積もって」）
目的：brane 2・image 2 の E（theory §12.7）の条件 (i)（t = −1）に要る、SO の 4 箱（冪和の p_(2,2) 方向）の大きな N の値。
- 実装：src/rcsknots/refined_kauffman/daha_fast_lp.py に D_N 版（DNLP、dn_jones）を追加。評価の汎関数の点の扱いを一般化した（点の鏡映は D.reflect、π_1^{−1} の q の冪はベクトル κ。GL の値は変わらず、既存のテストが通る）。scripts/refined_kauffman/rk_c_dn4_lp.py、tests/refined_kauffman/test_rk_daha_fast_lp.py に 1 件追加（計 5 件、6 秒）。
- 照合【計算】：新版の (Z_λ, K_λ) が旧版（daha_fast.daha_jones_fast、RF）と |λ| ≤ 4 の全 λ で全係数一致（N = 2〜6）。
- 実測（logs/refined_kauffman/c_dn4_timing.log。N ≤ 6 は |λ| ≤ 4 の全 λ、N ≥ 7 は |λ| = 4 の 5 個）：
  | N | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
  |---|---|---|---|---|---|---|---|
  | 旧版（RF） | 7.4 s | 24 s | 79 s | — | — | — | — |
  | 新版（LP） | 3.3 s | 12 s | 25 s | 75 s | 183 s | 510 s | 1925 s |
  | 新版の最大 RSS | 0.11 GB | 0.18 GB | 0.45 GB | 0.84 GB | 1.19 GB | 2.09 GB | 2.14 GB |
  - D_N での速さは旧版の 2〜3 倍にとどまる（GL の 18 倍と違う）。N = 6 の内訳（cProfile）：基底の値の B の作用 40%、評価の汎関数 46%（点の数は最後に 73 個 = GL の約 3 倍。D_N の Weyl 群の軌道が大きい）、K_λ の展開（RF）13%。
  - N = 10 の測定時はスワップが 3.2 GB 使われていた（機械は 8 GB、他のアプリを含む）。RSS が N = 9 とほぼ同じなのにメモリ圧縮の影響と思われる時間の伸び（3.8 倍）があるので、N = 10 の値は過大の可能性が高い。外挿には N = 7〜9 の伸び（時間 2.5〜2.8 倍、メモリ 1.4〜1.75 倍）を使う。
- 必要な N【見込み】：SO の補間の窓は 2 箱が U^{−10}〜U^{−4}（7 項）、3 箱が U^{−15}〜U^{−6}（10 項）なので、4 箱は U^{−20}〜U^{−8}（13 項）と見込む（外挿）。ℓ(λ) ≤ N から N ≥ 4。補間に N = 4〜16、照合に N = 17 程度が要る。
- 外挿（実行していない）：
  | N | 12 | 13 | 14 | 15 | 16 | 17 |
  |---|---|---|---|---|---|---|
  | 時間（|λ| = 4 の 5 個） | 2〜3 時間 | 6〜9 時間 | 14〜24 時間 | 1.5〜3 日 | 4〜8 日 | 9〜22 日 |
  | メモリ | 6〜11 GB | 8〜20 GB | 11〜34 GB | 16〜60 GB | 22〜105 GB | 31〜184 GB |
- **結論**：この方法では N = 17 までは届かない。メモリが N = 12〜13 で 8 GB を超え、時間も N ≥ 15 で日の単位になる。高速化の余地は GL のとき（約分が 70%）ほど大きくなく、主な重さは多項式の項数そのもの（D_N の Laurent 多項式は ±ε_i の両方向に伸びる）である。
- 代わりの方法の候補（未検討）：
  1. 条件 (i) に要るのは t = −1 の留数だけなので、q を有限体の数値にして t の一変数（nmod_poly）で計算し、U = t^N で補間した後に t = −1 を代入、q は多点で補間する（多項式の変数が減り、係数が小さい）。
  2. 窓を狭める：大きな N の 1 点から係数を読む方法（c_j の s の次数の幅が N より狭ければ分離する。rk_c_read_largeN.py の考え方）。N = 10 前後で分離するなら N の上限が下がる。
  3. 4 箱の SO 全体ではなく、p_(2,2) 方向に要る λ の一次結合だけを計算する（K_λ の H_ν 展開の段階で組み合わせを先に作る）。基底の値 V[ν] はどの λ にも共通なので、計算量はほとんど変わらない見込み。

## 2026-10-26（続き 71）大きな N の補間に要る N を減らす（SO の 3 箱で検証、4 箱に適用）（ユーザーの指示「窓を狭めるに進め」。以後、補間の文脈で「窓」の語は使わない）
目的：SO の 4 箱の R_λ(U) = Σ_j c_j(q, t) U^j（U = t^N）を、続き 70 の見込み（N = 4〜17）より少ない N から決める。
- 4 箱の保存：N = 2〜9 の |λ| = 4 の 5 個の Z, K を daha_fast_lp.dn_jones で計算し保存（scripts/refined_kauffman/rk_c_dn4_save.py、logs/refined_kauffman/c_dn4_save.log。N = 8 で 270 s、N = 9 は他の計算と同時に走らせて 1587 s、最大 RSS 1.41 GB）。
- U の冪の範囲【計算、N ≤ 9】：5 個の λ すべてで U^{−20}〜U^{−8}（13 項）。s の次数の最小・最大の傾きは N = 2〜4 から一定（続き 70 の外挿どおり）。
- 方法（scripts/refined_kauffman/rk_c_bounded_interp.py）：r の次数 α ごとに、未知数 x_{j,β}（c_j の r^α s^β の係数）の一次方程式 Σ_j x_{j,d−2Nj} = [R(N)]_{r^α s^d} を、使う N すべてについて flint の有理数行列の rref で解く。R(N) の約分は flint（sympy の約分より 3〜6 倍速く、全係数一致を確認）。
  - 仮定 A【仮定】：帯を取る N_s で、各項 c_j U^j が R(N_s) の s の次数の範囲に収まる（β + 2N_s j ∈ [m(N_s), M(N_s)]）。
  - 仮定 B【仮定】（mode "alpha"）：r^α の部分の U の冪の範囲を、使う N のうち大きい二つの間の r^α の部分の s の次数の傾きから読む。
  - 採用の条件：解が一意で、使わない N すべてで全係数一致。
  - 試して不十分だったもの：使う N すべての帯の共通部分（帯の傾きが U の冪の範囲と同じなので、N_s 一つの帯と同じ情報しか無い）。r の次数ごとの帯を使う N すべてで課す版（小さい N で端の次数が打ち消し、(1,1,1)・(2,1) で矛盾）。
- 3 箱での検証【計算、データは N = 2〜13】：
  | λ | 従来（U の冪 10 項の Vandermonde） | 仮定 A のみ（N_s = 使う最小の N） | 仮定 A＋B |
  |---|---|---|---|
  | (3) | N = 4〜13 | N = 4〜10 で一意でない | N = 2〜6 で一意、N = 7〜13 で一致 |
  | (2,1) | N = 4〜13 | N = 4〜11 で一意、N = 2, 3, 12, 13 で一致 | N = 2〜8 で一意（N_s = 4, 5）、N = 9〜13 で一致 |
  | (1,1,1) | N = 4〜13 | N = 4〜10 で一意、N = 3, 11〜13 で一致 | N = 3〜9 で一意（N_s = 4）、N = 10〜13 で一致 |
  - 一意になる条件は、ほぼ「使う N の個数 ≥ α ごとの U の冪の項数の最大」。(1,1,1) は r^{−6} が 7 項（U^{−12}〜U^{−6}）で、N = 4〜9 の 6 個では階数が 6 足りない。
  - N_s を ℓ(λ) まで下げると矛盾する（N = ℓ(λ) では端の次数が打ち消して s の次数の範囲が狭い）。N_s = 4 で 3 箱の全部が通る。
  - (1,1,1) の N = 2 は ℓ(λ) > N で Z = 0 なので使わない（R(U) の式は N ≥ ℓ(λ) で成り立つ。(2,1) の N = 2、(3) の N = 2, 3 では一致）。
- 4 箱への適用【計算、データは N ≤ 9】：
  | λ | 使う N | 結果 |
  |---|---|---|
  | (4) | 2〜8（N_s = 2） | 一意、N = 9 で全係数一致 |
  | (2,2) | 2〜8（N_s = 4, 5, 6） | 一意、N = 9 で全係数一致 |
  | (3,1) | 2〜9（N_s = 4） | 一意だが照合に使う N が残らない。2〜8 では α ごとの U の冪が 8 項の 6 個の α で階数が 24 足りない |
  | (2,1,1) | 3〜9（N_s = 5, 6） | 一意でない。α ごとの U の冪が最大 9 項（r^{−22}, r^{−18}, r^{−12}） |
  | (1,1,1,1) | 4〜9（N_s = 5, 6） | 一意でない。α ごとの U の冪が最大 8 項（r^{−10}） |
  - α ごとの傾きは N = 9 でもまだ動いている（(1,1,1,1) の r^{−16}, r^{−12}, r^{−8}, r^{−6}、(2,1,1) の r^{−24}, r^{−22}）。端の次数の打ち消しが無ければ上端の傾きは N について単調に増え、下端は減るが、増減が入れ替わっているので打ち消しがある。読み取った α ごとの範囲は真の範囲より狭い可能性があり、判定は使わない N での照合に頼る。
- SO の中の転置の対称【計算】：3 箱の (3) と (1,1,1) の c_j は、j ごとの項数が違う（例：U^{−14} で 15 項と 12 項）。λ ↔ λᵀ は SO と Sp を入れ替えるので、SO の中では対応しない。4 箱の (1,1,1,1)・(2,1,1) を (4)・(3,1) から得ることはできない。
- 必要な N の見込み【見込み】：(3,1) の照合に N = 10、(2,1,1) と (1,1,1,1) に N = 10, 11（項数 9 なら (2,1,1) は N = 3〜11 で一意、照合に N = 12）。続き 70 の N = 17 より 5〜6 少ない。
- N = 10, 11 の計算の準備（実行前に記録）：
  - 見積もり：N = 10 は約 0.5 時間（続き 70 の 1925 s はスワップの影響で過大の可能性）、RSS 約 2.1 GB。N = 11 は N = 7〜9 の伸び（時間 2.5〜2.8 倍、メモリ 1.4〜1.75 倍）から 1.3〜1.5 時間、RSS 3〜3.7 GB。N = 12 は 3.5〜4 時間、RSS 4〜6.5 GB で、機械（8 GB）では他のアプリを閉じないと危ない。N = 12 は N = 11 の結果を見てから決める。
  - 再開：rk_c_dn4_save.py は N ごとに保存し、保存済みの (N, λ) を飛ばす。
  - 実行：nohup で N = 10, 11 を順に一つずつ（同時に走らせるとメモリが足りない）。ログは logs/refined_kauffman/c_dn4_save.log に追記、完了の目印は最後の行の "done"。
- N = 10, 11 の結果（logs/refined_kauffman/c_dn4_save.log）：N = 10 は 2384 s・最大 RSS 1.91 GB、N = 11 は 5689 s・2.49 GB（N = 10→11 で時間 2.4 倍、メモリ 1.3 倍）。
- N ≤ 11 での 4 箱【計算】（仮定 A＋B）：
  | λ | 使う N（N_s） | 結果 |
  |---|---|---|
  | (4) | 2〜8（2） | 一意、N = 9, 10, 11 で全係数一致 |
  | (2,2) | 2〜8（4） | 一意、N = 9, 10, 11 で全係数一致 |
  | (3,1) | 2〜9（4） | 一意、N = 10, 11 で全係数一致（2〜10 でも同じ解、N = 11 で一致） |
  | (2,1,1) | 3〜11（5） | 一意だが照合に使う N が残らない。3〜10 では r^{−20} で解が無い（N = 9→10 の傾きで読んだ範囲 U^{−18}〜U^{−13} が狭すぎる。N = 10→11 では U^{−18}〜U^{−10}） |
  | (1,1,1,1) | 4〜11（5） | r^{−8} だけ一意でない（U^{−16}〜U^{−8} の 9 項、階数が 11 足りない）。他の α は一意 |
  - α ごとの傾きは N = 11 でもまだ動く：(2,1,1) の r^{−20} の上端（N = 9, 10, 11 で −12, −13, −10）と r^{−14}、(1,1,1,1) の r^{−8} の上端（−10 → −8）。
- N = 12, 13 の計算の準備（実行前に記録）：
  - 目的：N = 12 で (2,1,1) を照合し、(1,1,1,1) の r^{−8} を決める。N = 13 で (1,1,1,1) を照合する。
  - 見積もり（N = 10→11 の伸びから）：N = 12 は 3.5〜4 時間・RSS 3.2〜3.5 GB、N = 13 は 9〜10 時間・RSS 4.2〜4.6 GB。合わせて約 14 時間。
  - 実行：nohup で N = 12, 13 を順に一つずつ。ログは c_dn4_save.log に追記、完了の目印は最後の行の "done"。N ごとに保存するので、止めても N = 12 の結果は残る。
- t = q との照合：(4) と (1,1,1,1) は NRZ 1310.2240 の [r] 彩色の Kauffman 多項式（rk_c_check_tq.py の方法）と比べられる。(1,1,1,1) の r^{−8} の部分は、r の次数が一つなので t = q（r = s）でも (j, β) が混ざらず、他の α がわかれば t = q の値から一意に決まる（その場合 t = q は照合でなく入力になる）。
- 3 箱での最終確認【計算】：scripts/refined_kauffman/rk_c_4box_largeN.py 3box。仮定 A＋B の解（(3) は N = 2〜6、(2,1) は 2〜8、(1,1,1) は 3〜9）が、従来の方法（N = 4〜13 の Vandermonde）の R_U_m1_{λ}.pkl と全係数一致。
- 4 箱の R_λ(U) の保存（rk_c_4box_largeN.py、logs/refined_kauffman/c_4box_largeN.log）：
  - (4), (2,2), (3,1)【計算、照合は N ≤ 11】：上の表の N で解き、使わない N すべてで全係数一致したものを data/refined_kauffman/daha/R_U_m1_{4,22,31}.pkl に保存。U^{−20}〜U^{−8}、整数係数。
  - (4) の t = q【計算】：NRZ 1310.2240 の [4] 彩色の trefoil の Kauffman 多項式（rk_c_check_tq の対応）と比が 1。補間に使っていない独立な照合。
  - (1,1,1,1)【計算、r^{−8} の部分は t = q を入力にした】：r^{−8} 以外は N = 4〜11 の一次方程式で一意。r^{−8} の部分を [1^4] の t = q の値（F_{[1^4]}(ν, t) = F_{[4]}(−ν, t^{−1})）から作ると U^{−16}〜U^{−8}（データの傾きから読んだ 9 項と同じ）になり、N = 4〜11 の全データと全係数一致（ν → +ν では不一致）。R_U_m1_1111.pkl に保存。r^{−8} の部分の方程式は未知数 297・階数 286 なので、この一致は t = q の値が残り 11 次元の自由度の中にあることの確認になっている。ただし t = q は照合ではなく入力として使ったので、(1,1,1,1) の t = q での照合は無い。N = 13 の照合で補う。
  - (2,1,1)：N = 3〜11 で一意だが照合が無いので保存していない。N = 12 で照合する。t = q の式（NRZ）は [r] と [1^r] だけなので使えない。
  - 4 個とも整数係数。(4) の項数 621 は N ≥ 7 の R(N) の項数と同じ（大きな N では c_j U^j の s の次数が重ならない）。
- N = 12, 13 の計算は 2026-10-26 23:08 に開始（PID 14948）。
- N = 12 の結果（c_dn4_save.log）：26684 s（約 7.4 時間、見積もりの約 2 倍）、最大 RSS 2.25 GB。プロセスの CPU 時間は約 12600 s で、うちシステム時間が約 8300 s。機械（8.6 GB）の空きメモリが約 1.2 GB しかなく、スワップや圧縮に時間を取られたとみられる（夜間のスリープの影響も除けない）。
- N ≤ 12 での 4 箱【計算】（2026-10-27）：
  | λ | 使う N（N_s） | 結果 |
  |---|---|---|
  | (4), (2,2) | 2〜8（2, 4） | 保存済みの解が N = 9〜12 で全係数一致 |
  | (3,1) | 2〜9（4） | 保存済みの解が N = 10〜12 で全係数一致 |
  | (2,1,1) | 3〜11（5） | 一意、N = 12 で全係数一致。R_U_m1_211.pkl に保存（U^{−20}〜U^{−8}、574 項、整数係数） |
  | (1,1,1,1) | 4〜12（5, 6） | すべての α で一意。解は t = q から r^{−8} を作った保存済みの R_U_m1_1111.pkl と全係数一致（N = 4〜11 で作った解が N = 12 で一致、とも言える） |
  - これで (1,1,1,1) はデータだけで決まり、t = q での NRZ の [1^4]（ν → −ν）との一致は入力ではなく照合になった【計算】。ただし N = 4〜12 で決めた解には、使わない N での照合がまだ無い（N = 13 で行う）。
  - (2,1,1) には t = q での照合が無い（NRZ の式は [r] と [1^r] だけ）。照合は N = 12 の一つだけ。
  - rk_c_4box_largeN.py の CFG4 の (1,1,1,1) を N = 4〜12 に変更（tq1111 は N ≤ 11 のまま）。
- 5 個の λ ⊢ 4 の R_λ(U) が揃った（data/refined_kauffman/daha/R_U_m1_{4,31,22,211,1111}.pkl）。必要だった N は最大 12（続き 70 の見込み N = 17 より 5 少ない）。
- N = 13 の見積もりの見直し【見込み】：N = 11→12 の伸び（実時間 4.7 倍、CPU 時間でも 2 倍強）と、ピークが 4.2〜4.6 GB で空きメモリを大きく超えることから、実時間 18〜30 時間以上。24 時間を超えるおそれがあるので、続けるかどうかをユーザーに確認する。N = 13 の役割は (1,1,1,1) と (2,1,1) の使わない N での追加の照合だけで、R_λ(U) を決めるのには要らない。
- N = 13 の計算は、ユーザーの指示で止めた（2026-10-27。N = 13 は保存されていない）。

## 2026-10-27（続き 72）n = 2（brane 2・image 2）の E を条件 (i)（t = −1、SO の 4 箱）と (ii)（t = q）で決める試み（ユーザーの指示「5 個の R_λ(U) を使って n = 2 の E を決める計算に進んで」）
- SO の 4 箱の自由エネルギー（scripts/refined_kauffman/rk_c_4box_so_f.py、logs/refined_kauffman/c_4box_so_f.log、data/refined_kauffman/daha/largeN_4box_so_m1.pkl）：rk_c_3box_largeN の SO 側を次数 4 に広げた（D_N の K_λ（|λ| ≤ 4）の S_μ 展開の逆、W_λ = R_λ(U) K_λ(t^ρ)、全 d の Möbius 反転）。計算 7 s（保存の sympy 変換に 44 分）。
  - 【計算】次数 ≤ 3 の f が largeN_3box_m1 と全係数一致。Möbius 反転は tests/refined_kauffman/test_rk_4box_so_f.py（2 件）。
  - 【計算】次数 4 の f の分母は単項式 × (1 − t) に、p_4 と p_(2,2) の方向だけ (1 + t) が加わる。p_(3,1), p_(2,1,1), p_(1^4) には (1 + t) が無い。theory §2 の「不変な状態は brane と image の境界の対をもち、p₂ 型の方向（ψ²[s_λ]）にだけ現れる」のとおり（brane の表現 λ ⊢ 2 で ψ²[s_(2)] = (p_2² + p_4)/2、ψ²[s_(1,1)] = (p_2² − p_4)/2）。
- 条件 (i) の規格化【導出】：不変な状態（brane の表現 λ、電荷 m、符号 ε）は f^{cov} に m s_λ(Y)s_λ(Y′)、g に ½εm ψ²[s_λ](Y)/(1 + t) を与える。p_1(Y)²p_1(Y′)² の係数は ¼m（§12.6 の 4P の 4）、ψ²[s_λ] の p_2² の係数は ½ なので E = Σεm = 4(1 + t) g_{p_2²}。f^{cov} の分母は (1 − t) だけ（2 箱・3 箱の保存値で確認【計算】）なので、t = −1 では g の代わりに SO の f を使える：E|_{t=−1} = [4(1 + t) f_{(2,2)}]_{t=−1}。n = 1 では E = 2(1 + t) g_(2)（Schur）で、従来どおり。
  - 【計算】規格化の検算：(i) と (ii) の目標を s B_tw で割り r = s = i で比べると、n = 2 の ν^{−22}, −18, −14, −10, −6 で 0 でない値（−36, 132, −180, 108, −24）が一致（n = 1 でも 0 でない 4 ブロックで一致）。他のブロックは両辺 0 で検査にならない。ν の奇数のブロックは (i) の目標にだけ現れ、t = −1 で 0。
- 整数計画（scripts/refined_kauffman/rk_f_twisted22_i.py、logs/refined_kauffman/f_twisted22_i.log。設定は §12.7 と同じ、単項式ごとの符号、B_tw = z_q (z_q^+)²）：
  - n = 1【計算】：(i)+(ii) で E が一意（34 変数、§9 の再現）。data/refined_kauffman/daha/twisted_E_n1_m1.pkl に保存。
  - n = 2【計算】：ν^{−24}, −20, −16 で解なし。他の 9 ブロックは解があるがほとんど決まらない（ν^{−22}：128 中 125 が自由、ν^{−18}：181 中 177、ν^{−4}：11 中 11、ν^{−2}：2 中 0）。
  - ν^{−24}（最低のブロック）の矛盾は単項式 1 個で起きる【計算】：Pk = −s^{19} r^{−41} Σ_{i,j=0}^{11} q^{max(i,j)} t^{min(i,j)} で、最低の単項式 q⁰t⁰ の状態は 1 個（c = 1）。(i) は e = +c、(ii) は e = −c を要求する（(i) の目標は i r^{−41}(1 + q)(1 + 3q² + 5q⁴ + … + 11q^{10})、(ii) は −r^{−22}(q^{12} − 1)[6]_{q²}）。範囲・整数性を外した一次方程式でも、B_tw を z_q^a (z_q^+)^b（a + b = 2, 3）に替えても、捻りを q・t 平面の回転に替えても、Pk の台の上では解が無い。ν^{−20}, −16 は |e| ≤ |c| + 2h に緩めると Σh = 1 で解がある。
  - 解が無い 3 ブロックはどれも 2β の β が偶数（F = g_(2) + g_(1,1) の c = 2 の部分）で、r = s = i の検算が両辺 0 で符号を確かめられないブロックに当たる。診断として (ii) の β 偶の部分の符号を反転させた (1 + q)ψ²[F(−ν)] にすると全 12 ブロックで解があるが、ほとんど決まらない（ν^{−24}：78 中 72 が自由）。同じ反転を n = 1 に当てると反転したブロックで解が無くなるので、一様な規則ではない【計算】。係数の調整としては採用しない。
    - **訂正（2026-10-29、続き 82）**：n = 1 に当てたのは ν で書いた同じ符号の (1 + q)ψ²[F(−ν)] で、全体の符号まで変えていた。クロスキャップの数 c で書くと、この反転は「c = 2 の部分の符号の反転」で、n = 1 の c = 2 の部分は t = q で 0 なので n = 1 は変わらない。c について一様な規則 (ii′) として採用を再検討した（theory §12.18）。
- 結論：
  - 【計算】n = 2 では、条件 (i)（SO の 4 箱のデータ）と (ii)（t = q の Adams 型の仮説）を合わせても E は決まらない。(ii) をそのまま使うと 3 ブロックで矛盾し、とくに最低のブロックでは一つの状態に逆の符号を求める。
  - 【予想】破れているのは (ii) の c = 2 の部分の形（AgSc の強い整数性の t = q での形が、β 偶の凍結した状態について n = 2 では符号を含めて違う）か、E が Pk の台の外の単項式（打ち消した状態）を含むことのどちらか。データからは区別できない。
  - 【未解決】二つの特殊化（t = −1, t = q）だけでは、各ブロックの 100 前後の変数を決めるには情報が足りない。
- 次の候補：E を brane の表現 λ ∈ {(2), (1,1)} で分ける（E = E_(2) + E_(1,1)）。p_4 方向の (1 + t) の極が [2(1 + t) f]_{p_4}|_{t=−1} = ½(E_(2) − E_(1,1))|_{t=−1} を与え、被覆の (2,2) セクターの s_λ(Y)s_λ(Y′) 成分（composite [2;2], [11;11]）で P も λ ごとに分けられるので、変数ごとの台が狭まり、条件が倍になる。

## 2026-10-27（続き 73）n = 2 の E を brane の表現 λ ∈ {(2), (1,1)} で分けて解き直す（ユーザーの指示「E を λ ∈ {(2), (1,1)} で分け、p_4 方向の条件と被覆の P の λ ごとの分解を加えて n = 2 の E を解き直して」）
- スクリプト：scripts/refined_kauffman/rk_f_twisted22_lam.py（cover：被覆の分解、solve [--tr]：整数計画、diagnose：条件の組ごとの切り分け）。ログ：logs/refined_kauffman/f_twisted22_lam_cover.log、f_twisted22_lam.log、f_twisted22_lam_tr.log、f_twisted22_lam_diag.log。データ：data/refined_kauffman/daha/largeN_cover_fixed22_lam_m1.pkl、twisted_E_lam{,_tr}_n2_m1.pkl。ν のブロックのキャッシュ：data/cache/twisted22_lam_blocks.pkl（SO の f の約分に約 10 分）。
- 機能テスト：tests/refined_kauffman/test_rk_twisted22_lam.py（冪和 → Schur の分解、条件ごとの重み）。rk_b_twisted_multiplet.solve の条件に、変数ごとの重み w(key) を任意の第 3 要素として加えた（既存の呼び出しはそのまま。tests/refined_kauffman 全体 88 件が通る）。
- 被覆の (2,2) セクターの λ ごとの分解【計算】：log Z の p_ρ(Y)p_σ(Y′)（ρ, σ ⊢ 2）の係数 a_{ρσ} から Adams の ½ψ²[x₁₁] p_2p_2′ を除き、Schur の M_{λμ} にした。
  - Y ↔ Y′ の対称性 M_{(2),(1,1)} = M_{(1,1),(2)}、和 ΣM = 4a_{(1,1),(1,1)}、P_(2) + P_(1,1) + 2P_off = 4P（§12.6 の保存値）がすべて成り立つ。Adams の項を除くと (1 − t)M の分母は単項式だけ。
  - **Schur のまま λ で分けた P_λ = (1 − t)M_{λλ} は z_q（= (r − 1)(r + 1)/r）でしか割り切れない**。z_q³ は和にだけ現れる。冪和の方向ごとに境界の因子が違うため（KN の hatted の φ_(1,1) = z_q、φ_(2) = z_q^+）。
  - **(1 − t)a_{ρσ} / (z_q φ_ρ φ_σ) は 4 成分すべて Laurent**。hatted の P̂_λμ = Σ χ_λ(ρ)χ_μ(σ)(1 − t)a_{ρσ}/(z_q φ_ρ φ_σ) を使う（余分の z_q は brane と image を結ぶ境界のモードで、brane–image (1),(1) の P ∝ z_q と同じ。p_1²p_1′² では z_q³ を再現）。P̂_λλ/s は整数係数の Laurent。ν のブロックは P̂_(2)：−24〜−2、P̂_(1,1)：−24〜−6、P̂_off：−24〜−4。
  - 最低のブロック【計算】：P̂_(2)/s の ν^{−24} = −s^{19} r^{−41}(1 + qt) Σ_{0≤j≤i≤5} q^{2i} t^{2j}（42 項）、P̂_(1,1)/s は −s^{19} r^{−37}(1 + qt) Σ_{0≤j≤i≤4} q^{2i} t^{2j}（30 項）、P̂_off/s は 36 項。合計 42 + 30 + 2·36 = 144 が §12.6 の Σ_{i,j=0}^{11} q^{max}t^{min} の項数に一致する。
- 捻ったトレースの規格化【仮説】：不変な状態（brane λ・image λ、hatted の電荷 m、符号 ε）は g に ½εm s z_q^+ ψ²[Σ_ρ φ_ρ χ_λ(ρ) p_ρ/z_ρ]/(1 + t) を与える（ψ²φ_(1,1) = z_q z_q^+、ψ²φ_(2) = q + q^{−1}）。これから
  4(1 + t)g_{p_2²} = s z_q (z_q^+)² (E_(2) + E_(1,1))（§12.7 の B_tw を再現）、4(1 + t)g_{p_4} = s z_q^+ (q + q^{−1}) (E_(2) − E_(1,1))。
  - 【計算】t = −1 の SO の目標 4(1 + t)f_{(2,2)}、4(1 + t)f_{(4)} は、それぞれ s z_q (z_q^+)²、s z_q^+ (q + q^{−1}) で割ると Laurent（境界の因子の仮説と矛盾しない）。
- 整数計画（条件 (i22)：t = −1 の p_2²、(i4)：t = −1 の p_4、(ii)：t = q の和。E_λ の変数は P̂_λλ/s の単項式、|e| ≤ |c|、偶奇）【計算】：
  | ν | 変数 (2) / (1,1) | 結果 |
  |---|---|---|
  | −24 | 42 / 30 | 解なし |
  | −22 | 118 / 89 | 解なし |
  | −20 | 164 / 127 | 解なし |
  | −18 | 174 / 137 | 解あり、一意でない 157 / 129 |
  | −16 | 167 / 130 | 解なし |
  | −14 | 161 / 124 | 解あり、一意でない 160 / 123 |
  | −12 | 137 / 102 | 解あり、一意でない 133 / 100 |
  | −10 | 106 / 77 | 解あり、一意でない 103 / 76 |
  | −8 | 72 / 47 | 解あり、一意でない 69 / 45 |
  | −6 | 34 / 17 | 解あり、一意でない 30 / 17 |
  | −4 | 11 / 1 | 解あり、一意でない 11 / 0 |
  | −2 | 2 / 1 | 解あり、一意（0 / 0） |
  - 続き 72（λ で分けない）より解の無いブロックが 1 個増えた（ν^{−22}）。解のあるブロックもほとんど決まらない。
- 切り分け（diagnose、f_twisted22_lam_diag.log）【計算】：
  - ν^{−24}：(i22)+(ii) と (i4)+(ii) が、範囲・整数性を外した一次方程式でも解なし。t = −1 と t = q の矛盾で、λ の分け方によらない（続き 72 の単項式 1 個の矛盾と同じ）。P̂ の台の上では、P に見えない対を許しても解が無い。
  - ν^{−20}：**t = −1 の条件だけ（(i22) だけ、(i4) だけ）でも解が無い**。連続の緩和（|e| ≤ |c| だけ）でも解が無く、対角の P̂_λλ の状態数では、SO 側が求める E を収容できない。対角でない P̂_off の単項式を（比較のために）変数に加えると解が出る。P に見えない対を許すと最少 Σh = 2（(i22) だけ）、全条件で Σh = 6。
  - ν^{−22}, −16：各条件だけなら解があり、(i22)+(ii) で解が無い（連続の緩和でも）。最少の P に見えない対は Σh = 4、1。
- 別の仮説「Ĩ は brane の λ を image の λᵀ に移す」（不変な状態は P̂_off、--tr）【計算】：全ブロックで解なし。ν^{−2} の目標は t = −1 で 0 でないのに P̂_off に ν^{−2} のブロックが無い。この仮説は否定される。
- 結論：
  - 【計算】λ ごとに分けると、p_4 方向の条件と変数ごとの範囲が加わる分だけ条件が厳しくなり、解の無いブロックは 4 個に増えた。E は決まらない。
  - 【計算】最低のブロックの矛盾は t = −1 と t = q の間のもので、λ で分けても解消しない。続き 72 の予想のうち「(ii) の c = 2（β 偶）の部分が n = 2 では成り立たない」側を支持する（解の無い 4 ブロック中 3 ブロック −24, −20, −16 が β 偶、−22 は β 奇）。
  - 【予想】ν^{−20} では t = −1 だけで状態数が足りないので、E は（hatted の）被覆の固定セクターの指数に見えない状態（ボソンとフェルミオンの打ち消し）を含む。§3 の指摘（E は P の台に収まらない状態を含みうる）が、λ で分けたことで、(ii) の仮説によらない形で見えた。
  - 【未解決】打ち消した状態を入れると変数が無制限になり、t = −1, t = q の二つの特殊化では E は決まらない。E を決めるには、指数でなく状態の数（Poincaré 多項式、§3）か、Ĩ_int の作用の幾何的な決定が要る。

## 2026-10-27（続き 74）E を状態の数・Ĩ_int の幾何から決める方法の検討（ユーザーの指示「指数でなく状態の数（Poincaré 多項式）か Ĩ_int の幾何的な決定から E を決める方法を検討して」）
- スクリプト：scripts/refined_kauffman/rk_f_twisted_tensor.py（n1 / n2 [--all-inv] / n2scan / n2partial / n1partial / lower）。ログ：logs/refined_kauffman/f_twisted_tensor_n1.log、f_twisted_tensor_n2.log、f_twisted_tensor_lower.log。機能テスト：tests/refined_kauffman/test_rk_twisted_tensor.py（2 件）。
- 出発点の観察【計算】：§12.7 の trefoil の n = 1 の一意な E は、P と同じ形に因数分解する。
  P = s²(t − ν²)(ν² − 1)(q − 1)(ν²q − 1)·N₈/(ν¹²r¹⁰)、E = s²(t − ν²)(ν² + 1)(q + 1)(ν²q + 1)·Ñ₈/(ν¹²r¹⁰)、Ñ₈ は N₈ の 8 項の符号の付け替え（全 8 状態が不変）。
- 仮説（方法 A、テンソル積の模型）：固定セクターの状態の空間は 2 状態の因子（境界のモード q、ν²、ν²q、t/ν²）と残り N の積で、状態の数は N の係数の絶対値（打ち消しは N の中だけ）。Ĩ は因子ごとに働き、N の部分では |e_m| ≤ |c_m|・偶奇。
- n = 1【計算】：
  - 因子の捻り 16 通りのうち (i)+(ii) を満たすのは、(t − ν²) は不変・(ν² − 1)・(q − 1)・(ν²q − 1) は和、の 1 通りだけ（trefoil、T(2,5) で同じ）。(q − 1) → (q + 1) は §11.3 の境界のモード、(ν² − 1) → (ν² + 1) は §11.8 の framed unknot と同じ。
  - 変数は trefoil で 34 → 8、T(2,5) で 107 → 40（N の単項式。T(2,5) の N の状態数 Σ|c| = 46）。trefoil は一意（§12.7 の解と同じ）。**T(2,5) は初めて解をもつ**（指数の台では解が無かった、§9）が、40 中 11 が一意でない。全状態が不変（e = ±c）とすると T(2,5) は解なし。
  - (i) だけでは trefoil でも 8 中 4 が一意でない（(ii) が要る）。
- N の構造【計算】：T(2,3)・T(2,5) の N の ν のブロックは、T(2,5) の ν^{−14} を除いて全係数が同符号（ブロックごとの符号 × 正の多項式）。T(2,5) の下の 4 ブロックは **ν^{−(6+2k)} = Sym^k V / r^{4+2k}（k = 0..3、V = 1 + t + qt）** と全係数で一致する（trefoil は k = 0, 1）。ν^{−14} は +Sym⁴V と負の上の部分の和と読め、3 個の単項式（t³, t⁴, qt⁴ の位置）で打ち消しが起きている【予想：分け方】。
- 方法 B（Ĩ の幾何的な作用）の手がかり【計算】：trefoil の一意な E の N 部分の符号 e/c は、ブロック内で q の指数とともに反転し t の指数では変わらない：e = σ_β (−1)^{(q の指数)} c。下の系列では V の対合 diag(1, 1, −1)（qt のセルだけ反転）の Sym^k への誘導作用 × ブロックごとの符号 σ_k に当たる。
  - T(2,5) の下の 4 ブロックに V の対角な対合（4 通り）を課し、上のブロックを一般の e にすると、(i)+(ii) を満たすのは **diag(1, 1, −1) だけ**（trefoil と同じ）。σ_k は一意に決まり、残りの自由度は上のブロックの 20 変数中 7。
  - 同じ規則（q の偶奇）を全ブロックに当てると T(2,5) は解なし。上の系列の構造（Sym^k V の双対にあたるもの）は同定できていない。
- n = 2【計算】：hatted の P̂_λ（続き 73）も、境界のモードを除いたあと同じ 3 因子 (t − ν²)(ν² − 1)(ν²q − 1) をもつ（P̂_(2)、P̂_(1,1)、P̂_off とも）。N_λ の単項式は (2) 562、(1,1) 388（P̂ の単項式 1188、880）。
  - 因子の捻り 8 通りのうち、**t = −1 の条件（p_2² と p_4）だけで解があるのは n = 1 と同じ捻り (t − ν²)(ν² + 1)(ν²q + 1) だけ**（他は一次方程式でも解なし）。SO の 4 箱のデータが、n = 1 と独立に同じ捻りを選ぶ。続き 73 の ν^{−20} の状態数の不足も、この模型では起きない。
  - しかし t = −1 と (ii) を合わせると、どの捻りでも一次方程式から解が無い。(ii) の c = 1（β 奇）の部分だけを加えても解なし。t = −1 だけでは一次方程式の零空間の次元が 950 中 796 で、ほとんど決まらない。
  - N_λ のブロックは符号が混在し、下のブロックには (q − t) の因子が残る（例：N_(2) の ν^{−10} = −s³(q − t)(1 + q + qt)(1 + t + qt)/r^{11}、ν^{−8} = −s³(q − t)/r⁷）。n = 1 のような単純な Sym^k の形にはまだ分けられていない。
- 検討の結論：
  - 【計算】状態の数を「指数の単項式の係数」から「因数分解した N の係数」に替えると（方法 A）、n = 1 の T(2,5) と n = 2 の t = −1 が整合し、因子への Ĩ の作用は n = 1, 2 で同じ一つに決まる。指数で見えなかった状態（続き 73 の ν^{−20}）は、因子の積の展開での打ち消しとして説明できる。
  - 【計算】Ĩ の N 部分への作用は、下の系列では V の対合 diag(1, 1, −1) の Sym^k への誘導作用（Lefschetz の不動点公式の形、e^{iπJ₁} 型）で、trefoil 全体と T(2,5) の下の系列で確かめた。§11.8 の framed unknot の t の塔の (−1)^k（e^{iπJ₂} 型）と対になる。
  - 【予想】N の下の系列 Sym^k V は、重み 1, t, qt の 3 個のセルをもつ空間（偶数次のコホモロジー、例えば P²）の k 次の対称積の Poincaré 多項式で、被覆の Ĩ 不変な曲線のモジュライの点の部分に当たる。Ĩ は qt のセルを反転する。
  - 【未解決】(1) 上の系列（T(2,5) の ν^{−14}〜−20）の空間と Ĩ の作用、(2) n = 2 の N_λ の分解（(q − t) の因子と符号の混在）、(3) n = 2 で (ii) が t = −1 と矛盾する理由。(3) は、方法 A の模型で t = −1 だけが整合することから、(ii) の n = 2 での形の問題である可能性が高い【予想】。
- 次の手順の案：
  1. 上の系列の同定：trefoil と T(2,5) の N の上のブロックを、下の系列の双対（t の向きを反転した Sym^k）と打ち消しの組で書けるか調べ、同じ対合を当てて T(2,5) の残り 7 変数が決まるかを見る（小さい計算）。決まれば T(2,7) の 2 箱（データの有無を確認）で予言として照合する。
  2. n = 2 の N_λ を、n = 1 の V, V′ = 1 + q + qt と (q − t) の積・対称積で書けるかを、ブロックごとに調べる（下のブロックから）。書ければ同じ対合で E_λ を作り、t = −1 の条件で検証する（(ii) を使わない）。
  3. (ii) の n = 2 での形：2 で決めた E_λ の t = q の値と −(1 + q)ψ²[F] を比べて、どこが違うかを見る。

## 2026-10-27（続き 75）N の上のブロックの同定と、同じ対合での T(2,5) の E の決定（ユーザーの指示「trefoil と T(2,5) の N の上のブロックを同定し、同じ対合で T(2,5) の残りの 7 変数が決まるか調べて」）
- スクリプト：rk_f_twisted_tensor.py の knot_rule（rule）、full_rule（full）、lower_series（lower）。ログ：logs/refined_kauffman/f_twisted_tensor_rule.log、f_twisted_tensor_lower.log。データ：data/refined_kauffman/daha/twisted_E_n1_m2.pkl（T(2,5) の一意な E）。
- N のブロック n_k（ν^{−(6+2k)} r^{−(4+2k)} の係数、q, t）【計算】：
  | k | T(2,3) | T(2,5) |
  |---|---|---|
  | 0〜2m−1 | 1、V | Sym^k V（k = 0..3） |
  | 2m | −tV′（V′ = 1 + q + qt） | **t³V − t(1 + q + 3qt + qt² + 2q²t + 3q²t² + q³t² + q³t³)**（正の部分は t³V だけ） |
  | 2m+1.. | −t³ | −t³(1 + qt)V、−t⁵V、−t⁷ |
  - 上のブロックは下の Sym^k V の単純な双対（(q, t) → 単項式の写像）ではない：T(2,5) の k = 5 は qt⁴ に 2 状態をもち、3 個の異なる重みの Sym² からは出ない。T(2,3) の N だけは n_{3−k}(q, t) = −t³(q/t)^k n_k(1/q, 1/t) の双対性をもつ。
  - 生成関数 Σ n_k x^k に (1 − x)(1 − tx)(1 − qtx) を掛けると 1 − (…)² x^{2m} + …：x^{2m} の係数は T(2,3) で −(1 + t + qt)²、T(2,5) で −(1 + t + qt + qt² + q²t²)²（項数 3, 5 は各結び目の簡約 superpolynomial の項数と同じ）。それより先の係数には (q − t) の因子が出る。
- T(2,5) の残りの 7 変数（すべて k = 4 のブロック）【計算】：
  - 下の規則（V の対合 diag(1, 1, −1)、e = σ (−1)^{q の指数} c）を k ≤ 3 と **k = 4 の正の部分 t³V** に課すと、(i)+(ii) で E が**一意に決まる**（σ も、k = 4 の負の部分の 8 変数も）。(i) だけでは 17 中 9 が一意でない。k = 4 の負の部分の |c| = 1 の単項式に (−1)^q を加えても整合し、(−1)^t は解なし。
  - 一意な解の符号の構造：**k ≤ 2m のブロック（下の系列と移り目のブロックの正負の両方）は (−1)^q で、k > 2m のブロック（T(2,5) の k = 5, 6, 7）は (−1)^t で、ブロックと c の符号ごとに σ がそろう**。(−1)^t は V の対合 diag(1, −1, −1)（下の対合を双対で移したもの：基底 1 ↔ qt の反転と最高次の符号）に当たる。T(2,3) は全ブロックが (−1)^q（上の k = 2 は V′ 上の diag(1, −1, −1)、同じ型）。
  - Ĩ で入れ替わる対が 4 組ある（k = 5 の qt⁴ の 2 状態、k = 4 の (q, t) = (1, 2), (2, 2), (2, 3)）。k = 5 の qt⁴ の対は (1 + qt)⊗(1 + qt) の非対角の 2 状態で、入れ替えの対合（ψ² 型）の形をしている【予想】。
  - 規則（k ≤ 2m で (−1)^q、k > 2m で (−1)^t、多重度のある単項式では入れ替わる対を許す）を全ブロックに課した整数計画（full_rule）：T(2,3) は **(i) だけで一意**（(ii) が要らない）、T(2,5) は (i) だけでは 40 中 6 が一意でなく、(i)+(ii) で一意（上の解と同じ、対 4 組）。
- 結論：
  - 【計算】T(2,5) の n = 1 の E は、テンソル積の模型（続き 74）と「下の系列と移り目のブロックで V の対合 diag(1, 1, −1)、上の系列で双対の対合」の規則のもとで、(i)+(ii) から一意に決まる。指数の台の上では解が無かった（§9）ものが、状態の数と Ĩ の作用の構造を入れることで決まった。
  - 【予想】N は「下の系列（Sym^k V、Ĩ は V の 1 個のセルの反転）」と「その双対の上の系列」の差で、移り目のブロック（k = 2m）で両者が重なって打ち消す。入れ替わる対は上の系列の (1 + qt) のテンソル積の非対角に現れる。
  - 【未解決】上の系列と移り目の負の部分の空間の同定（生成関数の分子の幾何的な意味）、n = 2 への適用（N_λ のブロックの符号の混在と (q − t) の因子）、T(2,5) で (ii) がまだ要ること。
- 次の手順の案：(1) 同じ規則を n = 2 の N_λ に当て、t = −1 の条件だけでどこまで決まるかを見る（N_λ の下のブロックから、移り目の位置を探す）。(2) T(2,7) の 2 箱のデータがあれば、規則で E を作って (i)+(ii) を予言として照合する。

## 2026-10-28（続き 76）n = 1 の符号の規則を n = 2 の N_λ に当て、t = −1 だけでどこまで決まるか（ユーザーの指示「同じ符号の規則を n = 2 の N_λ に当て、t = −1 の条件だけでどこまで決まるか調べて」）
- スクリプト：rk_f_twisted_tensor.py の n2_rule（n2rule）。整数計画は full_rule から切り出した rule_milp（符号の規則・入れ替わる対・対の数の最小）。ログ：logs/refined_kauffman/f_twisted_tensor_n2rule.log。テスト：tests/refined_kauffman/test_rk_twisted_tensor.py::test_rule_milp_sign_rule_and_pairs。
- 設定：trefoil、続き 74 の捻り E_λ = s(t − ν²)(ν² + 1)(ν²q + 1) Σ e m、条件は t = −1 の i22 + i4 だけ。規則は部分（λ、ν のブロック、c の符号）ごとに e = σ (−1)^{q の指数} sign(c) w か σ (−1)^{t の指数} sign(c) w（w ≤ |c|、w ≡ c mod 2、|c| − w は入れ替わる対）。
- N_λ のブロック【計算】：N_(2) は ν^{−8}〜ν^{−24}、N_(1,1) は ν^{−12}〜ν^{−24}。上の端と下の端以外は符号が混在する。
  - 上の端の 2 ブロック（N_(2) の ν^{−8}、ν^{−10}）には (q − t) の因子がある。ν^{−10} は −(q − t)(1 + q + qt)(1 + t + qt) を単項式に掛けたもの。
  - 下の端の ν^{−24} は −t^{17/2} q^{−41/2}(1 + qt)·Sym^5 W（N_(1,1) では Sym^4 W）で、W = 1 + q² + q²t²。全部が負で、n = 1 の V = 1 + t + qt を (q, t) → (q², t²) とした形である。
- 結果【計算】（すべて t = −1 だけ）：
  - **ν の切れ目で (−1)^q と (−1)^t を分ける形では、全部分に規則を課すと解が無い**。λ に共通の切れ目 11 通りすべてで解なしになった。
  - **部分ごとの判定**：一つの部分だけに規則を課し、他は一般の e とした。
    - (−1)^q だけ可：下の系列側（N_(2) の ν^{−10}+、ν^{−12}−、ν^{−14}+、N_(1,1) の ν^{−12}、ν^{−14}+）。
    - (−1)^t だけ可：上の系列側（ν^{−22}−、両 λ）。
    - 中間の ν^{−16}〜ν^{−20} は部分ごとに交互になる。q だけ可が ν^{−16}−、ν^{−18}+、ν^{−20}−、t だけ可が ν^{−18}−、ν^{−20}+ で、λ = (2), (1,1) で同じ型である。
    - ほかに N_(2) の ν^{−10}− は t だけ可、q・t のどちらも可の部分が 9 個ある。
  - **両方可の 9 個の組 512 通りのうち、全部分に規則を課して解があるのは 128 通り**。そこで選択が決まるのは ν^{−16}+ の 2 個（両 λ とも t）だけである。
  - **入れ替わる対の数**：規則なしの t = −1 では最小 20（状態の数 15790）だが、規則を課すと最小 2845〜2893 に増える（128 通りのすべてで）。n = 1 の T(2,5) の対は 4 組だった。
  - **一意性**：対の数の最小が最も小さい組で、e が一意でない単項式は 950 中 743 である。一意なのは端のブロック（N_(2) の ν^{−8}、ν^{−10}、両 λ の ν^{−24}、計 84 単項式）だけで、中間の ν^{−16}〜ν^{−20} は 9 割以上が一意でない。
- 結論：
  - 【計算】n = 1 の符号の規則（ブロックと c の符号ごとの σ に (−1)^q か (−1)^t）は n = 2 の N_λ に**そのままの形では合わない**。課すこと自体はできるが、1/3 以上の状態が入れ替わる対になることを要求する（規則なしなら対 20 個で足りる）。t = −1 だけで決まるのは端のブロックだけである。
  - 【予想】n = 2 では Ĩ が働くセルの空間が V ではなく W = 1 + q² + q²t²（V の Adams 型 (q, t) → (q², t²)）で、下の端の ν^{−24} がその対称積になっている。符号は q² の指数（q の指数の半分）の偶奇で決まる可能性がある。中間ブロックの符号の混在と (q − t) の因子は、V 側の系列と W 側の系列の重なり（n = 1 の移り目のブロックの類似）と考えられる。
  - 【未解決】N_λ を V 側・W 側の系列に分ける同定と、それに合わせた n = 2 の規則。
- 次の手順の案：(1) N_λ の生成関数を n = 1 と同様に作り、(1 − x)(1 − tx)(1 − qtx) と (1 − x)(1 − q²x)(1 − q²t²x) のどちらで分子が簡単になるかを調べて系列を同定する。(2) W 側に q² の指数の偶奇の規則を課し、対の数がどこまで減るかを見る。

## 2026-10-28（続き 77）N_λ の系列の同定（ユーザーの選択「系列の同定に進む」、続き 76 の次の手順 (1)）
- スクリプト：rk_f_twisted_tensor.py の n2_series（n2series）、補助の times・greedy_cells・mirror_sign。ログ：logs/refined_kauffman/f_twisted_tensor_n2series.log。テスト：tests/refined_kauffman/test_rk_twisted_tensor.py::test_greedy_cells_and_mirror。
- 方法：(1 + ε m) を掛けて状態の数が減る因子を貪欲に探し、系列のセルを読んだ。ν 方向（端から 4 ブロック）と ν のブロックの中の二通りで調べた。対照の n = 1 の T(2,5) では、上の端から V のセルが正しく選ばれる。
- 結果【計算】：
  - **ν 方向の系列は無い**。n = 2 の N_λ では、どちらの端からも状態の数は 2〜7% しか減らない。続き 76 の案 (1)「(1 − x)(1 − tx)(1 − qtx) と (1 − x)(1 − q²x)(1 − q²t²x) のどちらか」は、どちらでもなかった。
  - **ブロックの中の塔**：ν^{−14}〜ν^{−22} では、両 λ とも貪欲法が (1 − qt)(1 − t)(1 − q) を選び、状態の数は 1/2〜1/7 に減る。端の ν^{−12}、ν^{−24} では (1 − qt) だけが選ばれる。全体では N_(2) が 10968 → 2844、N_(1,1) が 4822 → 1684 になる。
  - **塔の形**：単体の形の打ち切りなら貪欲法は q, t の 2 個しか選ばず、3 方向の弦の積なら 3 個を最も長い方向から選ぶ（対照の計算）。N_λ では qt が最初に選ばれるので、qt 方向が最も長い弦の積に近い（格子の対角線の縞）。
  - **鏡映対称**：全 16 ブロックが (q, t) → (t^{−1}, q^{−1}) で対称（符号 +）。
  - **端のブロック**：
    - N_(2) の ν^{−8} = (t − q)、ν^{−10} = (t − q)VV′（V′ = 1 + q + qt）。ν^{−12} は (t − q)Sym²V Sym²V′ でも (t − q)(VV′)² でも書けない。
    - ν^{−24} = −(1 + qt)h_K(W) = −(1 + qt)ψ²h_K(V′)、K = 5（λ = (2)）、4（λ = (1,1)）。
  - **λ の違い**：λ = (2) と q^{−2}·λ = (1,1) は各ブロックの q の指数の小さい側で一致し、差は大きい側に集まる。λ は塔の打ち切りの段数を変えるように見える。
  - **訂正**：続き 76 と theory §12.12 の「W は V の Adams 型 V(q², t²)」は誤りで、W = V′(q², t²) = ψ²V′（theory §12.12 に訂正を追記）。
- 結論：
  - 【予想】n = 2 の固定セクターの状態の空間は、ν のブロックごとに重み q, t, qt のセルの塔を打ち切ったもので、W = ψ²V′ は brane と image を入れ替える作用のトレースとして現れる。
  - 【予想】Ĩ の規則は、V のセルではなく塔のセルと W で書くべきである。続き 76 で V の規則が n = 2 に合わなかったことと整合する。
  - 【未解決】中間ブロックの分子、打ち切りの規則、塔のセルに合わせた符号の規則。
- 次の手順の案：§2′ の次の候補 (1)、(2)。

## 2026-10-28（続き 78）塔のセルに合わせた n = 2 の符号の規則（ユーザーの選択「符号の規則を試す」、続き 77 の次の候補 (2)）
- スクリプト：rk_f_twisted_tensor.py。n2_rule に規則の候補 opts と一様な規則の表示を加えた（n2tower）。部分ごとの規則の選択を変数にした整数計画 choice_milp と n2_choice（n2choice）を新設した。ログ：logs/refined_kauffman/f_twisted_tensor_n2tower.log。テスト：tests/refined_kauffman/test_rk_twisted_tensor.py::test_choice_milp_selects_pattern。
- 規則の候補【導出】：
  - 塔のセル (q, t, qt) の対角な対合が単項式の上で矛盾なく定まるのは ε_qt = ε_q ε_t のときだけで、符号は 0・q・t・q+t の 4 通り（(−1)^{αa+βb}）になる。
  - 鏡映で対称なのは 0 と q+t だけである。
  - 下の端の ψ²h_K(V′) の部分は、ψ² のトレースが g² = 1 のトレースなので不変になる。
- 結果【計算】（trefoil、t = −1 の i22 + i4 だけ。部分は λ・ブロック・c の符号）：
  - 0・q・t・q+t のどれを一様に課しても解なし。
  - 鏡映で対称な 0 と q+t だけでも解なし。両 λ の ν^{−18}−、ν^{−22}−、ν^{−24}− と λ = (2) の ν^{−10}− は t だけが可である（ν^{−24} では q ≡ t）。
  - 4 通りから部分ごとに選ぶと組は約 6.8 億通りになるので、選択を変数にした整数計画で入れ替わる対の数を直接最小にした。
    - q・t だけでは最小 2845（続き 76 を 1 秒で再現）、4 通りでは最小 2052（48 秒）。規則なしなら最小 20。
    - 4 通りの最適解でも、対は中間ブロックに残る（λ = (2) の ν^{−18}− で 472 など）。
  - qt 方向の縞ごとに 0 か q+t とする規則は解なし。4 通りを候補にした縞ごとの計算は、1200 秒で打ち切って未決定。
- 結論：
  - 【計算】部分ごとの対角な対合の形の規則は、V のセルでも塔のセルでも、t = −1 の条件と合わせると大量の入れ替わる対を要する。塔のセルにしても対は 2845 → 2052 に減るだけである。
  - 【予想】n = 2 の Ĩ は対角な対合では書けず、弦をまたぐ置換を含む。あるいは σ の群の分け方が合っていない。
  - 【未解決】n = 2 の E を決める規則。
- 次の手順の案：§2′ の次の候補 (1)。規則なしで対 20 の解の符号の型から規則を読み取る案が最も直接的である。ただし、その解の一意性の判定には時間がかかる（続き 76 の一意性の判定は約 800 秒）。

## 2026-10-28（続き 79）規則なしで入れ替わる対が最小の解の符号の型（ユーザーの選択「符号の型を読む」、続き 78 の次の候補 (1)）
- スクリプト：rk_f_twisted_tensor.py の n2_pattern（n2pattern）。rule_milp に max_pairs・ranges・samples・flip を加えた（対の数の上限、e の範囲、乱数の目的関数による解の標本、まとめて反転）。choice_milp に max_pairs・ranges を加えた。並列には multiprocessing の fork の Pool を使う（ProcessPoolExecutor はサンドボックスで os.sysconf を読めず止まる）。ログ：logs/refined_kauffman/f_twisted_tensor_n2pattern.log。データ：data/refined_kauffman/daha/twisted_n2_pattern.pkl（最小の解、標本、e の範囲。再開できる）。テスト：test_choice_milp_ranges、test_rule_milp_ranges、test_rule_milp_samples、test_rule_milp_flip。
- 計算時間と見積もり：
  - 最小の解は choice_milp で 502 秒、rule_milp で 108 秒。
  - 単項式 1 個の e の範囲は約 125 秒で、全 950 個では約 33 時間の見積もりになった（24 時間超なので、この方法では実行しなかった）。
  - 代わりに次の手順で、約 1 時間で全単項式を判定した（theory §12.15）：標本 16 個（8 並列、1000 秒。並列でもほとんど速くならない）、行空間の判定、まとめて反転 2 回（338 秒、654 秒）、残り 2 個を 1 個ずつ。
- 結果【計算】（trefoil、t = −1 の i22 + i4 だけ、規則なし、対の数 20）：
  - **e が一意な単項式は 950 中 308**（e = c が 161、e = −c が 145、対を含むものが 2）。
  - 下の端の ν^{−24} は全部一意で、符号は q の指数の偶奇（(1 + qt) の qt の項だけ反転）。λ = (2) の ν^{−8} も全部一意。
  - 中間ブロックでは、一意な単項式は q の指数の低い側に偏る。下半分で λ = (2) が 104/194、λ = (1,1) が 85/139、上半分で 75/368、44/249。
  - λ = (2) と q^{−2}·λ = (1,1) で c が等しく両方一意な単項式 92 個のうち、89 個で e も等しい（λ によらない）。
  - どの解でも対を含むのは、両 λ の ν^{−18} の対応する一つの単項式（c = −9、e = 5）だけ。残り 16 対は動く。
  - 一意な単項式だけで見ても、中間ブロックの部分ごとの e/c の符号は σ(−1)^d（d = 0、q、t、q+t）の形にならない。
- 結論：
  - 【計算】t = −1 と「入れ替わる対が最小」で、n = 2 の E の約 1/3 が決まる。下の端全体と、中間ブロックの低い q の領域である。
  - 【計算】決まった符号は対角な対合の形ではない。続き 78（対角な対合なら対は 2052 以上）と合わせて、n = 2 では二つの原理が両立しない。
  - 【未解決】どちらの原理が正しいか、決まった符号の型の意味。
- 次の手順の案：§2′ の次の候補 (1)、(2)。


## 2026-10-28（続き 80）決まった e の構造（ユーザーの選択「e の構造を読む」、続き 79 の次の候補 (1)）
- スクリプト：rk_f_twisted_tensor.py の row_min（1 行の対の最小化、動的計画法）と n2_rows（n2rows）。ログ：logs/refined_kauffman/f_twisted_tensor_n2rows.log（45 秒、ほぼ全部が条件の行列の構築）。テスト：test_row_min、test_row_min_matches_rule_milp（2 行の小さな例で rule_milp と対の最小・e の範囲が一致）。
- 経過：
  - 一意な行の低い q の部分に (1 ∓ q)(1 ∓ t)(1 ∓ qt) を掛けても E は簡単にならなかった（ν^{−22} で状態数 45 → 40、N は 45 → 18）。
  - そこで、条件の行列の階数 154 が行の数と同じことに気づいた。
- 結果（theory §12.16）：
  - 【導出＋計算】条件 (i) は行ごとの交代和 S = E_λ(q, −1, ν) だけを決める（rank A = rank R = rank [A; R] = 154）。
  - 【計算】対の最小化は行ごとに独立になる。行ごとの動的計画法で、対の最小 20、一意な単項式 308 とその値が、続き 79 の整数計画と完全に一致した。
  - 【計算】全単項式が一意な行は 63/154。内訳は飽和 |S| = Σ|c| が 37 行（e = sign(S)(−1)^{(s−1)/2}|c|）、最小の解が 1 通りの行が 26（うち対を含む 2）。一部だけ一意な行に 87 個ある。一意な割合は行の状態数が 1〜5 で 98/105、61 以上で 114/592。
  - 【計算】λ によらない 89 組のうち 81 組は、N の行と S がともに等しい行（各ブロックの低い q の 33 行）に入る。e が違う 3 組は行の問題が違う。
  - 【計算】E_λ(q, −1, ν) は因数分解せず、N(q, ±1) の符号の付け替えでもない（例：ν^{−22} の q の冪 1 で |S| = 3、|N(q, ±1)| = 1）。λ = (1,1) では S は 1 行を除いて 0 以下。
- 結論：
  - 【計算】t = −1 だけで決まる e には、それ自体の構造は無い。条件が与える S と、行の中の符号の組合せで尽きる。
  - 続き 79 の結果は、次のように説明できる。「低い q の領域」は行が小さいこと、「λ によらない」は行の問題が同じことによる。「対角な対合と両立しない」は、対の最小化が行ごとに別々に符号を組むことによる。
  - E の構造を得るには、行をまたぐ条件 (ii)（t = q）が要る。
- 次の手順の案：§2′ の次の候補 (2)（(ii) を加える）。(ii) は行をまたぐので行ごとの分解は使えず、整数計画になる。

## 2026-10-28（続き 81）n = 2 に (ii)（t = q）を加える（ユーザーの選択「(ii) を加える」、続き 80 の次の候補 (2)）
- スクリプト：rk_f_twisted_tensor.py の n2_quotient（割った形の目標、data/cache に保存）、quotient_rows（割った形の整数の条件行列）、hidden_milp（見えない対 h を許す整数計画）、pairs_onehot（入れ替わる対の最小を 0/1 変数で解く）、n2_ii（n2ii）。ログ：logs/refined_kauffman/f_twisted_tensor_n2ii.log（約 20 分、うち 600 秒は打ち切り）。テスト：test_hidden_milp、test_pairs_onehot_matches_row_min、test_n2_quotient_rows_equivalent。
- 前提：続き 74 で、この模型の (i)+(ii) は一次方程式から解が無いと分かっている。そこで、矛盾の場所と形を特定することを目標にした。
- 結果（theory §12.17）：
  - 【導出＋計算】条件を前因子と境界の因子で割ると、目標はどれも Laurent 多項式になり、条件は ν のブロックごとに閉じる（階数は元と同じ 154、142、285）。
  - 【計算】交点 q = t = −1 では、全ブロックで (i) と (ii) の値が一致する。矛盾の原因は N の台と |e| ≤ |c| である。
  - 【計算】N の台では、ν^{−8}〜ν^{−16} は解あり、ν^{−18}・ν^{−22} は見えない対が 1・3 個要り、ν^{−20} と ν^{−24} は一次方程式から解なし。
  - 【計算】ν^{−24} では、(i) と |e| ≤ |c| で e が全部決まり、(ii) の目標はちょうど H = −Σ_λ E(q, q) である（Σ_λ E(q, q) = ((1 − q²)/(1 + q²))Σ_λ N(q, q)）。(ii) はブロック全体の符号の反転を求め、その違いは因子 (1 − qt) のため交点で見えない。目標を −H にすると ν^{−24} だけは解けるが、ν^{−12}〜ν^{−22} は −H で解なし（ブロックごとの符号の付け替えでは直らない）。
  - 【計算】ν^{−20} は、λ = (2) の最も低い q の行に N の台の外の単項式を 1 個足すと一次方程式が解ける。
  - 【計算】見えない対を箱の台に許すと Σh の最小は 11（1 段・2 段広げると 10）で、置き場所は箱の縁に寄り、台で変わる。
  - 【計算】解のあるブロックでも、入れ替わる対の最小は ν^{−10} で 0 → 2、ν^{−12} で 4 → 61、ν^{−14} で 1 → [27, 131]、ν^{−16} で 4 → [8, 273]（[ ] は 300 秒で打ち切り）に増える。
  - rule_milp（McCormick）は ν^{−12} の対の最小化で 30 分以上かかったので止め、0/1 変数の定式化に替えた（3 秒）。
- 結論：
  - 【計算】(ii) を今の形で加えると、n = 2 の E は決まらない（4 ブロックで解なし、残りでも入れ替わる対の最小の原理と衝突する）。
  - 【予想】(ii) の n = 2 での形が破れている（続き 72・74 の予想を補強）。最もはっきりした証拠は ν^{−24} の符号である。
  - −H などの付け替えは診断にとどめ、採用しない。
  - 【未解決】n = 2 での (ii) の正しい形。
- 次の手順の案：§2′ の次の候補 (2′)（(ii) の n = 2 の形を理論から導く）。

## 2026-10-29（続き 82）n = 2 の (ii) の形を導く（ユーザーの選択「(ii) の形を導く」、続き 81 の次の候補 (2′)）
- スクリプト：rk_f_twisted_tensor.py の crosscap_split（F を c = 1, 2 に分ける。c = 1 は β ≡ n − 1 mod 2）、n2_quotient("c")（(ii′) の割った形の目標、data/cache/twisted_n2_quotient_c.pkl）、n1_crosscap、onehot_ranges（対が最小の解での e の範囲）、n2_iic（n2iic）。ログ：logs/refined_kauffman/f_twisted_tensor_n2iic.log（35 分、うち約 25 分は 300 秒の打ち切り）。テスト：test_crosscap_split、test_n1_crosscap_c2_vanishes_at_tq、test_n2_quotient_c_intersection。
- 導出（theory §12.18）：
  - 【導出】BFM の奇数 d の形と本研究の全 d の形から g = g^{odd} − ½ψ²[g^{odd}]。t = q で整数性から決まるのは E ≡ −(1 + q)ψ²[F] (mod 2(1 + q)) だけ（AgSc §4 の強い整数性）。
  - 【導出】等式 (ii) は「凍結した状態の二重被覆が同じ符号 ε で E に入る」という仮定による。一般には ε は状態ごと。
  - 【計算】n = 1 の F = g_(1) + 1 の c = 2 の部分（ν 奇、trefoil 4 項、T(2,5) 11 項）は t = q で 0。したがって n = 1 は ε の c 依存を見ない。(ii′) でも n = 1 の解（trefoil 一意、T(2,5) 自由 11/40）は完全に同じ。
  - 【仮説】(ii′)：E|_{t=q} = −(1 + q)ψ²[F^{c=1} − F^{c=2}]（ε = (−1)^c）。ν で書くと (−1)^n(1 + q)ψ²[F(−ν)]。
  - 探索の経過：まず F の c = 1, 2 の部分の相対符号の 2 通り（F_奇 − F_偶、F_偶 − F_奇）をブロックごとの一次方程式と見えない対の最小で調べた（−F は続き 81）。n = 2 で全ブロックが解けるのは F_奇 − F_偶（= F^{c=1} − F^{c=2}）だけで、F_偶 − F_奇 は ν^{−12} 以下が全部解なし。次に n = 1 で同じ 4 通りを調べ、F と F_偶 − F_奇（c = 2 の反転）が同じ解を与えることを見つけた。
- 結果（n = 2、trefoil、theory §12.18 の表）：
  - 【計算】交点 q = t = −1 では全ブロックで (i) と一致。目標が変わるのは ν^{−18}〜ν^{−24} だけ（ν^{−24} で H′ = −H）。
  - 【計算】N の台・|e| ≤ |c| で、9 ブロック中 8 ブロックが見えない対なしで解ける（(ii) では 5）。ν^{−22} は Σh ≥ 1（候補 4 単項式、必ず h > 0 のものは無い。(ii) では Σh ≥ 3）。
  - 【計算】ν^{−24} は (i)+(ii′) で対 0、72/72 が一意で、(i) だけの解と一致。
  - 【計算】入れ替わる対の最小（(i) だけ → (i)+(ii′)）：ν^{−10} 0 → 2、ν^{−12} 4 → 61、ν^{−14} 1 → [28, 131]、ν^{−16} 4 → [8, 273]、ν^{−18} 7 → [15, 209]、ν^{−20} 4 → [15, 53]。上の 4 ブロックは (ii) と同じ（目標が同じ）。
- 結論：
  - 【計算】c = 2 の符号を反転すると、§12.17 の矛盾のうち下の 4 ブロックのもの（ν^{−24} の全体の符号など）が、ν^{−22} の見えない対 1 個を除いて解消する。
  - 【仮説】(ii′) を n = 2 の (ii) の形とする。ε = (−1)^c は計算で選んだもので、物理的な導出（Ĩ の H¹ のゼロモードへの持ち上げの符号 (−1)^{b₁(N)}）は予想の段階。
  - 続き 72 の「同じ反転を n = 1 に当てると解なし、一様な規則ではない」を訂正：ν で書いた反転を全体の符号ごと当てていた。c で書けば一様。
  - 【未解決】上のブロックの対の数の増加（c = 2 とは無関係）と、ν^{−22} の見えない対 1 個。
- 次の手順の案：§2′ の次の候補 (2″)（上のブロックの対の増加の切り分け）。

## 2026-10-29（続き 83）上のブロックの対の増加を切り分ける（ユーザーの選択「対の増加を切り分ける」、続き 82 の次の候補 (2″)）
- スクリプト：rk_f_twisted_tensor.py の n2_split（n2split）、n2_hidden22（n2h22）。pairs_onehot に extra（台の外の単項式と |e| ≤ |c| + 2·extra を許し、見えない対を (|e| − |c|)/2 で数える）を加えた。
  - ログ：logs/refined_kauffman/f_twisted_tensor_n2split.log（83 分。うち大半は 600 秒の打ち切り）、f_twisted_tensor_n2h22.log（55 秒）。
  - テスト：test_n2_top_block_square_at_tq、test_pairs_onehot_extra。
  - 実行時の注意：`/usr/bin/time -l` はサンドボックス内で最大メモリの行を出さない（sysctl が拒否される）。終了の目印を "maximum resident" にした待機が終わらず、背景で約 17 時間待ち続けた。目印は "real" の行にする。
- 結果（theory §12.19 の表）：
  - 【計算】入れ替わる対の最小を、(ii′) だけ / (i) だけ / (i)+(ii′) / 箱の台で見えない対を許す、の順に並べる：
    - ν^{−10}：0 / 0 / 2 / 1
    - ν^{−12}：26 / 4 / 61 / 39
    - ν^{−14}：7 / 1 / [31, 115] / [16, 49]
    - ν^{−16}：0 / 4 / [8, 273] / [5, 253]
    - ν^{−18}：0 / 7 / [16, 209] / [11, 288]
    - ν^{−20}：2 / 4 / [16, 53] / [10, 37]
    - [ ] は 600 秒の打ち切りの [下界, 見つけた解]。
  - 【計算】どのブロックでも (i)+(ii′) は (ii′) だけと (i) だけの和より大きい。増加は 2 つの条件の両立から来る。台は部分的な原因にとどまる。
  - 【計算＋導出】ν^{−12} の構造：
    - t = q で N_λ = −X_λ²（X_(2) = (r² + 1)²(r⁴ + 1)²/r⁷、X_(11) = (r² + 1)²(r⁴ + 1)/r⁵）、H′ = Σ_λ ψ²[X_λ] である。
    - ψ² は冪を 2 倍するので、H′ は対角線の一つおきにしか値を持たない。これが (ii′) だけの対の出どころである。
    - λ ごとに E_λ(q, q) = ψ²[X_λ] を課すと (i) と両立し、対の最小は 94 である。入れ替えの描像での対は 120 + 28 = 148 個。
  - 【予想】Ĩ が二重被覆の 2 つの写しを入れ替える描像では、対が多いのは自然である。上のブロックで合わないのは、対の最小という原理のほうである。
  - 【計算】ν^{−22}：Σh = 1 の解で h > 0 になりうる単項式は 4 個で、すべて下端の角（下端から q の冪 0〜2、t の冪 1〜3、|c| ≤ 3）にある。そのうち λ = (2) の 2 個では e が c と逆の符号に限られる。
- 結論：
  - 【計算】対の増加は (i) と (ii′) を同時に満たすことによる。N の台は主因ではない。c = 1 の部分の形が誤っている可能性は、この計算では排除できていない（ν^{−12} の平方の構造は (ii′) の目標が自然な形であることを示唆する）。
  - 【予想】解を選ぶ原理を、対の最小から t = q の入れ替えの構造の持ち上げに変えるべきである。
  - 【未解決】平方の構造を ν^{−14} 以下と一般の t に広げること、ν^{−22} の見えない対の原因。
- 次の手順の案：§2′ の次の候補 (2″)（平方の構造の拡張と、それに基づく解の選び方）。

## 2026-10-29（続き 84）平方の構造を広げる（ユーザーの選択「平方の構造を広げる」、続き 83 の次の候補 (2″)）
- スクリプト：rk_f_twisted_tensor.py の n2_square（n2square）。ログ：logs/refined_kauffman/f_twisted_tensor_n2square.log（53 秒）。テスト：test_n2_top_block_square_at_tq に ν^{−14} の確認を追加、test_top_square_from_g。
- 結果（theory §12.20）：
  - 【計算】t = q では、ν^{−8}・ν^{−10} の N_λ は 0 である。ν^{−14}〜ν^{−24} の N_λ は、λ の和と差を含めて ± 平方ではなく、X_λ でも割り切れない。
  - 【計算】ν の級数全体の形式的な平方根は ν^{−8} から Laurent 多項式でなくなる（分母 (r² + 1)²）。
  - 【計算】一般の t の N_λ の ν^{−12}・ν^{−14} は既約で、積の構造は無い。
  - 【計算】X_11 = q^{−2}G_2 = −q G_11、X_2 = (q + q^{−1})X_11。G_λ は g_λ の t = q での ν^{−3} の係数である。
  - 【予想】ν^{−12} の平方は、最上の次数でだけ積に分かれる「brane ⊗ image」の部分で、合成表現の HOMFLY の最上の次数の分解と同じ型である。
- 結論：
  - 【計算】平方の構造は ν^{−12}・t = q に限られ、続き 83 の (2″) の案は ν^{−12} より下では使えない。
  - 【未解決】低いブロックで解を選ぶ原理。
- 次の手順の案：合成表現の skew 型の分解（n = 1 の状態の積の項）で、N の積の部分と残りを分けられるかを調べる。n = 1 と n = 2 の ν の次数と規格化の対応から始める。

## 2026-10-30（続き 85）skew 型の分解を試す（ユーザーの選択「skew 型の分解を試す」、続き 84 の次の候補）
- スクリプト：rk_f_twisted_tensor.py の n2_skew_parts、n2_skew（n2skew）。ログ：logs/refined_kauffman/f_twisted_tensor_n2skew.log（6 分）。テスト：test_n2_skew_parts。
- 前提：被覆のデータは log Z（連結な部分）なので、Koike の公式の積の項は Z の積として消えていて、そのままでは当てられない。代わりに、ν^{−12} の積の構造を全ブロックに延長して、N を積の部分と残りに分けた（theory §12.21）。
- 結果：
  - 【計算】H′ = ψ²[Φ] の Φ はどのブロックでも [3]_q で割り切れ、Y_(2) = (q + q^{−1})Φ/[3]_q、Y_(11) = Φ/[3]_q と λ ごとに分けられる（ν^{−12} では Y_λ = X_λ）。ただし、ν^{−13} 以下の Φ は q の整数冪と半整数冪が混ざる。
  - 【計算】積の部分 −Y_λ(ν)Y_λ(−ν)（image の側に符号 (−1)^c）を引いた残りの係数の絶対値の和：ν^{−12} で 0、ν^{−14} で 1392 → 144・456 → 80。ν^{−16} 以下では減らないか、かえって増える。
  - 【計算】(i) + λ ごとの E_λ(q, q) = ψ²[Y_λ] の対の最小：ν^{−8} 0、ν^{−10} 2、ν^{−12} 94、ν^{−14} [172, 257]。ν^{−16}〜ν^{−24} は解なし。
  - 【導出】ψ²[Y] ≡ Y² (mod 2) なので、λ ごとの目標が e ≡ c と両立するには、残り R_λ が対角線ごとに偶数であることが要る。R_λ は ν^{−16} で初めて奇数の係数を持ち、解なしのブロックと一致する。
- 結論：
  - 【計算＋導出】積の構造は上の 2 ブロック（ν^{−12}・ν^{−14}）までで、ν^{−16} で偶奇によって破れる。
  - 【未解決】ν^{−14} で対が増えること、低いブロックで解を選ぶ原理、Φ の半整数冪の意味。
- 次の手順の案：
  - (a) 低いブロックを brane と image をつなぐ状態として扱う。合成表現の相関関数 ⟨P_α P_β^∨⟩ = W_{[α;β]} + c₁W_adj + c₀ の c₁・c₀ の項（Koike の σ ≠ ∅ に当たる）から、どの状態がどのブロックに入るかを読む。
  - (b) 対の最小の原理と (ii′) の組み合わせ（続き 82）に戻り、T(2,5) の n = 2 で (ii′) を照合する。

## 2026-10-30（続き 86）c₁・c₀ の項を読む（ユーザーの選択「c₁・c₀ の項を読む」、続き 85 の次の候補 (a)）
- スクリプト：
  - rk_f_twisted22_lam.py の cover_parts（parts）：log Z の (2,2) 成分を comp・c1・c0・low の 4 項に分けて保存する（119 秒）。
  - rk_f_cover_fixed22.py の cover_Z に pair の引数を加えた。
  - rk_f_twisted_tensor.py の n2_parts_N・nu_series・n2_parts（n2parts、22 分。大半は 14 通りの部分和の約分）。
- ログ：logs/refined_kauffman/f_twisted22_parts.log、f_twisted_tensor_n2parts.log。データ：data/refined_kauffman/daha/largeN_cover_fixed22_lam_parts_m1.pkl。
- テスト：test_c1_c0_unrefined（既定）、test_cover_parts_sum（results の印、保存済みのデータの判定）。
- 結果（theory §12.22）：
  - 【計算】4 項の和は log Z の (2,2) 成分と N_λ に一致する。どの項も単独では Laurent でなく、真部分集合の和（14 通り）もすべて Laurent でない。
  - 【計算】一般の t で、ν^{−1} の級数の最高次は comp ν^{−10}、c1 ν^{−8}、c0・low ν^{−6}。N の最上部 ν^{−8} は打ち消しの結果で、合成表現の項が最上部を担うわけではない。
  - 【計算】ν を含む分母の因子（1 − q^a A/t^b 型）が項を鎖のようにつなぐ（comp–c1、comp–c0、c1–c0、c0–low）。
  - 【計算】t = q では c₁ = 1、c₀ = 1 または 0（unrefined の Koike）。hatted な各項は ν について Laurent になるが、q の分母（(1 − q)^k(1 + q)^k）が残る。N_λ のどのブロックでも、どの項の係数も Laurent でない。
- 結論：
  - 【計算】Koike の σ ごとの項は、単独では整数の状態の数にならない。
  - 【予想】Koike の分解は指標の分解で、BPS 状態の分解ではない。被覆の状態を積の部分とつなぐ状態に分ける原理は、相関関数の項からは読めない。
  - 【未解決】低いブロックで解を選ぶ原理。極を持つ部分を項の間で配り直す方法は、一意でなく場当たり的な調整になるので試していない。
- 次の手順の案：skew 型・積の構造の線（続き 83〜86）はここで区切る。続き 82 の (ii′) と対の最小の原理に戻り、T(2,5) の n = 2 で (ii′) を照合する（2 箱・4 箱のデータの作成から必要。実行前に見積もる）。

## 2026-10-30（続き 87）T(2,5) の n = 2 で (ii′) を照合するためのデータ作成の見積もり（ユーザーの選択「(2‴) 見積もりから」、続き 86 の次の候補 (2‴)）
目的：trefoil の n = 2 で (ii′)（ε = (−1)^c）が ν^{−18}〜ν^{−24} の矛盾を解消した（続き 82）ことを T(2,5) で照合する。計算に入る前に、要るデータと時間・メモリを小さい N で測って外挿する（CLAUDE.md の規則）。
- 要るデータ（trefoil の n = 2 の入力を m = 2 に置き換えたもの）：
  | データ | 役割 | m = 2 の状態 |
  |---|---|---|
  | SO の 2 箱 g（largeN_2box_m2.pkl）、随伴 Radj_A_m2.pkl | (ii′) の目標 F = g_(2) + g_(1,1) の c = 1, 2 の分離、被覆の (1;1) | 保存済み（続き 32） |
  | GL の composite (α;1)（\|α\| = 2、3 箱）Rcomp_A_m2_{α} | 被覆の log Z の (2;1) の項 | 無し。今回 N = 2〜16 を保存 |
  | GL の composite (α;β)（\|α\| = \|β\| = 2、4 箱）Rcomp22_A_m2_{α}_{β} | 被覆の固定セクター P̂_λ（E の未知数の台 N_λ） | 無し。今回 N = 4〜9 を保存 |
  | SO の 4 箱 R_U_m2_{λ}（\|λ\| = 4、D_N の DAHA）→ largeN_4box_so_m2 の f_(2,2), f_(4) | 条件 (i)（t = −1）の目標 R22・R4 | 無し。今回 N = 2〜6 を保存 |
- スクリプト：rk_c_composite22_run.py の run・window に m の引数を追加（既定 1）。rk_c_dn4_save.py と rk_c_fast_composite.py はもとから m を受ける。
- ログ：logs/refined_kauffman/c_dn4_save_m2.log、c_fast_composite_m2.log、c_composite22_run_m2.log、c_composite22_window_m2.log、c_fast_composite_m2_range.log、c_dn4_m2_range.log。
- U の冪の範囲【計算】（s の次数の最小・最大の N に対する傾き。どれも小さい N から一定）：
  | データ | m = 1 | m = 2 |
  |---|---|---|
  | composite (α;1)（x = A^{1/2} の冪） | −14〜−6（9 項） | −24〜−12（13 項） |
  | composite (2;2), (11;11)（U = t^N の冪、偶数冪だけ） | −10〜−4（7 項） | −18〜−8（11 項） |
  | composite (2;11) | −9〜−4（6 項） | −15〜−8（8 項） |
  | SO の 4 箱（U = t^N の冪、5 個の λ すべて） | −20〜−8（13 項） | −36〜−16（21 項） |
- 実測（1 コア）【計算】：
  | 計算 | m = 1 | m = 2 | 比（同じ N） | N ごとの伸び（m = 2） |
  |---|---|---|---|---|
  | composite (α;1)、N = 16、α 1 個 | 45 秒 | 164 秒 | 3.7 | 約 1.5 |
  | composite (2;2) 型 3 個、N = 9 | 22 秒・0.25 GB | 237 秒・0.36 GB | 11 | 時間 1.8〜1.9、メモリ 1.4〜1.5 |
  | SO の 4 箱（\|λ\| = 4 の 5 個）、N = 6 | 66 秒・0.43 GB | 905 秒・0.54 GB | 14 | 時間 2.4〜3.2（N = 3→4→5→6：45, 143, 382, 905 秒）、メモリ 1.35〜1.45 |
- 外挿【見込み】：
  - composite (α;1)：13 項なので N = 2〜18（ノード 6〜18、照合 2〜5）。残りの N = 17, 18 は合計 15 分以内。
  - composite (2;2) 型：11 項なので m = 1 と同じ組み方では N = 4〜16（ノード 6〜16、照合 4, 5。m = 1 のように上側の照合点を残すなら N = 17, 18 も要る）。N = 10 で 7 分、12 で 25 分、14 で 1.4 時間・2.3 GB、16 で 5 時間・5 GB、17 で 9 時間・7 GB、18 で 17 時間・10 GB。N ≤ 16 の合計は逐次で約 10 時間、N ごとに並列にしても律速の N = 16 で壁時計 5 時間（メモリの和が 8 GB を超えないよう N = 16 は単独で走らせる）。r の次数ごとに U の冪の範囲を絞る方法（続き 71、rk_c_bounded_interp）を composite に移せば N ≤ 12〜13 で済む見込み（合計 1〜1.5 時間。移植に 1〜2 時間）。
  - **SO の 4 箱（D_N）：届かない。** N = 6 の 905 秒から N ごとの伸び 2.4〜2.7 で外挿すると N = 8 で 1.5〜2 時間、9 で 4〜5 時間、10 で 9〜13 時間、11 で 1〜1.5 日、12 で 2〜4 日、13 で 5〜10 日、14 で 2〜4 週間。メモリは N = 12 で 3〜4 GB、14 で 6〜8 GB。21 項を素朴な補間で決めるには N = 4〜24 が要り、r の次数ごとに範囲を絞る方法でも（m = 1 では 13 項に N ≤ 11 が要った）N ≤ 16〜17 と見込まれ、月の単位になる。24 時間を超えるので実行しない（研究上のルール）。
- (ii′) だけの照合に識別力があるか【計算、trefoil の保存済みデータ】：n = 2 の各ブロックで、(ii)（続き 81）だけ・(ii′)（続き 82）だけの入れ替わる対の最小を対角線ごとの row_min で求めた。
  | ブロック | ν^{−8} | ν^{−10} | ν^{−12} | ν^{−14} | ν^{−16} | ν^{−18} | ν^{−20} | ν^{−22} | ν^{−24} |
  |---|---|---|---|---|---|---|---|---|---|
  | (ii) だけ | 0 | 0 | 26 | 7 | 0 | 0 | 0 | 0 | 0 |
  | (ii′) だけ | 0 | 0 | 26 | 7 | 0 | 0 | 2 | 0 | 0 |
  どちらも全ブロックで解があり、ほとんど同じである。続き 81・82 で (ii) と (ii′) を分けた ν^{−18}〜ν^{−24} の矛盾は、(i) と組み合わせて初めて現れる（ν^{−24} では (i) が e を一意に決め、(ii) はその符号違いを求めた）。
- 結論：
  - 【計算】T(2,5) の n = 2 で (ii′) を照合するには SO の 4 箱（条件 (i)）が不可欠で、GL の被覆の側（composite (α;1)、(2;2)）だけでは (ii) と (ii′) を区別できない。
  - 【見込み】SO の 4 箱の m = 2 は、現行の D_N の DAHA（約分なしの LP 版）では時間・メモリとも届かない。GL の側は 5〜10 時間（絞る方法なら 1〜2 時間）で作れる。
  - したがって (2‴) の T(2,5) の照合は、SO の 4 箱の計算法を変えない限り実行できない。候補は続き 70 の案 1（q を有限体の数値にし t の一変数で計算して U = t^N で補間し、後で t = −1 を代入する。実装 1〜2 日、速さの利得は未測定）。
- 停止点：見積もりの結果を報告し、次の方針の判断を求める。動いている計算は無い（D_N の m = 2 は N = 6 まで保存して止めた。composite (α;1) は N = 16 まで、(2;2) は N = 9 まで保存済み）。

## 2026-10-30（続き 88）T(2,7) の 2 箱のデータ作成と、n = 1 の符号の規則の照合（ユーザーの選択「候補 (4) へ」：T(2,5) で (ii) を要らなくする根拠、T(2,7) での照合。まず 2 箱のデータ作成の見積もりから）
- スクリプト：データ作成は既存の rk_c_fast_daha_dn.py（`3 2 20 1 2 11`）、rk_c_fast_adj.py（`3 2 26`）、rk_c_read_largeN.py（`3 1 2 11`）、rk_c_2box_largeN_m.py（`3`）。照合は rk_f_twisted_tensor.py に m の引数と新しい入口を加えた：rule・full・lower に m、blocks（N のブロックの表と生成関数の分子）、scan（k ≥ 2m のブロックごとに q / t を選ぶ全通り）、n1choice（ブロックと c の符号ごとに 0/q/t/qt と σ を選ぶ整数計画、n2_choice の n = 1 版）、n1upper（下は q の規則、上は一般の e で、対が最小の解の上のブロックの e と範囲）、n1freeset（一意でない単項式の値の組み合わせのうち実現できるもの）。rule_milp に fix（e の値の固定）を加えた。
- ログ：logs/refined_kauffman/c_fast_daha_dn_m3.log、c_fast_adj_m3.log、c_read_largeN_m3.log、c_2box_largeN_m3.log、f_twisted_tensor_rule_m3.log、f_twisted_tensor_lower_m3.log、f_twisted_tensor_full_m3.log、f_twisted_tensor_blocks.log、f_twisted_tensor_scan_m2.log、f_twisted_tensor_scan_m3.log、f_twisted_tensor_n1choice.log、f_twisted_tensor_n1choice_free_m3.log、f_twisted_tensor_n1upper.log、f_twisted_tensor_n1freeset_m3.log。
- データ：data/refined_kauffman/daha/Z_N{2..20}_m3_{1,2,11}.pkl（D_N の生の値）、Wadj_N{2..26}_m3.pkl、R_U_m3_{1,2,11}.pkl、Radj_A_m3.pkl、largeN_2box_m3.pkl（T(2,7) の 2 箱の f, f^cov, g）。
- 見積もりとデータ作成【計算】：
  - 続き 87 の見積もりと違い、2 箱（m = 3）は軽い。D_N の 2 箱は N = 2〜20 の全部で 447 秒（N = 20 で (2) が 64 秒、(1,1) が 12 秒、RSS 0.4 GB）、GL の随伴は N = 2〜26 で 363 秒。
  - U の冪の範囲（傾きから）：SO の [1] は U^{−13}〜U^{−6}（8 項、ノード N = 13〜20）、SO の (2)・(1,1) は U^{−26}〜U^{−12}（15 項、ノード N = 6〜20）、随伴は A^{1/2} = ν の冪 −26〜−12（15 項、ノード N = 12〜26）。
  - 補間の時間は sympy の約分が支配的で、[1] が約 1 分、(2) が約 2 時間 8 分、(1,1) が約 35 分。全体で約 3 時間。
  - 照合はすべて通った：補間に使わなかった N での一致（[1] は N = 2〜12、(2)・(1,1) は N = 2〜5、随伴は N = 2〜11）、g_(2) + g_(1,1) と (1+t)(g_(2) − g_(1,1)) が整数 Laurent、U = t², t³ で KN の T(2,7) の固定 N の g と一致。
- N のブロック（n_k = ν^{−(6+2k)} r^{−(4+2k)} の係数、blocks）【計算】：
  | k | T(2,7) | 状態の数 |
  |---|---|---|
  | 0〜5 | Sym^k V（V = 1 + t + qt） | 1, 3, 6, 10, 15, 21 |
  | 6（= 2m、移り目） | −t(…)、正の部分は **t³ Sym³ V**（10 状態）、負の部分 31 状態 | 41 |
  | 7（= 2m+1） | −t³(q⁴t⁴ + q³t⁴ + 2q³t³ + 3q²t³ + 3q²t² − qt⁴ + qt³ + 3qt² + 2qt − t⁴ − t³ + t + 1)、正の部分は t⁶V（3 状態） | 21 |
  | 8 | −t⁵(1 + qt)(Sym²V − qt) | 10 |
  | 9, 10, 11 | −t⁷ Sym²V、−t⁹ V、−t¹¹ | 6, 3, 1 |
  - 移り目の正の部分は T(2,5) で t³V、T(2,7) で t³ Sym³V（t³ Sym^{2m−3} V の形、T(2,3) では無い）。
  - 生成関数 Σ n_k x^k に (1 − x)(1 − tx)(1 − qtx) を掛けると、x^{2m} = x⁶ の係数は −(1 + t + qt + qt² + q²t² + q²t³ + q³t³)²（7 項は T(2,7) の簡約 superpolynomial の項数）、x⁹〜x^{13} には (q − t) の因子が出る（続き 75 の T(2,3)・T(2,5) と同じ型）。
- 符号の規則の照合【計算】：
  - lower（下の系列に V の対合、上は一般）：下のブロック k = 0〜5 は Sym^k V に一致。解があるのは diag(1, 1, −1) だけで、σ_k は 6 個とも一意、上のブロックの一意でない変数は 56 中 27。
  - rule（下の系列と移り目の正の部分に q の規則、残りは一般）：(i) だけでは σ の群 7 中 1、一般の変数 46 中 28 が一意でない。**(i)+(ii) でも一般の変数 46 中 7 が一意でない**（σ は 7 個とも一意）。7 個はすべて k = 7 のブロック（q⁰t⁶、q⁰t⁷、q¹t⁵、q¹t⁶、q¹t⁷、q²t⁵、q²t⁶）。T(2,5) では (i)+(ii) で一意だった。
  - full（続き 75 の規則を全ブロックに：k ≤ 2m で (−1)^q、k > 2m で (−1)^t、σ はブロックと c の符号ごと）：**(i) だけでも (i)+(ii) でも解なし**。
  - scan（k ≥ 2m のブロックごとに q / t を選ぶ）：T(2,7) は 64 通りすべて解なし（(i)、(i)+(ii) とも）。T(2,5) は 16 通り中 2 通り（続き 75 の規則。k = 7 は 1 状態で q と t が区別できない）。
  - n1choice（ブロックと c の符号ごとに 0/q/t/qt と σ を選ぶ）：(i) だけなら解あり（対の最小 2。k = 6− は 0、k = 7− は q、k = 7+ は t、k = 8〜10 は t）。**(i)+(ii) では解なし**。下の系列も選ばせても同じ。
  - n1upper（下は q の規則、上は一般の e）、(i)+(ii)：対の最小 15（T(2,5) は 4）、一意でないのは上の 7 個だけ。決まった e は次の構造をもつ。
    - k ≥ 8：すべて sign(e/c)(−1)^t = +1（t の規則、σ = +。対は k = 8 の q¹t⁷ の 1 個だけ）。
    - k = 6 の負の部分：e ≠ 0 の単項式はすべて sign(e/c)(−1)^q = −1 で |e| ≤ 1（q の規則、σ = −。|c| = 2, 3, 5 の単項式が e = 0, ±1 になり、対 10 個）。
    - k = 7：符号が混在し、7 個が一意でない。
  - n1freeset：7 個の値の組み合わせ 128 通りのうち、(i)+(ii) と両立するのは 3 通り（対の最小 15, 16, 16）。どれも k = 7 の正負の部分で σ(−1)^{0,q,t,q+t} の形にならない（例：対 15 の組は正の部分 (q⁰t⁶, q⁰t⁷, q¹t⁷) で e = (−1, −1, −1)、sign(e/c)(−1)^t = (−1, +1, +1)）。
- 結論：
  - 【計算】続き 75 の規則は T(2,7) で破れる。(i) だけでも解が無く、ブロックと符号ごとの対角な対合（0/q/t/qt）に広げても (i)+(ii) で解が無い。破れるのは上の系列の最初のブロック k = 2m+1 で、下の系列・移り目（q）と k ≥ 2m+2（t）の規則は T(2,7) でも保たれる。
  - 【計算】T(2,7) の n = 1 の E は、下の規則と (i)+(ii) では一意に決まらない（k = 7 の 7 変数、実現できる組 3 通り）。
  - 候補 (4) の目的「T(2,5) で (ii) を要らなくする根拠」は得られなかった。逆に、規則は m を上げると n = 1 でも破れるので、T(2,3)・T(2,5) で規則が合ったのは状態の数が少ないためと考えるべきである。n = 2 の trefoil で規則が合わなかった（続き 76・78）のと同じ型の破れが、n = 1 の T(2,7) で現れた。
  - 【予想】k = 2m+1 のブロックでは Ĩ は対角でない（弦をまたぐ置換を含む）。このブロックは T(2,5) で −t³(1 + qt)V（続き 75 で (1 + qt)⊗(1 + qt) の非対角の対を見た）、T(2,7) で 21 状態、正の部分 t⁶V を含む。
  - 【未解決】k = 2m+1 のブロックの Ĩ の作用、T(2,7) の実現できる 3 組のどれが正しいか（n = 2 の T(2,7) の条件は重く、計画に無い）。
- 停止点：候補 (4) の結果を報告し、次の方針の判断を求める。動いている計算は無い。

## 2026-10-30（続き 89）n = 1 の E を状態の空間の同定から作る（ユーザーの選択「候補 (5) へ」：k = 2m+1 のブロックの構造）
- スクリプト：rk_f_twisted_tensor.py の principle（principle。`principle 1 2 3`、`--save` で twisted_E_n1_m{m}.pkl に保存）。ログ：logs/refined_kauffman/f_twisted_tensor_principle.log、f_twisted_tensor_principle_save_m3.log。データ：data/refined_kauffman/daha/twisted_E_n1_m3.pkl（T(2,7) の E の予言）。テスト：tests/refined_kauffman/test_rk_twisted_tensor.py::test_principle_reproduces_unique_E（T(2,3)・T(2,5)、約 20 秒）。理論：theory §12.23。
- 経過：
  - 続き 88 のブロックの表から、上の系列を上から数えて（j = K − k、K = 4m−1）−t^{K−2j} で割ると Sym^j V になり、差が (q − t) で割り切れることに気づいた（T(2,5) の k = 5：Sym²V + t(q − t)、T(2,7) の k = 8：Sym³V + t²(q − t)）。移り目は Sym^{2m}V − S²（生成関数の x^{2m} の係数）。
  - 続き 75 の T(2,5) の一意な E をこの同定で読み直すと、移り目は「Sym⁴V に q の規則、S ⊗ S に入れ替え（トレース ψ²S）」、上は「Sym^j V に双対の対合、補正の対はもとのセルの符号」で全単項式が説明できた（移り目の対 3 個は S ⊗ S の非対角の対、k = 5 の対 1 個は補正の対）。
  - 同じ構成を T(2,7) に当てると、k ≠ 7 はすべて続き 88 の一意な値に一致し、k = 7（j = 4）だけが補正 R_4 = t(1 + 2t + 2qt + 2t² + 2qt² + q²t²) の読み方で変わった。セルを個別にずらす素朴な読み方では |e| ≤ |c| に反する単項式が出て (i) を満たさない。R_4 = t[(1 + t + qt)² + t²] = t(S′ ⊗ S′ ⊕ t²)（S′ = T(2,3) の簡約 superpolynomial）と読み、S′ ⊗ S′ に入れ替えを課すと、続き 88 の実現できる 3 組のうち対 16 の組になり、(i)・(ii) をともに満たした。
- 結果【計算】（m = 1, 2, 3）：
  | | 台の外 | |e| ≤ |c|・偶奇の違反 | 対 | (i) | (ii) | 保存した一意な E と一致 |
  |---|---|---|---|---|---|---|
  | T(2,3) | 0 | 0 | 0 | 満たす | 満たす | True |
  | T(2,5) | 0 | 0 | 4 | 満たす | 満たす | True |
  | T(2,7) | 0 | 0 | 16 | 満たす | 満たす | （予言。保存した） |
  - 補正：T(2,5) の j = 2 は R = t、R′ = −t。T(2,7) の j = 3 は R = t²、R′ = t²、j = 4 は R = t(S′² + t²)、R′ = −t(ψ²S′ + t²) = −t − 2t³ − q²t³。
- 結論：
  - 【計算】n = 1 の E は、状態の空間の同定（下：Sym^k V、移り目：Sym^{2m}V − S ⊗ S、上：双対の Sym^j V ＋ (q − t) の補正）と空間ごとの Ĩ（V の対合の誘導作用、テンソル積の入れ替え）から、条件 (i)・(ii) を入力にせずに作れる。T(2,3)・T(2,5) で一意な E を再現し、T(2,7) で (i)・(ii) を予言として満たす。
  - 【計算】候補 (4) の問い「T(2,5) で (ii) を要らなくする根拠」の答え：符号の規則（ブロックごとの σ(−1)^d）が中間的な記述だったためで、状態の空間で書けば (i) も (ii) も要らない。
  - 【予想】R_{m+1} = t^{m−2}(S′ ⊗ S′ + t²) の形と Ĩ の作用は T(2,7) の 1 例から読んだ仮説で、T(2,9) で検証が要る（j = m+2 の補正の構造も未同定）。
  - 【未解決】補正 (q − t)R_j の由来（下と上の双対性の破れ）、n = 2 への適用。
- 停止点：候補 (5) の結果を報告し、次の方針（§2′ の候補 (6)〜(8)）の判断を求める。動いている計算は無い。

## 2026-10-30（続き 90）n = 2 の N_λ に n = 1 の同定を当てる（ユーザーの選択「候補 (6) n = 2 へ」）
- スクリプト：rk_f_twisted_tensor.py の n2_fock（n2fock）。ログ：logs/refined_kauffman/f_twisted_tensor_n2fock.log。理論：theory §12.24。
- 方法：trefoil の n = 2 の N_λ（λ = (2), (1,1)、続き 76〜77 のデータ）のブロックを、上の端から r^5 s^{−3} で規格化して B_k（k = 0..8、ν^{−8−2k}）とし、n = 1 の下の系列と同じく「B_1/B_0 を 1 段のセル、B_2 を PE の 2 次」とみなして残りを測った。λ の組合せ（N_2 + N_11、N_2 − (q + 1/q)N_11、N_2 − q²N_11）、符号付きの PE、Sym²V·Sym²V′、(VV′)² も試した。別に、各ブロックの qt = −1 と t = −1 での値を求めた。
- 結果【計算】：
  - B_0 = (t − q)/q、B_1 = B_0·VV′/q²。最初の 1 段は V ⊗ V′（brane ⊗ image）の積になる。
  - B_2（300 状態）から Fock 型の項を引いた残りの状態数：B_0·PE_2(B_1/B_0) で 320、Sym²V·Sym²V′ で 316、(VV′)² で 340、符号付きの PE で 264。λ の組合せでも 324〜386。どれも B_2 と同程度で、ν 方向の Fock 系列は 2 段目で破れる。残りの t = q での値はどれも −(q+1)⁴(q²+1)⁴/q⁶ = −X_(2)²（続き 83）。
  - qt = −1 での各ブロックの状態数：N_(2) は 2, 2, 8, 4, 43, 21, 25, 11, 0（ν^{−8}〜ν^{−24}、もとは 2, 14, 300, 1686, 3313, 3099, 1901, 611, 42）、N_(11) は 2, 2, 15, 15, 13, 11, 0。ν^{−22} の残差は両 λ で同じ。t = −1 の値は 2〜61 状態で残る。
- 結論：
  - 【計算】n = 1 の ν 方向の同定（Sym^k V の系列）は n = 2 には直接当てはまらない。続き 77 の「ν 方向の系列は無い」を Fock 型の検査で確認した。
  - 【計算】n = 2 の状態の大部分は (1 + qt) を因子にもち、ブロックの中の塔（続き 77）が主である。
  - 【予想】n = 2 の空間は塔のセル（q, t, qt）と (1 + qt) の 2 状態の因子で書くべきで、入れ替えの構造は ν^{−12} の −X_λ² と (1 + qt) ⊗ (1 + qt) の非対角の型で入る。
  - 【未解決】塔の打ち切りの規則と (1 + qt) の商の同定。式からの導出（続き 86 の有理因子）が要る。
- 停止点：候補 (6) の結果（否定的）を報告し、次の方針の判断を求める。動いている計算は無い。

## 2026-10-30（続き 91）N を式から読む：unknot の参照、κ の Fock 系列、brane × image の積（ユーザーの選択「(6b) 式から導く」）
- スクリプト：rk_f_twisted_tensor.py の n1_terms（n1terms、m を指定）、新設の rk_f_n2_univ.py（univ：枠付き unknot の N^{univ}_λ の分母と形式級数、w1：枠なし unknot の確認）。ログ：logs/refined_kauffman/f_twisted_tensor_n1terms.log、f_n2_univ.log、f_n2_univ_w1.log。データ：data/refined_kauffman/daha/n2_N_univ_m1.pkl（枠付き unknot の N^{univ}_λ の分子・分母）。理論：theory §12.25。
- 方法：
  - n = 1：N = (1 − t)[R_adj K_adj + κ − W₁²]/D を 3 項（T1 = R_adj K_adj、T2 = κ、T3 = −W₁²）に分け、ν = 1/y の形式級数を多項式の割り算で求めてブロックごとに見た（sympy の series は 10 分で終わらず、割り算なら数秒）。
  - n = 2：rk_f_cover_fixed22 の pair_values で結び目のデータ（Rcomp22・Radj・Rcomp の読み込みと Rdata_nu）を 1 に置き換え、同じ hatted の手順で N^{univ}_λ を作った。枠 w = (t/q)^{1/2} はそのまま（枠付き unknot）と w = 1（枠なし）の 2 通り。
- 結果【計算】：
  - n = 1（m = 1, 2）：T2 = κ の項のブロックは k = 0〜9 すべてで Sym^k V（結び目によらない）。T1、T3 は単独では (q − 1)(t − 1) の極をもつ。T1 + T3 のブロックは k < 2m で 0、k = 2m で −S²、k > 2m で −Sym^k V − …。
  - n = 1 の unknot：w = 1 なら N = 0。枠付き（w = (t/q)^{1/2}）なら N = −(ν² − 1)(q − t)/(ν²q(t − ν²)(q − 1)(t − 1)(ν²q − 1)) で Laurent でない。
  - n = 2 の unknot（w = 1）：x₁₁ = 0、N_(2) = N_(11) = 0。
  - n = 2 の枠付き unknot：N^{univ}_λ は 0 でなく Laurent でない。分母 ν⁴ r (ν² − t)(q − 1)³(q + 1)³(q² + 1)²(t − 1)³(t + 1)³(t² + 1)²(ν²q − 1)(ν²q − t²)(ν²q² − t)（両 λ で同じ）、分子は 260 項・Σ|c| = 1234（λ = (2)）、1244（λ = (1,1)）で (√q − √t) を因子にもつ。ν^{−2} の形式級数は ν^{−2} から始まり、各ブロックの分母は (q∓1)³(q²+1)²(t∓1)³(t²+1)²、分子は (r − s) を因子にもつ。
  - 最初の patch の順序の誤り（関数を複製した名前空間に patch が効かず trefoil の R が混ざった）を直して再実行した。
- 結論：
  - 【計算】枠なしの unknot では n = 1 でも n = 2 でも N = 0。N は結び目のデータ R − 1 と枠 w の寄与だけからなる。
  - 【導出＋計算】n = 1 の §12.23 の構造は式の項に対応する：Sym^k V ← κ、−S ⊗ S ← W₁W₁′（積、入れ替え）、上の系列と (q − t) ← (R_adj − 1)K_adj と枠。
  - 【導出】n = 2 の塔（続き 77）は枠付き unknot の参照の q, t の分母（2 箱の Macdonald の次元・q_n²・κ）に由来し、R − 1 の項が打ち切る。ν のセルは {t, q⁻¹, t²q⁻¹, tq⁻²}。
  - 【予想】n = 2 の Ĩ は、参照の塔のセルには対角に、R − 1 の積の項には入れ替えとして働く。
  - 【未解決】R − 1 の項ごとの分解と打ち切りの仕組み、(i) での照合（候補 (6c)）。
- 停止点：(6b) の結果を報告し、次の方針の判断を求める。動いている計算は無い（ログの再生成だけが背景で動く。数分）。

## 2026-10-30（続き 92）n = 2 の N_λ を結び目のデータ R − 1 の項ごとに分ける（ユーザーの選択「(6c) R − 1 の項へ」）
- 方法（theory §12.26、scripts/refined_kauffman/rk_f_n2_rterms.py）：結び目のデータを 4 群 W（R_1, R_2, R_11）・adj（R_adj）・c（[α;1] の R_c）・b（[α;β] の R_b）に分け、群の部分集合 S ごとに「S だけ結び目の値、他は基準値」の N_λ(S) を計算し、包除で T(S) に分ける。基準は R = 1（unknot）と R = 0（Wilson loop なしの背景）。別に ε の次数ごとにも分けた。ν^{−2} の形式級数は flint の有理関数の多項式の割り算で取る（series_u。sympy の級数は 1 本 70 秒だったが数秒）。sympy への変換は Poly.from_dict で行う（to_sympy の単項式の和は数千項で 10 分以上かかった。最初の実行を途中で止めて直した）。
- 結果（logs/refined_kauffman/f_n2_rterms_{terms_b1,terms_b0,analyze_b1,analyze_b0,background,fock,hyptop,hyptop2}.log）：
  - 基準 R = 1：0 でない項は ∅, W, adj, c, b, W+adj, W+c。どれも（ブロックごとにも）Laurent でなく、塔の分母をもつ。N(S) は S が全部のときだけ Laurent。ε の次数ごとの N_d もすべて Laurent でない。
  - 基準 R = 0：0 でない項は ∅, W, adj, b, W+adj, W+c（T(c) = 0）。背景 T(∅) は閉じた形 −r⁷s³(q − t)ν⁶/[(A − t)(Aq − 1)(Aq − t)(Aq − t²)(Aq² − t)(A²q² − t²)]（λ = (2)。λ = (1,1) は分子 −r³s⁷(q − t)ν²）をもち、ブロックは Laurent で単項式の数 k(k + 1)、係数は三角数：真空 (q − t) と 6 セル（ν^{−2} あたり t, q⁻¹, tq⁻¹, t²q⁻¹, tq⁻²、ν^{−4} あたり t²q⁻²）の Fock 空間。n = 1 の κ の項 Sym^k V の n = 2 版。塔の分母は Wilson loop を含む項だけにある。
  - 結び目の部分 K = N − T(∅)：λ = (2) は ν^{−10} から K₁₀ = −r^{−11}s³(q − t)(1 + qt)²（背景の 5 セルと合わせて (q − t)VV′）、λ = (1,1) は ν^{−12} から −r^{−13}s³(1 + qt)²X₈、X₈ = (1 + q)²(1 + t)² + t(q − t)。λ = (2) の K₁₂ は −(1 + qt) × 23 項。k ≥ 26 で K = −T(∅)。結び目の 5 項のどの真部分集合の和もブロックごとに Laurent でない。
  - 上端の照合：λ = (2) の ν^{−8}, ν^{−10} で、仮説の族（真空の符号 ε_q, ε_t、背景のセルの対角な符号 ε₁..ε₅、積 (1 + qt)² に入れ替え σ(1 + q²t²) または対角 σ(1 + ε₆qt)²）768 通りのうち、(i) の行（24 本）を満たすのは 2 通りで、(ii) の行（合わせて 36 本）も満たす。二つは同じ E（E₈ = r^{−7}s³(q − t)、E₁₀ = r^{−11}s³(q − t)(q²t + qt² − qt − q − t − 1 − q²t²)）。入れ替えの方のセルの符号 (q⁻¹, tq⁻¹, t) = (−, −, +) は n = 1 の V の対合と同じ。
  - f_n2_rterms_background.log の「(1+qt) | True」の列は負の冪の Poly で判定した誤りで、背景のブロックは (1 + qt) で割り切れない（fock のログで正しく判定）。
- 結論：
  - 【計算】n = 2 の N_λ = Wilson loop なしの背景（閉じた形の 6 セルの Fock 空間、塔なし）＋ 結び目の部分（塔の分母は結び目の項にだけある）。
  - 【計算】結び目の部分は n = 1 の移り目 −S⊗S の位置（λ = (2) の ν^{−10}）で −(q − t)(1 + qt)⊗(1 + qt) として始まる。
  - 【計算＋予想】Ĩ の仮説「背景のセルに対角（n = 1 と同じ符号）、積に入れ替え」は上端 2 ブロックで (i)+(ii) を満たし E を一意に決める。n = 1 の原理（続き 89）が n = 2 の上端まで延びる最初の証拠。
  - 【未解決】結び目の部分の ν^{−12} 以降の分解（積と補正）、仮説の ν^{−12}（両 λ が結合）への延長、下の端の双対系列と背景の関係。
- テスト：tests/refined_kauffman/test_rk_n2_rterms.py（背景の閉じた形とブロックの Laurent 性・k(k+1)、T(c) = 0、上端 2 ブロックの K₁₀ と VV′）を追加。refined_kauffman のテストは 26 件通る（約 2 分）。
- 停止点：(6c) の結果を報告し、次の方針の判断を求める。動いている計算は無い。

## 2026-10-31（続き 93）Ĩ の仮説を ν^{−12} まで延ばす（ユーザーの選択「(6d) ν^-12 へ延ばす」）
- 方法（theory §12.27、rk_f_n2_rterms.py の hyp12・hyp12ii・hyp12b・hyp12bii・hyp12c・hyp12ci）：
  - ν^{−12} の結び目の部分 K₁₂ を、t = q の値 −X_λ² を再現する平方 Y⊗Y と (q − t) の補正に分ける（λ = (2) は候補の Y 18 通りを全部試す）。
  - ν^{−8}, ν^{−10} の e を続き 92 の上端の仮説に固定し、ν^{−8}..ν^{−12} の単項式だけを含む (i)・(i)+(ii) の行で ν^{−12} の e を整数計画で絞る。λ = (1,1) は仮説の族（Y の因子と補正に入れ替え/対角、符号、真空）187 通りを当て、λ = (2) は背景の誘導符号＋結び目の部分の平方（因子ごとに入れ替え/対角）＋補正（状態ごとに符号）で範囲を求める。
- 結果（logs/refined_kauffman/f_n2_rterms_hyp12.log・hyp12ii・hyp12b・hyp12bii・hyp12c・hyp12ci）：
  - λ = (1,1)：K₁₂ = −r^{−13}s³[Y⊗Y + (q − t)t(1 + qt)⊗(1 + qt)]、Y = (1 + q)(1 + t)(1 + qt)。族 187 通りのうち λ = (2) に解を残すのは 1 通り（全因子に入れ替え、σ = −1、真空 (−,−)）で、(i) だけで決まる。
  - λ = (2)：K₁₂ = −[Y⊗Y + (q − t)R] はどの Y でも書け、補正が最小なのは Y₂ = (1 + qt)Y（R は 17 状態、Σ|R| 82）。平方の Ĩ はどの Y でも「全因子に入れ替え、σ_Y = +1」の 1 通り（162 通り中）。
  - 補正 R の e は (i) と (i)+(ii) で範囲が完全に同じ：(ii) は (q − t) の対を見ない。Y₂ で一意なのは 17 中 3、対の最小（20）でも 12。構造を入れない整数計画では (i)+(ii) でも 49 中 4 個しか決まらない。
  - 下の端：ν^{−24} = −(1 + qt)h_K(W) を確認。ν^{−22} は t = q で h_{K−1}(W) で割り切れず、W の Fock の双対系列ではない。
- 結論：
  - 【計算】ν^{−12} の結び目の部分は両 λ で −Y⊗Y − (q − t) × 補正（n = 1 の移り目と上の系列の補正の型）。
  - 【計算＋予想】Ĩ の仮説（背景のセルに対角、平方に全因子の入れ替え）は ν^{−12} まで延びる。λ = (1,1) は一意、λ = (2) は Y によらず入れ替え。
  - 【導出】(q − t) の補正の符号は条件 (ii) に見えないので、条件からは決まらない。n = 1 の「対の符号はもとのセル t·m の ε」に当たる原理が n = 2 で要る。
  - 【未解決】λ = (2) の Y の選択、補正の符号の原理、ν^{−14} 以下への延長、下の端の双対系列の同定。
- 途中の失敗：hyp12c の最初の版は N に無い単項式を KeyError で飛ばしていて全候補が「台の外」になった。N に無い単項式では寄与の和が 0、ある単項式では |e| ≤ |c| という制約に直した。sympy の Poly で r², s² を生成元にすると失敗するので、Q = r², T = s² に写す _qt_poly を入れた。
- テスト：test_rk_n2_rterms.py に test_nu12_square_structure（両 λ の ν^{−12} の形）と test_nu12_hypothesis_rows（λ = (1,1) の E で λ = (2) に解があり、符号を変えた E では解が無い）を追加。set_groups が差し替える結び目のデータを各テストの後に戻す fixture も入れた（戻さないと後の test_principle_reproduces_unique_E が落ちる）。refined_kauffman の全テスト（results を除く 116 件）は、この修正の前の全体実行で 115 件が通り、落ちた 1 件は修正後に順序どおりの再実行で通った（約 8 分）。
- 停止点：(6d) の結果を報告し、次の方針の判断を求める。動いている計算は無い。

## 2026-10-31（続き 94）(q − t) の補正の符号の原理を探す（ユーザーの選択「(6e) 補正の符号の原理」。途中からユーザーの指示「これから 24 時間は自己判断で研究せよ」）
- 方法（theory §12.28、rk_f_n2_rterms.py の hyp12d・hyp12e・hyp12f）：λ = (2) の ν^{−12} の補正 R（Y₂ を固定、17 状態）の Ĩ の符号 e_m について、
  (a) n = 1 の読み「状態 t·m を q·m に移す」の多重度の検査（Y の 18 通り）、(b) 固定した符号の候補（一様、偶奇、主空間の状態の符号）の (i) の行と台での照合、
  (c) ν^{−14} 以下の e を自由にして深い段の行を加えたときの R の e の範囲、(d) 結び目の部分・補正を符号付きの平方の和 Σ z m Z⊗Z に分け各平方に入れ替えとする原理（_square_fit。Σ|z| の最小を制約なし／(i) を課して比較。基底 4 通り、背景の符号 4 通り）。
- 結果（logs/refined_kauffman/f_n2_rterms_hyp12d.log・hyp12e・hyp12f）：
  - (a) どの Y でも移せない m がある（Y₂ で 17 中 8）。入れ替えの下では状態の ε が定義できない。
  - (b) 候補はすべて (i) の行を 14〜18/44 本破る。
  - (c) (i) の行を ν^{−24} まで加えても一意は 3/17。(ii) を加えると ν^{−18} 以下で解なし（続き 81 の (ii) の破れ）。
  - (d) λ = (1,1)：基底 {1+q, 1+t, 1+qt} で最小 3 項、入れ替えの E は答え E₁₂(1,1) に全 49 単項式で一致（V, V′ の基底では 4 項で不一致）。λ = (2)：K₁₂ 全体の最小は 22（B1）・17（B2）・13（B3）、(i) を課すと 26・26・19 で最小の分解は両立しない。補正 R でも同じ（13 → 17 など）。(1 + qt) の因子に対角な Ĩ は 4 基底すべてで解なし。背景の符号 4 通りでも 22 → 26〜27。
- 結論：
  - 【導出＋計算】補正の符号は条件からは決まらない。
  - 【計算】n = 1 の読みは n = 2 に移せない。λ = (1,1) の答えは「(1 + q), (1 + t), (1 + qt) の積の平方の符号付きで最小の和に入れ替え」と言い換えられるが、この原理は λ = (2) の ν^{−12} で (i) と両立しない。
  - 【未解決】λ = (2) の補正の Ĩ。この線はいったん区切る。
- 自己判断で始めた次の候補 (7)：T(2,9) の 2 箱のデータ作成。
  - 見積もり：D_N の 2 箱は m = 4 で同じ N の m = 3 の 2.5〜3 倍（N = 12 で 33 秒、N = 20 で 139 秒、N = 26 で 387 秒、最大 RSS 0.57 GB）。U の冪の範囲は m = 1〜3 の規則（2 箱は幅 8m）から U^{−34}〜U^{−16}（19 項）と見込み、実測も同じ。補間は rk_c_read_largeN（sympy の約分と Newton 補間、m = 3 で 3 時間）では日の単位になる見込みなので、rk_c_bounded_interp（続き 71、flint の約分）を 2 箱に使う rk_c_bounded_save.py を書いた。
  - 検証：m = 3 で保存済みの R_U_m3_{1,2,11} と全係数一致（[1] は N = 2〜4、(2) は N = 3〜8、(1,1) は N = 4〜10 で一意。小さい N を含めると帯の仮定が N = 3, 4 で破れるので使う N の下限も動かす）。
  - 作成：R_U_m4_{1,2,11}.pkl（[1] は U^{−17}〜U^{−8} の 10 項を N = 2〜4 から、(2)・(1,1) は 19 項を N = 3〜8・4〜10 から一意に決め、使わない N ≤ 26 の全部で照合。約分 1264・923 秒）、Wadj_N{2..34}_m4.pkl（N = 34 で 409 秒）、Radj_A_m4.pkl（ν の冪 −34〜−16、ノード N_A = 16〜34、N_A = 2〜15 で照合）、largeN_2box_m4.pkl（g_(2) + g_(1,1) と (1+t)(g_(2) − g_(1,1)) は整数 Laurent）。T2_9_fixedN_3box_g.pkl は無いので rk_c_2box_largeN_m.py の U = t², t³ の照合は飛ばすようにした。
  - ログ：logs/refined_kauffman/c_fast_daha_dn_m4.log、c_fast_adj_m4.log、c_bounded_save_m3.log、c_bounded_save_m4.log、c_2box_largeN_m4.log。
- テスト：test_rk_n2_rterms.py に test_correction_sign_principles（λ = (1,1) の最小の分解が答えを再現、一様な符号の補正は (i) を破る）を追加。rk_f_twisted_tensor.principle に j ≥ m+2 の補正を N から読む枝（extra_R）を加え、m = 1, 2 のテストと T(2,7) の再現が変わらないことを確認。
- 途中の失敗：hyp12e の最初の版は背景の e と ψ²Y₂ の e を辞書の結合で重ねて上書きしていた（足すように直した）。λ = (1,1) の検査を (i) の行（λ = (2) の e の行）で行っていたので、答えとの比較に直した。hyp12f は 1 通り 30 分以上かかるので 32 通りから 4 通りに絞った。

## 2026-10-31（続き 95）T(2,9) の 2 箱で n = 1 の原理を検証する（自己判断で選んだ候補 (7)）
- 方法（theory §12.29）：続き 94 で作った largeN_2box_m4.pkl から N のブロック（blocks 4）を出し、principle 4（j ≥ m+2 の補正を N から読み、そのブロックの e を自由にして (i)・(i)+(ii) の解の有無と範囲）、principleR 4（Sym^j V の部分は原理で固定、補正 R_j の対の符号だけ自由。仮説 4 つの照合）、n1upper 4（下は q の規則、上は一般の e）を実行した。
- 結果（logs/refined_kauffman/f_twisted_tensor_{blocks,principle,principleR,n1upper}_m4.log）：
  - **R_{m+1} の予想が破れた**：R₅ = t²(V² + t²)、V = 1 + t + qt。続き 89 の S′ = T(2,2m−3) の superpolynomial（T(2,5) の 5 項）ではない。m ≤ 3 では S′ = V で区別がつかなかった。principle のコードを V 固定に直した（m = 1, 2 のテストと T(2,7) の再現は変わらない）。
  - j ≤ m+1 の原理は T(2,9) でも (i)+(ii) と整合する（k = 9 を自由にして解あり）。
  - j = m+2 の補正 R₆ = t·P は 35 状態（R_m, R_{m+1}, R_{m+2} で 1, 10, 35 = C(3,3), C(5,3), C(7,3)）。P は V³ を含まず、{1, V, V², V³} の非負結合で最小 7 項。Ĩ は (i) で 15 中 2 個しか決まらず、(i)+(ii) でも同じ（(ii) は対を見ない）。一様 ±、入れ替えの仮説は (i) を破る。
- 結論：【計算】原理は j ≤ m+1 で T(2,9) まで正しい（R_{m+1} は V 固定に訂正）。【計算】j = m+2 の補正の Ĩ は条件に決まらない。n = 2 の構図と同じ。【予想】補正の状態の数 C(2j′+3, 3)。T(2,11) で R_{m+3} の形を同定する。
- 途中の失敗：principleR の最初の版は対の単項式が N の台の外に出る（Sym の部分と打ち消す）ケースで assert に当たった。台の外では寄与の和が 0 という制約に直した（続き 93 の hyp12c と同じ）。hyp12e の基底 B4 の λ = (2) は 2 時間を超えたので打ち切った（結果はログに無い）。
- 次：(7b) T(2,11)（m = 5）の 2 箱のデータ。測定中（D_N N ≤ 12、随伴 N ≤ 20）。見積もりは次の続きに書く。

## 2026-10-31（続き 96）T(2,9) の R_{m+2} の Ĩ の候補と、T(2,11) のデータ作成（自己判断）
- R_{m+2} の Ĩ（theory §12.29 の追記、principleR 4 --neg / --pos、logs/refined_kauffman/f_twisted_tensor_principleR_m4{,_neg,_pos}.log）：
  - 補正の対の符号に一様な向きを課すと e ≥ 0 は解なし、e ≤ 0 は解あり（一意 6/15）。符号 −(−1)^m は R_m, R_{m+1} に続く交代 (−1)^{m+j′+1}。
  - 一意な値は「P = V⊗V·ψ²V + 2qt²V + t⁴ + q²t² の V⊗V に入れ替え、残りは対角」の対角の部分と単項式ごとに合う。符号の変種 32 通りのうち (i) の行と台を満たすのは 1 通り：R′_{m+2} = −(−1)^m t^{m−3}[ψ²V·ψ²V + t⁴ − q²t²]（2qt²V は対、対 13）。(i)+(ii) も満たす。
  - principle に j = m+2 の枝として入れ、T(2,9) の E 全体を作った：N の台の外 0、|e| ≤ |c|・偶奇の違反 0、入れ替わる対 40、**(i)・(ii) をともに満たす**。予言として twisted_E_n1_m4.pkl に保存（principle 4 --save、logs/refined_kauffman/f_twisted_tensor_principle_m4.log）。
  - 【予想】補正の対角の部分 D₀ = 1、D₁ = ψ²V + t²、D₂ = ψ²V² + t⁴ − q²t²。T(2,11) で検証する。
- T(2,11)（m = 5）の 2 箱のデータ作成の見積もり（候補 (7b)）：
- 目的：n = 1 の補正 F_{j′}（R_{m+j′} = t^{m−1−j′}F_{j′}）の 3 点目 F₃（84 状態の見込み）と、F₁ = V² + t²・F₂ = P が m によらないことの確認（theory §12.29）。
- 測定（logs/refined_kauffman/c_fast_daha_dn_m5.log、c_fast_adj_m5.log）：D_N の 2 箱は N = 12 で 52 秒（m = 4 の 1.6 倍、N あたり約 1.2 倍）、最大 RSS 0.4 GB。随伴は N = 16 で 56 秒（N あたり 1.1〜1.3 倍）。
- 見積もり：U の冪の範囲は幅 8m の規則から U^{−42}〜U^{−20}（23 項）。絞った補間は m = 4 で N = 3〜8・4〜10 から決まったので m = 5 は N ≤ 12 程度で決まる見込み。D_N は N ≤ 20 まで作り（残りは照合用、合計約 20 分）、随伴は ν の冪 −42〜−20 のノード N_A = 20〜42（N = 42 で約 25 分、合計 4〜5 時間）。24 時間を超えない。
- 再開：どちらも N ごとに保存し保存済みは飛ばす。実行は nohup の独立プロセス、完了の目印はログの「完了」。
- 随伴の補間も絞った方法に替えた（rk_c_bounded_save.py --adj、logs/refined_kauffman/c_bounded_save_adj_m{3,4}.log）：R_adj(N) = W/(qn² − κ) を flint で約分し（data/cache にキャッシュ）、x = A^{1/2} = s^N の冪ごとの係数を r の次数ごとに解く（rk_c_bounded_interp に UNIT・PAR を追加。SO は 2・2、随伴は 1・1）。m = 3, 4 とも **N_A = 3〜8 だけで一意**に決まり、使わない N_A（≤ 26・≤ 34）の全部で照合が通り、Newton 補間で作った保存済みの Radj_A_m{3,4} と全係数一致（最初の比較は係数 0 の j の保存の有無で False になっていたので直した）。したがって m = 5 の随伴は N_A ≤ 20 で足り、N_A = 21〜42 の計算（4〜5 時間の見込み）は N = 23 で止めた。rk_c_2box_largeN_m.py は Radj_A_m{m}.pkl があれば随伴の Newton 補間を飛ばす。

## 2026-10-31（続き 97）T(2,11) の 2 箱で原理と予想を検証し、補正の母関数を読む（自己判断で選んだ候補 (7b)）
- データ作成（theory §12.30）：D_N の 2 箱 N ≤ 20（合計約 40 分）、随伴 N_A ≤ 25。R_U_m5_{1,2,11}（絞った補間、(2) は N = 3〜8、(1,1) は N = 4〜10、約分 1078・615 秒）、Radj_A_m5（絞った補間 --adj、N_A = 3〜9、287 秒）、largeN_2box_m5（g の整数性 OK）。logs/refined_kauffman/c_bounded_save_m5.log・c_2box_largeN_m5.log。
- 結果（logs/refined_kauffman/f_twisted_tensor_{blocks,principle}_m5.log）：
  - R₅ = t⁴ = t^{m−1}、R₆ = t³(V² + t²)、R₇ = t²·P（T(2,9) の F₂ と全係数一致）：**F₁, F₂ は m によらない**。
  - R₈ = t·F₃、F₃ は 28 単項式・84 状態 = C(9,3)（予想どおり）。
  - k = 11 以外を原理（j = m+2 の Ĩ は T(2,9) で選んだ候補、符号 −(−1)^m = +）で固定し k = 11 を自由にすると (i)・(i)+(ii) に解あり：**R′_{m+2} の予言が T(2,11) で的中**。k = 11 は (i)+(ii) で 6/43 一意、対の最小 25。
  - 母関数の予想（correction_F、test_correction_generating_function）：Σ F_{j′}z^{j′} = [1 + (2t + t² + qt + 2qt²)z + qt³z²]/[(1 − z)(1 − qtz)(1 − t²z)(1 − q²t²z)]。F₀〜F₃ を再現し、F₄ = 165、F₅ = 286 状態（正係数、C(2j′+3, 3)）。
- 結論：【計算】原理は j ≤ m+2 で T(2,11) まで整合。【計算】F_{j′} は m によらない。【予想】F_{j′} の母関数（4 セル {1, qt, t², q²t²} の Fock 空間＋1・2 量子の補正）。【未解決】j′ ≥ 3 の Ĩ。
- テスト：test_correction_generating_function を追加（母関数が F₀〜F₂ を再現、状態数 C(2j′+3, 3)、j′ ≤ 5 で正係数）。principle の j ≥ m+3 の枝は母関数との一致を表示する。
- 次：(7c) T(2,13)（m = 6）で F₄ を検証（D_N・随伴 N ≤ 20 を作成中。D_N は m = 5 の N = 20 で 268 秒なので m = 6 で 1 時間程度の見込み）。principleR 5（R₈ の対の符号の範囲、e ≤ 0 / ≥ 0）を実行中。

## 2026-10-31（続き 98）n = 1 の補正の Ĩ の一般則（自己判断）
- 方法（theory §12.31）：T(2,11) の k = 11（j′ = 3）の補正 R₈（28 単項式・84 状態）の対の符号を principle_R（j′ = 2 の枝を追加）で (i) で絞り、e ≤ 0 / ≥ 0 の制約、最小 Σ|e| での一意性を見た。
- 結果（logs/refined_kauffman/f_twisted_tensor_principleR_m5{,_neg,_pos}.log）：e ≥ 0 は解なし、e ≤ 0 は解あり（一意 8/28、最小 Σ|e| = 16 = 4² で 23/28 が一意）。符号 (−1)^{m+j′+1}、対角の部分の状態数 (j′+1)²。D₀, D₁, D₂ が ψ² の像（X₀ = 1、X₁ = V + t、X₂ = Sym²V + tV）であることから X_{j′} = Sym^{j′}V + t·Sym^{j′−1}V（母関数 (1 + tz)/((1 − z)(1 − tz)(1 − qtz))）と読み、X₃ の仮説 e = −t²ψ²X₃ が (i)・(i)+(ii) の行と台を満たした（対 34）。
- principle を一般の j′ に拡張（correction_F・diag_X、--free=j で指定のブロックを自由に）。T(2,11) の E 全体：台の外 0、違反 0、対 80、(i)・(ii) をともに満たす（twisted_E_n1_m5.pkl に保存）。T(2,9) は保存済みと一致。
- 途中の失敗：principle_R に j′ = 2 の枝が無く m = 5 で解なしになっていた（直した）。diag_X の最初の式は c₀ = 2 になっていた（直した）。符号付き Fock の母関数（セルに符号）で D を作る試みは D₂ を再現せず破棄。
- テスト：test_diag_X を追加（X の形、母関数、(F − ψ²X)/2 が非負の整数）。
- 結論：【計算＋予想】n = 1 の E は T(2,3)〜T(2,11) の全ブロックで条件なしに作れる。【予想】F_{j′}・X_{j′} の母関数。【未解決】物理的な由来、n = 2 への移植。
- n = 2 への移植の必要条件【計算】：λ = (2) の ν^{−12} の補正 R（18 通りの Y のそれぞれ）について、「対角の部分 = 単項式 μ₀ × ψ² の像、残りは対」なら R の奇数係数の単項式が一つの偶奇の類に乗るはずだが、**どの Y でも 3〜4 個の類に分かれる**（Y₂ では t, q²t², q³t³, q⁴t⁵）。λ = (1,1) の補正 t(1 + qt)² は t·ψ²(1 + qt) で成り立つ。λ = (2) には n = 1 の型の Ĩ は（この平方の分け方では）載らない。§12.28 の否定的な結果と整合する。
- 次：T(2,13) の E 全体（連鎖を実行中）。終わったら n = 1 の原理を claims/RESULTS 向けにまとめる。n = 2 の λ = (2) は平方の分け方を変える（Y_a ⊗ Y_b + Y_b ⊗ Y_a など）必要がある。

## 2026-10-31（続き 99）T(2,13) で n = 1 の原理と母関数を検証（自己判断）
- データ（logs/refined_kauffman/c_fast_daha_dn_m6.log・c_fast_adj_m6.log・c_bounded_save_m6.log・c_2box_largeN_m6.log）：D_N の 2 箱 N ≤ 20（N = 20 で 630 秒）、随伴 N_A ≤ 20（N = 19 で 956 秒）。R_U_m6_{1,2,11}（U^{−50}〜U^{−24}、27 項。(2) は N = 3〜9、(1,1) は N = 4〜11 で一意、約分 1920・1060 秒）、Radj_A_m6（N_A = 3〜9、261 秒）、largeN_2box_m6（g の整数性 OK）。被覆の側の彩色データ T2_13_{1,2,11}.json が無くて g の作成で止まったので、scripts/framed_knot/stage2_torus.py の KNOTS に T2_13 を加え compute を直接呼んで作成（各 1 秒以下。参照表には無いので main は使えない）。
- 結果（theory §12.31 の追記、logs/refined_kauffman/f_twisted_tensor_{blocks,principle}_m6.log）：上の系列 j′ = 0..4 の補正がすべて母関数の予想 t^{m−1−j′}F_{j′} と全係数で一致（**F₄ = 165 状態の予言の的中**、F₂・F₃ は m = 4, 5, 6 で同じ）。原理で作った E 全体（対 140）は台の外 0・違反 0 で **(i)・(ii) をともに満たす**（D₄ = ψ²X₄ の予言の的中）。twisted_E_n1_m6.pkl に保存。
- n = 2 の次の段の観察【計算】：ν^{−14} の結び目の部分は t = q で正（λ = (2)：q(1 + q²)²[3]_q × 11 項、λ = (1,1)：q(1 + q²)[3]_q × 9 項）で平方でなく、背景は t = q で 0。ν^{−16} は負。「−Y⊗Y + (q − t)R」の形は ν^{−12} だけ（続き 84 の観察と整合）。N₁₄ は両 λ とも既約で、1 + qt、V、V′、1 + q、1 + t、q − t のどれでも割り切れない（t = q では q(1 + q²)^{2 or 1}[3]_q × 回文多項式、qt = −1 では小さい）。積による同定は効かず、n = 1 と同じく「状態の空間の生成関数」として読む必要がある。
- 補正の分解 F = ψ²X ⊕ 2P（theory §12.31 の追記、test_correction_pairs_generating_function）：対の片方 P の母関数は t z (V′ + q t² z)/((1 − z)(1 − qtz)(1 − t²z)(1 − q²t²z))、V′ = 1 + q + qt。F の母関数の分子はこの分解そのもの。
- 論文：paper/claims.md に P17〜P22（符号の規則の破れ、n = 1 の原理、補正の母関数、補正の Ĩ、(ii) が対を見ないこと、n = 2 の背景と結び目の部分）と D15（絞った補間）を追加。
- 結論：【計算】n = 1 の原理（状態の空間の同定と Ĩ）は T(2,3)〜T(2,13) の全ブロックで (i)・(ii) と整合する。【予想】F_{j′}・X_{j′} の母関数（j′ ≤ 4 で検証）。
- 次：(A) 母関数の物理的な由来と論文 §7 への反映、(B) n = 2 の ν^{−14} 以降の同定、(C) λ = (2) の補正の Ĩ。
