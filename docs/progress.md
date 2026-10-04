# 進捗記録

## 現在のステージ

ステージ 4 完了。ステージ 5（refined ケーブル公式と DAHA の照合）に着手。

## 完了した作業

### 2026-10-05 リポジトリと環境のセットアップ

- git を初期化した（リモートなし）。ディレクトリ構成、.gitignore、pyproject.toml、justfile を作った。
- `.venv`（Python 3.14.8）に sympy, numpy, mpmath, python-flint, pytest を入れ、rcsknots を editable インストールした。
- 外部ツールを確認した：SageMath 10.9 は使える（conda 環境 `sage`）、Wolfram Engine は未アクティベート。詳細は README.md「環境」節。
- `.venv/bin/pytest`：30 件すべてパス（ダミーテスト 1 件、下記の Macdonald 基盤テスト 29 件）。

### 2026-10-05 前セッションの先行作業（新しい構成に移行する前に行ったもの）

セットアップ指示の前に、指示書に沿って次の作業を進めていた。成果は捨てずに新しい構成へ移した。いずれもステージ 0 以降で再確認する。

- 論文 12 本の PDF と TeX ソース、1703.05408 v3 の ancillary ノートブック 2 つを `references/` に取得した（`just fetch-refs`）。
- ノートブックを `scripts/nb2txt.py` でテキスト化し、BPS 抽出の手順（逆 M 行列 → a 次数シフト → 変数変換 → g 展開）を読み取った。
  抽出メモは docs/literature/1703.05408.md の末尾に「暫定」として置いた。
- `src/rcsknots/` に Macdonald 基盤（冪和基底での P_λ、ノルム、主特殊化、指標）を実装した。
  テスト：q=t で Schur に戻る（n≤8）、直交性とノルムの閉公式（n≤9）、単項式基底での三角性（n≤8）、
  1703.05408 付録 A の明示式（n≤3）、主特殊化の閉公式と N=2,3,4 での直接代入の一致（n≤6）。（計算機で確認）
- `knots.py`（トーラス結び目とケーブルの評価、未検証）と `notebook_ref.py`（ノートブックの T(2,2p+1) 閉公式の移植）を書いた。
  ノートブックの式は、三葉結び目 [1] の reduced 値が DGR 規約の a²q⁻² + a²q²t² + a⁴t³ に一致することを確認した。
- 素朴な refined Rosso–Jones（P_λ[p_k→p_{mk}] を q,t のまま P_ν で展開し、T_ν^{n/m} を掛ける）は、ノートブックの値と
  単項式倍の違いでは一致しなかった（T(2,1), T(2,3), T(2,5), T(2,7) の [1],[2],[1,1]、T 固有値の符号 3 通りで確認）。
- 文献調査のサブエージェント 2 つに AS/Gorsky–Negut/Mironov–Morozov と DAHA 系の論文を読ませた。前者の報告：
  1106.4305・1201.3339 は refined の係数が素朴なプレシズムでは得られず、補正因子 γ が要ると述べている。
  これは上の不一致と整合する（docs/literature/notes_AS_GN.md）。

### 2026-10-05 ステージ 0：文献と規約の確定

- やったこと
  - docs/conventions.md の初版を作成。基準は KN の (a,q,t)。各項目を確定／暫定／未確定に分類し、出典の行番号を付けた。
  - 変数対応（KN ↔ DGR ↔ DAHA）を三葉結び目で自動テストにした：`tests/test_stage0_conventions.py`（2 件パス）。
  - novelty_report §4 の 6 項目を調べ、結果を §5 に追記した（被引用 16 件、書誌 27 件、DAHA の明示例、既存の refined ケーブル公式）。
  - DGR（math/0505662）を取得した。`just wolfram` レシピを追加した（ユーザー許可済み）。
- 結果
  - 研究計画 §2.2 の素朴なプレシズム仮説は誤り（文献：DMMSS・MMSS・Danilenko・Morozov、本研究の試算）。
    修正版のケーブル公式は「DAHA の定義を companion の段で線形に分解したもの」（conventions.md §7、仮説部分を明記）。
  - KN の 4_1 は [3],[2,1],[1,1,1] まで既知（付録 C）。新規性 (N2) を「ケーブルについての [3],[2,1],[1,1,1]」に絞った。
  - K1 の [1],[2] と K2〜K5 の [1] の refined 不変量は CD14 に明示例がある（単項式倍を除く）。ステージ 5 の照合対象。
  - DAHA は framing の単項式を捨てるので、ケーブルの単項式正規化は本研究で決める必要がある（conventions.md §5、未確定）。

### 2026-10-05 ステージ 2：トーラス結び目の再現（ユーザーがステージ 0 の 2 点を承認した後）

- やったこと
  - GN（1304.3328）の楕円 Hall 代数の式 (fix3) を実装した（`ehall.py`, `gn.py`）。多数の評価点をまとめて扱う mod p の numpy 演算（`field.py`）で計算する。
  - ノートブックのトーラス結び目の値を Wolfram で書き出した（`scripts/export_notebook_torus.wls` → `data/reference/notebook_torus.txt`、27 件）。
  - mod p 補間で厳密な多項式を復元し（`exact.py`）、正規化の原理 (N1)(N2) を実装した（`normalize.py`）。補間に使っていない点で、別の素数を使って検証している。
  - 非 refined の Rosso–Jones（Schur 関数、sympy で厳密）を実装した（`homfly.py`）。
- 結果（すべて計算機で確認）
  - ノートブックの全 27 件（T(2,3),T(2,5),T(2,7),T(2,9) の |λ|≤3 全 6 彩色、T(3,4) の [1],[2],[1,1]）について、
    GN の構成 + (N1)(N2) の正規化で**全係数が一致**した（`data/torus/*.json`、`scripts/stage2_torus.py`、ログ `logs/stage2_torus.log`）。
  - 全 27 件の t=q 極限が Rosso–Jones の framing 0 の彩色 HOMFLY-PT に一致した。三葉結び目 [1] は DGR の a²q⁻² + a²q²t² + a⁴t³（ステージ 0 のテスト）。
  - ミラー・転置対称性（[2]↔[1,1], [3]↔[1,1,1]）が全係数で成り立った。
  - EHA の内部整合：P̃_{k,0} = φ(p_k) 倍（k=1,2,3）、交換関係が 3 つの三角形で同じ定数 κ=(1−q)/t。
  - 正規化定数 c_k と基底の読み方（I_μ = H̃_μ/∇_μ）は GN 本文から直接読めず、数値で決めた（conventions.md §6）。
  - 実行時間：トーラス結び目 1 件あたり 0.1〜2.5 秒。
  - `.venv/bin/pytest`：115 件パス、2 件スキップ（ノートブックに無い彩色）。

### 2026-10-05 ステージ 3：refined BPS 抽出の検証（ユーザーがステージ 2 を承認した後）

- やったこと：f^q, f^t を生成関数の plethystic Log から |μ|≤3 で導出し（`bps.py`）、評価点での数値計算と補間で \hat f を厳密に復元し（`bpsnum.py`, `extract.py`）、
  ノートブックの `gfactor` と同じ g 展開で \hat N を取り出した。4_1 の入力も Wolfram でノートブックから書き出した。KN の表は TeX から読み取った（`scripts/parse_kn_tables.py`）。
- 結果（計算機で確認）
  - KN 付録 C の 13 表（T(2,3) の |ρ|≤3 全 6 表、4_1 の全 6 表、T(2,5) の [1]）を全項で再現した。
  - M 行列はノートブックの Mq2, Mq3 と一致。f^q は付録 B と一致（尺度 (t/q)^{|μ|/2} の違いは記録）。f^t は付録 B と一致、ただし f^t_{[1,1,1]} は付録 B の誤植でノートブックの式と一致。
  - f^t = M(q^{−1}) \hat f が T(2,3), T(2,5), 4_1 の全評価点で成立。
  - T(2,5) の |ρ|=2,3（KN 本文に表なし）もすべて非負、2J_r の偶奇はそろう。4_1 では偶奇がそろわない（KN の記述どおり）。
  - a 次数シフトの規約が KN の表の中で統一されていないことを見つけた（T(2,3) とそれ以外で逆）。本研究は nontorus 規約に統一。
  - 実行時間：T(2,3) の全 6 表で約 6 秒。
  - `.venv/bin/pytest`：142 件パス、2 件スキップ。

### 2026-10-05 ステージ 4：非 refined ケーブル公式の検証

- やったこと：Sage でケーブルのブレイド（companion の各ストランドを p 本の平行線にし、(σ₁⋯σ_{p−1})^{q−p·w} を追加）から HOMFLY-PT を計算した
  （`scripts/braid_homfly.sage` → `data/reference/braid_homfly.txt`）。ケーブル公式は sympy 版（`homfly.py`）と mod p 数値版（`homfly_num.py`、補間で厳密化）を用意した。
- 結果（計算機で確認、`tests/test_stage4_cable_homfly.py` 39 件）
  - framing 補正の向き：writhe 1 で描いた自明な結び目の (2,q₂)-ケーブル = T(2,q₂)（q₂ = 1,3,5,−3）。p=1 で companion に戻る。
  - Sage の変数対応 L = i𝐚, M = −i(𝐪−𝐪^{−1}) と鏡像（正のブレイド ↔ DGR の負の結び目）で、三葉結び目が DGR の式になる。
  - ケーブル公式（θ_λ = q^{κ_λ}a^{|λ|/2}、framing 補正 θ_λ^{−pq}、DGR の向き）の [1] が K1〜K6 の 9 結び目すべてでブレイドの HOMFLY-PT に一致。
  - |λ|≤3 の全彩色で、9 結び目すべてについて H_{λ^T}(a,q) = H_λ(a,q^{−1}) と H_λ(a,q=1) = H_□(a,1)^{|λ|} が厳密に成立。
  - sympy での厳密計算は |λ|=3 で遅すぎたので数値版に切り替えた（K1 の [3] で約 3 秒）。

## 次にやること

1. （済）docs/novelty_report.md がユーザーにより配置された（2026-10-05）。docs/research_plan.md は、元の指示書 CLAUDE_TASK_refined_cable_knots.md を移して配置した（2026-10-05）。
   研究計画 §3 のディレクトリ構成（refined_cables/…、REPORT.md、LOG.md）は現在の構成（docs/、src/rcsknots/、references/ など）と違う。現在の構成に読み替える。
2. （済）ステージ 0。ユーザーの確認を待つ。
3. （済）ステージ 2。ユーザーの確認を待つ。
4. （済）ステージ 3：BPS 抽出コード（M 行列、f^q、\hat f、g 展開）を実装し、T(2,3) の表と 4_1 の表を再現する。ノートブックの \hat f の出力も Wolfram で書き出して照合に使う。
5. （済）ステージ 4：非 refined ケーブル公式（`homfly.py` の `cable_unreduced`）を、Sage のブレイド HOMFLY と照合する。

## 未解決の問題

- （原因判明）素朴なプレシズムでトーラス結び目が再現できない件：素朴な持ち上げ自体が誤り（文献で既知）。GN/DAHA の構成に切り替える。
- ケーブルの単項式（framing）正規化の決め方（未確定、conventions.md §5）。
- T 行列の t 指数の符号：EHA で計算するので不要になった（未決着のまま、conventions.md §3）。
- （解決）Wolfram Engine はユーザーがアクティベートした（2026-10-05）。サンドボックス外でのみ動く。
