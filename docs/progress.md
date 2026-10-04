# 進捗記録

## 現在のステージ

セットアップ完了。ステージ 0 未着手（下記「前セッションの先行作業」は、ステージ 0 で再確認する）。

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

## 次にやること

1. docs/novelty_report.md を配置してもらう（未配置）。docs/research_plan.md は、元の指示書 CLAUDE_TASK_refined_cable_knots.md を移して配置した（2026-10-05）。
   研究計画 §3 のディレクトリ構成（refined_cables/…、REPORT.md、LOG.md）は現在の構成（docs/、src/rcsknots/、references/ など）と違う。現在の構成に読み替える。
2. ステージ 0：規約表を docs/conventions.md に確定する。停止点なので、確定後に報告して確認を求める。

## 未解決の問題

- 素朴なプレシズムでトーラス結び目が再現できない件（上記）。原因は未特定（γ 補正の要否をステージ 2 で調べる）。
- Wolfram Engine のアクティベート（ユーザー本人の操作が必要）。
