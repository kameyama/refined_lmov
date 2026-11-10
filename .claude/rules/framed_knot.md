---
paths:
  - "docs/framed_knot/**"
  - "scripts/framed_knot/**"
  - "tests/framed_knot/**"
  - "data/framed_knot/**"
  - "logs/framed_knot/**"
  - "src/rcsknots/*.py"
---
# framed_knot の規則

Kameyama–Nawata "Refined large N duality for knots"（arXiv:1703.05408）は、主にトーラス結び目・トーラス絡み目について refined Chern–Simons 不変量を resolved conifold 上の refined BPS 不変量 \hat N_{R,g,β,J_r} に書き換え、正値性予想を検証した。このテーマでは、これをトーラス結び目以外の結び目と framed knot に拡張し、G 展開の前後の正整数係数性を検証する。

## 文書
- 進捗：docs/framed_knot/progress.md。主要な結果：RESULTS.md（計算時の判定とログを一次の記録として整理する）。全体の整理：REPORT_overview.md。
- 計画と停止条件：research_plan.md（ステージ 0〜7）と research_plan_*.md。既知の結果と新規性の範囲：novelty_report.md。未解決：open_problems.md。
- 規約：conventions.md。変数 (a, q, t) の対応、Macdonald 多項式の正規化、S・T 行列と γ 因子、framing（ケーブルの q は companion の Seifert framing で測る）、reduced/unreduced と U(N)/SU(N)、プレシズム、refined BPS の抽出式を定めている。
- 論文の準備：paper/（構成案 outline.md、主張の格付け claims.md、記号 notation.md）。

## 正値性の判定
- 重視するのは、\hat f_ρ から得る refined LMOV 不変量の、G 展開の前（太字の係数 \hat𝐍_{ρ;i,j,k}）と後（\hat N_{ρ,g,β,J_r}）の両方での正整数係数性である。f^q, f^{\bar t} の係数は refined LMOV 不変量ではなく、補足にとどめる。
- rCS の候補を選ぶときに rCS の正値性を課さない（open_problems.md §F）。
- ヤング図 [1] は特別に扱う。\hat f_{[1]} = (a^{1/2} − a^{−1/2}) rCS_{[1]} なので、[1] に限り、homological な不変量（ホモロジーの Poincaré 多項式として解釈される量）はすべて正整数係数でなければならない。G 展開が閉じて正整数係数になるかも [1] で確認できなければならない。
- 新しい結び目では、まず [1] で G 展開の前と後の両方の正整数係数性を確認し、その結果を彩色付きの判定より先に報告する。

## 照合するもの（tests/framed_knot/）
1703.05408 の値の再現、t = q での colored HOMFLY との一致、ブレイドから計算した HOMFLY との一致、DAHA の明示例との一致、R ↔ Rᵀ の対称性。

## 計画に入れないもの
- 非代数的なケーブルについて、DAHA の値を入力とする正値性の計算（3 箱への拡張や破れ方の分類を含む）。理由は open_problems.md D5。
- novelty_report.md で既知とされた結果（トーラス結び目の refined 不変量、4_1 の [1]/[2]/[1,1] の BPS 不変量、DAHA による iterated torus の定義など）を、このテーマの成果として書かない。

## 停止点
ステージ 0 の後（規約表の確定）、ステージ 2 の後（トーラス結び目での 1703.05408 の再現）、ステージ 5 の後（refined ケーブル公式と DAHA の照合）、research_plan*.md の停止条件に当たったとき。
