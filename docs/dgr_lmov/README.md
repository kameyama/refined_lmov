# 別の研究：DGR・四重次数型の構造と refined LMOV 不変量（2026-10-11 開始）

このディレクトリは、本研究（KN 1703.05408 の refined LMOV 不変量の正値性を、トーラス結び目を超えて検証する研究。docs/framed_knot/）とは**別の研究**の文書である。
対象（refined CS 不変量と KN の refined LMOV 不変量）と道具（コード・データ）は共有するが、目的と規則が異なる。

## 目的
- DGR（math/0505662）や四重次数の色付きホモロジー（GGS 1304.3481）の考え方（微分 d_{n|m}、対称性、色付き微分、δ 次数、t_r/t_c）を、refined LMOV 不変量の spin 表示・genus 表示で、どう表せるかを調べる。
  - 同じ構造はそのままは成り立たないので、非 refined と refined の LMOV のデータを観察して 0 から考える（ユーザーの最初の指示）。
- 見つけた構造と非 refined のデータから refined LMOV と rCS を構成し、まだ計算されていない値を予言する（計画 D7〜D11）。特に、まだ数の少ない [2,1]・[3,1] が「unrefined の LMOV ＋ いくつかの仮定」から構成できるかを調べる。
- \hat f_ρ の側で自然に定義できる基底（SU(2)_L の既約表示）・特殊化の点（a = t^N q^{−M}）・種数 0 の部分に、DGR の対称性・微分・正値性・指数増大性がどう対応するかを整理する（計画 S0〜S6）。

## 本研究との違い（規則）
| 項目 | 本研究（docs/framed_knot/、CLAUDE.md） | この研究 |
|---|---|---|
| 目的 | KN の正値性予想の検証 | 新しい構造の発見と、構造からの構成・予言 |
| refined LMOV の正値性 | 検証の対象（予言として確かめる） | 構成の条件に使ってよい（2026-10-11、ユーザーの指示「この計画では refined LMOV に正値性は課して良い。今回は KN の検証をする計画ではなく新しい研究である」） |
| rCS の正値性 | rCS の候補を選ぶ際に課さない（open_problems.md §F） | §F は適用しない（2026-10-11、ユーザーの指示「今回は研究自体が異なるので衝突しない」）。Poincaré 多項式である rCS（対称表現など）には課してよい。非長方形の rCS は負の係数をもつので、数学的な理由で課さない |
| 停止点 | docs/framed_knot/research_plan*.md | research_plan.md の §3 |

共通の規則（CLAUDE.md のうち研究の進め方に関するもの）はこの研究にも適用する。
- 主張の分類（【導出】【計算（範囲）】【予想】【仮説】）
- 文献情報の確認
- 既知の結果を成果として書かないこと（本研究の novelty_report.md で既知とされたものを含む）
- 場当たり的な調整の禁止
- 自動テスト
- 外部への送信の事前許可

## ファイル
| ファイル | 内容 |
|---|---|
| research_plan.md | 計画（段階 D0〜D11 と S0〜S6、状態の表、停止点、記録の規則） |
| report.md | 本文（結果と導出。§0〜§10 が最初の観察、§11〜§14 が σ の辞書・既約表示・特殊化の点・種数 0、後半が計画の段階 D1〜D11 ごとの記録） |
| progress.md | 進捗（日付ごと、新しいものが上） |
| open_problems.md | 未解決（H1〜H18） |
| novelty.md | 新規性と既知の結果（文献調査） |

## コード・データ・記録
- コード：src/rcsknots/dgr_lmov/dgr.py（因子・特殊化の道具）、src/rcsknots/dgr_lmov/lmovstruct.py（既約表示・σ・特殊化の点）、src/rcsknots/dgr_lmov/large_n_inverse.py（GT1/GT2 の逆解き）、scripts/dgr_lmov/dgr_*.py（D 段階）、scripts/dgr_lmov/lmovstruct_*.py（S 段階）、scripts/dgr_lmov/lmov21_*.py, scripts/dgr_lmov/lmov31_linear.py（D11 の一次の解析）。
- テスト：tests/dgr_lmov/test_dgr_lmov.py, test_lmovstruct.py, test_lmov21_linear.py（既定の実行に含まれる）。保存値の再判定は tests/dgr_lmov/test_lmovstruct_results.py と results 印のテスト（既定の実行から外れる）。
- データ：data/dgr_lmov/（判定の要約、4 箱以上の判定 highbox/、構成の追跡 follow/、2 箱の予言 d10/、3 箱 box3/）、data/dgr_lmov/lmovstruct/（既約表示の調査、801 結び目、特殊化の点、D11 の零空間と多面体）。
- 記録：logs/dgr_lmov/dgr_*.log、logs/dgr_lmov/lmov*.log。
- 4 箱以上の入力（data/framed_knot/highbox/）は本研究のデータで、git の外にある（本体のチェックアウトのものを読むだけ）。

## 作業を再開するとき
1. progress.md の先頭「現在の状況（再開用）」で、確立したこと・撤回したこと・実行中の計算・再開の手順を確認する。
2. research_plan.md の状態の表と停止点を確認する。
3. 規約は本研究と共通（docs/framed_knot/conventions.md）。変数は太字 (𝐚, 𝐪, 𝐭)、用語（spin 表示・genus 表示）は docs/framed_knot/paper/notation.md。
