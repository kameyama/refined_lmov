# refined_lmov

結び目の refined Chern–Simons 不変量と refined BPS（LMOV）不変量についての研究リポジトリ。
目的と規則の異なる 3 つの研究を置いている。コード（src/rcsknots/）とキャッシュの一部は共有する。

## 本研究：iterated torus 結び目の refined CS 不変量（framed_knot/）

Kameyama–Nawata "Refined large N duality for knots"（arXiv:1703.05408）の refined BPS 不変量の計算を、
トーラス結び目からその**ケーブル**（iterated torus 結び目）に拡張する。
対象は |λ| ≤ 3 の全ヤング図形 [1], [2], [1,1], [3], [2,1], [1,1,1]。

- 研究計画：[docs/framed_knot/research_plan.md](docs/framed_knot/research_plan.md)
- 新規性の確認：[docs/framed_knot/novelty_report.md](docs/framed_knot/novelty_report.md)
- 規約：[docs/framed_knot/conventions.md](docs/framed_knot/conventions.md)
- 進捗：[docs/framed_knot/progress.md](docs/framed_knot/progress.md)
- 結果のレポート：[docs/framed_knot/REPORT.md](docs/framed_knot/REPORT.md)

## 第 2 の研究：DGR・四重次数型の構造と refined LMOV（dgr_lmov/、2026-10-11 開始）

DGR や四重次数の色付きホモロジーの考え方（微分、対称性、δ 次数など）が refined LMOV 不変量の上でどう表せるかを調べ、
見つけた構造から未計算の値を構成・予言する。入口は [docs/dgr_lmov/README.md](docs/dgr_lmov/README.md)。

## 第 3 の研究：refined Kauffman 多項式の large N duality（refined_kauffman/、2026-10-14 開始）

SO(N)/Sp(N) Chern–Simons の彩色 Kauffman 不変量の large N duality（Mariño 0904.1088：SO(N) の自由エネルギーから
U(N) の composite 表現の寄与の 1/2 を引いた残りが、crosscap c = 1, 2 の整数 BPS 不変量で書ける）を、
refined の設定（Macdonald・Koornwinder 多項式、Aganagic–Schaeffer 1202.4456 の orientifold の refined topological string）に持ち上げ、
refined Kauffman 不変量から refined BPS 不変量を抽出する式を作って、整数性と正値性を調べる。
2018 年の Mathematica による試行錯誤（Dropbox、git の外。読むだけで書き換えない）と同じ方向で進めている。

- 入口と規則：[docs/refined_kauffman/README.md](docs/refined_kauffman/README.md)
- 規約：[docs/refined_kauffman/conventions.md](docs/refined_kauffman/conventions.md)
- 進捗（先頭に再開用の「現在の状況」）：[docs/refined_kauffman/progress.md](docs/refined_kauffman/progress.md)
- 本文：[docs/refined_kauffman/report.md](docs/refined_kauffman/report.md)
- 論文の草稿：[docs/refined_kauffman/paper/ja_v2/](docs/refined_kauffman/paper/ja_v2/)（第 2 版が本体。第 1 版 ja/ は記録として残す）

## 構成

3 つの研究のファイルを、各ディレクトリの中で `framed_knot/`（本研究）、`dgr_lmov/`（第 2 の研究）、`refined_kauffman/`（第 3 の研究）に分けて置く。

```
docs/framed_knot/     本研究の研究計画・新規性報告・規約・進捗・結果・文献メモ（literature/）・論文（paper/）
docs/dgr_lmov/        別の研究の文書（入口は README.md）
docs/refined_kauffman/  第 3 の研究（refined Kauffman 多項式の large N duality）の文書・論文（入口は README.md）
references/           論文 PDF と TeX ソース（追跡しない。references/README.md に取得元）
src/rcsknots/         計算コード（Python パッケージ。直下は研究どうしで共有する）
src/rcsknots/dgr_lmov/  別の研究だけが使うモジュール（dgr, lmovstruct, large_n_inverse）
src/rcsknots/refined_kauffman/  第 3 の研究だけが使うモジュール（koornwinder, daha_dn, daha_cn, flint_rf*, kauffman_*）
tests/framed_knot/    本研究の照合テスト（pytest）
tests/dgr_lmov/       別の研究の照合テスト
tests/refined_kauffman/  第 3 の研究の照合テスト
scripts/framed_knot/  本研究の単発の計算・検証スクリプト
scripts/dgr_lmov/     別の研究のスクリプト（framed_knot/ のスクリプトも import する）
scripts/refined_kauffman/  第 3 の研究のスクリプト（rk_*.py）
data/framed_knot/     本研究の計算結果（JSON と人が読める表）
data/dgr_lmov/        別の研究の計算結果（lmovstruct/ を含む）
data/refined_kauffman/  第 3 の研究の計算結果（daha/ の pkl など）
data/cache/           途中結果のキャッシュ（追跡しない。研究どうしで共有する）
logs/framed_knot/     本研究の検証ログ・実行時間とメモリの記録
logs/dgr_lmov/        別の研究のログ
logs/refined_kauffman/  第 3 の研究のログ
notebooks/            試行用
```

スクリプトは `PYTHONPATH=src:scripts/framed_knot:scripts/dgr_lmov` で実行する。第 3 の研究は `PYTHONPATH=src:scripts/refined_kauffman` で実行する（スクリプトどうしの import はこのパスで解決する）。

## セットアップとテスト

```bash
just setup
```

```bash
just fetch-refs
```

```bash
.venv/bin/pytest
```

## 環境

2026-10-05 に確認した内容。

| ツール | 状態 | バージョン | 実行方法 |
|---|---|---|---|
| Python | 使える | 3.14.8（Homebrew） | `.venv/bin/python` |
| sympy / numpy / mpmath | 使える | 1.14.0 / 2.5.3 / 1.3.0 | .venv |
| python-flint | 使える | 0.9.0 | .venv（`nmod`, `nmod_mat`, `nmod_poly`, `fmpz_mpoly` など） |
| pytest | 使える | 9.1.1 | `.venv/bin/pytest` |
| khoca | 使える（sl(2) のみ。N ≥ 3 は自明な結び目で誤った値を返すので使わない） | 1.5 | .venv。Sage の conda 環境の pari を参照してビルド：`CFLAGS=-I$S/include LDFLAGS="-L$S/lib -Wl,-rpath,$S/lib" pip install khoca`（`S` は conda 環境 `sage` のパス）。`py` パッケージが必要 |
| SnapPy | 使える | 3.3.2（knot_floer_homology 1.2.2） | .venv |
| scipy | 使える | 1.18.1 | .venv（`milp`、整数計画） |
| SageMath | 使える | 10.9（conda-forge、環境 `sage`、Python 3.14.7） | `DOT_SAGE=$TMPDIR/dotsage /opt/homebrew/Caskroom/miniforge/base/envs/sage/bin/sage file.sage` |
| Wolfram Engine | 使える（サンドボックス外のみ） | 15.0.0 for Mac OS X ARM (64-bit) | `wolframscript -code ...` / `wolframscript -file x.wls` |

SageMath の最小例での動作確認（すべて動作した。実行時間は合計約 12 秒）：

- `SymmetricFunctions(QQ['q','t'].fraction_field()).macdonald().P()`：P[2,1] を冪和基底で表示できた。
- `macdonald().Ht()`：H̃[2,1] = s[3] + (q+t) s[2,1] + qt s[1,1,1]。
- プレシズム：`p(P[2]).plethysm(p[2])` を P 基底に展開できた。
- `Link(BraidGroup(2)([1,1,1])).homfly_polynomial()`：L^-2 M^2 − 2L^-2 − L^-4（三葉結び目）。

Wolfram Engine（2026-10-05 にユーザーがアクティベート）：
- サンドボックス内では設定とライセンスを読めず "A WolframKernel location could not be determined" で失敗する。サンドボックス外で実行する。
- 1703.05408 のノートブックは `NotebookImport[file, "Input" -> "HeldExpression"]` で入力セルを取り出して実行できる。
  例：non-torus ノートブックの定義セルを評価し、`FunctionExpand[torusknotsym[1,1,a,q,t]/torusknotsym[1,0,a,q,t]]` が
  a(1 + qt − aq)/√(qt)（三葉結び目の reduced 値）になることを確認した。Python に移植した式（`src/rcsknots/notebook_ref.py`）とも一致する。
