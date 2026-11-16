# refined CS 不変量と refined BPS 不変量の研究

[English](README.md) | 日本語

結び目の refined Chern–Simons 不変量と refined BPS（LMOV）不変量についての研究リポジトリ。下の表のテーマについて、コード、テスト、計算結果、論文の原稿を置いている。テーマはディレクトリ名で呼ぶ。

| テーマ | 内容 | 入口 |
|---|---|---|
| framed_knot | Kameyama–Nawata "Refined large N duality for knots"（arXiv:1703.05408）を、トーラス結び目以外の結び目と framed knot に拡張し、refined LMOV 不変量の正値性を検証する | [progress.md](docs/framed_knot/progress.md)、[REPORT.md](docs/framed_knot/reports/REPORT.md)、[conventions.md](docs/framed_knot/conventions.md) |
| dgr_lmov | DGR・四重次数型の構造と refined LMOV 不変量 | [README.md](docs/dgr_lmov/README.md) |
| refined_kauffman | refined Kauffman 多項式の large N duality | [README.md](docs/refined_kauffman/README.md) |

作業の規則は [CLAUDE.md](CLAUDE.md) と .claude/rules/。

## 論文

| テーマ | 題名 | ソース |
|---|---|---|
| framed_knot | Refined large N duality for framed knots and positivity of refined LMOV invariants beyond torus knots | [docs/framed_knot/paper/en/](docs/framed_knot/paper/en/)（[PDF](docs/framed_knot/paper/en/main.pdf)） |
| refined_kauffman | Refined crosscap invariants and the large N duality for colored Kauffman invariants | [docs/refined_kauffman/paper/en/](docs/refined_kauffman/paper/en/)（[PDF](docs/refined_kauffman/paper/en/main.pdf)） |

## 言語

この研究は AI（Claude、Anthropic）との対話で進め、対話は日本語で行った。そのため、docs/ 以下の作業文書（研究計画、進捗、規約、報告）、作業の規則（CLAUDE.md、.claude/）、commit メッセージ、コードのコメントの多くは日本語で書いてある。上の表の論文と英語版の README は英語である。

英語版の [README.md](README.md) が本体で、この文書はその日本語版である。2 つは同じ内容に保つ。

## 構成

各テーマのファイルを、各ディレクトリの中でテーマ名のサブディレクトリに分けて置く。

```
docs/framed_knot/     framed_knot の研究計画・新規性報告・規約・進捗・結果・文献メモ（literature/）・論文（paper/）
docs/dgr_lmov/        dgr_lmov の文書（入口は README.md）
docs/refined_kauffman/  refined_kauffman の文書・論文（入口は README.md）
references/           論文 PDF と TeX ソース（追跡しない。references/README.md に取得元）
src/rcsknots/         計算コード（Python パッケージ。全テーマで共有する）
src/rcsknots/dgr_lmov/  dgr_lmov だけが使うモジュール（dgr, lmovstruct, large_n_inverse）
src/rcsknots/refined_kauffman/  refined_kauffman だけが使うモジュール（koornwinder, daha_dn, daha_cn, flint_rf*, kauffman_*）
tests/framed_knot/    framed_knot の照合テスト（pytest）
tests/dgr_lmov/       dgr_lmov の照合テスト
tests/refined_kauffman/  refined_kauffman の照合テスト
scripts/framed_knot/  framed_knot の単発の計算・検証スクリプト
scripts/dgr_lmov/     dgr_lmov のスクリプト（framed_knot/ のスクリプトも import する）
scripts/refined_kauffman/  refined_kauffman のスクリプト（rk_*.py）
data/framed_knot/     framed_knot の計算結果（JSON と人が読める表）
data/dgr_lmov/        dgr_lmov の計算結果（lmovstruct/ を含む）
data/refined_kauffman/  refined_kauffman の計算結果（daha/ の pkl など）
data/cache/           途中結果のキャッシュ（追跡しない。全テーマで共有する）
logs/framed_knot/     framed_knot の検証ログ・実行時間とメモリの記録
logs/dgr_lmov/        dgr_lmov のログ
logs/refined_kauffman/  refined_kauffman のログ
notebooks/            試行用
```

スクリプトは `PYTHONPATH=src:scripts/framed_knot:scripts/dgr_lmov` で実行する（refined_kauffman は `PYTHONPATH=src:scripts/refined_kauffman`）（スクリプトどうしの import はこのパスで解決する）。

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
| Wolfram Engine | 使える（Claude Code のサンドボックス外のみ） | 15.0.0 for Mac OS X ARM (64-bit) | `wolframscript -code ...` / `wolframscript -file x.wls` |

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
