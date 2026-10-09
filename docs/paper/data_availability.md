# データとコードの公開の準備（2026-10-09、執筆前の準備。送信はしていない）

論文の結果を再現するためのデータとコードの一覧と、公開の方法の選択肢。**公開（arXiv の補助ファイル、GitHub の公開、Zenodo などへのアップロード）は外部への送信なので、実施の前に許可を得る。**

## 1. 公開の対象
| まとまり | 中身 | 大きさ（git の管理下） | 論文での使い道 |
|---|---|---|---|
| data/torus/, data/cables/ | rCS_λ（\|λ\| ≤ 3、正規化 (N1)(N2') 済み） | 約 30 MB | 入力（§4）。CD14 の明示例との照合 |
| data/bps/ | \hat f_ρ と \hat N_{ρ,g,β,J_r}（framing 0、\|ρ\| ≤ 3）。21 結び目 | .json/.txt 約 240 MB、.json.xz 約 330 MB | 主結果の表（§7、generated/summary_table.tex） |
| data/lmov/ | \hat f, f^q, f^{\bar t} の計算結果 | 約 275 MB | 再現の中間結果（data/bps と重なる） |
| data/knot_table/, data/homology/ | 結び目表 801 個の [1] の入力と判定 | 約 1 MB | §5 |
| data/v4/ | ツイスト・pretzel の LP 証明書 | 約 0.5 MB | §6 の証明の検証 |
| data/b4/ | thick な 10_139, 10_152 の入力と判定 | 小 | §6 |
| data/x0/ | t = q の [2,1] の照合 | 小 | §7.3 |
| data/paper_summary.json | 表の素材（項数と負の数） | 小 | §7 の表 |
| logs/ | 計算時の判定のログ（一次の記録） | — | 主張の根拠（claims.md） |
| src/rcsknots/, scripts/, tests/ | コードと機能テスト | — | 再現 |

git の管理下の data/ は全体で約 876 MB（496 ファイル）。data/cache/（途中結果、約 8.8 GB）は git の管理外で、公開しない（スクリプトで再生成できる）。

## 2. 公開の方法の選択肢（判断はユーザー）
1. **arXiv の補助ファイル（anc/）**：arXiv は大きな投稿を推奨していない。表に出る量の要約（data/paper_summary.json、小さい結び目の \hat N の .txt、LP 証明書）に限るのが現実的。
2. **GitHub のリポジトリを公開する**（現在の origin：https://github.com/kameyama/refined_lmov）：コードとデータを一緒に出せる。1 ファイル 100 MB の上限があるので、K10 は .json.xz のまま。
3. **Zenodo などのデータリポジトリ**：DOI が付き、論文から恒久的に引用できる。大きな data/bps と data/lmov に向く。

推す案：コードは 2（GitHub）、大きなデータは 3（Zenodo の DOI を論文に書く）、arXiv には 1 の要約だけ。いずれも実施前に許可を得る。

## 3. 公開の前にすること（許可が出てから）
- 公開する範囲を決め、README.md の「再現の手順」（`just setup` → `.venv/bin/pytest -q` → 表の生成 `PYTHONPATH=src python scripts/paper_summary.py`）を論文用に整える。
- ユーザーの以前の計算（2017-2018 年の Mathematica のノートブック、未発表）を照合に使った部分は、そのノートブックを公開するかを別に判断する（tests/test_framing.py は出力 9 件をテストの中に書き写している）。
- ライセンスを決める。
