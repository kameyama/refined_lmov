# 参照データの出典

このディレクトリのファイルは、計算結果を文献と照合するためのものである。文献に由来するファイルと、このリポジトリで計算したファイルを分けて書く。使った成果を発表するときは、下の文献を引用する。

## 文献に由来するもの

| ファイル | 出典 | 作り方 |
|---|---|---|
| kn_tables.json | M. Kameyama and S. Nawata, "Refined large N duality for knots", arXiv:1703.05408 の表 | scripts/framed_knot/parse_kn_tables.py が TeX ソースの表を読み取る |
| notebook_torus.txt | 同じ論文の補助ファイル refined-LMOV-non-torus.nb（arXiv で公開） | scripts/framed_knot/export_notebook_torus.wls がノートブックの定義を評価して書き出す |
| kntz_package/values.txt | M. Kameyama, S. Nawata, R. Tao and H. D. Zhang, "Cyclotomic expansions of HOMFLY-PT colored by rectangular Young diagrams", arXiv:1902.02275 の補助ファイル HOMFLYPTinvariants.m（arXiv で公開） | scripts/framed_knot/x1_kntz_package.wls がパッケージを実行して書き出す |
| daha_examples.txt | I. Cherednik and I. Danilenko, "DAHA and iterated torus knots", arXiv:1408.4348 の明示例 | scripts/framed_knot/extract_daha_examples.py が TeX ソースから読み取る |
| gs_1112.0030/s2_thick.py | S. Gukov and M. Stošić, "Homological algebra of knots and BPS states", arXiv:1112.0030 の付録 B | TeX ソースの式を書き写した（行番号はファイルの冒頭） |
| stosic_2402.03066/quivers.json | M. Stošić, "Generalized knots-quivers correspondence", arXiv:2402.03066 | TeX ソースから書き写した（行番号はファイルの中） |
| yui_kr/homology.json | K. Nakagane and T. Sano, "Computations of HOMFLY homology", arXiv:2111.00388 のプログラム yui-kr の計算結果（MIT、yui_kr/LICENSE） | yui_kr/SOURCE.md |
| yui_kr/braid.csv, yui_kr/homfly.csv（git では追跡しない） | C. Livingston and A. H. Moore, KnotInfo: Table of Knot Invariants, knotinfo.org | scripts/framed_knot/fetch_knot_table.py が yui-kr のリポジトリから取得する |

題名と著者は arXiv の各論文のページで確かめた（2026-11-16）。

## このリポジトリで計算したもの

| ファイル | 作り方 |
|---|---|
| braid_homfly.txt | scripts/framed_knot/braid_homfly.sage（ブレイドの閉包の HOMFLY-PT を Sage で計算） |
| cable_homfly_p3.txt | scripts/framed_knot/cable_homfly.sage（KnotInfo のブレイドのケーブルの HOMFLY-PT を Sage で計算） |
| x0_cable_p2.txt, x0_cable_p22.txt, x0_cable_p3.txt, x0_double_twist.txt | scripts/framed_knot/x0_double_twist.sage |
| two_cable_homfly.txt | scripts/framed_knot/two_cable_homfly.sage |
| yui_kr/sage_homfly.txt | scripts/framed_knot/knot_table_homfly.sage（KnotInfo のブレイドから Sage で計算） |
