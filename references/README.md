# references/ — 取得した文献

このディレクトリの中身は git で追跡しない（この README だけを追跡する）。`just fetch-refs` で再取得できる。

各論文は `references/<arXiv id>/` に置く。
- `<id>.pdf`：https://arxiv.org/pdf/<id> から取得
- `src.tar`：https://arxiv.org/e-print/<id> から取得した TeX ソース（ancillary ファイルを含む）
- `src/`：その展開先

取得日：2026-10-05。

| arXiv id | 取得版 | 内容（研究計画書の記載による。原典での確認はステージ 0 で行う） | 補助ファイル |
|---|---|---|---|
| 1703.05408 | v3 | Kameyama–Nawata, Refined large N duality for knots | `src/anc/refined-LMOV-torus.nb`, `src/anc/refined-LMOV-non-torus.nb` |
| 1105.5117 | 最新版 | Aganagic–Shakirov | — |
| 1202.2489 | 最新版 | Aganagic–Shakirov | — |
| 1106.4305 | 最新版 | Dunin-Barkowski–Mironov–Morozov–Sleptsov–Smirnov | — |
| 1201.3339 | 最新版 | Mironov–Morozov–Shakirov–Sleptsov | — |
| 1304.3328 | 最新版 | Gorsky–Negut | — |
| 1408.4348 | 最新版 | Cherednik–Danilenko | — |
| 1509.08351 | 最新版 | Cherednik–Danilenko | — |
| 2410.14703 | 最新版 | Cherednik | — |
| 1410.0859 | 最新版 | Morton–Samuelson | — |
| 1605.00978 | 最新版 | Cherednik–Philipp | — |
| 1211.4596 | 最新版 | An. Morozov | — |

| 1604.07456 | v1 | Mellit, Toric braids and (m,n)-parking functions（Duke Math. J. 170 (2021) 4123–4169、arXiv のメタデータで確認） | — |
| 1909.00418 | v1 | Hogancamp–Mellit, Torus link homology | — |
| 1404.4616 | 最新版 | Bergeron–Garsia–Leven–Xin（compositional (km,kn)-shuffle 予想） | — |
| 1405.0316 | 最新版 | Bergeron–Garsia–Leven–Xin, Some remarkable new plethystic operators in the theory of Macdonald polynomials（TeX で確認） | — |
| FranksWilliams1987/ | — | Franks–Williams, Trans. AMS 303 (1987) 97–108（AMS の PDF とテキスト抽出 fw.txt） | — |
| 0905.2555 | 最新版 | Schiffmann–Vasserot（楕円 Hall 代数の定義と中心元） | — |
| math_0505148 | 最新版 | Burban–Schiffmann（取得のみ、未使用） | — |
| 2103.01212 | v2 | Gorsky–Hogancamp–Mellit, Tautological classes and symmetry in Khovanov–Rozansky homology（hard Lefschetz と sl_2 の作用） | — |
| 1704.07630 | v1 | Mellit, Homology of torus knots（Geom. Topol. 26 (2022) 47–70、arXiv のメタデータで確認） | — |
| 1203.5878, 1403.1845, 1105.1151, 1204.5448, 1403.0303 | v1 など | D1 の調査で取得（Hikita、Armstrong–Loehr–Warrington、Gorsky–Mazin I・II、Gorsky–Mazin–Vazirani） | — |
| 1304.3481, 2402.03066, 1707.04017, 1602.02769, 1811.03110, 2108.12645, 1512.07906, 1209.1416, 2012.12008, 2007.00579, 1205.1515, 1303.3709, 2407.20010, 2601.22922, 1308.3838 | 最新版 | B3・B4 の調査で取得（GGS、Stošić、KRSS、Wedrich、EKL ほか。題名は arXiv API と TeX で確認） | — |
| 2505.02059, 2505.05668, 2408.01832, 2508.18417, 2405.10161, 2105.11806, 2312.14905 | 最新版 | B4 の新規性の確認で取得 | — |
| 1112.0030, 1209.1409, 1512.07883 | 最新版 | 第 2 期以降に取得（GS、NRZS、GNSSS） | — |

DOI・巻・頁は未確認（上の表で確認したと書いたものを除く）。

派生ファイル：`1703.05408/nb-torus.txt`, `1703.05408/nb-non-torus.txt` は ancillary ノートブックを
`scripts/framed_knot/nb2txt.py` でテキスト化したもの。再生成するには次を実行する：

```bash
.venv/bin/python scripts/framed_knot/nb2txt.py references/1703.05408/src/anc/refined-LMOV-torus.nb > references/1703.05408/nb-torus.txt
```
| 2205.12832 | v2 | Avetisyan–Mkrtchyan, On refined Chern-Simons / topological string duality for classical gauge groups（arXiv のメタデータで確認） | — |
| 2107.08679 | v1 | Avetisyan–Mkrtchyan, On partition functions of refined Chern-Simons theories on S³ | — |
| 2302.14319 | v1 | Avetisyan–Mkrtchyan, Two-fold refinement of non simply laced Chern-Simons theories | — |
| 2201.12264 | v2 | Kim–Wei, Refined topological vertex with ON-planes | — |
| 1405.6019 | v2 | Piazzalunga–Uranga, M-theory interpretation of the real topological string | — |
| 1101.3216 | v1 | Witten, Fivebranes and Knots（arXiv のメタデータで確認） | — |
| 1007.0263 | v1 | Krefl–Walcher, Extended holomorphic anomaly in gauge theory（Lett. Math. Phys. 95 (2011) 67–88、arXiv API で確認） | — |
| 1010.2635 | v1 | Krefl–Walcher, Shift versus extension in refined partition functions（掲載誌未確認） | — |
| 1709.01928 | v4 | Kim–Yagi, Topological vertex formalism with O5-plane（Phys. Rev. D 97 (2018) 026011） | — |
| 2012.13303 | v2 | Hayashi–Zhu, More on topological vertex formalism for 5-brane webs with O5-plane（JHEP 04 (2021) 292） | — |
| hep-th/0012136 | v1 | Sinha–Vafa, SO and Sp Chern-Simons at large N（掲載誌未確認） | — |
| 0902.0616 | v1 | Krefl–Walcher, The real topological string on a local Calabi-Yau（掲載誌未確認） | — |
