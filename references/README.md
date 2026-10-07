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
| 1112.0030, 1209.1409, 1512.07883 | 最新版 | 第 2 期以降に取得（GS、NRZS、GNSSS） | — |

DOI・巻・頁は未確認（上の表で確認したと書いたものを除く）。

派生ファイル：`1703.05408/nb-torus.txt`, `1703.05408/nb-non-torus.txt` は ancillary ノートブックを
`scripts/nb2txt.py` でテキスト化したもの。再生成するには次を実行する：

```bash
.venv/bin/python scripts/nb2txt.py references/1703.05408/src/anc/refined-LMOV-torus.nb > references/1703.05408/nb-torus.txt
```
