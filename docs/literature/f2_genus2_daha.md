# 種数 2 の refined CS・DAHA・skein 代数の文献（2026-10-08、ユーザーの指摘による追加調査、open_problems.md §F2）

ユーザーの指摘：Arthamonov–Shakirov と Hikami が種数 2 の曲面上の代数と結び目を計算している。以前の調査では Arthamonov–Shakirov 1504.02620 を 1 行挙げただけで（novelty_report.md §2.4、【調査報告】）、原典は確認していなかった。Hikami は調べていなかった。

主張の印：【原典】TeX 本文で確認（行番号は references/<id>/src の TeX）、【要旨】arXiv の要旨だけで確認、【計算】この調査での確認、【見通し】未検証の推論。
文献情報は arXiv API のメタデータで確認した（DOI・掲載誌は arXiv に記載のあるものだけ書く）。

## 1. 中心の 3 本【原典】

### Arthamonov–Shakirov, "Refined Chern-Simons Theory in Genus Two"（arXiv:1504.02620、J. Knot Theory Ramifications 29 (2020) 2050044）
- 内容：rank 1（U_q(sl_2)）の refined TQFT の写像類群の表現を種数 2 に拡張（highergenus.tex :53-55）。q-6j 記号そのものは Macdonald 変形を持たないように見えるが、その 2 乗は変形を持つ（:75-82）。
- 結び目：種数 2 の曲面上の結び目（3 指数の pretzel 結び目、4_1, 6_1, 6_2, 8_1, 8_4, 8_2, … を含む）の refined Jones 多項式を計算（:528-）。
- 不変性：振幅 Z そのものは結び目の不変量ではない（4_1 と鏡像で値が違う、:520-525）。|Z|² = Z(q,t)Z(1/q,1/t) が不変量であると予想（Conjecture II、:441）。
- HOMFLY-PT ホモロジーの Poincaré 多項式とは一致しない（要旨、:66-68）。著者は refined の値を index と解釈する。
- **高ランクへの障害**（Discussion「Higher rank」、:739-758）：種数 2 以上では状態の基底が 3 価のグラフになり、A_n 型では頂点に**テンソル積の多重度**の装飾が要る。多重度の扱いは技術的に非常に難しく、対称・反対称表現に限る通常の手も写像類群の方法では使えない。
  → 本研究の X2 で一意性を妨げたもの（[2,1] ⊗ [2,1]̄ の多重度 2 と gauge の自由度）と同じ壁。

### Arthamonov–Shakirov, "Genus Two Generalization of A_1 spherical DAHA"（arXiv:1704.02947、DOI 10.1007/s00029-019-0447-1）
- 3 変数の差分作用素と x_{ij} + x_{ij}^{-1} の積で生成される代数を作り、種数 2 の閉曲面の写像類群の全ての関係を外部自己同型として満たすことを証明（paper-v3.tex の要旨）。A_1 spherical DAHA の種数 2 版。
- 未解決として明記（Conclusion :643-692）：種数 g > 2、**高ランク N > 2**、非 spherical 版、種数 2 の結び目の q,t 不変量の系統的な計算（続編の予定）。

### Hikami, "DAHA and skein algebra of surfaces: double-torus knots"（arXiv:1901.02743、DOI 10.1007/s11005-019-01189-5）
- A_1 型と C^∨C_1 型の rank 1 の DAHA を「貼り合わせ」、種数 2 の曲面の Kauffman bracket skein 代数の表現を作る（glueDAHAnew02.tex の要旨）。
- double-torus 結び目（種数 2 の Heegaard 曲面上の単純閉曲線、ツイスト結び目・ダブルツイスト結び目を含む）に DAHA 多項式を提案し、主にダブルツイスト結び目で計算（:8264-8290）。
- パラメータ x_u = x_d = q_u で彩色 Jones 多項式に一致すると予想（Conjecture、:7366-7372）。C^∨C_1 の 4 つのパラメータは特定の値 𝐭_⋆ に固定される（:6155-6180、穴の貼り合わせ）。
- 導入（:355-372）：Cherednik の DAHA 多項式はトーラス結び目とその子孫（全て非双曲）に限られ、Arthamonov–Shakirov の種数 2 の試みは既知の量子不変量との関係が不明、と述べる。
- 将来の課題として、非 reduced の DAHA 多項式・高ランク・高い種数を挙げる（:8291-8296）。

## 2. 関連論文（上の 3 本を引用するもの、Semantic Scholar の引用一覧から選んだ）【要旨】
| arXiv | 著者 | 内容 | 本研究との関係 |
|---|---|---|---|
| 2008.12695 | Cooke–Samuelson | 種数 2 の skein 代数の handlebody の skein module への作用。Arthamonov–Shakirov の種数 2 spherical DAHA の t = q の特殊化と同型 | t = q での正しさの裏付け（sl_2） |
| 2309.01011 | Arthamonov | 種数 2 の q,t 変形代数（[AS19]）の平坦性、語の問題、古典極限 | 代数の基礎（sl_2） |
| 2402.09106 | Hikami | 特定のパラメータでの C^∨C_1 型 DAHA の一般化（Hecke 作用素の「Heegaard 双対」）、種数 2 の skein 代数との関係、Dehn twist の自己同型 | 1901.02743 の続編（sl_2） |
| 2408.17021 | Hikami | 2 点穴あきトーラスの skein 代数の DAHA、クラスタ代数による古典極限 | 同上 |
| 2507.01665 | Hikami | 種数 2 の skein 代数の q 差分作用素から、reduced Askey–Wilson 多項式と種数 2 の handlebody の skein module の対応 | 同上 |
| 2205.04414 | Cooke–Lacabanne | 「高ランク」の Askey–Wilson 代数 = (n+1) 点穴あき球面の Kauffman bracket skein 代数 | 「高ランク」は穴の数の意味で、sl_N（N ≥ 3）ではない |
| 2511.21014 | Matson–Samuelson | stated skein 代数と SL_2 の universal spherical DAHA | sl_2 |
| 1711.07605 | Cherednik–Philipp | 平面曲線特異点の加群（任意のランク）と、列で彩色した代数的結び目の DAHA superpolynomial | 代数的結び目のみ |
| 1905.00208 | Awata–Kanno–Mironov–Morozov（Nucl. Phys. B 949 (2019) 114816） | 位相的頂点で tangle 計算をする計画の Macdonald 変形。合成表現の Hopf 多項式と位相的頂点の畳み込みの対応が refined では壊れる | F1（refined の Racah 行列・tangle 計算）の困難の具体例 |
| 1801.02407 | A. Morozov（Phys. Lett. B 782 (2018) 104-111） | twist satellite の彩色 HOMFLY-PT の分解。superpolynomial への持ち上げは正と負のねじれで規則が違うように見える | W0 の調査 E で「arXiv 番号・本文未確認」とした論文。番号と掲載誌はここで確認。E が引いた「非長方形の hyperpolynomial は正でない」という記述は要旨には無い（本文は未確認） |

## 3. 本研究の値との比較【計算】
- Arthamonov–Shakirov の 4_1, 6_1, 8_1 の refined Jones 多項式（j = 1、1504.02620 の表、:540-570）と、本研究の rCS_{[1]}（KNTZ の refined 式を KN の変数に読み替えたもの、conventions.md §16）を a = t²（sl_2）に特殊化したものを比べた。
  4_1 はどちらも 5 項で係数の符号の並び（+, −, +, −, +）も同じだが、指数の並び方が違う。単純な変数の対応（q, t の ±1, ±2 乗、qt, q/t など 13 種の組）で単項式倍の一致は無かった（6_1, 8_1 も同じ）【計算、探した対応の範囲で】。
- 原典の「種数 2 の refined Jones は HOMFLY-PT ホモロジーの Poincaré 多項式と合わない」（要旨）と整合する。本研究の rCS_{[1]} は DGR 型（ホモロジー由来）の値なので、種数 2 の refined は別の変形である可能性が高い【見通し】。
  ただし Arthamonov–Shakirov の振幅は |Z|² だけが不変量（予想）なので、比べるなら |Z|² どうしで、正規化の自由度も含めて調べる必要がある（未実施）。

## 4. 結論（open_problems.md §F2 への反映）
- 種数 2 の DAHA・refined CS・skein 代数の系統（Arthamonov–Shakirov、Hikami、Cooke–Samuelson など）は、**双曲結び目（ツイスト結び目・ダブルツイスト結び目・pretzel 結び目）に q,t の値を与える既存の構成**で、F2 の最も具体的な拡張元である。
- ただし全て **rank 1（sl_2、Kauffman bracket）** で、HOMFLY-PT の a は無く、彩色は sl_2 の表現（KN の対称表現 [r] の a = t² の特殊化にあたる）に限られる。sl_2 では [2,1] は [1] と区別がつかない。
- [2,1] に届くには sl_N（N ≥ 3）への高ランク化が要る。障害は Arthamonov–Shakirov が明記したテンソル積の多重度で、本研究の X2 で一意性を妨げたものと同じ。どの論文も未解決として挙げている。
- 本研究から見た位置づけ：(i) sl_2 の部分は「双曲結び目の refined 不変量の定義に基づく候補」として、KN 型の値（a = t² の特殊化）と比べる価値がある。ただし 3 の予備的な比較では単純な対応は無い。
  (ii) [2,1] の問題は、高ランク化の多重度の問題に帰着し、文献でも未解決。
