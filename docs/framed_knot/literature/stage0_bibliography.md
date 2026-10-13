# Stage 0 文献チェック（2026-10-05 調査）

ステータスの意味は次のとおり。
- 確認済み = arXiv API / INSPIRE API / Crossref で実際にレコードを確認した
- 未確認 = 記憶・推測による情報で、今回は一次ソースを確認していない

使用した API:
- arXiv: `https://export.arxiv.org/api/query?id_list=...`
- INSPIRE: `https://inspirehep.net/api/literature?q=refersto:recid:1517805`
- Semantic Scholar の citations API（補助として使用）

---

## 1. arXiv:1703.05408（Kameyama–Nawata）の書誌（確認済み）

- 題名: Refined large N duality for knots。著者: Masaya Kameyama, Satoshi Nawata
- arXiv 投稿 2017-03-15（v3）。https://arxiv.org/abs/1703.05408
- DOI: 10.1142/S0218216520410011（確認済み。arXiv / INSPIRE / Crossref の三つで一致）
- 掲載誌: **J. Knot Theory Ramifications 32 (2023), no. 06, 2041001**
  - Crossref の情報: published-online 2020-07-08、published-print 2023-05
  - DOI の "20" は online 公開の年に対応する。引用するときは「JKTR 32 (2023) 2041001」と書くのが正確
- INSPIRE recid 1517805。INSPIRE 上の被引用数は 14
- 注意: `refersto:arxiv:1703.05408` で検索するとヒットしなかった（無関係な2.5万件が返る）。`refersto:recid:1517805` を使う必要がある

## 2. 1703.05408 を引用している論文

### 2a. INSPIRE に載っている 14 件（すべて確認済み）

| arXiv | 著者 | 題名 | 掲載 |
|---|---|---|---|
| 1704.07412 | Diaconescu, Donagi, Pantev | BPS states, torus links and wild character varieties | CMP 359 (2018) |
| 1705.03467 | Kimura, Mori, Sugimoto | Refined geometric transition and qq-characters | JHEP 01 (2018) |
| 1706.00761 | Mironov, Morozov, Morozov, Sleptsov | Gaussian distribution of LMOV numbers | NPB 924 (2017) |
| 1707.02991 | Kucharski, Reineke, Stošić, Sułkowski | BPS states, knots and quivers | PRD 96 (2017) |
| 1707.04017 | Kucharski, Reineke, Stošić, Sułkowski | Knots-quivers correspondence | ATMP 23 (2019) |
| 1711.09884 | Huang, Sun, Wang | Blowup Equations for Refined Topological Strings | JHEP 10 (2018) |
| 1805.00993 | Kozçaz, Shakirov, Vafa, Yan | Refined Topological Branes | CMP 385 (2021) |
| 1810.01885 | Kimura, Sugimoto | Quantum mirror curve of periodic chain geometry | JHEP 04 (2019) |
| 1811.03110 | Ekholm, Kucharski, Longhi | Physics and geometry of knots-quivers correspondence | CMP 379 (2020) |
| 2003.07836 | Mironov, Morozov | Hopf superpolynomial from topological vertices | NPB 960 (2020) |
| 2104.00713 | Cheng, Sułkowski | Refined open topological strings revisited | PRD 104 (2021) |
| 2108.03696 | Shi Cheng | 3d N=2 brane webs and quiver matrices | JHEP 07 (2022) |
| 2312.00624 | Mariño, Rella | On the Structure of Wave Functions in Complex Chern-Simons Theory | SIGMA 22 (2026) |
| （arXiv なし） | Shi Cheng | Refined Topological Vertex and Supersymmetric Gauge Theories（学位論文, 2022） | — |

### 2b. Semantic Scholar にだけ載っていた引用（確認済み）

- 2105.02037: Shengmao Zhu, "New structures for colored HOMFLY-PT invariants"
- 1801.09363: Bai, Jiang, Liang, Mironov, Morozov, Morozov, Sleptsov, "Quantum Racah matrices up to level 3 and multicolored link invariants"（J. Geom. Phys. 132 (2018)）

### 2c. 著者本人による後続論文（確認済み）

- 1902.02275: Kameyama, Nawata, Tao, Zhang, "Cyclotomic expansions of HOMFLY-PT colored by rectangular Young diagrams"
  - 掲載: Lett. Math. Phys. 110 (2020) 2573（DOI 10.1007/s11005-020-01318-5）
  - INSPIRE の引用リストには入っていない（参照文献が INSPIRE に登録されていないためと思われる）
  - 内容は**非 refined**。長方形の色の HOMFLY-PT について、double twist knot の閉公式と cyclotomic expansion を予想したもの。cable は扱っていない
- 1804.08364: Chuang, Diaconescu, Donagi, Nawata, Pantev, "Twisted spectral correspondence and torus knots"
  - (l, kl−1) 型の**トーラス結び目**の refined CS 不変量と、twisted wild character variety を扱う
- 1905.05664: Ito, Kameyama（Khovanov homology と Vassiliev 不変量）。refined BPS とは無関係
- arXiv で Kameyama / Nawata の論文を全件確認した。cable / iterated torus knot の refined BPS / refined LMOV を扱った論文は**見つからなかった**

### 2d. 非トーラス結び目の refined BPS / LMOV に関係する引用論文（abstract は確認済み）

- **1707.02991 / 1707.04017（KRSS, knots–quivers）**
  - 6 交点までの全結び目、いくつかの無限族（twist knot など）、thick な結び目に quiver を対応させる
  - その結果、**対称表現 S^r** の LMOV 整数性（refined を含む motivic DT）が自動的に示される
  - 新しい superpolynomial の公式も導いている
  - 扱うのは対称表現だけ。cable と [2,1] は扱っていない
- **1811.03110（Ekholm–Kucharski–Longhi）**
  - knots–quivers 対応の幾何的・物理的な基礎づけ（3d N=2 理論、basic disks）
  - 「(doubly) refined LMOV 不変量」も論じている。対称表現の範囲
- **2104.00713（Cheng–Sułkowski）**
  - toric CY 上の refined open BPS 不変量が非負整数であることを、strip 型幾何などで示す
  - 結び目としてはトーラス・unknot 系の幾何に限られ、cable は扱っていない
- **1706.00761（MMMS）**
  - 多数の結び目と表現（非対称表現を含む）で、非 refined の LMOV 数の分布を調べている
  - refined 版ではない
- **2003.07836（Mironov–Morozov）**
  - Hopf link の色付き superpolynomial を、DIM の4点関数から正の形で構成する。結び目ではなく link
- **2105.02037（Zhu）**
  - 色付き HOMFLY-PT の strong integrality を skein 理論で証明。LMOV 構造の精密化（非 refined）

**結論（確認済みの範囲で）**
- 1703.05408 を引用する論文の中に、cable / iterated torus knot、または [2,1]・[3]・[1,1,1] の色について、**refined BPS（refined LMOV）不変量やその正値性を計算したものはない**
- 非トーラス結び目の refined な結果は、knots–quivers 系（対称表現・twist knot など）に限られている

## 3. arXiv 番号の検証（すべて確認済み）

全件で、想定していた著者と実際の著者が一致した。

| arXiv | 実際の第一著者（共著者） | 題名 | 年 | 掲載 / DOI |
|---|---|---|---|---|
| 1111.6195 | Cherednik | Jones polynomials of torus knots via DAHA | 2011 | arXiv に journal-ref なし（IMRN 2013 と記憶しているが**未確認**） |
| 1704.07630 | Mellit | Homology of torus knots | 2017 | Geom. Topol. 26 (2022) 47–70, 10.2140/gt.2022.26.47 |
| 1603.00407 | Elias (Hogancamp) | On the computation of torus link homology | 2016 | Compositio 155 (2019) 164–205 |
| 1909.00418 | Hogancamp (Mellit) | Torus link homology | 2019 | journal-ref なし |
| 1410.0859 | Morton (Samuelson) | The HOMFLYPT skein algebra of the torus and the elliptic Hall algebra | 2014 | Duke 166 (2017) 801–854 |
| 1201.3339 | Mironov (Morozov, Shakirov, Sleptsov) | Interplay between MacDonald and Hall-Littlewood expansions of extended torus superpolynomials | 2012 | JHEP 05 (2012) 070 |
| 1003.1568 | Oblomkov (Shende) | The Hilbert scheme of a plane curve singularity and the HOMFLY polynomial of its link | 2010 | Duke 161 (2012) 1277 |
| 1201.2115 | Oblomkov (Rasmussen, Shende) | The Hilbert scheme of a plane curve singularity and the HOMFLY homology of its link | 2012 | Geom. Topol. 22 (2018) 645 |
| 2307.16794 | Galashin (Lam) | Monotone links in DAHA and EHA | 2023 | journal-ref なし |
| 2407.18123 | Caprau (González, Hogancamp, Mazin) | Khovanov-Rozansky homology of Coxeter knots and Schröder polynomials for paths under any line | 2024 | journal-ref なし |
| 1504.02620 | Arthamonov (Shakirov) | Refined Chern-Simons Theory in Genus Two | 2015 | JKTR 29 (2020) 2050044 |
| 1105.5117 | Aganagic (Shakirov) | Knot Homology from Refined Chern-Simons Theory | 2011 | CMP 333 (2015) 187（INSPIRE） |
| 1202.2489 | Aganagic (Shakirov) | Refined Chern-Simons Theory and Knot Homology | 2012 | Proc. Symp. Pure Math. 85 (2012) 3（INSPIRE） |
| 1106.4305 | Dunin-Barkowski (Mironov, Morozov, Sleptsov, Smirnov) | Superpolynomials for toric knots from evolution induced by cut-and-join operators | 2011 | JHEP 03 (2013) 021 |
| 1304.3328 | Gorsky (Neguţ) | Refined knot invariants and Hilbert schemes | 2013 | J. Math. Pures Appl. 104 (2015) 403 |
| 1408.4348 | Cherednik (Danilenko) | DAHA and iterated torus knots | 2014 | AGT 16 (2016) 843 |
| 1509.08351 | Cherednik (Danilenko) | DAHA approach to iterated torus links | 2015 | journal-ref なし |
| 2410.14703 | Cherednik | Superpolynomials of algebraic links | 2024 | DOI 10.1093/qmath/haaf035（Quart. J. Math.） |
| 1605.00978 | Cherednik (Philipp) | DAHA and plane curve singularities | 2016 | AGT 18 (2018) 333 |
| 1211.4596 | Anton Morozov（単著） | The first-order deviation of superpolynomial in an arbitrary representation from the special polynomial | 2012 | JETP Lett. 97 (2013) 171 |
| 1211.6375 | Anokhina (Mironov, Morozov, Morozov) | Knot polynomials in the first non-symmetric representation | 2012 | NPB 882 (2014) 171 |
| 1203.5978 | Itoyama (Mironov, Morozov, Morozov) | HOMFLY and superpolynomials for figure eight knot in all symmetric and antisymmetric representations | 2012 | JHEP 07 (2012) 131 |
| 1205.1515 | Fuji (Gukov, Sułkowski) | Super-A-polynomial for knots and BPS states | 2012 | NPB 867 (2013) 506 |
| 1209.1409 | Nawata (Ramadevi, Zodinmawia, Sun) | Super-A-polynomials for Twist Knots | 2012 | JHEP 11 (2012) 157 |
| 1609.00143 | Kononov (Morozov) | Rectangular superpolynomials for the figure-eight knot | 2016 | TMP 193 (2017) 1630 |
| 1310.2240 | Nawata (Ramadevi, Zodinmawia) | Colored Kauffman Homology and Super-A-polynomials | 2013 | JHEP 01 (2014) 126 |
| math/0505662 | Dunfield (Gukov, Rasmussen) | The Superpolynomial for Knot Homologies | 2005 | Exp. Math. 15 (2006) 129 |

補足: 1211.6375 は [2,1] 色の HOMFLY を扱う論文、1211.4596 は任意の表現での superpolynomial の一次のずれを扱う論文で、どちらも [2,1] 関連の文献として使える。

## 4. (2,1) 色の三葉結び目ホモロジーを提案した論文（確認済み）

- **arXiv:1304.3481**: "Quadruply-graded colored homology of knots"
- 著者は **Gorsky, Gukov, Stošić の3人**。想定の「Gukov–Stošić」とは食い違う（Gorsky が加わっている）
- 掲載: Fundamenta Math. 243 (2018), no. 3, 209–299（DOI 10.4064/fm30-11-2017, INSPIRE）
- PDF 本文で次を確認した
  - **§4.6 "(2,1) ('hook') representation"** で、(2,1) 色の reduced HOMFLY を示している
  - そのうえで41生成元の四重次数付きホモロジーを提案している（式 (4.14)）。dcol 型の微分（d^col_{2|0}, d^col_{1|1}, d^col_{0|2}）でほぼ決まる
  - 長方形の場合と違い、exponential growth が成り立たないことも指摘している
  - 長方形と非長方形の違いは §3.11 で論じている
  - 目次: 1. Introduction, 2. Bosonic and Fermionic generators, 3. Gradings, Differentials, Mirror Symmetry, 4. Examples, 5. Colored HOMFLY homology for torus knots, 6. HOMFLY homology from A-model and B-model, 7. …
- 参考: Gukov–Stošić の2人だけの論文は 1112.0030 "Homological algebra of knots and BPS states"（GTM 18 (2012) 309）。こちらは対称表現と反対称表現の colored differential を扱う論文

## 5. arXiv:2610.00303（確認済み）

- 題名: "Rectangular Cyclotomic Expansions of Colored HOMFLY-PT Polynomials"
- 著者: Honghuai Fang, Tian Zhou
- 投稿: 2026-09-28（math.QA / math.GT, 32 pp）。https://arxiv.org/abs/2610.00303
- 投稿日は 9/28。「Oct 2026」というのは番号が 2610 であることによる
- 要旨
  - Kameyama–Nawata–Tao–Zhang（1902.02275）の長方形 cyclotomic expansion 予想を、任意の zero-framed 結び目の reduced colored HOMFLY-PT について証明した
  - double twist knot の明示公式（任意の長方形色、任意の非零ツイストパラメータの組）も証明した
  - 各結び目に Z[A^±, q^±] 係数の double Newton expansion を構成した。係数は長方形の幅と高さに依存しない
  - 証明では four-point Casimir の exterior-power 比較を使う
- **refined（t 変形）不変量や cable には触れていない**（abstract から判断。本文は未確認）

## 6. 「refined Rosso–Jones / cabling formula」に関する既存文献

| arXiv | 著者 | 内容（一行） | 状態 |
|---|---|---|---|
| 1408.4348 | Cherednik–Danilenko | DAHA-Jones / DAHA-superpolynomial をトーラス結び目から**任意の iterated torus knot**（任意のルート系・任意の重み λ）に拡張。多項式性・duality を示し、Puiseux 展開が自然に現れる。KR / ORS との一致は予想 | 確認済み |
| 1509.08351 | Cherednik–Danilenko | 上の構成を iterated torus **links**・代数的 link に拡張（Macdonald 多項式の積と splice diagram を使う）。Hopf link では「DAHA-vertex」が得られる | 確認済み |
| 1605.00978 | Cherednik–Philipp | 非トーラスの代数的結び目について、非色付き DAHA superpolynomial を flagged Jacobian factor で幾何的に記述する予想 | 確認済み |
| 2410.14703 | Cherednik | motivic superpolynomial の理論。代数的 link（行で色付け）、iterated torus link の漸化式、「DAHA = motivic」予想 | 確認済み |
| 1408.0483 | Samuelson | 球面 DAHA SH_{q,t}（A_1）の位相的実現。Cherednik の A_1 構成を全ての iterated cable of the unknot に一般化し、colored Jones への特殊化を証明。t=−q² で CD14 と比較 | 確認済み（IMRN, DOI 10.1093/imrn/rnx198） |
| 1410.0859 | Morton–Samuelson | 2次元トーラスの HOMFLYPT skein 代数 ≅ EHA の t=q 特殊化。iterated cable of the unknot に対して EHA から3変数多項式を構成し、λ 色 HOMFLY と CD14 の両方に特殊化することを示す（CD の Connection Conjecture の一つを証明） | 確認済み |
| 2307.16794 | Galashin–Lam | monotone link（Coxeter link と同じクラス。代数的 link を含む）に DAHA / EHA の元を対応させる。凸の場合は Schur 正値性を予想 | 確認済み |
| 1307.2216 | Anokhina–An. Morozov | **非 refined** の cabling procedure（色付き HOMFLY を cabled knot の基本表現 HOMFLY と射影子から得る） | 確認済み（TMF 178 (2014)） |
| 1405.0884 | Danilenko | 上の cabling 関係を Khovanov–Rozansky に持ち上げようとした論文。unknot の [1,1] と [2] で、**素朴な t 変形は存在しない**ことを報告 | 確認済み |
| 1801.02407 | A. Morozov | twist satellite（cabling の一般化）の色付き HOMFLY の分解則。分解則を **superpolynomial に持ち上げる際の問題**（正と負のツイストで規則が異なる）を指摘 | 確認済み（PLB 782 (2018)） |
| 1304.3328 | Gorsky–Neguţ | Aganagic–Shakirov の refined トーラス不変量を DAHA / Hilbert scheme の観点で扱う。**トーラスのみ** | 確認済み |
| 1111.6533 | Diaconescu–Shende–Vafa | 代数的結び目の Lagrangian cycle と large N 双対性（非 refined の HOMFLY / stable pair） | 確認済み |

### わかっていること（確認済みの abstract に基づく）

- iterated torus knot / link と代数的 link の refined 不変量は、**Cherednik–Danilenko の DAHA superpolynomial として定義されている**
  - 任意の色 λ（Macdonald 多項式の重み）に対して、DAHA の自己同型と評価を反復する**アルゴリズム**で与えられる
  - t=q で HOMFLY-PT に一致する（予想。iterated cable of the unknot については Morton–Samuelson が証明した）
  - 正値性と KR homology との一致は予想
- Morton–Samuelson の EHA による方法は、t=q の skein 代数上で cabling を扱う。一般の t では、CD と一致する3変数多項式を与える
- 非 refined の Rosso–Jones / cabling 公式（Anokhina–Morozov など）は確立している
- refined への持ち上げでは、素朴な t 変形がうまくいかないという報告がある（Danilenko 1405.0884、Morozov 1801.02407）

### わかっていないこと（今回の検索では見つからなかった。網羅的とは言えない）

- 「cable の refined 不変量を、companion の色付き refined 不変量（superpolynomial）の Macdonald 係数付き一次結合として書く」という Rosso–Jones 型の**明示公式**
  - DAHA / EHA の構成は内部で Macdonald 展開を使うので、実質的にはこれに近い
  - しかし、companion の「色付き superpolynomial」だけを入力とする閉公式として述べた論文は見つからなかった
- cable / iterated torus knot について、refined BPS（refined LMOV）不変量を抽出し、その正値性・整数性を調べた論文（1703.05408 の枠組みで）
- [2,1]・[3]・[1,1,1] 色の、cable の refined 不変量
- 注意: arXiv の全文検索はキーワード依存なので取りこぼしがありうる
  - Gorsky–Neguţ 系の EHA 論文の後続、Mellit / Hogancamp 系の cable に関する KR 計算は個別に確認していない（**未確認**）

## 参考 URL

- INSPIRE 被引用一覧: https://inspirehep.net/api/literature?q=refersto%3Arecid%3A1517805&size=100
- arXiv API（番号検証）: https://export.arxiv.org/api/query?id_list=1111.6195,1704.07630,...
- GGS §4.6: https://arxiv.org/pdf/1304.3481
- 2610.00303: https://arxiv.org/abs/2610.00303
