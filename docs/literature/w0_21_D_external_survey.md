> 出所：ユーザーが別の AI に調査させた結果（2026-10-08 にユーザーが docs/[2,1]文献調査.md として置いたもの）を、名前を変えてここに移した。本文は元のまま。
> 文献情報の確認状況は本文末尾の一覧のとおりで、本研究のセッションでは原典を確認していない。W0 の他の調査（A・B・C）との統合は docs/literature/w0_21_survey.md。

# 双曲結び目の [2,1] 彩色 refined 不変量：構成の文献調査（2026年10月時点）

結論を先に書くと、双曲結び目の [2,1] について、a, q, t の 3 変数がそろった不変量を定義どおりに与える構成は見つからなかった。最も近いのは次の 2 つで、どちらも t を含まない。
- Gruen–San Martín Suárez（arXiv:2508.15171）の sl_3 の F_K。sl_3 のすべての既約表現（随伴 (1,1) を含む）を扱い、正ブレイド結び目について定義されている。\[1\]
- Beliakova–Gorsky（arXiv:2101.08243）の gl_N 普遍不変量の円分展開。任意の結び目と任意の彩色について定義されている。

t 依存をもつ構成は、対称表現に限られるか、トーラス型（およびその links）に限られる。

## TL;DR
- **定義どおりの構成は「t なし」に限られる**：sl_3 の F_K（2508.15171）は、正ブレイド結び目のすべての sl_3 既約表現を扱う。ただし著者は §1 で「三葉結び目と cinquefoil（5_1）だけを明示計算した」と書いており、変数は x, y, q である（a も t もない）。Beliakova–Gorsky は任意の結び目に対して gl_N（N ごと）の円分展開を与え、補間 Macdonald 多項式の t は基底にだけ現れる。t をホモロジー次数とみなす解釈は、論文の中で「期待」として述べられているにとどまる。
- **t を含む構成は対称表現かトーラス型に限られる**：(a,q,t) 変形 F_K（2005.13349、2110.13768）は対称彩色のみで、著者自身が t 変形が一意でないことを認めている。Macdonald 変形した随伴不変量（2505.16569）はトーラス links のみである。補助のサブエージェントの探索でも、2014〜2026年10月の範囲で、非トーラス結び目の [2,1] を t 依存で与えた論文は 1211.6375 以外に見つからなかった。
- **実務上の推奨**：双曲結び目の [2,1] を検証するには、次の 3 つを拘束条件として組み合わせるのが現実的である。(1) Gu–Jockers（1407.5643）などが与える、定義どおりの [2,1] HOMFLY-PT（t=−1 の極限）。(2) sl_3 の F_K から得られる N=3 の随伴不変量（10_139 は正ブレイド結び目なので適用範囲に入る。ただしこれは筆者の推論）。(3) KN が BPS の段で予想した正値性。ただし (3) を ansatz を選ぶ条件に使うと、正値性予想の検証としては循環論法になる。

## Key Findings

### 系統 1：3d–3d 対応、Ẑ / F_K とその高ランク版・精密化

| 構成 | (i) 双曲結び目での定義と計算例 | (ii) 彩色、a/N 依存 | (iii) 量の型 | (iv) 正値性・整数性の主張 | (v) KN に入力できるか |
|---|---|---|---|---|---|
| Gukov–Manolescu F_K（1904.06057） | 予想として一般の結び目。4_1 は数値計算 | sl_2 のみ（x, q）。a なし | q 級数（係数は符号つきの整数） | 係数が整数であること（BPS 状態の数え上げとして） | 不可（a, t なし。彩色は sl_2 の大色極限）\[2\] |
| Park の inverted state sum（2106.03942） | 4_1、β5、β6 などの双曲結び目で明示計算\[3\] | sl_2 のみ | q 級数（index 型） | 整数性 | 不可 |
| Park 高ランク Ẑ, F_K（1909.13002） | F_K^G を予想として定義。明示式はトーラス結び目のみ\[4\] | 任意の半単純 g（sl_3 なら (1,1)=[2,1] 方向を含む）。N は固定 | q 級数（Weyl 反対称化した形） | 係数が整数の Laurent 級数（予想） | 不可（a, t なし）\[1\] |
| Gruen–San Martín Suárez F_K^{sl_3}（2508.15171） | **正ブレイド結び目のすべて**で定理として定義。計算例は三葉結び目と 5_1 のみ\[5\] | sl_3 の全既約表現（x=q^n, y=q^m で V_(n,m) の RT 不変量を復元）。N=3 に固定 | Z[q,q^{-1}][[x,y]] の級数（負の係数あり） | 整数係数（定理） | 不可（a, t なし）。N=3 での照合用\[1\] |
| Gruen の sl_N 対称大色 R 行列（2212.05222） | 正ブレイド結び目など | sl_N の対称表現のみ | q 級数 | 整数性 | 不可\[1\] |
| EGGKPS「Ẑ at large N」（2005.13349） | 4_1 の t 変形 F_K を super-A 多項式の漸化式から項ごとに計算\[6\] | 対称表現のみ。a 依存あり | 級数。(a,t) 変形は superpolynomial に由来 | a=−t^{-1} で t 変形 Alexander 多項式\[6\] | 変数 a, q, t はそろうが対称彩色のみ\[7\] |
| EGGKPSS「Branches, quivers, and ideals」（2110.13768） | 4_1、5_2 の分枝ごとの F_K^{(α)}（sl_N、任意の N）。quiver 形式\[8\]\[9\] | 対称表現のみ。a あり。概要に「t 変形版も与える」とあり、§6.4「Refinement」で扱われている | 級数（quiver 生成関数） | quiver 形式から整数性 | 対称彩色のみ。t 変形は一意でない（後述） |
| Kucharski「Quivers for 3-manifolds」（2005.13394） | 結び目補空間の quiver 対応 | sl_2（対称） | 級数 | DT 不変量の正値性からの整数性 | 不可\[10\] |
| 超群 sl(2\|1) の F_K（2508.10279） | plumbed 結び目補空間。計算例はトーラス結び目のみ\[11\] | sl(2\|1) | 3 変数級数 | — | 不可 |

要点：
- **2508.15171 は、本調査で見つかった中で唯一、「双曲結び目を含む族」について非対称・非 minuscule の彩色（sl_3 の随伴 (1,1)）を定理として扱う構成**である。主定理は「正ブレイド β_K について、縮約量子トレース T_{V_x^1⊗V_y^2}(β_K) が収束し、F_K^{sl_3}(x,y,q) は well-defined な結び目不変量になる」というもので、補題 2.3 は「V_(n,0)⊗V_(0,m) で縮約トレースを取り、重み (n,m) の元で評価すると、V_(n,m) の sl_3 Reshetikhin–Turaev 不変量になる」と述べている。ただし著者は「三葉結び目と 5_1 だけを明示計算した」と書いている（三葉結び目の結果は §4.3 で Park のトーラス結び目の式と照合済み）。計算量については、§3.4 で「交点 1 つあたりの内部ラベルの数が 3 から 11 に増える」と、§5 で sl_N の場合は「交点 1 つあたり ∼N^4 個のパラメータになる」と述べている。sl_N への一般化は「予想としての枠組み」にとどまる。\[5\]
- **推論（論文には書かれていない）**：10_139 は正ブレイド結び目で、双曲結び目でもある。したがって定理の適用範囲に入り、x=q^1, y=q^1（規約によるずれを含む）に特殊化すれば sl_3 の随伴彩色の RT 不変量が得られるはずである。これは HOMFLY-PT_{[2,1]}(a=q^3) と比べる検証点になる。一方、4_1、5_2、6_1 は fibered でないか正ブレイドでないので、この定理の範囲外である（4_1 は fibered だが正ブレイドではない）。
- **2110.13768 の注意書き**（ScienceDirect 版より）：著者は「superpolynomial と super-A 多項式から受け継いだ t 変形は、ループ数に基づく t 変形と整合しない。さらに、knot complement quiver の形が違えば t 変形も違ってくる」と明記している。つまり\[12\] F_K 系の t は、現時点では定義どおりに一意に決まるものではない。
- 2512.23700（Quantum invariants and fiberedness, 2025）は、「fibered でない結び目の F_K 不変量には定義がない」と述べ、これまでの例は Park の拡張か、ツイスト結び目の色つき Jones 多項式からの補間による孤立した例だけだとしている。fibered\[13\] でない 5_2、6_1 では、sl_2 の段階でもまだ定義が確立していない。
- GPV（1602.05302）、GPPV（1701.06567）の流れで、Ẑ を圏化するホモロジーの存在は予想されている。しかし今回の調査範囲では、結び目について「非対称彩色」と「t」が同時に入った拡張は見つからなかった。\[14\]

### 系統 2：円分展開・補間 Macdonald 多項式（ホモロジーを作らない、a/N 一様化の候補）

| 構成 | (i) | (ii) | (iii) | (iv) | (v) |
|---|---|---|---|---|---|
| Beliakova–Gorsky（2101.08243、Selecta Math. 30 (2024) 101） | **任意の結び目**（定理）。gl_N の RT 不変量と普遍 gl_N 不変量の円分展開\[15\] | **任意の彩色**（[2,1] を含む）。N ごと（a に一様な形ではない） | Laurent 多項式（index 型） | 双対係数の Laurent 整数性（定理）\[16\] | 不可。ただし t は補間 Macdonald 多項式の基底パラメータとして自然に現れる |
| KNTZ（1902.02275） | ダブルツイスト結び目（予想）、任意の結び目の円分展開（予想）\[17\] | 長方形のみ。a あり | HOMFLY-PT と Poincaré 多項式（Mathematica 付属）\[18\] | — | 長方形のみ |
| Chen–Liu–Zhu（2110.03616） | ダブルツイスト結び目 | 概要未確認 | — | — | 未確認 |
| Cyclotomic Newton expansions（2608.29585, 2026） | Beliakova–Gorsky の補間公式を rank 一様化 | gl_n | 整数値 | 整数性 | 未確認（概要のみ） |
| SU(n) 円分展開と Bailey 変換（2609.13766, 2026） | 2 本ブレイドのトーラス結び目のみ | — | — | — | 不可 |

要点：Beliakova–Gorsky は論文中で、「圏化された補間多項式は、q を量子次数、t をホモロジー次数とする補間 Macdonald 多項式に（変数変換の後で）対応すると期待する」と書いている。**これは期待の表明であり、\[19\]定理ではない**（障害として「斉次でない」ことが挙げられている）。それでも、任意の結び目と任意の彩色について「t が入りうる基底」を定義どおりに与えている点で、検証できる原理に最も近い。基底の側では、命題 11.5 で A=q^N に依存する補間多項式の HOMFLY-PT 版がすでに与えられている。未解決なのは、展開係数を a について一様化すること（N ごとではなく HOMFLY-PT 型にすること）と、展開係数の t 変形を決めることの 2 つである。

### 系統 3：双曲結び目の [2,1] HOMFLY-PT（t なしだが、定義どおりで検証用の基礎になる）

| 構成 | (i) | (ii) | (iii) | (iv) | (v) |
|---|---|---|---|---|---|
| Gu–Jockers（1407.5643、CMP 338 (2015)） | **8 交点までの 2 橋双曲結び目**（4_1、5_2、6_1 などを含む）\[20\] | {2,1}。a 依存あり（量子 6j 記号を一般の N で）\[20\] | HOMFLY-PT（index 型） | — | t がないので t=−1 の極限としてのみ使える |
| Mironov–Morozov–Morozov（1408.3076） | ツイスト結び目の族全体（evolution 法）。k\[21\] に依存しない 6 個の未知係数 u_α(A,q) を、既知の 7 個のツイスト結び目の H_[21] から決めている（Anokhina–Morozov 1307.2216 の 4 個と Gu–Jockers の 3 個） | [21]。a あり | HOMFLY-PT | — | 同上 |
| Mironov–Morozov–Morozov–Ramadevi–Singh（1504.00371） | double fat 図式で表せる結び目（pretzel 型を含む）。予想式を 3 本ブレイドの値と照合\[22\] | [21] を含む | HOMFLY-PT | — | 同上\[23\] |
| Mironov–Morozov 普遍 Racah（1511.09077） | arborescent 結び目 | 随伴（普遍形） | 随伴多項式 | — | 同上\[24\] |
| Qin（2101.00695） | pretzel (3,3,3)、(3,3,−3)\[25\] | [2,1] | HOMFLY-PT | — | 同上 |

要点：t=−1 での [2,1] HOMFLY-PT は、4_1、5_2、6_1、ツイスト結び目、pretzel について定義に基づく値が得られている。arXiv:1407.5643v3 の概要ページには、「8 交点までの 2 橋双曲結び目について、表現 {2,1} で彩色した HOMFLY 多項式を明示的に与える」と書かれている。ただし MMM（1408.3076）自身が「R=[21] では、特定の結び目や特定の N でさえ計算が既存の計算機の能力の境界にある」と書いているとおり、計算は重い。\[21\]

### 系統 4：t 変形（Macdonald 変形・evolution）の新しい試み

| 構成 | (i) | (ii) | (iii) | (iv) | (v) |
|---|---|---|---|---|---|
| Bishler–Mironov–Morozov「Macdonald deformation of Vogel's universality and link hyperpolynomials」（2505.16569, 2025） | **トーラス links T[2,2n] のみ** | 随伴（随伴の 2 乗を 6 つの既約表現に分解）。ADE 型に限る | hyperpolynomial（t 依存） | 普遍性の予想 | 双曲結び目は不可\[26\] |
| Mironov–Singh「Torus Knots in Adjoint Representation」（2512.23095, Phys. Rev. D 114 (2026) 046004） | トーラス結び目のみ | 随伴は N 依存の composite [2,1^{N-2}] であり、固定した [2,1] ではない\[27\] | HOMFLY-PT（t なし） | — | 不可\[28\] |
| Anokhina–Morozov–Popolitov（2104.14491）、Anokhina–Lanina–Morozov（2308.13095） | ツイスト結び目など | 「adjoint」は sl_2 の [2] を指し、[2,1] ではない\[29\] | Khovanov 多項式（t あり） | — | 不可 |

要点：Morozov 系以外では、ツイスト結び目や pretzel の [2,1] に t 依存の式を与えた論文は**見つからなかった**。A. Morozov（1208.3544, 2012）の「R=[2,1] と R=[2,2] という最も単純なフック図式と非フック図式について、少なくとも一部の superpolynomial を見つけることが極めて重要」という問題提起は、2026年10月の時点でも未解決のままである。\[30\]

### 系統 5：ホモロジー側（sl_3 foam、Webster ほか）

| 構成 | (i) | (ii) | (iii) | (iv) | (v) |
|---|---|---|---|---|---|
| Webster（1001.2020、1309.3796） | 定義上は任意の結び目。明示計算は少数 | 任意の g、任意の表現（[2,1] を含む）。N 固定 | ホモロジー。**非 minuscule 表現では、彩色 unknot のホモロジーが無限次元** | Poincaré 級数は非負だが多項式ではない | 不可（a なし。有限の Poincaré 多項式にならない）\[31\] |
| 圏化 sl_3 射影子と sl_3 RT 不変量（1109.1745） | 定義上は任意のタングル。計算例は簡単なもの | sl_3 の任意の既約表現（(1,1) を含む） | 複体（K^∠ の中）。無限複体 | — | 不可\[32\] |
| 彩色 sl_3 不変量の圏化（1503.08451） | 定義のみ（関手性は未解決）\[33\] | V_(m,n)（随伴を含む）\[33\] | ホモロジー | — | 不可 |
| Lie 理論的な彩色 Jones の圏化（2109.12889） | unknot の計算で、無限に多くの次数に非零ホモロジー | sl_2（sl_k にも拡張） | 非有界導来圏 | — | 不可\[34\] |
| J. Wang（2211.08409） | 三葉結び目と Hopf link（「非自明な結び目について初の完全計算」）\[35\] | 外積冪（minuscule）のみ\[35\] | ホモロジー | SU(N) 表現空間のコホモロジーと同型\[35\] | 不可 |
| Guicardi（2609.10941, 2026） | Verma 加群の KLRW 圏化（ブレイド補空間上の Khovanov）\[36\]\[37\] | sl_2 のみ | 圏化した R 行列\[37\] | — | 不可\[38\] |

要点：sl_3 foam や Webster の構成で、[2,1] 型（sl_3 の (1,1)）の彩色ホモロジーを双曲結び目について計算した例は**見つからなかった**。Webster のホモロジーは、彩色 unknot ですら非 minuscule の色では無限次元になる。この事実は、[2,1]\[31\] の Poincaré「多項式」が直接は得られないことと整合する。

**推測**：トーラス結び目の rCS_{[2,1]} が負の係数をもつことと、この無限次元性（Euler 標数をとる際の打ち消し）とは関係している可能性がある。ただし、これを示した文献は見つからなかった。

### 系統 6：計算ソフト

| ソフト | 対象 | [2,1] に対応するか |
|---|---|---|
| khoca（Lewark–Lobb） | 任意の links の sl_2 Khovanov ホモロジー。N>2 の sl(N) は 2 部結び目のみ。いずれも無彩色\[39\] | 非対応 |
| FoamHo（Lewark, 1212.2553） | sl_3 foam ホモロジー（無彩色、Morrison–Nieh のアルゴリズム）\[40\] | 非対応 |
| fkcompute（2607.12155, 2026） | Park の inverted state sum による sl_2 の F_K\[3\]\[41\] | 非対応（sl_2 のみ） |
| 彩色 Khovanov のデータベース（2505.03916, 2025） | sl_2 の対称彩色（khoca を利用）\[42\] | 非対応 |
| KNTZ の Mathematica 付属ファイル（1902.02275） | ダブルツイスト結び目の長方形彩色 HOMFLY-PT と Poincaré 多項式\[18\] | 非対応（長方形のみ） |

非対称表現で彩色した sl(N) や HOMFLY-PT のホモロジーを計算する公開ソフトは**見つからなかった**。2508.15171 の sl_3 F_K についても、公開コードは確認できなかった。

## 「双曲結び目の [2,1] に最も近い構成」の順位

1. **Gruen–San Martín Suárez の F_K^{sl_3}（2508.15171）**：定義どおりの構成で、定理であり、随伴 (1,1) を含み、双曲結び目を含む正ブレイド結び目の族を扱う。欠点は、t も a もないこと、量が級数であること、計算例がトーラス結び目だけであること。10_139\[1\] が適用範囲に入るのは筆者の推論で、計算した例はない。
2. **Beliakova–Gorsky（2101.08243）**：任意の結び目と任意の彩色について定義どおりで、t が入る基底（補間 Macdonald）をもつ。欠点は、N ごとの構成であること、t 解釈が「期待」にとどまること、展開係数 a_λ(K) を明示計算した例が §11.5 の 8 の字結び目（gl_2）しかなく、双曲結び目の [2,1] 係数は示されていないこと。
3. **[2,1] HOMFLY-PT（Gu–Jockers 1407.5643、MMM 1408.3076、1504.00371）**：a, q について定義どおりで、4_1、5_2、6_1、pretzel に具体的な値がある。t がないので、refined 版の t=−1 極限の拘束条件としてしか使えない。
4. **(a,q,t) 変形 F_K（2005.13349、2110.13768）**：3 変数はそろうが対称彩色のみで、著者自身が t 変形の非一意性を認めている。
5. **Webster 系のホモロジー**：定義上は [2,1] を含むが、非 minuscule では無限次元になり、計算例もない。
6. **Macdonald 変形した随伴不変量（2505.16569）**：t はあるがトーラス links のみ。\[26\]

順位づけの基準は「定義どおりであること」と「双曲結び目の族を実際に含むこと」を優先し、t をもつことは次点とした。KN の rLMOV 不変量にそのまま入力できる構成は、どの順位にも存在しない。

## Recommendations
- **短期（検証データの構築）**：Gu–Jockers の {2,1} HOMFLY-PT（4_1、5_2、6_1）を t=−1 の拘束条件として使い、KN の BPS 展開を t ansatz の段で試す。そのうえで、正値性を「選ぶ条件」として使った部分と「予言」として確かめた部分を分けて報告すれば、循環論法を避けられる。
- **N=3 の独立な検証**：10_139（正ブレイド）について 2508.15171 の F_K^{sl_3} を計算し、(n,m)=(1,1) に特殊化して HOMFLY-PT_{[2,1]}(a=q^3) と照合する。計算量が問題なので、著者が示唆する stratified state sum や inverted state sum の sl_3 版の開発を待つか、著者に問い合わせるのが現実的である。\[1\]
- **中期（原理の候補）**：Beliakova–Gorsky の補間 Macdonald 基底で、展開係数の t 変形と a 一様化を試みる。KNTZ が長方形でやったことの非長方形版にあたる。これは文献にない新しい研究課題になる。
- **避けるべきこと**：F_K 系の t 変形を [2,1] に流用すること。著者自身が非一意性を明記している。

## Caveats
- **調査済みリストへの影響**：2609.09437（概要のみ確認）は、「非長方形の色 [2,1] では、縮約した evolution に 5 つの部分 Young 図式以外の追加チャネルが必要」と述べているが、HOMFLY-PT に限られ t はない。1609.00143\[43\] は「非長方形表現や他の結び目への一般化は今後」としている。**調査済みリストの結論を変える新しい情報は見つからなかった**。\[44\]
- 2508.15171 は PDF の全文で確認したが、arXiv の概要ページは開いていない。10_139 への適用と、(1,1) への特殊化の規約（x,y の q シフト）は筆者の推論である。
- 検索の制約：web 検索の予算（18 回）を使い切った。サブエージェントが追加で 11 回検索した。そのため、F_K の t 変形の 2024〜2026年の続編、Webster 系の 2025〜2026年の計算例、3d mirror や Coulomb branch による非対称彩色（今回は該当なし）については、見落としの可能性が残る。
- **検索した語（「見つからなかった」項目について）**：refined F_K t-deformation; higher rank F_K arbitrary representation; knot complement quiver non-symmetric; sl3 foam adjoint colored computation; colored HOMFLY-PT homology non-symmetric figure-eight; khoca colored; superpolynomial [21] non-torus; refined Chern-Simons hyperbolic [2,1]; Gukov Putrov Vafa refined colored; Webster explicit computation; Macdonald deformation Vogel universality（サブエージェント分：superpolynomial figure-eight [2,1]; Khovanov–Rozansky [2,1] figure-eight; [21] superpolynomial twist knots differential expansion t-deformation など）。範囲は arXiv を中心に 2010〜2026年10月。

## 文献一覧

### arXiv の概要ページで確認できたもの
- 1001.2020 — B. Webster, "Knot invariants and higher representation theory I: diagrammatic and geometric categorification of tensor products", 2010\[45\]
- 1902.02275 — M. Kameyama, S. Nawata, R. Tao, H. D. Zhang, "Cyclotomic expansions of HOMFLY-PT colored by rectangular Young diagrams", 2019（v3 2020、Lett. Math. Phys.）
- 2211.08409 — J. Wang, "Colored sl(N) homology and SU(N) representations", 2022
- 1407.5643 — J. Gu, H. Jockers, "A note on colored HOMFLY polynomials for hyperbolic knots from WZW models", 2014（サブエージェントが概要ページで確認）
- 2609.10941 — P. Guicardi, "On Diagrammatic Categorification of Verma Modules I: Braiding", 2026（サブエージェント。一部はミラーサイトから取得）\[38\]

### arXiv の本文（PDF/HTML）または検索結果の抜粋でのみ確認したもの（概要ページは未確認）
- 2508.15171 — A. Gruen, L. San Martín Suárez, "A large color R-matrix for sl_3", 2025（全文を確認）\[1\]
- 2101.08243 — A. Beliakova, E. Gorsky, "Cyclotomic expansions for gl_N knot invariants via interpolation Macdonald polynomials", 2021
- 1904.06057 — S. Gukov, C. Manolescu, "A two-variable series for knot complements", 2019\[46\]
- 1909.13002 — S. Park, "Higher rank Ẑ and F_K", 2019\[1\]
- 2106.03942 — S. Park, "Inverted state sums, inverted Habiro series, and indefinite theta functions", 2021\[1\]
- 2005.13349 — T. Ekholm, A. Gruen, S. Gukov, P. Kucharski, S. Park, P. Sułkowski, "Ẑ at large N: from curve counts to quantum modularity", 2020\[7\]
- 2110.13768 — T. Ekholm, A. Gruen, S. Gukov, P. Kucharski, S. Park, M. Stošić, P. Sułkowski, "Branches, quivers, and ideals for knot complements", 2021
- 1408.3076 — A. Mironov, A. Morozov, A. Morozov, "On colored HOMFLY polynomials for twist knots", 2014
- 1504.00371 — A. Mironov, A. Morozov, A. Morozov, P. Ramadevi, V. K. Singh, "Colored HOMFLY polynomials of knots presented as double fat diagrams", 2015\[23\]
- 2101.00695 — W. Qin, "HOMFLY Polynomials of Pretzel Knots", 2021
- 2505.16569 — L. Bishler, A. Mironov, A. Morozov, "Macdonald deformation of Vogel's universality and link hyperpolynomials", 2025\[47\]
- 2512.23095 — A. Mironov, V. K. Singh, "Torus Knots in Adjoint Representation", 2025\[48\]
- 1309.3796 — B. Webster, "Knot invariants and higher representation theory", 2013\[49\]
- 1212.2553 — L. Lewark, "sl3-foam homology calculations", 2012
- 2505.03916 — K. Ritter von Merkl ほか, "Computing colored Khovanov homology", 2025
- 2512.23700 — "Quantum invariants and fiberedness", 2025（著者は未確認）
- 2508.10279 — "A supergroup series for knot complements", 2025（著者は未確認）
- 2607.12155 — "fkcompute: an efficient F_K invariant calculator", 2026（著者は未確認）
- 2601.18042 — "Skeins, q-series, and modularity", 2026（著者は未確認）
- 2608.29585 — "Cyclotomic Newton Expansions and a Rank-Uniform Integer-Valued Newton Completion", 2026（著者は未確認）
- 2609.13766 — "Cyclotomic expansions of colored SU(n) invariants of two-strand torus knots and Bailey transforms", 2026（著者は未確認）
- 1109.1745 — "A Categorification of Quantum sl_3 Projectors and the sl_3 Reshetikhin-Turaev Invariant of Tangles", 2011（著者は未確認）\[32\]
- 1503.08451 — "Categorification of the colored sl_3-invariant", 2015（著者は未確認）
- 2109.12889 — "A Lie theoretic categorification of the coloured Jones polynomial", 2021（著者は未確認）\[34\]
- 1208.3544 — A. Morozov, "Special colored superpolynomials and their representation-dependence", 2012（サブエージェント）
- 2104.14491 — A. Anokhina, A. Morozov, A. Popolitov, "Khovanov polynomials for satellites and asymptotic adjoint polynomials", 2021（サブエージェント）
- 2308.13095 — A. Anokhina, E. Lanina, A. Morozov, "Towards tangle calculus for Khovanov polynomials", 2023（サブエージェント）
- 2410.10468 — "BPS Spectra of Complex Knots", 2024（refined BPS 構造は |R|≤2 まで。著者は未確認）\[50\]

### arXiv では確認できず、学術誌サイトや引用文献リストでのみ確認したもの
- 1511.09077 — A. Mironov, A. Morozov, "Universal Racah matrices and adjoint knot polynomials: Arborescent knots", Phys. Lett. B 755 (2016)\[24\]
- 2212.05222 — A. Gruen, "The sl_N symmetrically large coloured R matrix", 2022（2508.15171 の参考文献から）\[1\]
- 2110.03616 — Q. Chen, K. Liu, S. Zhu, "Cyclotomic expansions for the colored HOMFLY-PT invariants of double twist knots", 2021（参考文献から。内容は未確認）
- 2005.13394 — P. Kucharski, "Quivers for 3-manifolds: the correspondence, BPS states, and 3d N=2 theories", JHEP 2020\[10\]
- 1602.05302 — S. Gukov, P. Putrov, C. Vafa, "Fivebranes and 3-manifold homology", JHEP 2017\[14\]
- 1701.06567 — S. Gukov, D. Pei, P. Putrov, C. Vafa, "BPS spectra and 3-manifold invariants", JKTR 2020\[51\]
- khoca — L. Lewark, A. Lobb（ETH の Web ページと GitHub。arXiv 論文ではない）

## 出典

1. <https://arxiv.org/pdf/2508.15171>
2. [A TWO-VARIABLE SERIES FOR KNOT COMPLEMENTS SERGEI GUKOV AND CIPRIAN MANOLESCU](https://web.stanford.edu/~cm5/surgeries.pdf)
3. [Gukov–Manolescu Series: Knot and 3-Manifold Invariants](https://www.emergentmind.com/topics/gukov-manolescu-series)
4. [Higher rank \$\\hat{Z}\$ and \$F\_K\$](https://arxiv.org/pdf/1909.13002)
5. [A large color 𝑅-matrix for 𝔰⁢𝔩\_3](https://arxiv.org/html/2508.15171)
6. [\$\\widehat{Z}\$ at large \$N\$: from curve counts to quantum modularity](https://arxiv.org/pdf/2005.13349)
7. [\$\${\\widehat{Z}}\$\$ at Large N: From Curve Counts to Quantum Modularity](https://link.springer.com/article/10.1007/s00220-022-04469-9)
8. [Branches, quivers, and ideals for knot complements](https://arxiv.org/html/2110.13768v1)
9. [Branches, quivers, and ideals for knot complements (Journal Article)](https://www.osti.gov/pages/biblio/1977165-branches-quivers-ideals-knot-complements)
10. [Quivers for 3-manifolds: the correspondence, BPS states, and 3d \$\$ \\mathcal{N} \$\$ = 2 theories](<https://link.springer.com/article/10.1007/JHEP09(2020)075>)
11. [A supergroup series for knot complements](https://arxiv.org/html/2508.10279v2)
12. [Branches, quivers, and ideals for knot complements - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0393044022000705)
13. [Quantum invariants and fiberedness](https://arxiv.org/html/2512.23700)
14. [Fivebranes and 3-manifold homology](https://www.researchgate.net/publication/298809263_Fivebranes_and_3-manifold_homology)
15. [Cyclotomic expansions for \$\\mathfrak{gl}\_N\$ knot invariants via interpolation Macdonald polynomials](https://arxiv.org/pdf/2101.08243)
16. [Cyclotomic Newton Expansions and a Rank-Uniform Integer-Valued Newton Completion](https://arxiv.org/html/2608.29585)
17. [Cyclotomic expansions of HOMFLY-PT colored by rectangular Young diagrams](https://link.springer.com/article/10.1007/s11005-020-01318-5)
18. [Cyclotomic expansions of HOMFLY-PT colored by rectangular Young diagrams](https://arxiv.org/html/1902.02275)
19. [Cyclotomic expansions for \$\$\\mathfrak {gl}\_N\$\$ link invariants via interpolation Macdonald polynomials](https://link.springer.com/article/10.1007/s00029-024-00990-y)
20. [A note on colored HOMFLY polynomials for hyperbolic knots from WZW models](https://arxiv.org/abs/1407.5643)
21. [On colored HOMFLY polynomials for twist knots](https://arxiv.org/pdf/1408.3076)
22. [Colored HOMFLY polynomials of knots presented as double fat diagrams - Journal of High Energy Physics](<https://link.springer.com/article/10.1007/JHEP07(2015)109>)
23. [The KNTZ trick from arborescent calculus and the structure of the differential expansion](https://link.springer.com/article/10.1134/S0040577920080036)
24. [Universal Racah matrices and adjoint knot polynomials: Arborescent knots - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0370269316000824)
25. [HOMFLY Polynomials of Pretzel Knots William Qin∗ Abstract](https://arxiv.org/pdf/2101.00695)
26. [Andrei D. Mironov](https://www.sciencedirect.com/author/7201531779/andrei-d-mironov)
27. [Torus Knots in Adjoint Representation](https://arxiv.org/html/2512.23095v2)
28. [Torus knots in adjoint representation](https://journals.aps.org/prd/abstract/10.1103/ncsr-z788)
29. [Towards tangle calculus for Khovanov polynomials](https://arxiv.org/pdf/2308.13095)
30. [1Introduction](https://arxiv.org/html/1208.3544)
31. [(PDF) Knot invariants and higher representation theory II: the categorification of quantum knot invariants](https://www.researchgate.net/publication/45919679_Knot_invariants_and_higher_representation_theory_II_the_categorification_of_quantum_knot_invariants)
32. [A Categorification of Quantum sl\_3 Projectors and the sl\_3 Reshetikhin-Turaev Invariant of Tangles](https://arxiv.org/pdf/1109.1745)
33. [Categorification of the colored 𝔰⁢𝔩\_3-invariant](https://arxiv.org/html/1503.08451)
34. [A Lie theoretic categorification of the coloured Jones polynomial](https://arxiv.org/html/2109.12889)
35. [\[2211.08409\] Colored sl(N) homology and SU(N) representations](https://arxiv.org/abs/2211.08409)
36. [\[2609.10941\] On Diagrammatic Categorification of Verma Modules I: Braiding](https://arxiv.org/abs/2609.10941)
37. [베르마 가군의 도식적 범주화 I부: 꼬임](https://www.alphaxiv.org/ko/abs/2609.10941)
38. [On Diagrammatic Categorification of Verma Modules I: Braiding](https://arxiv.org/pdf/2609.10941)
39. [Khoca - Lukas Lewark](https://people.math.ethz.ch/~llewark/khoca.php)
40. [sl3-FOAM HOMOLOGY CALCULATIONS LUKAS LEWARK](https://arxiv.org/pdf/1212.2553)
41. [fkcompute: an efficient \$F\_K\$ invariant calculator](https://arxiv.org/pdf/2607.12155)
42. [\[2505.03916\] Computing colored Khovanov homology](https://ar5iv.labs.arxiv.org/html/2505.03916)
43. [1Introduction](https://arxiv.org/html/2609.09437)
44. [Rectangular superpolynomials for the figure-eight knot](https://arxiv.org/pdf/1609.00143)
45. [\[1001.2020\] Knot invariants and higher representation theory I: diagrammatic and geometric categorification of tensor products](https://arxiv.org/abs/1001.2020)
46. [A two-variable series for knot complements](https://arxiv.org/pdf/1904.06057)
47. [High Energy Physics - Theory May 2025](http://arxiv.org/list/hep-th/2025-05?skip=225&show=500)
48. [Torus Knots in Adjoint Representation](https://arxiv.org/pdf/2512.23095)
49. [arXiv:1309.3796v3 \[math.GT\] 12 Aug 2015](https://arxiv.org/pdf/1309.3796)
50. [BPS Spectra of complex knots](https://arxiv.org/pdf/2410.10468)
51. [BPS spectra and 3-manifold invariants](https://www.worldscientific.com/doi/abs/10.1142/S0218216520400039)
