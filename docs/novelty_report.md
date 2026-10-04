# 新規性確認報告：iterated torus 結び目への refined large N 双対の拡張

最終更新：2026-10-05（チャットでの文献調査に基づく）

## 0. この文書の使い方

- 本研究の成果として書いてよい範囲と、既知の結果を区別するための文書である。§2 に挙げた結果を、本研究の成果として書かないこと。
- 文献情報には確認状況の印を付けている。
  - 【確認】調査時に arXiv の要旨・本文・出版情報で確認した。
  - 【調査報告】文献調査レポートの段階で参照したもの。番号・内容は原典で再確認していない。
  - 【未確認】番号や内容が確認できていない。
- §4 の「未確認事項」は、ステージ 0 で必ず解消すること。

---

## 1. 研究目標

Kameyama–Nawata "Refined large N duality for knots"（arXiv:1703.05408【確認】、J. Knot Theory Ramifications 掲載【確認】）の計算を、トーラス結び目以外の結び目クラスへ拡張する。具体的には、iterated torus（ケーブル）結び目について次の三つを行う。

1. |λ| ≤ 3 の全ヤング図形 [1], [2], [1,1], [3], [2,1], [1,1,1] で refined CS 不変量を計算する。
2. それを refined BPS 不変量 N̂_{R,g,β,J_r} に展開する。
3. 正値性・整数性を判定する。

---

## 2. 既知の結果（本研究の成果として書かないもの）

### 2.1 refined CS 理論とトーラス結び目

- **Aganagic–Shakirov**（arXiv:1105.5117, 1202.2489【確認】）
  - refined CS 理論は semi-free な円作用をもつ3次元多様体上で定義され、解は Macdonald 多項式に関係した S,T 行列の1パラメータ変形で与えられる。
  - Seifert 多様体とその中のトーラス結び目の不変量を得た。基本表現で彩色した多数のトーラス結び目について、Khovanov–Rozansky（KR）ホモロジーの Poincaré 多項式との一致を確認している。
  - S³ 上の refined CS の large N 双対が resolved conifold 上の refined topological string であることも示した。
- **Cherednik の DAHA-Jones 多項式**（arXiv:1111.6195【調査報告】）、**Gorsky–Negut**（arXiv:1304.3328【確認】）
  - AS の構成を DAHA の観点から捉え直し、Hilb^n(C²) と結びつけた。
  - Gorsky–Negut は Cherednik の superpolynomial 安定化予想を証明した。
- **KR ホモロジー側の証明**：
  - Mellit（arXiv:1704.07630【調査報告】）：任意のトーラス結び目の三重次数 KR ホモロジー
  - Elias–Hogancamp（arXiv:1603.00407【調査報告】）
  - Hogancamp–Mellit（arXiv:1909.00418【調査報告】）：正トーラス絡み目と Sym^l 彩色トーラス結び目

### 2.2 refined large N 双対（1703.05408 本体）

- **トーラス結び目・トーラス絡み目**【確認】：
  - U(N) refined CS の large N 双対を定式化した。
  - refined CS 不変量を resolved conifold の refined BPS 不変量に関係づけ、正値性予想として解釈した。
- **非トーラス結び目は第6節（v2 で追加）**【確認】：
  - 確認できた表は 4_1 の R=[1], [2], [1,1]（Table 16–18）だけである。
  - 4_1 の入力は refined CS の S,T 行列から導いたものではなく、既存の superpolynomial と考えられる【未確認：本文で要確認】。
- **補助ファイル**：Mathematica ノートブックが2つ添付されている【確認】。

### 2.3 iterated torus 結び目の refined 不変量（定義と一般理論は既知）

- **Cherednik–Danilenko**（arXiv:1408.4348、Algebr. Geom. Topol. 16 (2016) 843–898【確認】）
  - DAHA-Jones 多項式を、任意の約化ルート系と重みについて任意の iterated torus 結び目へ拡張した。
  - 多項式性・双対性を示した。非負係数の場合に reduced stable KR 多項式と一致すると推測している。
- **Cherednik–Danilenko**（arXiv:1509.08351【確認】）
  - 全 algebraic links を含む iterated torus links へ拡張した（Macdonald 多項式の積、splice diagram）。
  - Hopf link から DAHA-vertex を構成した。
- **Cherednik の講演要旨**【確認】
  - DAHA superpolynomial は現状 iterated torus links に対してのみ定義されるが、任意の彩色（Young 図形）で問題なく扱える。
- **Cherednik**（arXiv:2410.14703【確認】）
  - algebraic cable の DAHA superpolynomial と彩色トーラス絡み目のものを結ぶ漸化式を与えた。
- **Morton–Samuelson**（arXiv:1410.0859【調査報告】）
  - HOMFLYPT skein algebra of T² と楕円 Hall 代数の対応から、iterated cable の3変数多項式を構成した。Cherednik–Danilenko の多項式に特殊化する。
- **Cherednik–Philipp**（arXiv:1605.00978【確認】）
  - 非トーラスの代数的結び目に焦点を当てた研究。正値性について、長方形 Young 図形が例外的だと述べている。

> **注意**：ケーブル結び目の refined 不変量は「定義・計算法」としては既知である。本研究の新規性は、定義そのものではなく §3 の範囲にある。

### 2.4 非トーラス結び目の彩色不変量（[2,1] を含む）

- **Anokhina–Mironov–Morozov–Morozov**（arXiv:1211.6375【確認】）
  - 4_1 と 3-strand 結び目の [21] HOMFLY を厳密に計算した。トーラス以外で (反)対称表現を超えた最初の結果である。
  - 4_1 の [21] superpolynomial も議論しているが、special polynomial 近傍の展開についての予想に依存している。
- **Anton Morozov**（arXiv:1211.4596【確認】）
  - 任意表現の superpolynomial の一次ずれ公式を推測した。
  - 2012 年時点では [2,1] の superpolynomial の例が一つも知られていないと述べ、4_1・[2,1] への適用結果を「やや議論の余地がある」としている。
- **4_1 の (反)対称・長方形表現**：
  - Itoyama–Mironov–Morozov–Morozov（arXiv:1203.5978【確認】）
  - Fuji–Gukov–Sułkowski（arXiv:1205.1515【確認】）
  - Nawata–Ramadevi–Zodinmawia–Sun（arXiv:1209.1409【確認】）
  - Kononov–Morozov（arXiv:1609.00143【確認】）：長方形でない表現での彩色 superpolynomial の構成はまだ不明瞭と述べている。
  - Nawata–Ramadevi–Zodinmawia（arXiv:1310.2240【確認】）：[r,r] 彩色
- **genus 2 の refined CS**：Arthamonov–Shakirov（arXiv:1504.02620【調査報告】）
  - sl₂ のみ。4_1 で負係数が出て、Poincaré 多項式と一致しない。
  - q-6j 記号は単体では Macdonald 変形できず、その2乗だけが変形できる。
- **[2,1] 彩色のホモロジー**：Gukov–Stošić の quadruply-graded colored homology（arXiv 番号【未確認】）に「(2,1)-colored trefoil の提案」の節がある【確認：目次のみ】。

### 2.5 その他の関連（直接の競合ではない）

- **refined Rosso–Jones・evolution 法**：
  - Dunin-Barkowski–Mironov–Morozov–Sleptsov–Smirnov（arXiv:1106.4305【確認】）
  - Mironov–Morozov–Shakirov–Sleptsov（arXiv:1201.3339【調査報告】）
- **algebraic links と Hilbert scheme**：
  - Oblomkov–Shende（arXiv:1003.1568【調査報告】）
  - Oblomkov–Rasmussen–Shende（arXiv:1201.2115【調査報告】）
- **monotone/Coxeter links**：
  - Galashin–Lam（arXiv:2307.16794【調査報告】）
  - Caprau–González–Hogancamp–Mazin（arXiv:2407.18123【調査報告】）

---

## 3. 本研究の新規性の範囲（暫定）

以下は、§4 の未確認事項が解消されるまで**暫定**である。

| 項目 | 状態 | 根拠 |
|---|---|---|
| (N1) ケーブル結び目の refined BPS 不変量 N̂_{R,g,β,J_r}（1703.05408 の枠組み）の計算 | 新規と思われる | 1703.05408 の非トーラス結び目は 4_1 のみ。ケーブルを扱った後続研究は未確認 |
| (N2) 非トーラス結び目での [3], [2,1], [1,1,1] の refined BPS 不変量 | 新規と思われる | 1703.05408 の 4_1 は [1],[2],[1,1] のみ |
| (N3) 具体的なケーブル（[(2,3),(2,13)], [(2,3),(2,1)] など）の |λ|≤3 全彩色での refined 不変量の明示値 | 一部既知の可能性あり | DAHA の論文に明示例がある可能性。§4 で確認 |
| (N4) refined ケーブル公式（Rosso–Jones 型の Schur→Macdonald 持ち上げ）と DAHA の同値性 | 不明 | Morton–Samuelson が skein/EHA と DAHA を関係づけている。プレシズム形の公式との同値性が既知かは未確認 |
| (N5) 正値性予想の検証：代数的 vs 非代数的ケーブル、長方形 vs [2,1] | 新規と思われる | Cherednik–Philipp は DAHA 側で長方形の特殊性に言及しているが、BPS 不変量としての判定は未確認 |

**書いてはいけないこと**：
- ケーブル結び目の refined 不変量を「初めて定義した」「初めて計算可能にした」
- トーラス結び目の refined BPS 不変量
- 4_1 の [1],[2],[1,1] の結果

---

## 4. 未確認事項（ステージ 0 で解消すること）

1. 1703.05408 の第6節の全文を読み、4_1 以外の非トーラス結び目、[2,1] 彩色、ケーブルへの言及の有無を確認する。
2. 1703.05408 を引用している論文（INSPIRE / Google Scholar）を確認し、ケーブル結び目や非トーラス結び目の refined BPS 不変量を計算した後続研究がないか調べる。特に Kameyama、Nawata、その共著者の後続論文。
3. 1408.4348、1509.08351、2410.14703、1605.00978 に載っている明示例（結び目と彩色）を一覧にする。§3 の (N3) の範囲を確定させ、ステージ 5 の照合対象にする。
4. refined Rosso–Jones 型のケーブル公式が文献に既にあるかを確認する（Mironov–Morozov グループ、Morton–Samuelson、Cherednik）。§3 の (N4) を確定させる。
5. §2 の【調査報告】と【未確認】の文献情報を原典で確認し、印を更新する。
6. arXiv:2610.00303（double twist 結び目の rectangular cyclotomic 展開。2026年10月の新しいプレプリント、査読前）は本研究の対象外だが、関連文献として扱いを決める。
