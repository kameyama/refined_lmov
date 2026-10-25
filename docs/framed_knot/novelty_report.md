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

---

## 5. ステージ 0 での確認結果（2026-10-05 追記）

§4 の 6 項目を調べた結果。詳細と URL は `docs/framed_knot/literature/stage0_bibliography.md`（書誌・被引用）、
`docs/framed_knot/literature/1703.05408.md`（KN 本文）、`docs/framed_knot/literature/notes_DAHA.md`（DAHA の明示例）。§2・§3 の本文は書き換えていない。

1. **1703.05408 の §6 と付録 C**【確認：TeX ソース】
   - 付録 C の表：4_1 は [1],[2],[1,1] だけでなく **[3],[2,1],[1,1,1] もある**。5_2・6_1・6_3 は [1],[2],[1,1]。
     ノートブックのデータは 4_1, 5_2, 6_1, 6_2, 6_3, 7_2, 8_1。
   - 4_1 の [1],[2],[1,1] の入力は、彩色 HOMFLY ホモロジーの Poincaré 多項式（t_c 次数）を (cov) で書き換えたもの（§2.2 の【未確認】は解消）。
   - 4_1 の [2,1] は既存の superpolynomial ではない。KN がミラー対称性・指数増大性・両手型性を課して求めた候補で、負係数を含む（式 (fig8-21)）。
   - homologically thick な 9_42 と非トーラス絡み目では、形 (hat-f-2) に乗らないと述べている。
   - ケーブル・iterated torus・代数的結び目への言及はない。
   - **→ §3 (N2) の修正**：「非トーラス結び目での [3],[2,1],[1,1,1]」のうち 4_1 は既知（KN 付録 C）。
     新規と言えるのは「**ケーブル**（S,T 由来・DAHA 由来の refined 不変量をもつ非トーラス結び目）での [3],[2,1],[1,1,1]」に限る。
2. **1703.05408 を引用する論文**【確認：INSPIRE 14 件 + Semantic Scholar 2 件】
   - ケーブル・iterated torus 結び目、または [2,1]・[3]・[1,1,1] 彩色について refined BPS 不変量やその正値性を計算した論文は無かった。
   - 近いものは対称表現の LMOV（knots–quivers 1707.02991, 1707.04017、Ekholm–Kucharski–Longhi 1811.03110）と、toric CY 上の refined open BPS の非負性（Cheng–Sułkowski 2104.00713）。
   - Kameyama・Nawata の後続（1902.02275）は長方形色の cyclotomic expansion を扱い、refined ではなくケーブルも扱わない。
   - 1703.05408 の書誌：J. Knot Theory Ramifications 32 (2023) no.06, 2041001、DOI 10.1142/S0218216520410011【確認】。
3. **DAHA 論文の明示例**【確認：TeX ソース】（一覧は `notes_DAHA.md` §2.3 の表と §4）
   - K1 = (2,13)-ケーブル [(3,2),(2,1)]：[1] と [2]（2ω₁）が全係数で載っている（CD14 (1),(5)）。[1,1] は q=1, a=0 の値のみ（super-duality で [2] から単項式倍を除いて決まる）。
   - K3 = (2,15)-ケーブル：[1]（CD14 (2)）。
   - K2 = (2,1)-ケーブル、K4 = (2,3)-ケーブル：[1]（CD14 (d),(c)）。(2,5),(2,7),(2,−1) の [1] もある（K5）。
   - K6a = T(2,5) の (2,21)-ケーブル：DAHA-Betti 多項式のみ。
   - [3],[2,1],[1,1,1] の明示例は、どのケーブルについても無い。
   - **→ §3 (N3) の確定**：K1 の [1],[2] と K2・K3・K4・K5 の [1] は**既知**（単項式倍を除く）。ステージ 5 の照合対象にする。
     K1 の [1,1] は super-duality により実質既知。K1〜K5 の [3],[2,1],[1,1,1]、および K2〜K5 の [2],[1,1] は明示例が見つからなかった。
4. **refined のケーブル公式の既存文献**【確認：キーワード検索。網羅的ではない】
   - refined のケーブル不変量は DAHA の反復（CD14, CD15）で定義済みで、EHA 側（Morton–Samuelson 1410.0859）との一致も示されている。
   - 非 refined の cabling 公式は確立している（Anokhina–Morozov 1307.2216 など）。
   - Rosso–Jones の素朴な t 変形が失敗するという報告がある（Danilenko 1405.0884、Morozov 1801.02407、DMMSS 1106.4305 の γ 因子）。本研究の試算でも、素朴なプレシズムはトーラス結び目で失敗した（conventions.md §6）。
   - 「companion の彩色 refined 不変量だけを入力とする Rosso–Jones 型の明示公式」は見つからなかった。
   - **→ §3 (N4)**：素朴な形の公式は誤り（既知の失敗）。正しい形は DAHA の定義を companion の段で線形に分解したもの（conventions.md §7）で、定義から直ちに従う。
     したがって、これを本研究の新しい結果として強く主張するのは難しい。「計算法として明示し、検証した」程度の位置づけにとどめる。
5. **【調査報告】の文献情報**【確認：arXiv API】
   - 27 件の arXiv 番号すべてで、想定した著者・題名と一致した。
   - 例外：「Gukov–Stošić の quadruply-graded colored homology」は **Gorsky–Gukov–Stošić、arXiv:1304.3481**（Fundam. Math. 243 (2018) 209–299）。§4.6 で (2,1) 彩色の三葉結び目ホモロジーを提案している。
   - 1111.6195（Cherednik）の掲載誌は未確認。
6. **arXiv:2610.00303**【確認：要旨のみ】
   - Fang–Zhou, "Rectangular Cyclotomic Expansions of Colored HOMFLY-PT Polynomials"（2026-09-28 投稿）。
   - 長方形 cyclotomic expansion 予想と double twist 結び目の公式の証明で、要旨では refined 不変量・ケーブルに触れていない。
   - 扱い：関連文献として引用するにとどめ、照合対象にはしない。

### 修正後の新規性の範囲（暫定）

| 項目 | ステージ 0 後の判断 |
|---|---|
| (N1) ケーブルの refined BPS 不変量（KN の枠組み） | 新規と思われる（被引用論文に無い） |
| (N2) [3],[2,1],[1,1,1] の refined BPS 不変量 | 4_1 は既知。ケーブルについては新規と思われる |
| (N3) ケーブルの refined 不変量の明示値 | K1 の [1],[2]（と [1,1]）、K2〜K5 の [1] は既知。それ以外の彩色は明示例なし |
| (N4) refined ケーブル公式 | 素朴な形は誤り（既知）。正しい形は DAHA の定義から従うので、新規性の主張は弱い |
| (N5) 正値性の判定（代数的 vs 非代数的、長方形 vs [2,1]） | 新規と思われる |

---

## 6. 第 2 期（結び目表の方針）の新規性（2026-10-06 追記）

詳細は docs/framed_knot/REPORT_knot_table.md §6。文献確認は KN の引用論文（INSPIRE 14 件、Semantic Scholar 9 件）と主要論文の本文で行った。

| 項目 | 判断 |
|---|---|
| 命題 1：[1] の正値性が非負性と DGR の対称性（GHM 2103.01212）から従う | 新規と思われる（構成要素は既知） |
| 命題 2：2 箱の Laurent 性 ⇒ dim H^{[2]} = (dim H^{[1]})²、9_42 の破れの原因 | 導出は新規と思われる。9_42 の S² の次元 401（GGS 1304.3481）と、9_42 が形に乗らないこと（KN）は既知 |
| 結び目表 801 個の [1] の refined BPS | 新規と思われる |
| [1] の 2J_r の偶奇 ⇔ parity knot | 条件は既知（Chandler–Gorsky 2209.13058）。BPS との関係は新規と思われる |
| ツイスト結び目・pretzel の 2 箱（族として）と、t_c 次数の選び方への依存 | 新規と思われる。knots–quivers の結果は別の量（対称彩色の refined LMOV の整数性） |
| GS 1112.0030 の誤植 | 明示的な指摘は見つからなかった（GGS で暗黙に修正済み） |

書いてはいけないこと（追加）：9_42 の S² の次元 401、9_42 が KN の形に乗らないこと、parity knot の定義、入力の式（NRZS 1209.1409、GNSSS 1512.07883、Nakagane–Sano 2111.00388）。

---

## 7. 第 3 期（3 箱）の新規性（2026-10-06 追記）

詳細は docs/framed_knot/REPORT_phase3.md §6、文献は docs/framed_knot/literature/phase3_21_survey.md。

| 項目 | 判断 |
|---|---|
| KN の表に無いトーラス結び目 5 個（T(3,4) の 3 箱、T(3,5), T(2,11), T(3,7), T(4,5)）の 3 箱の正値性 | 新しい表だが、KN がトーラス結び目全般について予想しているので「予想の確認」にあたる |
| 代数的なケーブル K9 の 3 箱の正値性（第 1 期の K1, K3, K6a に続く） | 新規と思われる |
| 非トーラス結び目の [2,1] が構造的な条件で一意に決まらないこと（整数計画で明示） | 新規と思われる。曖昧さ自体は AMMM 1211.6375 も認めている |
| a = t² での [2,1] の U(2) への簡約の恒等式（refined） | 文献での扱いは未確認 |

書いてはいけないこと（追加）：4_1 の [2,1] の候補（KN、AMMM）、[2,1] 彩色 HOMFLY-PT の個別の値（Morozov ら）、Cherednik 2607.25666 の Coxeter 結び目の superpolynomial。

---

## 8. 第 4 期の新規性（2026-10-06 追記、文献確認はサブエージェントによる）

| 項目 | 判断 |
|---|---|
| U(2) への簡約の恒等式（a = t² で rCS_{[2,1]} = rCS_{[1]}·rCS_{[1,1]}） | **既知**（Cherednik–Danilenko 1408.4348 の安定化定理、AMMM 1211.6375 の要請 (D)）。成果として書かない |
| 非トーラス結び目の [2,1] が構造的な条件で一意に決まらないこと、解の差の因数分解 | 部分的に既知（AMMM が曖昧さを認める）。解の空間の構造の明示は新規と思われる |
| 2 段の反復ケーブル K10 の |λ| ≤ 3 の refined 不変量と 3 箱までの refined BPS の正値性（2026-10-08） | 新規と思われる（CD14 は無彩色の DAHA-Betti 多項式のみ。文献の再確認はしていない） |
| ツイスト結び目の \hat f の p についての漸化式と特性根 | superpolynomial の段階では実質的に既知（NRZS 1209.1409、Mironov–Morozov–Morozov 1306.3197）。\hat f の段階は新規と思われる |
| [1] の 2J_r の偶奇 ⇔ parity knot（Chandler–Gorsky 2209.13058）、11 交点以下の parity knot の一覧（非トーラスは 10_125, 10_140） | 条件の同値性と一覧は新規と思われる（原典に一覧は無い） |
| K8（T(3,4) の (2,25)-ケーブル）の 3 箱の正値性 | 新規と思われる（第 3 期の K9 に続く。CD14 は K8 の [1] のみ） |
| thin な双曲結び目の 3 箱で正値性を満たす [2,1] が存在すること（6_1 の整数解、5_2 の有理数解） | 新規と思われる。ただし [2,1] は一意でなく、予想の検証ではなく「矛盾しない」ことの確認 |

書いてはいけないこと（追加）：U(2) への簡約の恒等式そのもの、parity knot の定義、superpolynomial の p についての漸化式（evolution）。

## 9. B4（thick な非トーラス結び目の 2 箱）の新規性（2026-10-07 追記、文献確認はサブエージェントによる。根拠の 2 箇所は親セッションで原典を確認）

| 項目 | 判断 |
|---|---|
| 10_139, 10_152 の S²（t_c 次数）の Poincaré 多項式（次元 441 = 21², 841 = 29²） | 文献の範囲では見当たらない。ただし予想式（Stošić 2402.03066 の quiver 型の式 ＋ KRSS 1707.04017 の (quad-C)）に基づく。Stošić は 10_139, 10_152 について彩色 HOMFLY-PT 多項式の式のみを与え（749, 838 行）、ホモロジーは 9_42 のみ。Chen 1512.07906 は無彩色のみ |
| 10_139, 10_152 の 2 箱の \hat f_{[2]}, \hat f_{[1,1]} が G 展開の前後とも非負 | 新規と考えられる（予想式に基づく）。knots–quivers の DT 不変量の正値性（Efimov）からは従わない：(i) EKL 1811.03110（kq-draft.tex 1566-1568 行、親セッションで確認）は KN の t_c の精密化が分母 (q²;q²t²)_r をもち motivic 級数の形と両立しないと明言、(ii) quiver 由来の refined LMOV は Ω の符号付き整数結合で、trefoil の N_2 は全項負（EKL）、(iii) KN の \hat f_{[2]}, \hat f_{[1,1]} は [1,1] の入力と M(t)^{−1} を要し、対称表現だけの quiver 級数からは出ない。quiver から自動的に従うのは別の量（対称表現の quiver 型 refined LMOV）の整数性のみ |
| 8_19 = T(3,4) の S², S³ の t_c 次数が refined CS と一致すること | 既知（GNSSS 1512.07883 の BPS-SS.tex 4359 行：Shakirov 1308.3838 と r ≤ 3 で一致、T(3,5) も。親セッションで確認）。成果として書かない |
| Stošić の 8_19 の quiver に (quad-C) を当てはめた結果が EHA の [2]〜[1,1,1] と全項一致すること | 部分的に新しい（KRSS は T(3,4) の最低 a 次数の行の quiver のみ）。予想式の検証例という位置づけ |
| thick な結び目で t_c 次数が t_r 次数のデータから決まらない原因が C_ij にあること | (quad-C) から直ちに読める帰結で、新規性は弱い。t_c が t_r から決まらないこと自体は GGS 1304.3481 も同旨 |

書いてはいけないこと（追加）：T(3,4), T(3,5) の四重次数の S^r の式と refined CS との一致（GNSSS）、(quad-C) そのもの（KRSS）、10_139, 10_152, 8_19 の quiver（Stošić）、quiver 由来の refined LMOV の整数性（KRSS, EKL）。
調べた範囲：取得済み・新規取得の文献の TeX の全文検索、KN の被引用（INSPIRE 14 件、Semantic Scholar 9 件）、Stošić の被引用（INSPIRE 4 件、うち 1 件は本文未確認）、KRSS の被引用（2023 年以降、題名のみ）。Google Scholar は未使用。

## 10. framed knot の refined LMOV の新規性（2026-10-09 追記、文献確認はサブエージェントによる。原典は TeX を grep して該当箇所を読んだもの）

調べた範囲：INSPIRE の被引用（KN 14 件、Mariño–Vafa hep-th/0108064 144 件、Cheng–Sułkowski 2104.00713 12 件の全題名）、INSPIRE の全文検索、arXiv API（2026-01-01〜2026-10-09）、主要論文の TeX。

| 項目 | 判断 |
|---|---|
| refined の framing の因子 framing[λ] そのもの | **既知**。正しい形は Taki 0710.1776 の式 (modification) f_μ(t,q) = (−1)^{|μ|}t^{‖μ^t‖²/2}q^{−‖μ‖²/2}（IKV の因子を、Nekrasov の分配関数との一致のために修正したもの。ユーザーの指摘で IKV は誤り・誤植と判断）。‖λ‖² = 2n(λ') + |λ| から **framing[λ]^{−1} にちょうど一致**（親セッションで式を確認）。Cheng–Sułkowski 2104.00713 の f^p_λ も同じ形。IKV hep-th/0701156 の f_ν = (−1)^{|ν|}t^{n(ν)}q^{−n(ν^t)} は (t/q)^{|ν|/2} が欠けている。Diaconescu–Donagi–Pantev 1704.07412 は符号なしの (t/s)^{p\|λ\|/2} f_λ^{−p} |
| 単一ブレーン（対称・反対称表現）の refined open BPS の非負性 | 既知。Cheng–Sułkowski（framing は固定、resolved conifold は §3.3）、Banerjee–Ishtiaque–Jeong 2601.07662（任意の framing f を CS 準位 −f として入れるが、framing は q_1 の冪だけで (q/t) を含む形ではない。系統的な正値性の表は無い） |
| quiver での framing（C → C + f）と対称表現の refined LMOV の整数性 | 既知（KRSS 1707.04017、EKL 1811.03110。EKL は整数性を定義から、正値性は主張しない。KN の refinement とは t_c/t_r の取り方が違うと明記） |
| 非 refined の framed LMOV | 既知（Mariño–Vafa、Zhu 2105.02037・2510.15692、Wang–Zhu 2502.16609、Luo–Zhu 1611.06506・2106.02882、Chen–Liu–Peng–Zhu 1402.3571）。**用語の注意**：Zhu の "refined LMOV conjecture for framed links" は非 refined の整数性の強化で、KN の refined とは別物（Zhu 自身が Remark で明記） |
| framing 付きの refined torus link のデータから正の量（局所曲線幾何、closed の GV 展開） | 関連するが別の対象（DDP 1704.07412、CDDNP 1804.08364） |
| KN の \hat f_ρ を framing x 付きで作り、一般の ρ（\|ρ\| ≤ 3、[2,1] を含む）で G 展開の前後の正整数性を検証すること | **新規と思われる**（前例なし） |
| 𝐭^{−nx}\hat f^{(x)} が x の指数和になる構造、全ての x での 2 箱の正値性 | 新規と思われる |
| 符号 (−1)^{\|λ\|} を含む refined の framing が正値性に必要かの判定 | 新規と思われる（符号の役割を refined で論じた文献は見つからない） |
| 2026 年の新着 | refined LMOV の正値性、iterated torus の refined BPS、ツイスト・pretzel の refined 2 箱、双曲結び目の refined [2,1] を扱う論文は見つからない。周辺：Cherednik 2607.25666、Fang–Zhou 2610.00303（非 refined） |
| K7・K10 の \|λ\| = 3 の refined の値と refined BPS | 新規と思われる。CD14（1408.4348）は K1 の [1], [2]（[1,1] は super-dual）、K7 の [1]、K10 の DAHA-Betti 多項式のみ。Cherednik 1709.07589 も 2 段のケーブルは無彩色のみ |

書いてはいけないこと（追加）：refined の framing の因子の形（Taki・Cheng–Sułkowski）、単一ブレーンの refined open BPS の非負性、quiver の framing のシフト、非 refined の framed の整数性。
未確認：Cheng–Sułkowski の掲載頁。KN の書誌（INSPIRE は J. Knot Theory Ramifications 32 (2023) 2041001、DOI 10.1142/S0218216520410011 は 2020 年の巻を示唆。§5 の記載とあわせて要再確認）。

## 11. DGR 型の構造（別の研究、2026-10-11）

docs/dgr_lmov/novelty.md に移した。
