# 結び目不変量の圏化：アプローチ別の「構成物・計算範囲・破れる場所」の文献地図（2026年10月時点）

結論から言うと、HOMFLY-PT 圏化の主なアプローチのうち、厳密に定義されたホモロジーで全トーラス結び目まで計算が済んでいるのは Soergel 双加群版の三重次数ホモロジー（Mellit、Hogancamp–Mellit）だけです。非トーラス結び目の直接計算は「11 交点以下の素な結び目」（Nakagane–Sano）が最新の上限です。\[1\]iterated torus 結び目（ケーブル）については、DAHA 超多項式という代数的な構成でしか体系的な計算例がありません。その DAHA 超多項式は、三葉結び目の非代数的 2-ケーブル Cab(5,2)T(3,2)（Cherednik–Danilenko の規約）で初めて負の係数をもちます。\[2\]

## TL;DR
- **ホモロジー側**：sl(2) Khovanov ホモロジーは 16 交点以下の全素結び目まで計算されている。三重次数（HOMFLY-PT）ホモロジーは 11 交点以下の全素結び目、全トーラス結び目、正のトーラス絡み目、Sym^l 彩色トーラス結び目まで。DGR 予想の d_0（knot Floer への微分）は、「二つのスペクトル系列」という修正版の形で証明済み（Beliakova–Putyra–Robert–Wagner）。
- **代数側（refined CS / DAHA / 楕円 Hall 代数）**：ホモロジーではない。三葉結び目でも非長方形彩色（3-hook）では負の係数が必ず出る（Cherednik）。ケーブルへの拡張（Cherednik–Danilenko、Morton–Samuelson）では、代数的なケーブルは正値で、計算された範囲では Khovanov 多項式と整合する。非代数的なケーブルでは Cab(7,2)T(3,2) まで正値で、Cab(5,2)T(3,2) 以下では負の係数が出る。Cherednik–Danilenko の Cab(a,b) は b 本の撚り糸のケーブルで、a は Seifert framing で測った傾き（文献から導出）。
- **物理側（refined BPS / LMOV）**：Kameyama–Nawata の正値性予想は、トーラス結び目・トーラス絡み目、thin な非トーラス結び目（(2s−1,1,2t−1)-pretzel、6_2、6_3、三箱の 4_1）で成立が観察されている。thick な 9_42 と、Whitehead 絡み目を含む非トーラス絡み目では、要求される正値の形に書けない。ケーブルについての refined BPS の検証は、今回の調査範囲では見つからなかった。

---

**凡例**：【明記】文献に明記（節・式を添える）／【導出】文献の記述から私が導いたもの／【推測】私の推測。arXiv 番号は、arXiv の概要ページか PDF 本文で確認したものだけを書いた。それ以外は「未確認」とした。

## 1. アプローチの一覧表

| アプローチ | 構成するもの | 厳密か予想か | 代表論文（確認済み） | 計算された結び目 | うまくいかない結び目・状況 | 他との関係 |
|---|---|---|---|---|---|---|
| Khovanov（sl(2)）と変種 | 二重次数ホモロジー。Euler 標数は Jones 多項式 | 厳密 | Bar-Natan math/0201043, math/0410495, math/0606318；Shumakovitch math/0405474 | 16 交点以下の全素結び目、14 交点以下の全素絡み目、T(8,7)（48 交点）など\[3\] | p,q≥3 のトーラス結び目 T(p,q) は全て H-thick（Stošić, AGT 7）。計算時間は指数的に増える | Lee のスペクトル系列、Rasmussen の s 不変量。三重次数ホモロジーから d_2 で得られる |
| sl(N) KR（行列因子化・foam） | 二重次数ホモロジー。Euler 標数は sl(N) 多項式 | 厳密 | Rasmussen math/0607544；Mackaay–Stošić–Vaz 0708.2228；Wu 0907.0695（彩色） | 9 交点以下の全結び目（Rasmussen）。\[4\]Nakagane–Sano のデータから N>2 で直接\[5\] | N=2 では高次の微分が残りうる。\[5\]直接計算は困難 | HOMFLY→sl(N) のスペクトル系列（各 N>0）は定理\[4\] |
| HOMFLY-PT 三重次数（行列因子化／Soergel） | 三重次数ホモロジー。Euler 標数は HOMFLY-PT | 厳密 | Khovanov math/0510265；Mellit 1704.07630；Hogancamp–Mellit 1909.00418；Nakagane–Sano 2111.00388；Gorsky–Hogancamp 1712.03938（y-ified） | 全トーラス結び目、正のトーラス絡み目、(n,n) トーラス絡み目、11 交点以下の全素結び目 | 代数的絡み目の計算は未解決問題（Hogancamp–Mellit Problem 1.10）。ケーブルの直接計算例は見つからず | q↔t 対称性と sl(2) 作用（Gorsky–Hogancamp–Mellit）。ORS 予想はトーラス結び目で証明 |
| 彩色ホモロジー（対称・反対称・長方形、t_r/t_c） | 彩色版の三重〜四重次数ホモロジー | 対称彩色は厳密な構成あり（Webster–Williamson 0905.0486 など）。四重次数は予想\[6\] | Gukov–Stošić 1112.0030；Gorsky–Gukov–Stošić 1304.3481\[7\]\[8\] | Sym^l 彩色トーラス結び目（定理）。三葉結び目・ツイスト結び目・pretzel（予想・物理的） | 非長方形彩色では正値性・構造が保証されない（§2.6） | colored differential、mirror symmetry（行と列の交換）\[7\] |
| superpolynomial（DGR 型） | Poincaré 多項式を微分 d_N と対称性から組み合わせ的に決めたもの | 予想（組み合わせ的な仮定） | Dunfield–Gukov–Rasmussen math/0505662 | 小さい結び目、トーラス結び目 | thick な結び目、高次の微分が非自明な場合 | d_N（N>0）の部分は Rasmussen が定理化、d_0 は BPRW が修正版で証明\[5\]\[9\] |
| refined CS／DAHA／楕円 Hall 代数 | 三変数の多項式（ホモロジーではない） | 代数的には厳密に定義。ホモロジーとの一致は予想 | Aganagic–Shakirov 1105.5117；Cherednik 1111.6195；Cherednik–Danilenko 1408.4348, 1509.08351；Morton–Samuelson 1410.0859；Samuelson 1408.0483 | 全トーラス結び目、任意の iterated torus 結び目・絡み目（形式的にはいくらでも計算可能） | 非長方形彩色（三葉結び目でも）、非代数的ケーブル（Cab(5,2)T(3,2) 以下）で負の係数 | 無彩色のトーラス結び目では Mellit が三重次数ホモロジーとの一致を証明（AS・Cherednik・GN 予想） |
| 幾何（ORS Hilbert スキーム、OR、affine Springer） | 曲線特異点の Hilbert スキームの重み多項式、Hilbert スキーム上の層 | Euler 標数版は定理（Maulik 1210.6323）。ホモロジー版は予想 | Oblomkov–Shende 1003.1568；ORS 1201.2115；Maulik 1210.6323 | トーラス結び目、三葉結び目の (2,13)-ケーブル（HOMFLY 版） | 代数的絡み目に限る（非代数的ケーブルには定義がない） | DAHA 超多項式と一致すると予想（CD Conj. 2.4(iii)） |
| knot Floer と d_0 | 二重次数ホモロジー（Euler 標数は Alexander 多項式） | 厳密 | Beliakova–Putyra–Robert–Wagner 2210.00878 | — | DGR の元の形（単一の微分）ではなく、gl(0) ホモロジーを経由する二段階の形で証明\[9\] | 三重次数ホモロジー →（二つのスペクトル系列）→ HFK-hat\[9\] |
| 物理（refined LMOV、knots–quivers、Ẑ） | BPS 状態の数（整数・正値性の予想） | 予想（knots–quivers は LMOV 整数性の証明を与える） | Kameyama–Nawata 1703.05408；Kucharski–Reineke–Stošić–Sułkowski 1707.02991, 1707.04017 | トーラス結び目・絡み目、6 交点以下の全結び目（quiver）、thin な pretzel など | 9_42（thick）、非トーラス絡み目（Whitehead 絡み目） | 長方形彩色のトーラス結び目では refined CS ＝ ホモロジーの Poincaré 多項式という予想に依存 |
| 計算機 | プログラム | — | JavaKh、KhoHo、khoca、KnotJob、kr-calc\[10\]\[11\]\[12\] | §2.10 | 交点数が増えると指数的に重くなる | — |

## 2. アプローチごとの詳しい節

### 2.1 Khovanov ホモロジー（sl(2)）と変種
1. **構成**：二重次数のホモロジーで、次数付き Euler 標数が Jones 多項式になる。厳密に定義されている。Bar-Natan による tangle と cobordism の圏での定式化（math/0410495）があり、局所的な簡約にもとづく高速アルゴリズム（math/0606318）がある【明記】。\[13\]\[14\]
2. **代表論文**：Bar-Natan math/0201043（初期の計算と予想）、math/0410495、math/0606318。\[15\]Shumakovitch math/0405474（捩れ）。\[3\]Stošić「Homological thickness and stability of torus knots」AGT 7 (2007)（arXiv 番号は未確認）。Khovanov の原論文の arXiv 番号は未確認。
3. **計算された結び目**：
   - Bar-Natan math/0606318：35 交点の T(7,6) と 48 交点の T(8,7)。T(8,7) は JavaKh で数分【明記】。\[13\]
   - Shumakovitch math/0405474：整係数の Khovanov ホモロジーは 16 交点以下の全素結び目、14 交点以下の全素絡み目で既知。KhoHo は 17–19 交点まで実用的【明記】。\[3\]
   - arXiv:2205.11430（Chernov–Maguire）：鏡像を数えずに 20 交点以下の全素結び目 2,199,471,680 個を、トーラス結び目・ツイスト結び目と比較した。Khovanov 多項式がトーラス結び目かツイスト結び目と一致した結び目は、全てそれ自身がトーラス結び目かツイスト結び目だった。Jones 多項式が一致した 4 個については JavaKh で Khovanov ホモロジーを比較している【明記】。
   - arXiv:2304.13613（Kegel–Manikandan–Mousseau–Silvero「Khovanov homology of positive links and of L-space knots」）：正の絡み目では、第 1 Khovanov ホモロジーが消えることとファイバー絡み目であることが同値である。三葉結び目の (2,1)-ケーブルは、ファイバー結び目かつ強擬正なのに第 1 Khovanov ホモロジーが消えない最も簡単な例。T(2,3) の (2,3)-ケーブル K に対し、q≥2 の 2-ケーブル K_{2,q} の低次の Z_2 係数 Khovanov ホモロジーを決定（Thm 5.10）【明記】。
4. **うまくいかない状況**：p,q≥3 の全てのトーラス結び目 T(p,q) は H-thick（Stošić, Thm 3.3）【明記】。\[16\]thick な結び目では Jones 多項式から Khovanov ホモロジーを復元できない【導出】。
5. **変種**：Lee 変形、Bar-Natan 変形、odd Khovanov（Ozsváth–Rasmussen–Szabó）、Lipshitz–Sarkar の安定ホモトピー型（s 不変量の精密化 1206.3532）がある。今回これらの原典の概要ページは確認していないので、arXiv 番号は「未確認」とする。odd Khovanov は KnotJob（Schütz）で計算可能【明記：Schütz の論文】。\[12\]y-ification 版は Sano 2602.17435【明記】。\[17\]\[18\]
6. **他との関係**：三重次数ホモロジーに d_2 を取ると Khovanov ホモロジーが得られる（Rasmussen のスペクトル系列は一般には E_2 で退化しない）【明記：Rasmussen math/0607544】。\[19\]

### 2.2 sl(N) Khovanov–Rozansky ホモロジー
1. **構成**：行列因子化（Khovanov–Rozansky I）または foam（Mackaay–Stošić–Vaz 0708.2228）による二重次数ホモロジー。\[20\]Euler 標数は sl(N) 多項式。\[21\]厳密に定義されている。
2. **代表論文**：Rasmussen math/0607544「Some differentials on Khovanov-Rozansky homology」（Geom. Topol. 19 (2015)）。各 N>0 について HOMFLY→sl(N) のスペクトル系列があり、9 交点以下の全ての結び目の KR ホモロジーが決まる【明記：要旨】。\[4\]\[22\]\[23\]彩色版は Wu 0907.0695。\[20\]
3. **計算された結び目**：9 交点以下の全結び目（Rasmussen）。2 橋結び目は N>4 で「N-thin」（Rasmussen の Duke 論文を 2211.03743 が引用。原論文の arXiv 番号は未確認）。\[24\]Nakagane–Sano の 11 交点以下のデータでは、N>2 の微分 d_m(N) が次数の理由で全て自明なので、sl(N) ホモロジーが三重次数ホモロジーから直接得られる（2111.00388 §5.1）【明記】。\[5\]khoca（Lewark–Lobb）は sl(N) ホモロジーを計算でき、T(7,8) の例がある【明記：README】。\[11\]
4. **うまくいかない状況**：N=2 では高次の微分が非自明になりうる。Nakagane–Sano は一部の結び目について d_m(2)=0（m>1）を手で確認しただけである【明記】。\[5\]
5. **関係**：スペクトル系列の存在は定理。DGR が予想した「単一の微分 d_N」で済むかどうかは、結び目ごとに違う【導出】。\[5\]

### 2.3 HOMFLY-PT（三重次数）ホモロジー
1. **構成**：Khovanov–Rozansky II（行列因子化）と、Khovanov の Soergel 双加群の Hochschild ホモロジー（math/0510265）。三重次数で、Euler\[20\] 標数は HOMFLY-PT。reduced 版と unreduced 版は互いを決める（Rasmussen、Chandler–Gorsky 2209.13058 §2）【明記】。\[25\]
2. **代表論文と計算例**：
   - Elias–Hogancamp 1603.00407：(n,n) トーラス絡み目の全ての n【明記】。\[26\]
   - Mellit 1704.07630：任意のトーラス結び目の明示公式。Aganagic–Shakirov・Cherednik・Gorsky–Negut・ORS の予想の一部を証明【明記：要旨】。\[27\]
   - Hogancamp–Mellit 1909.00418：正のトーラス絡み目 T(m,n) と Sym^l 彩色トーラス結び目【明記】。\[28\]
   - Nakagane–Sano 2111.00388（JKTR 34(12) 2550062, 2025）：11 交点以下の全素結び目。プログラム\[18\] kr-calc【明記】。v1（2021）の時点を紹介した\[1\]\[25\]\[29\] Gorsky の講演スライドでは「10 交点以下と 11 交点の大部分」とされており、v2（2024）で全 11 交点に拡大したと読める【導出】。\[30\]\[31\]
   - Chandler–Gorsky 2209.13058（Exp. Math. 2024）：上のデータ全体について S 不変量を計算【明記】。\[18\]\[31\]\[32\]\[33\]
   - Turner 2410.03068：4 本の撚り糸の Coxeter braid【明記：題名】。\[34\]
3. **y-ified 版と link splitting**：Gorsky–Hogancamp 1712.03938（Geom. Topol. 26 (2022)）。結び目では\[25\] HY(K)=H(K)⊗C[x,y] となる（Gorsky のスライド）【明記】。Gorsky–Hogancamp–Mellit\[35\] 2103.01212 は、全ての結び目について sl(2) 作用と q↔t 対称性を示した【明記】。彩色版の\[30\]\[36\]\[37\] link splitting は文献リストに題名が出てくるだけで、arXiv 番号は未確認。\[38\]
4. **うまくいかない状況**：代数的絡み目（すなわちトーラス絡み目の iterated cable）の三重次数ホモロジーは未解決問題として明示されている（Hogancamp–Mellit 1909.00418, Problem 1.10）【明記】。今回の検索範囲（「triply\[39\] graded homology cables of torus knots」など）では、三重次数ホモロジーをケーブルに対して直接計算した例は見つからなかった。Nakagane–Sano の範囲（≤11 交点）には非自明なケーブルは入らない【導出：最小のケーブルである三葉結び目の (2,3)-ケーブルは 15 交点（2304.13613 が Dunfield の表を引用）】。\[40\]
5. **関係**：ORS 予想はトーラス結び目で証明済み（Gorsky–Mazin–Vazirani、Hogancamp–Mellit。2004.15024 と Gorsky のスライドによる）【明記】。\[35\]\[41\]

### 2.4 彩色ホモロジー
1. **構成**：Gukov–Stošić 1112.0030 は、彩色 HOMFLY ホモロジーに「colored differential」を含む一連の性質を定式化し、対称表現と反対称表現の間の mirror symmetry を見出した【明記：要旨】。\[7\]Gorsky–Gukov–Stošić 1304.3481（Fund. Math. 243 (2018)）は四つの独立な次数を予想し、トーラス結び目の長方形彩色ホモロジーを明示的に記述した【明記：要旨】。\[8\]\[42\]
2. **t_r と t_c**：四重次数予想では、t 次数が t_r と t_c の二成分に分かれる。図版の説明文には Q = q + t_r − t_c という tilde 版の次数が使われている（S^2 彩色の三葉結び目）【明記：図版】。\[43\]Kameyama–Nawata §6 は t_c 次数付きの Poincaré 多項式を入力に使う【明記：サブエージェントによる本文確認】。\[44\]
3. **計算例**：Sym^l 彩色トーラス結び目は定理（Hogancamp–Mellit）。三葉結び目の S^2 彩色（GGS）、ツイスト結び目 K_p の彩色超多項式（Nawata–Ramadevi–Zodinmawia–Sun 1209.1409。微分と特殊極限で検証）、8 の字結び目の全対称・反対称表現（arXiv:1203.5978、著者は未確認）、絡み目の彩色ホモロジーと「sliding」性（Gukov–Nawata–Saberi–Stošić–Sułkowski 1512.07883）【明記：要旨】。\[6\]\[45\]
4. **うまくいかない状況**：非長方形彩色では、トーラス結び目でも DAHA／refined CS 側に負の係数が出る（§2.6）。\[44\]\[46\]\[47\]そのため Poincaré 多項式との対応が崩れる【明記：Cherednik 1111.6195 §3.4】。ホモロジー側で非長方形彩色の三重次数ホモロジーを一般に構成したものは、今回の範囲では確認できなかった。

### 2.5 superpolynomial（DGR 型）
1. **構成**：Dunfield–Gukov–Rasmussen math/0505662（Exp. Math. 15 (2006)）。reduced\[48\] 超多項式 P(a,q,t) は非負の整係数をもち、十分大きな N について KhR_N(q,t)=P(a=q^N,q,t) を満たすことが要求される【明記】。微分\[49\] d_N（N∈Z）と対称性から組み合わせ的に決める。予想に基づく構成である。\[5\]
2. **規約**：DGR では標準のトーラス結び目 T_{a,b} は**負の交点**をもつ【明記】。HOMFLY\[49\] 多項式も P_DGR(q,a)=P(q,a^{-1})（Chandler–Gorsky §2.1）と規約が違う【明記】。比較では鏡像の取り違えに注意が必要。
3. **用語の区別**：本報告では「superpolynomial」を DGR の原義（三重次数ホモロジーの Poincaré 多項式で、正整数係数）に限る。Cherednik 系の「DAHA superpolynomial」、Kameyama–Nawata が言う「stable limit rCS（DAHA-superpolynomial とも呼ばれる）」【明記：1703.05408 §2】、Cherednik\[44\] 2607.25666 の「motivic superpolynomial」は、どれも負の係数を許す代数的・幾何的な多項式であり、区別して扱う。\[44\]
4. **うまくいかない状況**：Morozov 系の「cut-and-join 発展」による超多項式は、5_2 から 10_139 へ続く族では機能するが、4_1 では負の係数が出る（arXiv:1106.4305）【明記：要旨】。DAHA'KhR_2（DGR\[50\]\[51\] 流の「できるだけ全射」仮定による簡約）は、トーラス結び目では T(12,7) で初めて失敗する（Cherednik の先行論文による。CD 1408.4348 §2.5 の記述）【明記】。\[2\]\[52\]

### 2.6 refined Chern–Simons、DAHA、楕円 Hall 代数（ホモロジーではない）
1. **構成**：
   - Aganagic–Shakirov 1105.5117：Seifert 多様体上の refined CS 理論。\[53\]半自由な円作用を使い、トーラス結び目の不変量を定義する。sl(n) ホモロジーの Poincaré 多項式と一致すると予想【明記】。\[54\]
   - 続編の 1202.2489：指数（index）を計算しているので、Poincaré 多項式と一致する必要はないと明言している【明記】。\[55\]\[56\]
   - CMP 333 (2015) 版：基本表現のトーラス結び目について、指数が HOMFLY ホモロジーの Poincaré 多項式に等しいと予想【明記：要旨】。\[57\]
   - Cherednik 1111.6195：DAHA の PBW 定理と coinvariant による DAHA-Jones 多項式と DAHA 超多項式【明記】。\[53\]\[58\]\[59\]\[60\]
   - Gorsky–Negut（J. Math. Pures Appl. 104 (2015)、arXiv 番号は未確認）：安定化を証明し、AS と Cherednik の定義が同値であることを示した（Cor. 2.9）【明記】。\[61\]\[62\]
2. **負の係数（三葉結び目）**：Cherednik 1111.6195 §3.4 によると、T(3,2) の 3-hook 彩色 ω_1+ω_2 が負の係数が現れる最も簡単な例である。さらに、条件 (i)–(iii) を満たす正係数の多項式 H_{3,2}(b) は存在しえないことが示されている。Gorsky の観察では、HD/(1−t) を t の無限級数に展開すると正係数になる【明記】。\[46\]Kameyama–Nawata §4 も、非長方形彩色の rCS は変数変換 (4.3) の後でも一般に正負の係数をもつと明記している【明記】。\[44\]
3. **iterated torus への拡張**：
   - Cherednik–Danilenko 1408.4348：任意の既約ルート系と重みについて、任意の iterated torus 結び目に拡張。多項式性・双対性を示し、超多項式を構成（Thm 2.1–2.3）【明記】。\[2\]\[63\]t=q で HOMFLY-PT に帰着することは全ての iterated 結び目と重みについて予想され、sl_2 では証明されている（Conj. 2.4(i)、(2.29)）【明記】。\[2\]\[52\]
   - Morton–Samuelson 1410.0859：HOMFLYPT スケイン代数 ≅ 楕円 Hall 代数の特殊化を示した。iterated cable について λ-彩色 HOMFLYPT に特殊化する三変数多項式を構成し、それが CD の gl_N DAHA 多項式にも特殊化することを証明（CD の予想の一つを証明）【明記】。\[64\]\[65\]\[66\]同論文は、代数的な cabling 操作は標準的な位相的ケーブル構成と同じではないと注意している【明記】。\[64\]
   - Samuelson 1408.0483：sl_2 DAHA で iterated cable に拡張。\[61\]例として三葉結び目の (2,±5)-ケーブル【明記】。\[67\]CD はこの J_2 多項式が自分たちのものと大きく違い、q,t の分母をもつことを指摘している【明記：CD 序論】。\[2\]\[52\]
   - 絡み目への拡張は CD 1509.08351。代数的絡み目と Hopf 絡み目（DAHA-vertex）を扱う。(1−t)\[68\] のべきで割るだけでは、彩色付きの代数的絡み目で正値性を保証できない例を示している（§7）【明記：本文抜粋】。\[68\]\[69\]
4. **最新**：Cherednik 2607.25666「Instanton slices and their superpolynomials」（v1、84 頁）は、motivic 超多項式とインスタントン和を結びつけた。双曲結び目 K12n242、K12n725 をはじめ多くの非代数的結び目の超多項式に、motivic な解釈（スタックの F_q 点の数え上げ）を与えたと要旨にある【明記：要旨のみ、本文未読】。\[70\]\[71\]

### 2.7 幾何的な構成
1. **ORS**：Oblomkov–Rasmussen–Shende 1201.2115 は、平面曲線特異点の絡み目の HOMFLY ホモロジーの次元を、特異点に台をもつ点の Hilbert スキームの重み多項式で表す予想。t=−1 で Oblomkov–Shende の予想に帰着する【明記：要旨】。最小 a 次数の部分は、compactified Jacobian のコホモロジーの perverse filtration から復元される【明記】。\[2\]\[72\]\[73\]
2. **Euler 標数レベル**：Oblomkov–Shende 1003.1568（Duke 161 (2012)）は、トーラス結び目 y^k=x^n と、特異点 y^4=x^7−x^6+4x^5y+2x^3y^2（絡み目は三葉結び目の (2,13)-ケーブル）で予想を検証した【明記】。\[74\]\[75\]Maulik 1210.6323 が一般の場合に証明し、Diaconescu–Hua–Soibelman の彩色版も証明した【明記】。\[76\]
3. **Hilbert スキーム上の行列因子化・層**：Oblomkov–Rozansky の一連の論文（「HOMFLYPT homology of Coxeter links」1706.00124 など。他の arXiv 番号は未確認）。\[36\]Gorsky のスライドによれば、Hilb^n(C^2) 上の層 F_β の構成は Gorsky–Negut–Rasmussen、Gorsky–Hogancamp、Oblomkov–Rozansky によって提案されていて、これらは互いに同値と予想されている【明記：スライド】。\[30\]mirror symmetry 予想は Oblomkov–Rozansky と Gorsky–Hogancamp–Mellit が解決したとされる（2512.20597 の記述）【明記】。\[9\]
4. **適用範囲**：ORS は代数的絡み目（平面曲線特異点の絡み目）に限られる。非代数的ケーブルには、そもそも左辺の幾何的対象がない【導出】。CD は、Jacobian factor の Betti 数が知られているのはごく少数（Piontkowski の (4,2u,v) 系列など）だと述べている【明記：CD §2.2】。\[2\]\[52\]
5. **affine Springer fiber、character variety**：Jacobian factor は affine Springer fiber の一種である（CD §2.2 が Laumon を引用）【明記】。\[2\]\[52\]トーラス絡み目と wild character variety の関係は arXiv:1704.07412「BPS states, torus links and wild character varieties」（Diaconescu–Donagi–Pantev）が扱い、wild character variety のコホモロジーを refined stable pair 理論とトーラス絡み目の不変量に結びつけている【明記】。

### 2.8 knot Floer ホモロジーと d_0
- DGR は、三重次数ホモロジーから knot Floer への微分 d_0 を予想した。\[18\]Beliakova–Putyra–Robert–Wagner 2210.00878（JEMS 2025）は、その「適切に更新した版」を証明した。具体的には、reduced 三重次数ホモロジーと knot Floer ホモロジーが、中間の gl(0) ホモロジーを経由する**二つのスペクトル系列**で結ばれる【明記】。\[9\]\[77\]\[78\]したがって「単一の微分」としての元の形は証明されていない【導出】。
- 彩色版の knot Floer は arXiv:2508.21776（Alishahi–Gorsky–Liu、Gorsky のホームページで確認）【明記】。\[79\]
- Whitehead double の knot Floer（Hedden math/0606094）と thin な結び目のケーブル（Hanselman–Watson）は存在を引用で確認しただけで、本文は確認していない。

### 2.9 物理的な構成
1. **refined LMOV（Kameyama–Nawata 1703.05408）**：
   - **構成**：refined CS の large N 双対性を定式化し、トーラス結び目・絡み目の refined CS 不変量を、resolved conifold 上の M2–M5′ BPS 状態の次数付き次元と結びつけた。追加の U(1)_R が自明に作用するという仮定の下での**正値性予想**である【明記：要旨・§4】。\[44\]\[80\]\[81\]
   - **観察**：変数変換 (4.3)（a=−𝐚²𝐭、q^{1/2}=−𝐪𝐭、t^{1/2}=𝐪）の後、reformulated 不変量 f̂_ρ(T_{m,n}) の係数は全て非負整数になった。BPS 縮退度の表は T_{2,3} と T_{2,5}（Appendix C）【明記】。\[44\]
   - **非トーラス結び目（§6）**：入力は t_c 次数付きの彩色 HOMFLY ホモロジーの Poincaré 多項式である。(2s−1,1,2t−1)-pretzel 結び目（1512.07883 §5.3 の閉公式）と 6_2、6_3（Nawata–Oblomkov 1510.01795 §2.3）について、二箱までで正値の形 (4.5)/(4.6) に書けた。8 の字結び目 4_1 については、mirror symmetry (2.2)・指数増大 (2.3)・amphichirality (6.3) から三箱の refined 不変量 (6.4) を構成した。この (6.4) 自体は負の係数をもつが、reformulated 不変量はやはり正値の形になった【明記：サブエージェントによる本文確認】。\[44\]
   - **破れ**：thick な 9_42 では □ と二箱のデータから正値の形に書けない。Whitehead 絡み目を含む非トーラス絡み目では、基本表現でも拡張が成り立たない（§6、§7）【明記：サブエージェント確認】。\[44\]
2. **knots–quivers**：Kucharski–Reineke–Stošić–Sułkowski 1707.02991 は、LMOV 不変量を quiver の motivic DT 不変量で表し、LMOV 整数性の証明を与えた【明記：要旨】。1707.04017\[82\]\[83\] は無限族、6 交点以下の全結び目、いくつかの thick な結び目について quiver を同定した【明記：要旨の転載】。refined\[84\] 版では非負性が本質的になる（PRD 104, 106012 (2021)「Refined open topological strings revisited」。著者は未確認）【明記】。\[85\]\[86\]
3. **整数性のみの検証**：arXiv:2410.10468（Singh–Chbili、サブエージェントによる確認）は、[8_20]_{2p+1}、[8_21]_{2p+1}、[3_1]_{2p+1}、ツイスト結び目について |R|≤2 で整数性を確認し、破れはなかった。ただし\[87\] t 次数のない LMOV 整数なので、正値性の検証ではない【明記】。\[88\]
4. **3d–3d、Ẑ**：今回は検索予算の都合で、3d–3d 対応の原典は確認できなかった。Ekholm–Gruen–Gukov–Kucharski–Park–Sułkowski「Ẑ at large N: from curve counts to quantum modularity」2005.13349 は、(2,2p+1)-トーラス結び目について a-変形した F_K の閉じた式を与え、さらに超多項式にもとづく t-変形を提案している【明記】。Kucharski「Quivers for 3-manifolds」2005.13394 は題名のみ確認で、内容は未確認。KN は、L(p,1) の解析接続 WRT 不変量の refinement で同じ正値性現象があると述べている（GPV16、GPPV17 を引用）【明記】。\[44\]

### 2.10 計算機による直接計算
| プログラム | 対象 | 範囲（確認済み） |
|---|---|---|
| JavaKh（Bar-Natan–Green、後に Morrison が更新）\[89\] | Khovanov | T(8,7)（48 交点）を数分（math/0606318）\[13\] |
| KhoHo（Shumakovitch） | Khovanov（整係数、odd、reduced）\[89\] | 17–19 交点まで実用的。16 交点以下の全素結び目で既知（math/0405474）\[3\] |
| khoca（Lewark–Lobb） | sl(N)、同変版 | T(7,8) を F_7 で（README）。彩色 Khovanov データベースにも使用（2505.03916）\[11\]\[90\] |
| KnotJob（Schütz） | odd Khovanov、Sq^1 精密化 | 11 交点以下の計算表など（Schütz の論文）\[12\] |
| kr-calc（Nakagane–Sano） | HOMFLY-PT 三重次数 | 11 交点以下の全素結び目（2111.00388） |
| Colored-Torus-Knot-Homology（Conners） | 彩色トーラス結び目 | 彩色された正のトーラス結び目と一部の正のトーラス絡み目の reduced HOMFLY ホモロジーの再帰公式。彩色について正確に指数増大することを示し、Gorsky–Gukov–Stošić の 2013 年の予想を解決（Conners 2512.20597） |
| Artamonov のソフト | 彩色 HOMFLY 多項式 | CD 1408.4348 §3.1 で使用 |

## 3. 結び目のクラス × アプローチ の表

| クラス | Khovanov | sl(N) KR | HOMFLY-PT ホモロジー | 彩色ホモロジー | DAHA／refined CS | ORS 幾何 | refined BPS |
|---|---|---|---|---|---|---|---|
| トーラス結び目 | 計算あり（T(8,7) まで）。p,q≥3 は thick | 安定版で E_2 を計算（Gorsky–Hogancamp 系）\[19\] | **全て決定**（Mellit） | Sym^l は定理（HM）、長方形は予想（GGS） | 全て計算可。長方形は正（予想・観察）、非長方形は負あり | 証明済み（GMV、HM） | 正値性を観察（KN） |
| 2 橋結び目 | ≤16 交点で計算済み | N>4 で KR-thin（Rasmussen）\[24\] | ≤11 交点で計算済み | ツイスト結び目で予想的な公式（1209.1409）\[6\] | 定義なし（DAHA は iterated torus 専用） | 代数的でないので対象外 | 4_1（三箱まで）で正の形 |
| pretzel 結び目 | ≤16 交点で計算済み | — | ≤11 交点。HOMFLY ホモロジーは 5_2 と各 P(−3,3,2n+1) を判別する（Baldwin–Sivek 2211.03743、新しい計算表ではない） | (2s−1,1,2t−1) 族で閉公式（1512.07883） | 定義なし | 対象外 | (2s−1,1,2t−1) 族で正の形（KN §6） |
| 代数的な iterated torus | CD が KnotTheory で比較：Cab(13,2)、Cab(15,2)T(3,2) は DAHA'Kh＝Kh。Cab(19,3)T(3,2)、Cab(25,2)T(4,3) は差あり（予想の方向） | 直接計算なし（今回の範囲） | 直接計算なし。未解決問題（HM Problem 1.10） | なし | 正値（予想・計算例多数）。Betti 数の予測 | Euler 標数版は証明済み（Maulik）。(2,13)-ケーブルで検証（OS） | 検証例は見つからず |
| 非代数的な iterated torus | CD：Cab(7,2)T(3,2) で一致。Cab(5,2)、Cab(1,2)、Cab(−1,2) では −1↦1/t の置換で実際の Kh が得られる\[2\] | なし | なし | なし | Cab(7,2) まで正、Cab(5,2)T(3,2) 以下で負 | 定義なし | 見つからず |
| その他のサテライト（Whitehead double など） | Anokhina–Morozov–Popolitov 2104.14491 が、最も簡単な二つのサテライト族（ひねった Whitehead double と 2-ケーブル）の Khovanov 多項式を katlas.org のプログラムで明示的に計算し、パターンと companion への分解を予想 | なし | なし | なし | 定義なし | 対象外 | Whitehead **絡み目**で破れ（KN） |
| thick な結び目 | 多数 | N=2 で高次微分の可能性 | ≤11 交点で計算済み（9_42、10_125、11n135 など）\[25\] | — | — | — | 9_42 で破れ（KN） |
| 非トーラス絡み目 | 計算可能 | — | y-ified 版あり | 1512.07883 | CD 1509.08351（代数的絡み目） | 予想は代数的絡み目に適用 | Whitehead 絡み目で破れ |

## 4. 「特に詳しく調べてほしい点」への回答

### 4.1 iterated torus 結び目とサテライト：代数的か非代数的か
- **ケーブルの規約（重要）**：CD 1408.4348 §2.1 では、Cab(a,b)(K) は K の周りの小さなトーラス上に T(a,b) を K の向きに沿って置いたもの。iterated 結び目 𝒯(r⃗,s⃗) は、a_1=s_1、a_i=a_{i−1}r_{i−1}r_i+s_i を使って Cab(a_ℓ,r_ℓ)⋯Cab(a_2,r_2)(T(a_1,r_1)) と定義される（(2.2)–(2.3)）【明記】。\[2\]\[52\]この式と Cab(13,2)T(3,2)＝𝒯({3,2},{2,1}) から、Cab(a,2) は「2 本の撚り糸、Seifert framing で測った傾き a」のケーブルである。これは Oblomkov–Shende の「三葉結び目の (2,13)-ケーブル」が代数的であることとも整合する【導出】。代数的であるための条件は全ての r_i,s_i>0 であり、族 Cab(2m+1,2)T(3,2) では 2m+1≥13 にあたる。CD は Cab(13,2) をこの族で最初の代数的結び目と明記している【明記】。鏡像については Cab(a,b)(K*)＝(Cab(−a,b)(K))* が成り立つ（(2.9)）【明記】。\[2\]\[52\]
- **ユーザー指定の例**（三葉結び目の (2,q)-ケーブル、q は Seifert framing での傾き）：

| q | 代数的か | DAHA 超多項式（CD） | Khovanov との比較（CD） | 出典 |
|---|---|---|---|---|
| 13 | 代数的 | 正、a 次数 3。ℋ(q=1,t,a=0)=1+3t+4t²+4t³+4t⁴+3t⁵+2t⁶+t⁷+t⁸ が Piontkowski の (4,6,7) の Betti 数と一致\[52\] | DAHA'Kh＝実際の Kh | 1408.4348 §3.1【明記】 |
| 7 | 非代数的（pseudo-algebraic） | 正。ただし a 次数は s_1r_2−1=3 より小さい\[2\] | DAHA'Kh＝実際の Kh | §4.1【明記：抜粋】 |
| 5 | 非代数的 | **負の係数**（族で最初）。DAHA'Kh=…−q¹²t⁹−q¹⁴t¹⁰−q¹⁸t¹²−q²⁰t¹³\[2\] | −1↦1/t の置換で実際の Kh が得られる | §4.1【明記：抜粋】 |
| 3 | 非代数的 | ℋ(q=0,t,a=0)=1 とする記述がある（「2m+1≥3」）。個別の式は未確認 | 未確認 | §4.1（未確認） |
| 1 | 非代数的 | 負。ℋ(q=0,t,a=0)=1+a/t が初めて現れる\[2\] | tilde 正規化に問題あり\[2\] | §4.1【明記：抜粋】 |
| −1 | 非代数的 | 負。DAHA'Kh=1−q⁴+q²t+q⁴t²\[2\] | 置換規則が適用される | §4.1【明記：抜粋】 |

  CD の言葉では、Cab(5,2)T(3,2) が「the first with non-positive」DAHA 超多項式をもつ。\[2\]
- **T(5,2) の (2,1)-ケーブル**（CD の記法では Cab(1,2)T(5,2) にあたると導出）：今回確認できた CD の本文と抜粋には、この結び目の計算はなかった。CD が扱う T(5,2) のケーブルは代数的な Cab(21,2)T(5,2)（Puiseux 指数 (4,10,11)）と Cab(32,3)T(5,2)（Betti 多項式 (3.6)、Euler 数 3031）で、どちらも正である【明記】。未確認事項とする。\[2\]\[52\]
- **その他の代数的ケーブル（CD §3.2）**：Cab(19,3)T(3,2)（(6,9,10)、δ=21、Euler 数 365）、Cab(25,2)T(4,3)（(6,8,9)、Euler 数 227 が Piontkowski の公式 229/2+25v/2 の v=9 と一致）、Cab(20,3)T(3,2)、Cab(25,4)T(3,2)、二重 iterated の Cab(53,2)Cab(13,2)T(3,2)。Cab(19,3)T(3,2) と Cab(25,2)T(4,3) では DAHA'Kh と実際の Kh の差が q^iт^j(1−t²) 型の項になり、CD はこれを予想の確認とみなしている【明記】。さらに、−1 を t^{−1} に置き換えるだけでは Kh が得られない最も簡単な反例は Cab(7,2)T(4,3) だと明記されている【明記】。\[2\]\[52\]
- **アプローチごとの違い**：代数的ケーブルには、DAHA（正値）、ORS（Euler 標数版は定理）、Morton–Samuelson が全て適用できる。非代数的ケーブルでは、DAHA と楕円 Hall 代数の多項式は定義されるが正値性が崩れ、ORS は定義されない。ホモロジー側では、Khovanov ホモロジーだけが（KnotTheory などで）計算されている【導出】。

### 4.2 DAHA／refined CS が非代数的 iterated torus でホモロジーと一致しないという報告
- CD 1408.4348 の要旨は、KR 多項式との一致を「非負係数の場合」に限っている。Conj. 2.4(ii) は無彩色の pseudo-algebraic な結び目に限定されている【明記】。一般の iterated 結び目での DAHA と KhR の関係は「challenge」と書かれている【明記：序論】。\[2\]\[52\]\[63\]
- 非代数的ケーブルの具体例は 4.1 の表の通り（Cab(5,2)、Cab(1,2)、Cab(−1,2)T(3,2)）。ただし、ここで比較されているのは sl(2) Khovanov との比較（DAHA'Kh）であり、三重次数ホモロジーそのものとの比較ではない。三重次数ホモロジーはこれらの結び目では計算されていない【導出】。
- CD 1509.08351 は、代数的**絡み目**の彩色版でも (1−t) のべきで割る正規化だけでは正値性が足りない例を示している【明記：抜粋】。\[69\]
- 後続：Morton–Samuelson（楕円 Hall 代数で同じ多項式に特殊化）、Samuelson（sl_2、分母あり）、Cherednik–Elliot 1503.01441（トーラス結び目の composite 版）、Cherednik\[91\] 2607.25666（非代数的な双曲結び目の motivic 超多項式）。今回の検索範囲では、非代数的ケーブルの\[70\] DAHA 超多項式を三重次数ホモロジーの直接計算と比較した論文は見つからなかった。
- **ユーザーの研究（AS の非トーラス結び目への拡張）への含意【推測】**：非代数的ケーブルについて refined CS 型の不変量が負の係数をもつことは、それ自体が構成の誤りを示すわけではない。CD でも HOMFLY への特殊化は成り立っている。妥当性の検証基準は、(a) t=q での彩色 HOMFLY-PT との一致、(b) 代数的ケーブルでの正値性と ORS の Betti 数との一致、(c) 非代数的ケーブルでの CD 型の負の係数のパターン（たとえば Cab(5,2)T(3,2)）の再現、の三段階に分けるのがよい。

### 4.3 refined BPS（refined LMOV）不変量の正値性
- **成立が観察されたもの**：トーラス結び目・トーラス絡み目の任意の彩色。T_{2,3} と T_{2,5} の表（1703.05408 §4–5、Appendix C）。\[44\]thin な非トーラス結び目：(2s−1,1,2t−1)-pretzel、6_2、6_3（二箱まで）、4_1（三箱、(6.4) を構成して）（§6）【明記】。\[44\]
- **破れたもの**：9_42（thick）、Whitehead 絡み目を含む非トーラス絡み目（§6、§7）【明記：サブエージェント確認】。\[44\]
- **後続**：2410.10468 は整数性のみを検証（正値性は検証せず）。1902.02275（Kameyama–Nawata–Tao–Zhang）は要旨を見る限り正値性の検証ではない。\[88\]\[92\]\[93\]refined な knots–quivers 系（PRD 104, 106012）は、strip 幾何などで refined open BPS の非負性を示している【明記】。\[86\]\[94\]ケーブルや iterated torus 結び目での refined BPS の正値性検証は、今回の検索範囲（「refined BPS invariants positivity iterated torus knots cables」「Kameyama Nawata … follow-up」）では見つからなかった。

### 4.4 HOMFLY-PT ホモロジーの直接計算の最新範囲
- **非トーラス結び目**：11 交点以下の全素結び目（Nakagane–Sano、kr-calc、JKTR 2025）【明記】。\[1\]\[18\]\[29\]
- **族**：全トーラス結び目（Mellit）、正のトーラス絡み目と Sym^l 彩色トーラス結び目（Hogancamp–Mellit）、4 本の撚り糸の Coxeter braid（Turner）【明記】。\[27\]\[28\]\[34\]
- これより新しい結び目表の拡張（12 交点など）は、今回の検索範囲では見つからなかった。

## 5. 未確認の事項と、さらに調べるべき文献
- CD 1408.4348 §4.1 の全文（pp. 35–40）：Cab(3,2)T(3,2) の個別の式、式 (3.11)–(3.12) が対応する結び目（Figure 1 の Cab(11,2)T(3,2) か）、(4.16) の Samuelson との比較。
- T(5,2) の (2,1)-ケーブル：DAHA、Khovanov、HOMFLY のいずれでも計算例を確認できなかった。
- Lee、Bar-Natan 変形、odd Khovanov、Lipshitz–Sarkar 安定ホモトピー型、foam 評価（Robert–Wagner、Queffelec–Rose）の原典の arXiv 番号。
- 3d–3d 対応の原典と、Ẑ・F_K のケーブルでの計算例の有無。
- Ekholm–Kucharski–Longhi 1811.03110 と refined knots–quivers 系での正値性の破れの有無。
- Dowlin の Khovanov→knot Floer スペクトル系列（JAMS 2024）の arXiv 番号。
- Gorsky–Negut（J. Math. Pures Appl. 2015）、Rasmussen の 2 橋 KR 論文（Duke）、Wedrich の彩色 HOMFLY-PT の指数増大の arXiv 番号。
- 2303.16271、2004.15024、PRD 104 106012 の著者。

## 6. 文献一覧（確認状況つき）
| arXiv | 著者 | 題名 | 確認状況 |
|---|---|---|---|
| 1105.5117 | M. Aganagic, S. Shakirov | Knot Homology from Refined Chern-Simons Theory | 概要ページで確認 |
| 1202.2489 | M. Aganagic, S. Shakirov | Refined Chern-Simons Theory and Knot Homology | PDF で確認 |
| 1111.6195 | I. Cherednik | Jones polynomials of torus knots via DAHA | 概要ページと PDF で確認 |
| 1408.4348 | I. Cherednik, I. Danilenko | DAHA and iterated torus knots | 概要ページと PDF v2 で確認（§4.1 は抜粋のみ） |
| 1509.08351 | I. Cherednik, I. Danilenko | DAHA approach to iterated torus links | 概要ページと PDF 抜粋で確認 |
| 1503.01441 | I. Cherednik, R. Elliot | Refined composite invariants of torus knots via DAHA | 概要ページで確認 |
| 2607.25666 | I. Cherednik | Instanton slices and their superpolynomials | 概要ページで確認（本文未読） |
| 1408.0483 | P. Samuelson | Iterated torus knots and double affine Hecke algebras | PDF 抜粋で確認 |
| 1410.0859 | H. Morton, P. Samuelson | The Homflypt skein algebra of the torus and the elliptic Hall algebra | arXiv の HTML 版で確認 |
| 1703.05408 | M. Kameyama, S. Nawata | Refined large N duality for knots | 概要ページと PDF で確認 |
| 1902.02275 | M. Kameyama, S. Nawata, R. Tao, H. D. Zhang | Cyclotomic expansions of HOMFLY-PT colored by rectangular Young diagrams | ar5iv で確認 |
| 2410.10468 | Singh, Chbili（サブエージェント確認） | BPS Spectra of complex knots | PDF で確認 |
| 1707.02991 | P. Kucharski, M. Reineke, M. Stošić, P. Sułkowski | BPS states, knots and quivers | 概要ページで確認 |
| 1707.04017 | 同上 | Knots-quivers correspondence | PDF で確認 |
| 1811.03110 | Ekholm, Kucharski, Longhi | Physics and geometry of knots-quivers correspondence | 題名のみ確認 |
| math/0505662 | N. Dunfield, S. Gukov, J. Rasmussen | The superpolynomial for knot homologies\[45\] | 掲載誌の PDF で確認 |
| math/0607544 | J. Rasmussen | Some differentials on Khovanov-Rozansky homology\[4\] | 概要ページで確認 |
| math/0510265 | M. Khovanov | Triply-graded link homology and Hochschild homology of Soergel bimodules\[6\] | 文献リストで確認 |
| 1603.00407 | B. Elias, M. Hogancamp | On the computation of torus link homology | PDF で確認 |
| 1704.07630 | A. Mellit | Homology of torus knots | 概要ページで確認 |
| 1909.00418 | M. Hogancamp, A. Mellit | Torus link homology | 概要ページと PDF で確認 |
| 1712.03938 | E. Gorsky, M. Hogancamp | Hilbert schemes and y-ification of Khovanov-Rozansky homology\[36\] | 文献リストで確認 |
| 2103.01212 | E. Gorsky, M. Hogancamp, A. Mellit | Tautological classes and symmetry in Khovanov-Rozansky homology\[36\] | 文献リストで確認 |
| 2111.00388 | K. Nakagane, T. Sano | Computations of HOMFLY homology | 概要ページと PDF で確認 |
| 2209.13058 | A. Chandler, E. Gorsky | Structures in HOMFLY-PT homology\[18\] | PDF で確認 |
| 2410.03068 | J. P. Turner | Triply Graded Link Homology for Coxeter Braids on 4 Strands | PDF で確認 |
| 2602.17435 | T. Sano | A y-ification of Khovanov homology | PDF で確認 |
| 1112.0030 | S. Gukov, M. Stošić | Homological algebra of knots and BPS states\[7\] | 概要ページで確認 |
| 1304.3481 | E. Gorsky, S. Gukov, M. Stošić | Quadruply-graded colored homology of knots\[8\] | 概要ページで確認 |
| 1209.1409 | S. Nawata, P. Ramadevi, Zodinmawia, X. Sun | Super-A-polynomials for twist knots\[95\] | 掲載誌の概要で確認 |
| 1512.07883 | S. Gukov, S. Nawata, I. Saberi, M. Stošić, P. Sułkowski | Sequencing BPS spectra\[42\] | 掲載誌の概要で確認 |
| 1106.4305 | 未確認 | Superpolynomials for torus knots from evolution induced by cut-and-join operators | HTML で確認 |
| 1003.1568 | A. Oblomkov, V. Shende | The Hilbert scheme of a plane curve singularity and the HOMFLY polynomial of its link | 概要ページで確認 |
| 1201.2115 | A. Oblomkov, J. Rasmussen, V. Shende | The Hilbert scheme of a plane curve singularity and the HOMFLY homology of its link | 概要ページで確認 |
| 1210.6323 | D. Maulik | Stable pairs and the HOMFLY polynomial | PDF で確認 |
| 1706.00124 | A. Oblomkov, L. Rozansky | HOMFLYPT homology of Coxeter links | 文献リストで確認 |
| 2210.00878 | A. Beliakova, K. Putyra, L.-H. Robert, E. Wagner | A proof of Dunfield-Gukov-Rasmussen conjecture\[96\] | 著者ページと引用で確認 |
| 2508.21776 | A. Alishahi, E. Gorsky, B. Liu | Colored knot Floer homology: structures and examples\[79\] | PDF とホームページで確認 |
| math/0606318 | D. Bar-Natan | Fast Khovanov homology computations\[15\] | 著者ページの PDF で確認 |
| math/0405474 | A. Shumakovitch | Torsion of the Khovanov homology\[3\] | PDF で確認 |
| 0708.2228 | M. Mackaay, M. Stošić, P. Vaz | sl(N) link homology using foams and the Kapustin-Li formula\[20\] | 文献リストで確認 |
| 0905.0486 | B. Webster, G. Williamson | A geometric construction of colored HOMFLYPT homology\[6\] | 文献リストで確認 |
| 0907.0695 | H. Wu | A colored sl(N)-homology for links in S^3\[20\] | 文献リストで確認 |
| 2304.13613 | M. Kegel, N. Manikandan, L. Mousseau, M. Silvero | Khovanov homology of positive links and of L-space knots | PDF で確認 |
| 2104.14491 | L. Anokhina, A. Morozov, A. Popolitov | Khovanov polynomials for satellites and asymptotic adjoint polynomials | HTML で確認 |
| 2205.11430 | V. Chernov, M. Maguire | Conjectures on the Khovanov Homology of Torus Knots, Twist Knots, and Legendrian Simple Knots | PDF で確認 |
| 2211.03743 | J. Baldwin, S. Sivek | （題名未確認：HOMFLY ホモロジーによる 5_2 と P(−3,3,2n+1) の判別） | PDF で確認 |
| 1704.07412 | W.-Y. Diaconescu, R. Donagi, T. Pantev | BPS states, torus links and wild character varieties | PDF で確認 |
| 2505.03916 | K. Ritter von Merkl ほか | Computing colored Khovanov homology\[90\] | PDF で確認 |
| 2512.20597 | L. Conners | Doubly-graded exponential growth of colored torus knot homology | PDF で確認 |
| 1510.01795 | S. Nawata, A. Oblomkov | Lectures on knot homology | 文献リストで確認 |
| 2005.13349 | T. Ekholm, A. Gruen, S. Gukov, P. Kucharski, S. Park, P. Sułkowski | Ẑ at large N: from curve counts to quantum modularity | arXiv で確認 |

## 出典

1. [\[2111.00388\] Computations of HOMFLY homology](https://arxiv.org/abs/2111.00388)
2. [DAHA and iterated torus knots](https://arxiv.org/pdf/1408.4348)
3. [arXiv:math/0405474v2 \[math.GT\] 19 Jun 2018 TORSION OF THE KHOVANOV HOMOLOGY](https://arxiv.org/pdf/math/0405474)
4. [\[math/0607544\] Some differentials on Khovanov-Rozansky homology](https://arxiv.org/abs/math/0607544)
5. [Computations of HOMFLY homology](https://arxiv.org/pdf/2111.00388)
6. [Super-A-polynomials for twist knots](<https://link.springer.com/article/10.1007/JHEP11(2012)157>)
7. [\[1112.0030\] Homological algebra of knots and BPS states](https://arxiv.org/abs/1112.0030)
8. [\[1304.3481\] Quadruply-graded colored homology of knots](https://arxiv.org/abs/1304.3481)
9. [Emmanuel Wagner's research works](https://www.researchgate.net/scientific-contributions/Emmanuel-Wagner-2029980674)
10. [BEHAVIOR OF KNOT INVARIANTS UNDER GENUS 2 MUTATION](https://people.mpim-bonn.mpg.de/stavros/publications/mutation.pdf)
11. [GitHub - LLewark/khoca: Khoca is a C++/python program for research in knot theory, a branch of pure mathematics. It computes certain Khovanov-Rozansky homologies of knots. · GitHub](https://github.com/LLewark/khoca)
12. [A SCANNING ALGORITHM FOR ODD KHOVANOV HOMOLOGY DIRK SCH¨UTZ](https://maths.dur.ac.uk/users/dirk.schuetz/d37.pdf)
13. [FAST KHOVANOV HOMOLOGY COMPUTATIONS DROR BAR-NATAN](https://www.math.toronto.edu/~drorbn/papers/FastKh/FastKh.pdf)
14. [FAST KHOVANOV HOMOLOGY COMPUTATIONS](https://www.researchgate.net/publication/2128713_Fast_Khovanov_Homology_Computations)
15. [Khovanov homology theories and their applications](https://arxiv.org/pdf/1101.5614)
16. [Homological thickness and stability of torus knots](https://projecteuclid.org/journals/algebraic-and-geometric-topology/volume-7/issue-1/Homological-thickness-and-stability-of-torus-knots/10.2140/agt.2007.7.261.pdf)
17. [A y-ification of Khovanov homology Taketo Sano February 20, 2026 Abstract](https://arxiv.org/pdf/2602.17435)
18. [Khovanov-Rozansky homology over \$\\mathbb{F}\[U,V\]\$](https://arxiv.org/pdf/2608.29986)
19. [(PDF) On Stable Khovanov Homology of Torus Knots](https://www.researchgate.net/publication/225291404_On_Stable_Khovanov_Homology_of_Torus_Knots)
20. [Super-A-polynomials for Twist Knots](https://arxiv.org/pdf/1209.1409)
21. [Some differentials on Khovanov-Rozansky homology](https://arxiv.org/pdf/math/0607544)
22. [Some differentials on Khovanov–Rozansky homology](https://projecteuclid.org/journals/geometry-and-topology/volume-19/issue-6/Some-differentials-on-KhovanovRozansky-homology/10.2140/gt.2015.19.3031.full)
23. [The action of full twist on the superpolynomial for torus knots](https://www.researchgate.net/publication/324982202_The_action_of_full_twist_on_the_superpolynomial_for_torus_knots?_share=1)
24. [arXiv:2211.03743v2 \[math.GT\] 7 Jun 2023](https://arxiv.org/pdf/2211.03743)
25. [STRUCTURES IN HOMFLY-PT HOMOLOGY ALEX CHANDLER AND EUGENE GORSKY](https://arxiv.org/pdf/2209.13058)
26. [On the computation of torus link homology](https://arxiv.org/pdf/1603.00407)
27. [\[1704.07630\] Homology of torus knots](https://arxiv.org/abs/1704.07630)
28. [\[1909.00418\] Torus link homology](https://arxiv.org/abs/1909.00418)
29. [Taketo Sano's Website](https://www.taketo-sano.me/)
30. [Algebra, geometry and combinatorics of link homology Eugene Gorsky](https://www.math.ucdavis.edu/~egorskiy/Presentations/Montreal_slides.pdf)
31. [(PDF) Structures in HOMFLY-PT homology](https://www.researchgate.net/publication/363888773_Structures_in_HOMFLY-PT_homology)
32. [Khovanov-Rozansky homology of two-bridge knots and links](https://www.researchgate.net/publication/2121885_Khovanov-Rozansky_homology_of_two-bridge_knots_and_links)
33. [Structures in HOMFLY-PT Homology: Experimental Mathematics: Vol 33, No 4](https://www.tandfonline.com/doi/abs/10.1080/10586458.2023.2300412)
34. [Triply Graded Link Homology for Coxeter Braids on 4 Strands](https://arxiv.org/pdf/2410.03068)
35. [Triply graded link homology: results and structures Eugene Gorsky](https://www.math.ucdavis.edu/~egorskiy/Presentations/Georgia_slides.pdf)
36. [Row-Column Mirror Symmetry for Colored Torus Knot Homology](https://arxiv.org/pdf/2303.16271)
37. [research statement - eugene gorsky](https://www.math.ucdavis.edu/~egorskiy/Research_statement_Gorsky_2022.pdf)
38. [Doubly-graded exponential growth of colored torus knot homology](https://arxiv.org/pdf/2512.20597)
39. [TORUS LINK HOMOLOGY MATTHEW HOGANCAMP AND ANTON MELLIT](https://arxiv.org/pdf/1909.00418)
40. [KHOVANOV HOMOLOGY OF POSITIVE LINKS AND OF L-SPACE KNOTS](https://arxiv.org/pdf/2304.13613)
41. [Generalized affine Springer theory and Hilbert schemes on planar curves](https://arxiv.org/pdf/2004.15024)
42. [Cyclotomic expansions of HOMFLY-PT colored by rectangular Young diagrams](https://link.springer.com/article/10.1007/s11005-020-01318-5)
43. [The tilde-version of the S 2-colored homology of the trefoil knot (Q =...](https://www.researchgate.net/figure/The-tilde-version-of-the-S-2-colored-homology-of-the-trefoil-knot-Q-q-t-r-t-c_fig9_258841578)
44. [Refined large N duality for knots](https://arxiv.org/pdf/1703.05408)
45. [Sequencing BPS spectra](<https://link.springer.com/article/10.1007/JHEP03(2016)004>)
46. [Jones polynomials of torus knots via DAHA](https://arxiv.org/pdf/1111.6195)
47. [ANNALES DE LA FACULTÉ DES SCIENCES Mathématiques IVAN CHEREDNIK AND ROSS ELLIOT](https://www.numdam.org/item/AFST_2016_6_25_2-3_433_0.pdf)
48. [Knots, BPS States, and Algebraic Curves](https://link.springer.com/article/10.1007/s00220-016-2682-z)
49. [The Superpolynomial for Knot Homologies](https://projecteuclid.org/journals/experimental-mathematics/volume-15/issue-2/The-Superpolynomial-for-Knot-Homologies/em/1175789736.pdf)
50. [Evolution properties of the knot's defect](https://arxiv.org/pdf/2204.05977)
51. [Superpolynomials for torus knotsfrom evolution induced by cut-and-join operators](https://arxiv.org/html/1106.4305)
52. <https://arxiv.org/pdf/1408.4348v2>
53. [Colored Kauffman Homology and Super-A-polynomials](https://arxiv.org/pdf/1310.2240)
54. [\[1105.5117\] Knot Homology from Refined Chern-Simons Theory](https://arxiv.org/abs/1105.5117)
55. [Untitled Document](https://arxiv.org/html/1202.2489)
56. [arXiv:1202.2489v1 \[hep-th\] 12 Feb 2012 Reﬁned Chern-Simons Theory and](https://arxiv.org/pdf/1202.2489)
57. [Knot Homology and Refined Chern-Simons Index - ADS](https://ui.adsabs.harvard.edu/abs/2015CMaPh.333..187A/abstract)
58. [arXiv:1111.6195v2 \[math.QA\] 29 Nov 2011 JONES POLYNOMIALS OF TORUS KNOTS VIA](https://www.arxiv.org/pdf/1111.6195v2)
59. [\[1111.6195\] Jones polynomials of torus knots via DAHA](https://arxiv.org/abs/1111.6195)
60. [From Cherednik algebras to knot homology via cuspidal D- ...](https://arxiv.org/pdf/2407.00971)
61. [DAHA and iterated torus knots](https://www.researchgate.net/publication/264899405_DAHA_and_iterated_torus_knots)
62. [Refined knot invariants and Hilbert schemes - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0021782415000367)
63. [\[1408.4348\] DAHA and iterated torus knots](https://arxiv.org/abs/1408.4348)
64. [The Homflypt skein algebra of the torus and the elliptic Hall algebra](https://arxiv.org/html/1410.0859)
65. [THE HOMFLYPT SKEIN ALGEBRA OF THE TORUS AND THE ELLIPTIC HALL ALGEBRA](https://www.liverpool.ac.uk/~su14/papers/MortonSamuelsonDuke.pdf)
66. [THE HOMFLYPT SKEIN ALGEBRA OF THE TORUS AND THE ELLIPTIC HALL ALGEBRA](https://livrepository.liverpool.ac.uk/3001915/1/homfly_morton_samuelson_arxivV2_may2015.pdf)
67. [Iterated torus knots and double affine Hecke algebras](https://arxiv.org/pdf/1408.0483)
68. [\[1509.08351\] DAHA approach to iterated torus links](https://arxiv.org/abs/1509.08351)
69. [DAHA approach to iterated torus links](https://arxiv.org/pdf/1509.08351)
70. [\[2607.25666\] Instanton slices and their superpolynomials](https://arxiv.org/abs/2607.25666)
71. [Algebraic Geometry](https://arxiv.org/list/math.AG/pastweek?skip=14&show=1000)
72. [The Hilbert scheme of a plane curve singularity and the HOMFLY homology of its link](https://escholarship.org/uc/item/31f9w49x)
73. [\[1201.2115\] The Hilbert scheme of a plane curve singularity and the HOMFLY homology of its link](https://arxiv.org/abs/1201.2115)
74. [The Hilbert scheme of a plane curve singularity and the HOMFLY polynomial of its link](https://projecteuclid.org/euclid.dmj/1336142076)
75. [\[1003.1568\] The Hilbert scheme of a plane curve singularity and the HOMFLY polynomial of its link](https://arxiv.org/abs/1003.1568)
76. [Stable pairs and the HOMFLY polynomial](https://link.springer.com/article/10.1007/s00222-015-0624-6)
77. [Anna BELIAKOVA](https://www.researchgate.net/profile/Anna-Beliakova)
78. [COLORED KNOT FLOER HOMOLOGY: STRUCTURES AND EXAMPLES](https://arxiv.org/pdf/2508.21776)
79. [Homepage of Eugene Gorsky](https://www.math.ucdavis.edu/~egorskiy/publications.html)
80. [\[1703.05408\] Refined large N duality for knots](https://arxiv.org/abs/1703.05408)
81. [Refined large N duality for knots](https://www.worldscientific.com/doi/10.1142/S0218216520410011)
82. [\[1707.02991\] BPS states, knots and quivers](https://arxiv.org/abs/1707.02991)
83. [arXiv:2505.02059v2 \[math.GT\] 12 Jun 2025 A survey of knots and quivers](https://arxiv.org/pdf/2505.02059)
84. [(PDF) Knots-quivers correspondence](https://www.researchgate.net/publication/386914868_Knots-quivers_correspondence)
85. [Refined open topological strings revisited](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.104.106012)
86. [Refined open topological strings revisited](https://doi.org/10.1103/physrevd.104.106012)
87. [BPS Spectra of complex knots](https://arxiv.org/pdf/2410.10468)
88. <https://arxiv.org/html/2410.10468>
89. [Search for torsion in Khovanov homology](https://arxiv.org/pdf/1701.04924)
90. [COMPUTING COLORED KHOVANOV HOMOLOGY KARIM RITTER VON MERKL](https://arxiv.org/pdf/2505.03916)
91. [\[1503.01441\] Refined composite invariants of torus knots via DAHA](https://arxiv.org/abs/1503.01441)
92. [\[1902.02275\] Cyclotomic expansions of HOMFLY-PT colored by rectangular Young diagrams](https://ar5iv.labs.arxiv.org/html/1902.02275)
93. <https://arxiv.org/abs/1902.02275>
94. [References For: Phys. Rev. D 104, 106012 (2021) - Refined open topological strings revisited](https://journals.aps.org/prd/references/10.1103/PhysRevD.104.106012)
95. [BPS counting for knots and combinatorics on words](<https://link.springer.com/article/10.1007/JHEP11(2016)120>)
96. [University of Zurich - Institut für Mathematik - Universität Zürich](https://math.uzh.ch/en/people?L=1&key1=578)
