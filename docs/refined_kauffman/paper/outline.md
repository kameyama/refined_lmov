# 論文の構成案：refined Kauffman 多項式の large N 双対（2026-10-20、執筆前の準備）

もとにした文書：docs/refined_kauffman/report.md、tex/report.tex、theory_orientifold_index.md、progress.md（続き 1〜40）、literature.md、conventions.md。
印は【導出】【計算】【予想】【仮説】（CLAUDE.md の分類）。主張の根拠と格付けは claims.md に置く。

**決定事項（案。ユーザーの確認待ち）**
- 論文は 1 本。本研究（docs/paper/、framed knot の KN）とは別の論文にする。
- 準備の段階では日本語で書き、最終版は英語にする（docs/paper/ と同じ）。
- 著者は 亀山昌也 の単著（docs/paper/ と同じ扱い）。2018 年のノートブックの試行は「著者が以前に行った計算」と書く。
- 正値性は「未解決」として書き、障害を定式化する。正値性が成り立つと主張しない。
- 投稿先は決めない。

---

## 0. 論文の軸

**主張の骨格**
1. **定式化**：SO(2N) の refined CS（AgSc 1202.4456）の結び目不変量から、refined のクロスキャップ部分 g = f − ½f^{cov} を取り出す枠組みを作る。
   - 生成関数は双対 Cauchy 核（t-brane）。
   - 規格化は「一様な規格化」（行列式と荷電共役で不変）。
   - reformulation は Mariño の M(q^{−1})。
2. **構造**：全 d の構造 g = ½E/(1+t) + C を、orientifold の開弦の M 理論の指数（対・不変・凍結の分類）から導く。
   - 閉弦の厳密な結果（AgSc §6、AM 2205.12832）でも成り立つことを確かめる。
   - AM の「Klein 瓶」の項は、この構造の並べ替えであることを示す。
3. **計算**：DAHA-Jones 多項式（D_N と GL_N）で、refined SO(2N) の多箱を一般の N で求める。
   - 対象は T(2,3) の 2・3 箱と T(2,5) の 2 箱。
   - 整数性、(1+t) の極の位置（Adams）、t = q での Mariño との一致を示す。
4. **普遍因子**：refined だけの部分は (𝐭² − 1)(𝐚² − 𝐪²) で割り切れる。正確には ĝ − ĝ|_{t=q}(ν → ν(q/t)^{1/2}) が割り切れる。
   - 成り立つ範囲：[1]、2 箱、3 箱、Hopf のすべて。
   - ずれの向き (q/t)^{1/2} は一意に決まる。
5. **正値性の障害**：[1]（トーラス結び目）では判定式の形で符号がそろう。多箱ではそろわない。
   - 障害は E と C の分離と、凍結した状態の多重項の基底にある。
   - これは指数からは決まらない【導出】。
   - 試した方法と、その失敗の範囲を示す。

**超える対象（既知の結果）**
- Mariño 0904.1088 の unrefined の構造。trefoil・4_1 の数値は Mariño、NRZ 1310.2240 にある。
- AgSc 1202.4456 の refined SO(2N) の値（[1] の T(2,2m+1) の一般式、Hopf）と、大きな N の閉弦。
- AM 2205.12832 の閉弦の GV 型展開。
- Cherednik の DAHA-Jones 多項式の定義（1111.6195）。D 型の t = q での Kauffman との一致（Cherednik–Danilenko 1509.08351）。
- これらの再現は規約の固定と検証のためで、成果としては書かない。

**書き方の方針**
- 定理・命題は【導出】に限る。計算の結果は「計算結果」として範囲つきで書く。
- 仮説（h/g 分解の規格化、DAHA の値を refined CS の値とみなすこと、一様な規格化が framing 0 であること）は明記する。
- 正値性について、unrefined で符号が混在することを、refined で正値性を期待しない理由にしない（README の規則）。
- 失敗した試み（Mariño 型の素朴な refined 化、G 展開、基底の取り替え、固定セクターの整数計画など）は、障害の定式化に要るものだけを §7 に短く書き、詳細は付録または記録に回す。

---

## 1. 章立て

仮題（日本語の作業題）：refined Kauffman 不変量の large N 双対：orientifold の refined 指数とクロスキャップの構造
英題案：*Refined large N duality for Kauffman invariants: orientifold indices and the structure of refined crosscap amplitudes*

| 節 | 内容 | もとの文書 |
|---|---|---|
| 要旨 | 枠組み、構造、計算、普遍因子、正値性の障害 | report.tex 要旨 |
| §1 序論 | Mariño の Kauffman の large N 双対、AgSc の refined orientifold、KN（U(N) の refined large N）。本論文の問い（refined の SO で BPS の整数性・正値性はどう見えるか）と結果の一覧（表） | literature.md §1・§4 |
| §2 準備 | 変数、refined SO(2N) の不変量（AgSc、Koornwinder）、双対 Cauchy 核、被覆と composite、reformulation M(q^{−1})、一様な規格化 | report.tex §2・§5・§9 |
| §3 orientifold の refined 指数 | 対・不変・凍結の分類と g = ½E/(1+t) + C（命題）、E と C の分離の不定性（命題）、E を含まない方向（系）、閉弦での確認（AgSc §6 と AM の Klein 瓶） | theory_orientifold_index.md §1〜§4、§8 |
| §4 unknot と Hopf 絡み目 | 双対 Cauchy 型の整数性、unknot = Mariño 型 + (−1)、Hopf の ν~、(1,b) と hatted の hook の式、E_22・E_24 の決定と Adams | report.tex §6・§7 |
| §5 トーラス結び目の [1] | h/g 分解、一様な規格化、(N)(P)(c=1)(c=2)、Π_m の閉じた形 | report.tex §8・§9 |
| §6 一般の N の多箱 | DAHA-Jones（D_N と GL_N）の計算法と検証、補間、大きな N の 2 箱・3 箱、整数性と極、Mariño との辞書、普遍因子とずれの向き | report.tex §11・§12、theory §8.3 |
| §7 正値性の現状と障害 | [1] の判定式、多箱の符号の混在、E と C の分離、捻ったトレースの整数計画（T(2,5) で失敗）、composite の DAHA の非正値性、凍結した状態の多重項（未解決） | report.tex §11・§12、theory §5〜§10 |
| §8 議論 | 未解決の問題、他の結び目（入力の問題）、閉弦の Kähler の三重化との関係 | report.tex §14 |
| 付録 | 生成した式（[1]、2 箱の大きな N、Hopf）、DAHA の実装の要点、固定 N の検査 | tex/gen/*.tex |

---

## 2. 図表

- 表 1（序論）：対象ごとの整数性・構造・正値性の状況（report.tex §13 の表を更新）。
- 表 2（§6）：一般の N の補間の範囲（彩色、N のノード、照合点）。
- 表 3（§7）：正値性の判定の試みと結果（基底、二次量子化、固定セクター、捻ったトレース、composite）。

## 3. 残りの作業

1. claims.md の各主張について、テストと記録を確かめる。
2. 日本語の TeX 草稿（ja/main.tex、sections/）を書き、コンパイルする。
3. references.bib を作る。原典で確認した文献だけを載せ、ほかは「未確認」と注記する。
4. ユーザーの確認：決定事項（案）、題、英語版への移行の時期。
