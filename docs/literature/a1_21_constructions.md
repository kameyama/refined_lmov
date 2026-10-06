# open_problems A1：双曲結び目の [2,1] の refined 不変量を定義に基づいて得る枠組みの調査（2026-10-07）

調査はサブエージェントによる（arXiv の API・概要・本文の検索）。引用・転載はしていない。

## 結論
- 双曲結び目（4_1 など）の [2,1] 彩色の refined 不変量について、「定義があり、計算例があり、a 依存性もある」枠組みは見つからなかった。
  各候補は、彩色が対称・反対称・長方形・列（minuscule）に限られるか、階数 1（sl₂）しかないか、定義だけで計算例が無いかのいずれか。
- [2,1] で定義があるのはホモロジー側の Cautis（1611.09924、HOMFLY clasp による彩色三重次数ホモロジー）と Elias–Hogancamp（1801.00191、full twist の圏論的対角化による Young 冪等元）だけ。
  どちらも双曲結び目での計算例は無い（本文で確認）。得られるのは Poincaré 多項式（正値）で、KN の rCS（index）との関係は別に規約を決める必要がある。
- Arthamonov–Shakirov（1504.02620、genus 2 の refined CS）は rank 1（sl₂）のみで、a 依存性が無く [2,1] は扱えない。高階数化は未解決と本文に明記。
  twist 結び目（4_1, 6_1, 8_1 など）の値は超多項式（KhR の Poincaré 多項式）と一致しないと明記（refined CS は index であり Poincaré 多項式とは限らない）。
- その他（Mironov–Morozov の Racah・微分展開は t = q のみ、Hikami・Cooke–Samuelson は rank 1、Galashin–Lam は無彩色、Cherednik 2607.25666 は minuscule のみ、
  knots–quivers は対称表現のみ、Aganagic のミラー対称性による構成は minuscule）は、現時点で双曲結び目の [2,1] に届かない。

## 推奨（現実的な順）
1. Elias–Hogancamp／Cautis の [2,1] 射影子で彩色した三重次数ホモロジーを 4_1 で直接計算する。
   計算基盤は Nakagane–Sano の kr-calc（yui-kr）。3 本ケーブル（9 本・36 交点）に無限射影子を挿入するので計算量が大きく、Rust の導入（ユーザーの許可）も要る。
2. Arthamonov–Shakirov の genus 2 の refined CS（rank 1）を twist 結び目で再現し、a = t² の U(2) への簡約などの検証に使う。[2,1] そのものは得られない。

## 主な文献（確認方法：A = arXiv API、B = 概要、C = 本文検索）
1611.09924 Cautis（A,B,C）、1801.00191 Elias–Hogancamp（A,B）、1504.02620 Arthamonov–Shakirov（A,B,C）、1704.02947 Arthamonov–Shakirov（A,B）、
1211.6375 AMMM（A,B,C）、1612.00422 A. Morozov（A,B,C）、2609.09437 Kononov–Morozov（A,B）、2607.25666 Cherednik（A,B,C）、2405.00875 Conners（A,B,C）、
2505.05668 Kucharski–Noshchenko（A,B,C）、2207.14104 Aganagic（A）、1901.02743・2402.09106 Hikami（A,B）、2307.16794 Galashin–Lam（A）。ほか一覧はサブエージェントの報告による。
