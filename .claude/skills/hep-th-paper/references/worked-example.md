# 検出から照合までの通し例

tests/fixture_ai.tex（架空の段落）を、検出 → 報告 → 最小限の修正 → 照合の順に処理した記録である。スクリプトの出力は、スキルのディレクトリで 2026-11-16 に実際に実行したもの。

## 1. 検出

```
$ python3 scripts/check_tells.py tests/fixture_ai.tex --bib tests/fixture.bib
tests/fixture_ai.tex:4: [soft 構造] 内容のない導入: In recent years
tests/fixture_ai.tex:4: [soft 構造] 内容のない導入: garnered significant attention
tests/fixture_ai.tex:4: [soft 語彙] 同じ段落に 12 語: garnered, delve, intricate, novel, robust, showcasing, elegant, delve, leveraging, pivotal, underscoring, pave the way
tests/fixture_ai.tex:6: [soft 構造] 意義づけで終わる文: showcasing
tests/fixture_ai.tex:7: [soft 誠実性] 出典のない「よく知られている」: It is well known that
tests/fixture_ai.tex:14: [soft 誠実性] .bib のエントリに DOI も eprint もない（実在を原典で確認）: Macdonald95
tests/fixture_ai.tex:14: [hard 誠実性] \cite のキーが .bib にない: TODO-missing
tests/fixture_ai.tex:14: [soft 構造] 意義づけで終わる文: underscoring
tests/fixture_ai.tex:15: [soft 構造] not only ... but also: not only
tests/fixture_ai.tex:15: [soft 構造] まとめ文: Overall
tests/fixture_ai.tex:15: [soft 構造] ヘッジの重なり: may potentially
[soft 構造] 引用している文献が 3 本（20 本未満）
本文: 92 語、語彙 12 回（1,000 語あたり 130.43）、文長の平均 20.2 語・変動係数 0.39
（語ごとの回数は省略）
層別: 誠実性 hard 1・soft 2、構造 hard 0・soft 8、語彙 hard 0・soft 1、表層 hard 0・soft 0
hard 1 件、soft 11 件、引用 3 本
```

スクリプトが拾わなかったものを読んで足す。この段落では、"simple, elegant, and general" の三点列挙（1 回だけなのでスクリプトは出さない）と、ディスプレイ式の末尾に句読点がないこと。

## 2. 報告（ここでは編集しない）

| 層 | 確度 | 位置 | 該当箇所 | 根拠 | 直さなくてよい可能性 | 直し方 |
|---|---|---|---|---|---|---|
| 誠実性 | hard | fixture_ai.tex:14 | `\cite{..., TODO-missing}` | キーが .bib にない | なし | 著者が文献を特定する。inspire-bibtex で作る |
| 誠実性 | soft | :7 | "It is well known that the invariants are non-negative ..." | 文に引用がない | 直前の文献で示されているなら引用を移すだけでよい | 直さない。要出典として著者に確認する |
| 誠実性 | soft | :14 | `Macdonald95` | DOI も eprint もない | 書籍なので正当 | 書誌を原典で確認する |
| 構造 | soft | :4 | "In recent years, ... has garnered significant attention" | 内容のない導入 | なし | 先行研究が何を示したかから始める |
| 構造 | soft | :6, :14 | ", showcasing ...", ", underscoring the pivotal role of ..." | 意義づけで終わる文が 2 つ | なし | 削る |
| 構造 | soft | :15 | "Overall, this is not only ... but also may potentially pave the way ..." | まとめ文、対比、ヘッジの重なりが 1 文にある | なし | 引用との整合だけを述べる |
| 語彙 | soft | :4–16 | delve ×2, intricate, novel, robust, ... | 1 段落に 12 語（1,000 語あたり 130） | robust が技術的な意味なら残す（ここでは違う） | 平易な語にするか削る |
| 表層 | soft | :9–12 | ディスプレイ式 | 末尾に句読点がない | なし | ピリオドを付ける |

集計：誠実性 3、構造 4、語彙 1、表層 1。引用 3 本のうち未定義 1 本。直すのは構造・語彙・表層で、誠実性の 3 件は著者の確認事項として残す。

## 3. 最小限の修正

結果は tests/fixture_fixed.tex。数式・引用・参照・数値には触れていない。"It is well known that ..." の文は、引用を足すことも主張を弱めることもせず、そのまま残した。

## 4. 照合

```
$ python3 scripts/check_invariants.py --before tests/fixture_ai.tex tests/fixture_fixed.tex
数式・表: 一致（3 個）
引用: 一致（4 個）
参照: 一致（1 個）
ラベル: 一致（1 個）
本文の数値: 一致（2 個）
定理環境: 一致（0 個）
見出し: 一致（1 個）
不変量はすべて一致
```

式の末尾に付けたピリオドは差分にならない。検出をやり直すと、語彙は 12 回から 0 回、構造の soft は 8 件から 1 件（引用の本数）になり、誠実性の 3 件は残る。

## 5. まとめの報告

直したもの（構造 4、語彙 1、表層 1）、残したもの（なし）、著者が確かめるもの（未定義の引用 1、出典のない "well known" 1、書誌の確認 1）を挙げて終える。
