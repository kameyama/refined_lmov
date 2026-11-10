---
name: inspire-bibtex
description: 引用文献の書誌情報を INSPIRE-HEP で調べ、INSPIRE が出す BibTeX から refs.bib のエントリを作る。論文に文献を足すとき、refs.bib を直すとき、arXiv 番号・掲載誌・巻・頁を確認するとき、「文献を調べて」「bibtex を作って」「引用を追加して」と言われたときに使う。書誌情報を記憶から書かない。
---

# INSPIRE で文献を調べて BibTeX を作る

書誌情報（著者、題名、掲載誌、巻、頁、年、DOI、arXiv 番号）は INSPIRE-HEP（https://inspirehep.net）の記載を使う。記憶から書かない。

## 手順

1. 文献ごとに「引用キー」と「識別子」を決める。識別子は arXiv 番号（`1703.05408`、`hep-th/0405083`）か DOI。どちらも分からなければ、先に INSPIRE で探す（下の「識別子が分からないとき」）。
2. スクリプトで BibTeX を取る。サンドボックスのネットワーク許可に `inspirehep.net` が要る。

   ```bash
   python3 .claude/skills/inspire-bibtex/inspire_bibtex.py KN17=1703.05408 CD14=10.2140/agt.2016.16.843 > $TMPDIR/new.bib
   ```

   - 出力は INSPIRE の BibTeX そのままで、キーだけを指定した引用キーに替え、`shorthand = {キー}` を足してある（このリポジトリの論文は biblatex の alphabetic で、ラベルを bib のキーと同じにする）。
   - 標準エラーに、文献ごとの結果が出る：`ok`（掲載誌あり）、`no-journal`（INSPIRE に掲載誌の情報が無い）、`not-found`（INSPIRE に記録が無い）、`ambiguous`（複数の記録に当たった）。
3. 出力を読み、著者と題名が意図した文献かを確かめてから、refs.bib の該当エントリを置き換える。フィールドを手で書き換えない（INSPIRE の題名が arXiv 版と違う、巻が予想と違う、という場合もそのまま使い、気付いた点はユーザーに伝える）。
4. `no-journal` と `not-found` の文献だけ、次の順で補う。補ったエントリには、エントリの直前に `% source: Crossref（INSPIRE に掲載誌の情報なし）` のように出典をコメントで書く。
   - DOI が分かっていれば Crossref（`https://api.crossref.org/works/<DOI>`）の container-title・volume・page・published で掲載誌の情報を足す。
   - Crossref でも確認できなければ arXiv のみのエントリ（`@misc`）にし、「掲載誌は未確認」と報告する。
5. コンパイルして文献表を見る（手順は .claude/rules/ の論文の規則）。

## 識別子が分からないとき

INSPIRE の検索 API を題名や著者で引く。

```bash
curl -s -G "https://inspirehep.net/api/literature" --data-urlencode 'q=a kameyama and t "refined large N"' --data-urlencode "format=bibtex" --data-urlencode "size=5"
```

- 著者は `a 姓`、題名の語句は `t "..."`、arXiv 番号は `arxiv:番号`、DOI は `doi:番号`。
- ブラウザで見るなら `https://inspirehep.net/literature?q=...`。各記録の「cite」から同じ BibTeX が取れる。
- 複数当たったら、著者・題名・年で 1 件に絞ってから、その arXiv 番号か DOI をスクリプトに渡す。

## 注意

- INSPIRE は hep 系の文献が中心で、数学の文献（Macdonald の本、Koornwinder 多項式の論文など）は記録が無いことがある。無い文献は手順 4 で扱い、INSPIRE に無いことを報告する。
- 外部に送るのは arXiv 番号・DOI・題名の検索語だけにする。ユーザーの名前や連絡先を User-Agent などに入れない。
