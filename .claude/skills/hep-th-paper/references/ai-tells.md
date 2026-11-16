# 「AI が書いた」と読まれるパターン

目次：読み方 / 誠実性 / 構造 / 語彙 / 表層 / 直さなくてよいもの / 直し方の例

## 読み方

ここにあるパターンは、見直す箇所を教える手がかりであって、判定ではない。人が書いた論文にも現れる。1 件ずつ該当する文を読み、「直さなくてよいもの」に当たらないかを確かめてから指摘にする。該当なしも正当な結果である。

パターンは重い順に 4 つの層に分ける。語彙は世代ごとに変わるが、構造の癖は残る。専門家が強く反応するのも上の 2 層である。

| 層 | 中身 | 扱い |
|---|---|---|
| 誠実性 | 引用の実在と正確さ、主張と根拠の対応、記号の整合、本文と式の一致 | 文章の修正としては直さない。著者が確かめる項目として報告する |
| 構造 | 箇条書き、太字ラベル、囲み、検証ログ、決まり型の文と段落 | 直す |
| 語彙 | LLM が過剰に使う語、自己評価の語 | 1 語では指摘しない。同じ段落での重なりと全体の密度で判断する |
| 表層 | em ダッシュ、Unicode、Markdown、見出しの書式 | 直すが、優先度は低い |

確度は 2 段階である。hard は見つけたら必ず直す（誠実性の hard は報告して解消する）。soft は読んで判断する。scripts/check_tells.py が機械的に拾えるものには ✔ を付けた。✔ のないものは、読まないと見つからない。行末に `% tells-ok` と書いた行は、スクリプトが hard から外す（理由のある例外にだけ使い、理由を同じコメントに書く）。

スクリプトは、コメント・数式・表・`\cite` や `\ref` の引数を伏せた本文だけで語彙と文の形を見る。

## 誠実性

- ✔ hard：`\cite` のキーが .bib にない。
- ✔ soft：引用している .bib のエントリに DOI も eprint も ISBN もない。実在と書誌を原典で確かめる（skill `inspire-bibtex`）。書籍は正当にこうなることがある。
- ✔ hard：LLM のメタコメント（"Here is ..."、"Let me ..."、"As requested"、"I hope this helps"、"Certainly"）。arXiv が未確認の LLM 出力の証拠として名指ししている。
- ✔ soft："verified"、"confirmed"、"rigorous"、"established"。証明か、範囲を明記した計算が原稿の中に対応しているかを読む。
- ✔ soft："It is well known that"、"It is known that"、"It is widely believed" の文に引用がない。引用を足すのではなく「要出典」として報告する。
- 証明のない主張が Theorem / Proposition になっている。
- 推測を導出のように書いている（"it follows that"、"one can show" の先に導出がない）。"always"、"universally"、"for all knots" が、確かめた範囲を超えている。
- 既知の結果を、この論文の成果のように書いている。
- 引用の内容が本文の説明と合わない（著者、年、何を示した論文か）。
- 記号の不整合：同じ量に 2 つの記号、定義より前の使用、マクロの定義と本文のずれ。規約の揺れ（(q, t) の対応、正規化、framing、reduced / unreduced）。編集中に教科書の標準形へ戻してしまうのが典型である。
- 本文と式の不一致：本文が "positive"、"vanishes"、"symmetric"、"polynomial" と書いている量が、式ではそうなっていない。
- 計算の範囲が書かれていない（どの結び目、どの彩色、どの次数まで）。

## 構造

原稿の組み立て：

- ✔ hard：本文の `itemize` / `enumerate` / `description`。散文に直す。
- ✔ hard：`\textbf{Key insight:}`、`\textbf{Main result.}` のような太字ラベルで始まる段落。`\paragraph{...}` を小見出しとして並べること。
- ✔ hard：tcolorbox、mdframed、framed、`\fbox` による「まとめ箱」「要点の囲み」。
- ✔ hard："Key insight"、"Key idea"、"Main takeaway"、"Takeaway"、"Summary:"、"TL;DR"。
- ✔ hard：検証ログ。PASS / FAIL、スクリプトのファイル名（`.py`、`.sage`、`.wls`）、コードの出力の貼り付け、実行時間の羅列、"Verified numerically for N=2,...,8 ✓"。
- ✔ hard：abstract の中の `\cite`。文献に触れるときは著者名だけを書く。
- ✔ soft：abstract が 1,920 文字を超える。arXiv の abstract 欄はこれより長いものを受け付けない（https://info.arxiv.org/help/prep.html で 2026-11-16 に確認）。文字数は TeX の記法ごと数えている。
- ✔ soft：内容のない導入（"In recent years, X has attracted significant attention"）。未解決の具体的な問題から書き始める。
- ✔ soft：寄与の列挙（"Our main contributions are ..."、番号つきの寄与）。結果を散文で述べる。
- ✔ soft：節の冒頭の "In this section, we ..." が 3 回以上。
- ✔ soft：Remark が 5 個以上。
- ✔ soft：引用が少ない（目安は 20 本未満）。分野の基本文献に触れていない。予想や公式に出典がない。
- abstract・introduction・節の冒頭・節末・discussion が同じ文の言い直しになっている。結論には限界と未解決点を書く。
- 標準的な定義の長い再説明と、固有の記法・規約の説明不足が同居している。
- 途中で使われなくなった節や記号、同じ内容の重複。
- 物理的な動機づけがなく、式の列挙になっている。逆に、動機づけの文が一般論で、この論文の設定に触れていない。

文と段落の決まり型：

- ✔ soft：段落末・節末のまとめ文（"In summary,"、"To summarize,"、"Overall,"、"Taken together,"）。情報を足さないなら削る。
- ✔ soft：意義づけで終わる文（", highlighting ..."、", underscoring ..."、", paving the way for ..."、"This highlights ..."）。削るか、具体的な帰結を書く。
- ✔ soft："not only ... but also"、"not just X but Y"、"This is not merely ..." の対比。Y を直接述べる。
- ✔ soft：1 つの段落に三点列挙（"A, B, and C"）が 2 回以上。項目の数は内容で決める。
- ✔ soft：ヘッジの重なり（"may potentially"、"could possibly"）。ヘッジは 1 つにし、根拠を添える。主張の強さは変えない。
- ✔ soft：Moreover / Furthermore / Additionally / In addition で始まる段落が 3 個以上。論理の関係（therefore、in contrast）が要るところだけ残す。
- どの段落も同じ長さ・同じ組み立て（主張 → 説明 → まとめ）。スクリプトは文長の平均と変動係数を出すが、閾値は置いていない。`--voice` で著者の過去論文の値と並べて見る。
- 短い文の連続と、同じことの言い直し。

## 語彙

✔ delve, crucial, pivotal, notably, importantly, interestingly, remarkably, it is worth noting, it is important to note, it should be noted, comprehensive, novel, robust, seamless, intricate, elegant, striking, powerful, groundbreaking, landscape, tapestry, realm, multifaceted, meticulous, garner, shed light, pave the way, underscore, showcase, leverage, a testament to, plays a key role, deep connection。

- スクリプトは、同じ段落に 2 語以上あるときだけ指摘し、1 語だけの出現は末尾の集計（語ごとの回数、1,000 語あたりの密度）に出す。`--all` で 1 語の出現も 1 件ずつ出る。
- 直すときは、平易な語にするか削る。自己評価の語（novel、powerful、elegant、groundbreaking）は、具体的な事実に置き換える（"a novel and powerful framework" ではなく、何を何に拡張した方法かを書く）。
- 文頭の "Importantly," "Notably," "Interestingly," は、ほとんどの場合そのまま削れる。
- abstract と title は 1 語でも読む。ここの自己評価は 1 回で目立つ。
- 語の一覧は、LLM の語彙の偏りを調べた研究（Kobak ほか、Liang ほか、Geng–Trotta、Juzek–Ward）と Wikipedia の "Signs of AI writing" に挙がる語から取った。ユーザーの調査メモ（2026-11）からの孫引きで、原典は確認していない。モデルの世代が変わると古くなるので、語の一覧より構造の層を優先する。

## 表層

- ✔ soft：1 つの段落に em ダッシュ（`---`）が 2 つ以上。カンマ、括弧、コロン、文の分割に直す。
- ✔ hard：Markdown と Unicode の混入。`**bold**`、バッククォートで囲んだ語、行頭の `#`、✓ ✗ → ⇒、Unicode の引用符（“ ” ‘ ’）、絵文字。
- ✔ hard：見出しのコロン構文（"Main Result: A Refined Duality"）。
- ✔ soft：`eqnarray`。
- 見出しの乱発（`\subsubsection`）、見出しのタイトルケースの揺れ。

## 直さなくてよいもの

指摘にする前に、ここと模範論文（exemplars.md）に照らす。当たるものは指摘から外すか、「直さなくてよい可能性」として報告に書く。

著者の語彙：著者自身の arXiv:1703.05408（2017 年。LLM 以前）をこのスクリプトで測ると、本文 8,468 語に語彙の一覧の語が 7 回（1,000 語あたり 0.83）現れる。remarkable / remarkably が 3 回、striking が 2 回、"worth mentioning" と groundbreaking が 1 回ずつで、同じ段落に 2 語が重なる箇所が 2 つある。"In this section, we" は 5 回、"not only" は 1 回、出典のない "It is known that" は 2 回ある（2026-11-16 に測定）。つまり、これらは著者の普段の書き方の範囲にある。`--voice <著者の過去論文>` を付けると、原稿での頻度が過去論文の 2 倍以下の語は指摘から外れる。

分野の用語と言い回し：

- landscape（string landscape、landscape of vacua）、swampland、threefold（Calabi–Yau threefold）、enhance（symmetry enhancement。一覧には入れていない）。
- robust（変形に対して robust という技術的な意味）、crucial（"crucial for anomaly cancellation" のように論理的に必須の場合）、novel（既存の名称の一部）、rich（"rich structure"）。
- 人名と重なる語（Garner's algorithm）。
- overall が形容詞の場合（overall sign、overall factor、overall normalization）。
- "it is easy to see"、"one finds"、"we note that"、"let us"、"we shall"、"To this end"、"It is known that" に引用が付いたもの、"This completes the proof."。
- 単著でも "we"。計算と手順の説明の受動態。
- 定理環境の中の定型表現。引用符の中の他人の文。

構造の例外（模範論文にもある）：

- Discussion の未解決問題の短い一覧 1 か所。対称性や性質を式つきで並べ、本文がラベルで参照する短い一覧。
- 主要な予想 1 つを囲む枠（arXiv:1703.05408 にある）。
- 実際に 3 つある対象の列挙（著者名、結び目、変数）。
- 表の空欄の `---`。

## 直し方の例

箇条書きの寄与：

```latex
% 直す前
Our main contributions are:
\begin{itemize}
  \item \textbf{Framed duality:} We extend the refined large $N$ duality to framed knots.
  \item \textbf{Positivity:} We verify positivity for 12 knots.
\end{itemize}

% 直した後
In this paper we extend the refined large $N$ duality of \cite{Kameyama:2017ryw} to framed
knots. The framing enters through ..., and the reformulated invariants are given by
\eqref{eq:framed-f}. We then compute the refined LMOV invariants of the 12 knots listed in
table~\ref{tab:knots} for all colors with up to three boxes, and find that they are
non-negative integers in every case.
```

検証ログ：

```latex
% 直す前
\textbf{Verification.} We verified numerically that (3.5) holds for $N=2,\dots,8$ $\checkmark$.

% 直した後
We have checked \eqref{eq:main} for $2 \le N \le 8$ by comparing all coefficients of both
sides; the results are collected in appendix~\ref{app:data}.
```

格上げされた主張（誠実性。claims.md の格付けか著者の判断に従って直す）：

```latex
% 直す前
\begin{proposition} The invariants $\widehat N_{\rho,g,\beta,J_r}$ are non-negative integers. \end{proposition}

% 直した後
\begin{conjecture} For every knot $K$ and every $\rho$, the invariants
$\widehat N_{\rho,g,\beta,J_r}(K)$ are non-negative integers. \end{conjecture}
We have confirmed this for the knots in table~\ref{tab:knots} with $|\rho| \le 3$.
```

語彙と決まり型（文体だけの修正。式・引用・数値は同じ）：

```latex
% 直す前
We delve into the intricate interplay between $f_\lambda$ and the BPS spectrum, showcasing a
novel and robust framework. Overall, this is not only consistent with \cite{KN17} but also may
potentially pave the way for future studies.

% 直した後
We study how $f_\lambda$ is related to the BPS spectrum. The result is consistent with \cite{KN17}.
```

検出から照合までを通した例は [worked-example.md](worked-example.md) にある。
