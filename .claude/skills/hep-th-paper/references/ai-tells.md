# 「AI が書いた」と読まれるパターン

hard は見つけたら必ず直す。soft は 1 件ずつ読み、不自然なものを直す（"robust" も "remarkable" も、人が書いた論文に普通に現れる。問題は密度と、中身のない使い方である）。content は機械では判定しきれないので、読んで判断する。

scripts/check_tells.py が機械的に拾えるものには ✔ を付けた。行末に `% tells-ok` と書いた行は、スクリプトが hard から外す（理由のある例外にだけ使う）。

## hard：構造

- ✔ 本文の `itemize` / `enumerate` / `description`。散文に直す。例外は Discussion の未解決問題の短い一覧 1 か所と、対称性や性質を式つきで並べる短い一覧（模範論文にもある）。
- ✔ `\textbf{Key insight:}`、`\textbf{Main result.}` のような太字ラベルで始まる段落。`\paragraph{...}` を小見出しとして並べること。
- ✔ tcolorbox、mdframed、framed、`\fbox` による「まとめ箱」「要点の囲み」。例外は主要な予想 1 つを囲む枠で、これは模範論文（arXiv:1703.05408）にもある。
- ✔ "Key insight"、"Key idea"、"Main takeaway"、"Takeaway"、"Summary:"、"TL;DR"。
- ✔ Markdown と Unicode の混入：`**bold**`、バッククォートで囲んだ語、行頭の `#`、✓ ✗ → ⇒、絵文字。
- ✔ 検証ログ：PASS / FAIL、スクリプトのファイル名（`.py`、`.sage`、`.wls`）、コードの出力の貼り付け、実行時間の羅列、"Verified numerically for N=2,...,8 ✓"。
- ✔ LLM のメタコメント："Here is ..."、"Let me ..."、"As requested"、"I hope this helps"、"Certainly"。arXiv が未確認の LLM 出力の証拠として名指ししている。
- ✔ 見出しのコロン構文（"Main Result: A Refined Duality"）。
- ✔ `\cite` のキーが .bib にない。

## hard：主張

- 証明のない主張が Theorem / Proposition になっている。
- "verified"、"confirmed"、"proved"、"rigorous"、"established" に、証明も範囲つきの計算も対応していない（✔ 語の出現だけを soft として出す。対応の有無は読んで判断する）。
- 推測を導出のように書いている（"it follows that"、"one can show" の先に導出がない）。
- 既知の結果を、この論文の成果のように書いている。

## soft：語彙

✔ delve, crucial, pivotal, notably, importantly, interestingly, remarkably, it is worth noting, it is important to note, it should be noted, comprehensive, novel, robust, seamless, intricate, elegant, striking, powerful, rich, landscape, tapestry, shed light, pave the way, underscore, showcase, leverage, a testament to, plays a key role, deep connection。

文頭の "Importantly," "Notably," "Interestingly," は、ほとんどの場合そのまま削れる。

## soft：文と段落の形

- ✔ 1 つの段落に em ダッシュ（`---`）が 2 つ以上。カンマ、括弧、文の分割に直す。
- ✔ "not only ... but also"、"not X but Y" の対比の繰り返し。
- 三点列挙の連続（"A, B, and C" が続く）。項目の数は内容で決める。
- ✔ 節の冒頭の "In this section, we ..." が 3 回以上。
- ✔ 節末・段落末のまとめ文（"In summary,"、"To summarize,"、"Overall,"、"Taken together,"、"This completes ..."）。
- どの段落も同じ長さ・同じ組み立て（主張 → 説明 → まとめ）。
- 短い文の連続と、同じことの言い直し。同じ事実を abstract・introduction・節の冒頭・節末・discussion で 5 回述べない。
- ✔ Remark が 5 個以上。
- ✔ `eqnarray`。

## content：専門家が強く反応するもの

- ✔ 引用が少ない（目安は 20 本未満）。分野の基本文献に触れていない。
- 予想や公式に出典がない。
- 標準的な定義の長い再説明と、固有の記法・規約の説明不足が同居している。
- 規約の揺れ（(q, t) の対応、正規化、framing、reduced / unreduced）。編集中に教科書の標準形へ戻してしまうのが典型である。
- 計算の範囲が書かれていない（どの結び目、どの彩色、どの次数まで）。
- 途中で使われなくなった節や記号、同じ内容の重複（"zombie sections"）。
- 物理的な動機づけがなく、式の列挙になっている。逆に、動機づけの文が一般論で、この論文の設定に触れていない。

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
\eqref{eq:framed-f}. We then compute the refined LMOV invariants of the twelve knots listed in
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

格上げされた主張：

```latex
% 直す前
\begin{proposition} The invariants $\widehat N_{\rho,g,\beta,J_r}$ are non-negative integers. \end{proposition}

% 直した後
\begin{conjecture} For every knot $K$ and every $\rho$, the invariants
$\widehat N_{\rho,g,\beta,J_r}(K)$ are non-negative integers. \end{conjecture}
We have confirmed this for the knots in table~\ref{tab:knots} with $|\rho| \le 3$.
```
