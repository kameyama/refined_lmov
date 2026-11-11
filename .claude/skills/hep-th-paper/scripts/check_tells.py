"""LaTeX 原稿から「AI が書いた」と読まれるパターンを機械的に拾う。

使い方：python3 check_tells.py <.tex ファイルかディレクトリ> ... [--bib refs.bib ...]
ディレクトリは下の *.tex をすべて読む。--bib を省くと、原稿の \\bibliography / \\addbibresource から探す。
出力は「ファイル:行: [hard|soft] 内容」。hard が 1 件でもあれば終了コード 1。
行末に % tells-ok と書いた行は hard から外す。パターンの説明は references/ai-tells.md。
"""
import re
import sys
from pathlib import Path

HARD = [
    (r"\\begin\{(itemize|enumerate|description)\}", "箇条書き（散文に直す。例外は ai-tells.md）"),
    (r"\\begin\{(tcolorbox|mdframed|framed)\}|\\fbox\{", "囲み（主要な予想 1 つの枠以外は外す）"),
    (r"^\s*(\\noindent\s*)?(\\textbf\{|\{\\bf\s)[^}]{0,60}[:.]\s*\}", "太字ラベルで始まる段落"),
    (r"\\paragraph\{", "\\paragraph の小見出し"),
    (r"\b(Key (insight|idea|point|takeaway)s?|Main takeaway|Takeaways?:|Summary:|TL;DR)", "要点ラベル"),
    (r"\*\*[^*\n]+\*\*|(?<!`)`[^`'\n]+`(?!`)|^#{1,4} ", "Markdown 記法"),
    (r"[✓✔✗✘→⇒←↔≈≤≥×·•…]|[\U0001F300-\U0001FAFF\u2600-\u27BF]", "Unicode の記号・絵文字（LaTeX のコマンドで書く）"),
    (r"\b(PASS(ED)?|FAIL(ED)?)\b|\b[\w/-]+\.(py|sage|wls|ipynb)\b|\\checkmark", "検証ログ（PASS/FAIL、スクリプト名、チェックマーク）"),
    (r"\b(Here is|Here's|Let me|As requested|I hope this|Certainly[,!]|Sure[,!])", "LLM のメタコメント"),
    (r"\\(sub)*section\*?\{[^}]*\w:\s[^}]*\}", "見出しのコロン構文"),
]
SOFT_WORDS = (
    r"delv\w+|crucial\w*|pivotal|notably|importantly|interestingly|remarkabl\w+|"
    r"it is worth noting|it is important to note|it should be noted|comprehensive|novel|robust\w*|"
    r"seamless\w*|intricate|elegant\w*|striking\w*|powerful|landscape|tapestry|shed(s|ding)? light|"
    r"pav(e|es|ing) the way|underscor\w+|showcas\w+|leverag\w+|testament|plays? a key role|deep connection"
)
SOFT = [
    (rf"\b({SOFT_WORDS})\b", "語彙"),
    (r"\bnot only\b", "not only ... but also"),
    (r"\b(In summary|To summarize|In conclusion|Overall|Taken together|This completes)\b", "まとめ文"),
    (r"\b(verif\w+|confirm\w+|rigorous\w*|establish\w*)\b", "証明か範囲つきの計算に対応しているか確認"),
    (r"\\begin\{eqnarray\*?\}", "eqnarray（align を使う）"),
]


def strip_comment(line):
    return re.sub(r"(?<!\\)%.*", "", line)


def tex_files(paths):
    for p in map(Path, paths):
        yield from sorted(p.rglob("*.tex")) if p.is_dir() else [p]


def main(argv):
    args, bibs = [], []
    it = iter(argv)
    for a in it:
        if a == "--bib":
            bibs.append(Path(next(it)))
        else:
            args.append(a)
    if not args:
        print(__doc__)
        return 2

    hard = soft = 0
    cited, remarks, in_this_section = {}, 0, 0
    for path in tex_files(args):
        lines = path.read_text(encoding="utf-8").splitlines()
        dashes = 0
        for n, raw in enumerate(lines, 1):
            line = strip_comment(raw)
            for m in re.finditer(r"\\(?:bibliography|addbibresource)\{([^}]+)\}", line):
                for name in m.group(1).split(","):
                    bib = path.parent / name.strip()
                    bibs.append(bib if bib.suffix == ".bib" else bib.with_name(bib.name + ".bib"))
            for m in re.finditer(r"\\[a-zA-Z]*cite[a-zA-Z]*\*?(?:\[[^\]]*\])*\{([^}]+)\}", line):
                for key in m.group(1).split(","):
                    cited.setdefault(key.strip(), f"{path}:{n}")
            remarks += len(re.findall(r"\\begin\{remark\}", line))
            in_this_section += len(re.findall(r"\bIn this (sub)?section,? we\b", line))
            if not line.strip():
                dashes = 0
            dashes += len(re.findall(r"(?<!-)---(?!-)|—", line))
            if dashes >= 2:
                print(f"{path}:{n}: [soft] 1 つの段落に em ダッシュが 2 つ以上")
                soft += 1
                dashes = 0
            for pat, msg in HARD:
                m = re.search(pat, line)
                if m and "tells-ok" not in raw:
                    print(f"{path}:{n}: [hard] {msg}: {m.group(0).strip()[:50]}")
                    hard += 1
            for pat, msg in SOFT:
                for m in re.finditer(pat, line, flags=re.I):
                    print(f"{path}:{n}: [soft] {msg}: {m.group(0)}")
                    soft += 1

    if remarks >= 5:
        print(f"[soft] Remark が {remarks} 個（5 個以上）")
        soft += 1
    if in_this_section >= 3:
        print(f"[soft] \"In this section, we\" が {in_this_section} 回")
        soft += 1
    if cited:
        if len(cited) < 20:
            print(f"[soft] 引用している文献が {len(cited)} 本（20 本未満）")
            soft += 1
        bib_keys = set()
        for bib in bibs:
            if bib.exists():
                bib_keys |= set(re.findall(r"^\s*@\w+\{\s*([^,\s]+)\s*,", bib.read_text(encoding="utf-8"), flags=re.M))
        if bib_keys:
            for key, where in cited.items():
                if key not in bib_keys:
                    print(f"{where}: [hard] \\cite のキーが .bib にない: {key}")
                    hard += 1
        else:
            print("[soft] .bib が見つからないので引用キーを照合していない（--bib で指定する）")
            soft += 1

    print(f"hard {hard} 件、soft {soft} 件、引用 {len(cited)} 本")
    return 1 if hard else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
