"""文章を直す前後で、数式・表・引用・参照・ラベル・数値・定理環境・見出しが変わっていないことを確かめる。

使い方：python3 check_invariants.py <.tex かディレクトリ> ...            # git の HEAD と作業ツリーを比べる
        python3 check_invariants.py <.tex かディレクトリ> ... --rev <rev>   # 指定したコミットと比べる
        python3 check_invariants.py --before <直す前> <直した後>           # 2 つのファイルかディレクトリを比べる

渡したファイル全体を 1 つの原稿として比べるので、文をファイル間で移しても差分にならない。
数式は空白と \\label、末尾の句読点を無視して比べる（式に句読点を付ける修正は差分にならない）。
差分があれば「- 直す前にだけあるもの」「+ 直した後にだけあるもの」を出し、終了コード 1。
文体だけを直したのなら差分は 0 件のはずである。構成や引用を意図して変えた場合は、差分の各行を説明できることを確かめる。
"""
import re
import subprocess
import sys
from collections import Counter

from latex_prose import CITE, LABEL, MATH, REF, mask, strip_comments, tex_files

THEOREMS = "theorem|proposition|lemma|corollary|conjecture|definition|remark|hypothesis|result|claim|example|proof"


def normalize_math(segment):
    s = re.sub(r"^(\\begin\{[^}]*\}|\\\[|\$\$|\\\(|\$)|(\\end\{[^}]*\}|\\\]|\$\$|\\\)|\$)$", "", segment)
    s = re.sub(r"\s+", "", LABEL.sub("", s))
    return re.sub(r"(?:\\[,;:!]|\\q?quad|~)*[.,;]((?:\\end\{[a-z*]+\})*)$", r"\1", s)


def invariants(text):
    text = strip_comments(text)
    return {
        "数式・表": Counter(normalize_math(m.group(0)) for m in MATH.finditer(text)),
        "引用": Counter(key.strip() for m in CITE.finditer(text) for key in m.group(1).split(",")),
        "参照": Counter(key.strip() for m in REF.finditer(text) for key in m.group(2).split(",")),
        "ラベル": Counter(LABEL.findall(text)),
        "本文の数値": Counter(re.findall(r"(?<![\w.\\])\d+(?:[.,]\d+)*(?!\w)", mask(text))),
        "定理環境": re.findall(r"\\begin\{(" + THEOREMS + r")\*?\}", text),
        "見出し": [re.sub(r"\s+", " ", m.group(0)) for m in re.finditer(r"\\(?:sub)*section\*?\{[^}]*\}", MATH.sub("[MATH]", text))],
    }


def read_tree(paths):
    return "\n".join(p.read_text(encoding="utf-8") for p in tex_files(paths))


def read_git(rev, paths):
    names = subprocess.run(["git", "ls-tree", "-r", "--name-only", rev, "--", *paths],
                           check=True, capture_output=True, text=True).stdout.split("\n")
    return "\n".join(
        subprocess.run(["git", "show", f"{rev}:./{name}"], check=True, capture_output=True, text=True).stdout
        for name in sorted(names) if name.endswith(".tex")
    )


def main(argv):
    paths, rev, before = [], "HEAD", None
    it = iter(argv)
    for a in it:
        if a == "--rev":
            rev = next(it)
        elif a == "--before":
            before = next(it)
        else:
            paths.append(a)
    if not paths:
        print(__doc__)
        return 2

    old = invariants(read_tree([before]) if before else read_git(rev, paths))
    new = invariants(read_tree(paths))
    differs = 0
    for name in old:
        a, b = old[name], new[name]
        if a == b:
            print(f"{name}: 一致（{sum(a.values()) if isinstance(a, Counter) else len(a)} 個）")
            continue
        differs += 1
        removed, added = Counter(a) - Counter(b), Counter(b) - Counter(a)
        print(f"{name}: 不一致" + ("" if removed or added else "（個数は同じで順序が違う）"))
        for sign, items in (("-", removed), ("+", added)):
            for item, count in sorted(items.items()):
                print(f"  {sign} {item[:120]}" + (f"  ×{count}" if count > 1 else ""))
    print("不変量はすべて一致" if not differs else f"{differs} 種類の不変量が変わっている")
    return 1 if differs else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
