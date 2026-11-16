"""LaTeX 原稿を本文と「変えてはいけないもの」に分ける共通部品。

check_tells.py（検出）と check_invariants.py（修正前後の照合）が使う。
"""
import re
from pathlib import Path

# 本文として読まない環境。数式のほか、表・図・コードも中身ごと伏せる。
ENVS = "equation|align|alignat|flalign|gather|multline|eqnarray|displaymath|math|verbatim|lstlisting|tikzpicture|tabular|longtable"
MATH = re.compile(
    r"\\begin\{(" + ENVS + r")(\*?)\}.*?\\end\{\1\2\}"
    r"|(?<!\\)\\\[.*?(?<!\\)\\\]"
    r"|(?<!\\)\$\$.*?(?<!\\)\$\$"
    r"|(?<!\\)\\\(.*?(?<!\\)\\\)"
    r"|(?<!\\)\$.*?(?<!\\)\$",
    re.S,
)
CITE = re.compile(r"\\[a-zA-Z]*cite[a-zA-Z]*\*?(?:\[[^\]]*\])*\{([^}]*)\}")
REF = re.compile(r"\\(eqref|[cC]ref|autoref|pageref|ref)\*?\{([^}]*)\}")
LABEL = re.compile(r"\\label\{([^}]*)\}")
OTHER_ARGS = re.compile(
    r"\\(?:url|href|input|include|includegraphics|bibliography|addbibresource|bibliographystyle)\*?(?:\[[^\]]*\])*\{[^}]*\}"
)
FILL = "\x00"  # 伏せた行の目印。空行と区別して、段落が切れないようにする。


def tex_files(paths):
    for p in map(Path, paths):
        yield from sorted(p.rglob("*.tex")) if p.is_dir() else [p]


def strip_comments(text):
    return re.sub(r"(?<!\\)%.*", "", text)


def _placeholder(tag):
    return lambda m: tag + ("\n" + FILL) * m.group(0).count("\n")


def mask(text):
    """コメント・プリアンブル・数式・引用と参照の引数を伏せた本文を返す。行数は変えない。"""
    text = re.sub(r"(?m)^[ \t]*%.*$", FILL, text)
    text = strip_comments(text)
    head, sep, body = text.partition(r"\begin{document}")
    if sep:
        text = "\n" * head.count("\n") + body
    text = MATH.sub(_placeholder("[MATH]"), text)
    text = CITE.sub(_placeholder("[CIT]"), text)
    for pat in (REF, LABEL, OTHER_ARGS):
        text = pat.sub(_placeholder("[REF]"), text)
    return re.sub(r"\\(?:begin|end)\{[^}]*\}", "", text)


def paragraphs(masked):
    """(行番号を返す関数, 1 行にまとめた本文) を段落ごとに返す。"""
    start, buf = 1, []
    for n, line in enumerate(masked.split("\n") + [""], 1):
        if line.strip():
            if not buf:
                start = n
            buf.append(re.sub(r"\s+", " ", line.replace(FILL, "")).strip())
        elif buf:
            ends, pos = [], 0
            for piece in buf:
                pos += len(piece) + 1
                ends.append(pos)
            flat = " ".join(buf)
            if flat.strip():
                yield (lambda p, s=start, e=ends: s + sum(x <= p for x in e)), flat
            buf = []
