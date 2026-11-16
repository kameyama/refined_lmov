"""LaTeX 原稿から「AI が書いた」と読まれるパターンの候補を機械的に拾う。

使い方：python3 check_tells.py <.tex ファイルかディレクトリ> ... [--bib refs.bib ...] [--voice <著者の過去論文>] [--all]
ディレクトリは下の *.tex をすべて読む。--bib を省くと、原稿の \\bibliography / \\addbibresource から探す。
--voice に著者が自分で書いた論文の .tex を渡すと、著者が普段使う語を語彙の指摘から外し、文長を比べる。
--all は、段落に 1 語だけ現れた語彙も 1 件ずつ出す（既定では末尾の集計にだけ出す）。

出力は「ファイル:行: [hard|soft 層] 内容」。層は 誠実性 > 構造 > 語彙 > 表層 の順に重い。
hard が 1 件でもあれば終了コード 1。行末に % tells-ok と書いた行は hard から外す。
拾うのは候補であって判定ではない。パターンの説明と、直さなくてよい場合は references/ai-tells.md。
語彙・文の形は、数式・引用・参照・表・コメントを伏せた本文だけを見る。
"""
import re
import statistics
import sys
from collections import Counter
from pathlib import Path

from latex_prose import CITE, mask, paragraphs, strip_comments, tex_files

LAYERS = ["誠実性", "構造", "語彙", "表層"]

# 行ごとに見るもの（コメントだけを除いた行）。
HARD = [
    ("構造", r"\\begin\{(itemize|enumerate|description)\}", "箇条書き（散文に直す。例外は ai-tells.md）"),
    ("構造", r"\\begin\{(tcolorbox|mdframed|framed)\}|\\fbox\{", "囲み（主要な予想 1 つの枠以外は外す）"),
    ("構造", r"^\s*(\\noindent\s*)?(\\textbf\{|\{\\bf\s)[^}]{0,60}[:.]\s*\}", "太字ラベルで始まる段落"),
    ("構造", r"\\paragraph\{", "\\paragraph の小見出し"),
    ("構造", r"\b(Key (insight|idea|point|takeaway)s?|Main takeaway|Takeaways?:|Summary:|TL;DR)", "要点ラベル"),
    ("構造", r"\b(PASS(ED)?|FAIL(ED)?)\b|\b[\w/-]+\.(py|sage|wls|ipynb)\b|\\checkmark", "検証ログ（PASS/FAIL、スクリプト名、チェックマーク）"),
    ("誠実性", r"\b(Here is|Here's|Let me|As requested|I hope this|Certainly[,!]|Sure[,!])", "LLM のメタコメント"),
    ("表層", r"\*\*[^*\n]+\*\*|(?<!`)`[^`'\n]+`(?!`)|^#{1,4} ", "Markdown 記法"),
    ("表層", r"[✓✔✗✘→⇒←↔≈≤≥×·•…“”‘’]|[\U0001F300-\U0001FAFF\u2600-\u27BF]", "Unicode の記号・引用符・絵文字（LaTeX のコマンドで書く）"),
    ("表層", r"\\(sub)*section\*?\{[^}]*\w:\s[^}]*\}", "見出しのコロン構文"),
]

# 段落ごとに見るもの（伏せた本文）。
SOFT = [
    ("誠実性", r"(?i)\b(verif\w+|confirm\w+|rigorous\w*|establish\w*)\b", "証明か範囲つきの計算に対応しているか確認"),
    ("構造", r"(?i)\bnot only\b", "not only ... but also"),
    ("構造", r"(?:^|(?<=[.!?] ))(In summary|To summarize|In conclusion|Overall|Taken together|All in all),", "まとめ文"),
    ("構造", r", (highlighting|underscoring|showcasing|emphasizing|paving the way|opening (the door|new avenues)|shedding light|offering (new )?insights?)\b"
            r"|\bThis (highlights|underscores|showcases)\b", "意義づけで終わる文"),
    ("構造", r"(?i)\b(may|might|could) (potentially|possibly|perhaps)\b|\b(potentially|possibly) (may|might|could)\b", "ヘッジの重なり"),
    ("構造", r"\bIn recent years\b|\b(attracted|garnered|received|drawn) (much|significant|considerable|great|a lot of|increasing) (attention|interest)\b", "内容のない導入"),
    ("構造", r"(?i)\b(main|key|our) contributions? (are|is|can be|of this)\b", "寄与の列挙（散文で結果を述べる）"),
]
WELL_KNOWN = re.compile(r"(?i)\b(it is (well[- ])?known that|as is well known|it is widely (believed|accepted|known))[^.!?]*")
TRIAD = re.compile(r"\b[a-z-]+, [a-z-]+(?: [a-z-]+)?,? and [a-z-]+")  # 小文字だけ：著者名の並びを拾わない
CONNECTIVE = re.compile(r"^(Moreover|Furthermore|Additionally|In addition),")

# 語彙。1 語では指摘しない。同じ段落に 2 語以上あるときだけ出す。
WORDS = [
    r"delv(e|es|ed|ing)", r"crucial\w*", r"pivotal", r"notably", r"importantly", r"interestingly", r"remarkabl\w+",
    r"it is worth noting|it is important to note|it should be noted|worth mentioning", r"comprehensive\w*", r"novel",
    r"robust\w*", r"seamless\w*", r"intricate\w*", r"elegant\w*", r"striking\w*", r"powerful", r"groundbreaking",
    r"(?<!string )landscape(?! of (string |flux )?vacua)", r"tapestry", r"realm", r"multifaceted", r"meticulous\w*",
    r"garner(s|ed|ing)?(?!')", r"shed(s|ding)? light", r"pav(e|es|ed|ing) the way", r"underscor\w+", r"showcas\w+", r"leverag\w+",
    r"testament", r"plays? an? (key|crucial|pivotal|central) role", r"deep connections?",
]
WORD_RES = [re.compile(rf"(?i)\b({w})\b") for w in WORDS]


def abstracts(text):
    """abstract 環境と \\abstract{...} の (開始位置, 中身) を返す。"""
    for m in re.finditer(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", text, re.S):
        yield m.start(1), m.group(1)
    for m in re.finditer(r"\\abstract\s*\{", text):
        depth, i = 1, m.end()
        while i < len(text) and depth:
            depth += {"{": 1, "}": -1}.get(text[i], 0)
            i += 1
        yield m.end(), text[m.end():i - 1]


def read(paths):
    """原稿を読み、指摘の候補と集計を返す。"""
    r = {"found": [], "bibs": [], "cited": {}, "words": 0, "sentences": [], "vocab": [],
         "remarks": 0, "in_this_section": 0, "connectives": 0}
    add = lambda path, n, kind, layer, msg: r["found"].append((str(path), n, kind, layer, msg))
    for path in tex_files(paths):
        raw = path.read_text(encoding="utf-8")
        text = strip_comments(raw)
        for m in re.finditer(r"\\(?:bibliography|addbibresource)\{([^}]+)\}", text):
            for name in m.group(1).split(","):
                bib = path.parent / name.strip()
                r["bibs"].append(bib if bib.suffix == ".bib" else bib.with_name(bib.name + ".bib"))
        for m in CITE.finditer(text):
            for key in m.group(1).split(","):
                r["cited"].setdefault(key.strip(), (path, text.count("\n", 0, m.start()) + 1))
        r["remarks"] += len(re.findall(r"\\begin\{remark\}", text))
        for start, body in abstracts(text):
            for m in CITE.finditer(body):
                add(path, text.count("\n", 0, start + m.start()) + 1, "hard", "構造", f"abstract の中の \\cite（著者名だけを書く）: {m.group(1)}")
            size = len(re.sub(r"\s+", " ", body).strip())
            if size > 1920:
                add(path, text.count("\n", 0, start) + 1, "soft", "構造", f"abstract が {size} 文字（arXiv の abstract 欄は 1,920 文字まで）")
        for n, (raw_line, line) in enumerate(zip(raw.splitlines(), text.splitlines()), 1):
            for layer, pat, msg in HARD:
                m = re.search(pat, line)
                if m and "tells-ok" not in raw_line:
                    add(path, n, "hard", layer, f"{msg}: {m.group(0).strip()[:50]}")
            if re.search(r"\\begin\{eqnarray\*?\}", line):
                add(path, n, "soft", "表層", "eqnarray（align を使う）")

        for para_index, (line_of, flat) in enumerate(paragraphs(mask(raw))):
            plain = re.sub(r"\\[a-zA-Z]+\*?", " ", flat)
            r["words"] += len(re.findall(r"[A-Za-z][A-Za-z'-]*", re.sub(r"\[(MATH|CIT|REF)\]", " ", plain)))
            r["sentences"] += [n for n in (len(s.split()) for s in re.split(r"(?<=[.!?]) +(?=[A-Z])", plain)) if n >= 3]
            r["in_this_section"] += len(re.findall(r"\bIn this (sub)?section,? we\b", flat))
            r["connectives"] += bool(CONNECTIVE.match(flat))
            for layer, pat, msg in SOFT:
                for m in re.finditer(pat, flat):
                    add(path, line_of(m.start()), "soft", layer, f"{msg}: {m.group(0).strip(' ,')}")
            for m in WELL_KNOWN.finditer(flat):
                if "[CIT]" not in m.group(0):
                    add(path, line_of(m.start()), "soft", "誠実性", f"出典のない「よく知られている」: {m.group(1)}")
            if len(re.findall(r"(?<!-)---(?!-)|—", flat)) >= 2:
                add(path, line_of(0), "soft", "表層", "1 つの段落に em ダッシュが 2 つ以上")
            triads = TRIAD.findall(flat)
            if len(triads) >= 2:
                add(path, line_of(0), "soft", "構造", f"1 つの段落に三点列挙が {len(triads)} 回: {triads[0]}")
            for i, pat in enumerate(WORD_RES):
                for m in pat.finditer(flat):
                    r["vocab"].append((i, str(path), para_index, line_of(m.start()), m.group(0).lower()))
    return r


def bib_entries(bibs):
    entries = {}
    for bib in bibs:
        if bib.exists():
            text = bib.read_text(encoding="utf-8")
            for m in re.finditer(r"^\s*@(\w+)\s*\{\s*([^,\s]+)\s*,(.*?)(?=^\s*@\w+\s*\{|\Z)", text, flags=re.M | re.S):
                entries[m.group(2)] = m.group(3)
    return entries


def per_thousand(count, words):
    return 1000 * count / max(words, 1)


def main(argv):
    args, bibs, voice_paths, show_all = [], [], [], False
    it = iter(argv)
    for a in it:
        if a == "--bib":
            bibs.append(Path(next(it)))
        elif a == "--voice":
            voice_paths.append(next(it))
        elif a == "--all":
            show_all = True
        else:
            args.append(a)
    if not args:
        print(__doc__)
        return 2

    r = read(args)
    voice = read(voice_paths) if voice_paths else None
    found = r["found"]
    add = lambda layer, msg, where=("", 0): found.append((str(where[0]), where[1], "soft", layer, msg))

    # 語彙：著者が過去の論文で同じくらい使っている語は外し、残りを段落ごとに数える。
    counts = Counter(i for i, *_ in r["vocab"])
    voice_counts = Counter(i for i, *_ in voice["vocab"]) if voice else Counter()
    usual = {i for i in counts
             if voice and per_thousand(counts[i], r["words"]) <= 2 * per_thousand(voice_counts[i], voice["words"])}
    by_paragraph = {}
    for i, path, para, line, word in r["vocab"]:
        if i not in usual:
            by_paragraph.setdefault((path, para), []).append((line, word))
    for (path, _), hits in by_paragraph.items():
        hits.sort()
        if len(hits) >= 2:
            found.append((path, hits[0][0], "soft", "語彙", f"同じ段落に {len(hits)} 語: " + ", ".join(w for _, w in hits)))
        elif show_all:
            found.append((path, hits[0][0], "soft", "語彙", f"単発: {hits[0][1]}"))

    if r["remarks"] >= 5:
        add("構造", f"Remark が {r['remarks']} 個（5 個以上）")
    if r["in_this_section"] >= 3:
        add("構造", f"\"In this section, we\" が {r['in_this_section']} 回")
    if r["connectives"] >= 3:
        add("構造", f"Moreover / Furthermore / Additionally / In addition で始まる段落が {r['connectives']} 個")
    cited = r["cited"]
    if cited:
        if len(cited) < 20:
            add("構造", f"引用している文献が {len(cited)} 本（20 本未満）")
        entries = bib_entries(bibs + r["bibs"])
        if not entries:
            add("誠実性", ".bib が見つからないので引用キーを照合していない（--bib で指定する）")
        for key, where in cited.items() if entries else []:
            if key not in entries:
                found.append((str(where[0]), where[1], "hard", "誠実性", f"\\cite のキーが .bib にない: {key}"))
            elif not re.search(r"(?im)^\s*(doi|eprint|isbn|url)\s*=", entries[key]):
                add("誠実性", f".bib のエントリに DOI も eprint もない（実在を原典で確認）: {key}", where)

    for path, n, kind, layer, msg in sorted(found, key=lambda f: (f[0] == "", f[0], f[1], LAYERS.index(f[3]))):
        print(f"{path}:{n}: [{kind} {layer}] {msg}" if path else f"[{kind} {layer}] {msg}")

    def describe(x):
        s = x["sentences"]
        rate = per_thousand(len(x["vocab"]), x["words"])
        cv = statistics.pstdev(s) / statistics.mean(s) if s else 0
        return f"{x['words']} 語、語彙 {len(x['vocab'])} 回（1,000 語あたり {rate:.2f}）、文長の平均 {statistics.mean(s) if s else 0:.1f} 語・変動係数 {cv:.2f}"

    print(f"本文: {describe(r)}")
    if voice:
        print(f"著者の過去論文: {describe(voice)}")
    for i, c in counts.most_common():
        note = f"、著者の過去論文 {voice_counts[i]} 回" + ("（普段の範囲なので指摘から外した）" if i in usual else "") if voice else ""
        print(f"  語彙 {WORDS[i][:40]}: {c} 回{note}")
    tally = Counter((kind, layer) for _, _, kind, layer, _ in found)
    print("層別: " + "、".join(f"{layer} hard {tally['hard', layer]}・soft {tally['soft', layer]}" for layer in LAYERS))
    hard = sum(kind == "hard" for _, _, kind, _, _ in found)
    print(f"hard {hard} 件、soft {len(found) - hard} 件、引用 {len(cited)} 本")
    return 1 if hard else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
