"""INSPIRE-HEP から BibTeX を取り、引用キーを差し替えて shorthand を足す。

使い方：python3 inspire_bibtex.py KEY=識別子 [KEY=識別子 ...]
識別子は arXiv 番号（1703.05408、hep-th/0405083）か DOI（10. で始まる）。
標準出力に BibTeX、標準エラーに文献ごとの結果（ok / no-journal / not-found / ambiguous）。
"""
import re
import sys
import time
import urllib.parse
import urllib.request

API = "https://inspirehep.net/api/literature?format=bibtex&size=3&q="


def fetch(ident):
    query = ("doi:" if ident.startswith("10.") else "arxiv:") + ident
    with urllib.request.urlopen(API + urllib.parse.quote(query), timeout=40) as res:
        return res.read().decode("utf-8")


def rekey(entry, key):
    entry = re.sub(r"^(@\w+)\{[^,]+,", lambda m: f"{m.group(1)}{{{key},\n    shorthand = {{{key}}},", entry.strip(), count=1)
    return entry + "\n"


def main(args):
    failed = False
    for arg in args:
        key, ident = arg.split("=", 1)
        text = fetch(ident)
        count = len(re.findall(r"^@", text, flags=re.M))
        if count == 0:
            print(f"{key}: not-found ({ident})", file=sys.stderr)
            failed = True
        elif count > 1:
            print(f"{key}: ambiguous ({ident}, {count} 件)", file=sys.stderr)
            failed = True
        else:
            has_journal = re.search(r"^\s*(journal|booktitle|publisher)\s*=", text, flags=re.M)
            print(f"{key}: {'ok' if has_journal else 'no-journal'} ({ident})", file=sys.stderr)
            print(rekey(text, key))
        time.sleep(1)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
