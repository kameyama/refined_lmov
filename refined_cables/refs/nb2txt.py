"""Mathematica .nb の Input/Output/Text セルを線形テキストに変換する（読解用の簡易ツール）。"""
import re
import sys


def tokenize(s):
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c == '"':
            j = i + 1
            buf = []
            while s[j] != '"':
                if s[j] == '\\':
                    buf.append(s[j:j + 2])
                    j += 2
                else:
                    buf.append(s[j])
                    j += 1
            yield ('str', ''.join(buf))
            i = j + 1
        elif c in '[]{},':
            yield (c, c)
            i += 1
        elif c.isspace():
            i += 1
        elif s.startswith('(*', i):
            i = s.index('*)', i) + 2
        else:
            j = i
            while j < n and s[j] not in '[]{},"' and not s[j].isspace():
                j += 1
            yield ('sym', s[i:j])
            i = j


def parse(tokens):
    toks = list(tokens)
    pos = 0

    def expr():
        nonlocal pos
        kind, val = toks[pos]
        pos += 1
        if kind == '{':
            items = []
            while toks[pos][0] != '}':
                items.append(expr())
                if toks[pos][0] == ',':
                    pos += 1
            pos += 1
            return ('list', items)
        node = (kind, val)
        while pos < len(toks) and toks[pos][0] == '[':
            pos += 1
            args = []
            while toks[pos][0] != ']':
                args.append(expr())
                if toks[pos][0] == ',':
                    pos += 1
            pos += 1
            node = ('call', node, args)
        return node

    out = []
    while pos < len(toks):
        out.append(expr())
        if pos < len(toks) and toks[pos][0] == ',':
            pos += 1
    return out


def head(node):
    return node[1][1] if node[0] == 'call' and node[1][0] == 'sym' else None


def render(node):
    k = node[0]
    if k == 'str':
        s = node[1]
        s = s.replace('\\[Equal]', '==').replace('\\[Rule]', '->').replace('\\n', '\n')
        s = re.sub(r'\\\[(\w+)\]', r'\\\1', s)
        return s.replace('\\<', '').replace('\\>', '').replace('\\"', '"')
    if k == 'sym':
        return node[1]
    if k == 'list':
        return ''.join(render(x) for x in node[1])
    h = head(node)
    a = node[2]
    if h == 'RowBox':
        return render(a[0])
    if h == 'SuperscriptBox':
        return f"{render(a[0])}^({render(a[1])})"
    if h == 'SubscriptBox':
        return f"{render(a[0])}_{render(a[1])}"
    if h == 'FractionBox':
        return f"(({render(a[0])})/({render(a[1])}))"
    if h == 'SqrtBox':
        return f"Sqrt[{render(a[0])}]"
    if h in ('UnderoverscriptBox', 'SubsuperscriptBox'):
        return f"{render(a[0])}_({render(a[1])})^({render(a[2])})"
    if h in ('UnderscriptBox', 'OverscriptBox'):
        return f"{render(a[0])}_({render(a[1])})"
    if h in ('StyleBox', 'TagBox', 'InterpretationBox', 'FormBox', 'TooltipBox', 'DynamicBox', 'ButtonBox', 'BoxData', 'TextData', 'Cell', 'GridBox'):
        if h == 'GridBox':
            return '\n'.join(' | '.join(render(c) for c in row[1]) for row in a[0][1])
        return render(a[0]) if a else ''
    return ''


def main(path):
    text = open(path, encoding='utf-8', errors='replace').read()
    start = text.find('Notebook[')
    tree = parse(tokenize(text[start:]))
    nb = tree[0]

    def walk(node):
        if node[0] == 'call':
            h = head(node)
            if h == 'Cell' and len(node[2]) >= 2 and node[2][1][0] == 'str':
                style = node[2][1][1]
                if style in ('Input', 'Output', 'Text', 'Section', 'Subsection', 'Subsubsection', 'Title', 'Code'):
                    body = render(node[2][0])
                    print(f"=== [{style}]\n{body}\n")
                    return
            for x in node[2]:
                walk(x)
        elif node[0] == 'list':
            for x in node[1]:
                walk(x)

    walk(nb)


if __name__ == '__main__':
    sys.setrecursionlimit(100000)
    main(sys.argv[1])
